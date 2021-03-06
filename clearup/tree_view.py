from __future__ import print_function
import json

import os
import re
import subprocess
import time
import hashlib

import sys
from os.path import join, dirname, isfile
from collections import defaultdict
from os.path import abspath, join, dirname, splitext, basename

import six
from Bio import SeqIO, Phylo
from clearup.ultrafast.test_ultrafast_fp import plot_heatmap
from flask import Flask, render_template, abort, request, url_for
from manage import compare_pairwise

from ngs_utils.bed_utils import Region
from ngs_utils.file_utils import safe_mkdir, file_transaction, can_reuse, verify_file, which
from ngs_utils.file_utils import can_reuse, safe_mkdir
from ngs_utils import logger as log

import clearup
from clearup.genotype import build_tree
from clearup.model import Project, db, Sample, Run, get_or_create_run
from clearup.utils import read_fasta, FASTA_ID_PROJECT_SEPARATOR
from clearup import app, DATA_DIR

PROJ_COLORS = [
    '#000000',
    '#1f78b4',
    '#b2df8a',
    '#33a02c',
    '#fb9a99',
    '#e31a1c',
    '#fdbf6f',
    '#ff7f00',
    '#cab2d6',
    '#6a3d9a',
    '#ffff99',
    '#b15928',
]


def run_analysis_socket_handler(project_names_line):
    log.debug('Recieved request to start analysis for ' + project_names_line)
    ws = request.environ.get('wsgi.websocket', None)
    if not ws:
        raise RuntimeError('Environment lacks WSGI WebSocket support')

    def _run_cmd(cmdl):
        log.debug(cmdl)
        proc = subprocess.Popen(cmdl.split(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, env=os.environ)
        for stdout_line in iter(proc.stdout.readline, None):
            if not stdout_line:
                break
            if not six.PY2:
                stdout_line = stdout_line.decode()
            if '#(' not in stdout_line.strip():
                _send_line(ws, stdout_line)
        log.debug('Exit from the subprocess')

    manage_py = abspath(join(dirname(__file__), '..', 'manage.py'))
    _run_cmd(sys.executable + ' ' + manage_py + ' analyse_projects ' + project_names_line)
    run = Run.find_by_project_names_line(project_names_line)
    if not run:
        _send_line(ws, 'Run ' + str(run.id) + ' for projects ' + project_names_line + ' cannot be found. Has genotyping been failed?', error=True)

    ws.send(json.dumps({'finished': True}))
    return ''

def _send_line(ws, line, error=False):
    if error:
        log.err(line.rstrip())
    else:
        log.debug(line.rstrip())
    ws.send(json.dumps({
        'line': line.rstrip(),
        'error': error
    }))


def run_processing(project_names_line, redirect_to, email=None):
    pnames = project_names_line.split('--')

    log.debug(f'Recieved request to start analysis for {project_names_line}')
    run_id = hashlib.sha256(str(project_names_line).encode()).hexdigest()
    run_log = join(DATA_DIR, f'run_{run_id}.log')
    if isfile(run_log):
        msg = f'''<p>Run for projects {project_names_line} already started. Please, wait until it finished. 
                  <br>
                  <p>Follow the log at:</p>
                  <pre>{run_log}</pre>
                  <p>And reload the page when it\'s finished.</p>'''
    else:
        manage_py = abspath(join(dirname(__file__), '..', 'manage.py'))
        vardict = which('vardict')
        assert vardict, 'vardict is not in PATH. Are you running from "clearup" environment?'
        back_url = f"http://{clearup.HOST_IP}:{clearup.PORT}{redirect_to}"
        cmdl = f'{sys.executable} {manage_py} analyse_projects {project_names_line} --back_url={back_url}'
        if email:
            cmdl += f' --email={email}'
        log.debug(cmdl)
        process = subprocess.Popen(cmdl, stderr=subprocess.STDOUT, stdout=open(run_log, 'w'),
                         env=os.environ, close_fds=True, shell=True)

        msg = f'''<p>Starting analysis with a command line:</p>
                  <pre>{cmdl}</pre>
                  <p>Process is running under ID={process.pid}. Follow the log at:</p>
                  <pre>{run_log}</pre>
                  <p>And reload the page when it\'s finished.</p>'''

    if email:
        msg += f'<br>When the run finished, you will be notified by an email sent to {email}'

    return render_template(
        'submitted.html',
        projects=pnames,
        title='Comparing projects ' + ', '.join(pnames),
        project_names_line=project_names_line,
        redirect_to=redirect_to,
        message=msg
    )

def run_processing_printing_console(project_names_line, redirect_to):
    pnames = project_names_line.split('--')
    return render_template(
        'processing.html',
        projects=pnames,
        title='Comparing projects ' + ', '.join(pnames),
        project_names_line=project_names_line,
        redirect_to=redirect_to
    )


def render_phylo_tree_page(project_names_line, email=None):
    run = Run.find_by_project_names_line(project_names_line)

    if not Run.is_ready(run) or run.rerun_on_usercall:
        return run_processing(project_names_line,
            redirect_to=url_for('phylo_tree_page', project_names_line=project_names_line),
            email=email)

    # log.info('Runing ultrafast')
    # pairwise_dict = compare_pairwise(run)
    # plot_heatmap(pairwise_dict, run.work_dir_path(), ' '.join(p.name for p in run.projects))

    log.debug('Prank results found, rendering a tree for run ' + str(run.id))
    fasta_file = verify_file(run.fasta_file_path())
    if not fasta_file:
        raise RuntimeError('Run ' + project_names_line + ' does not contain ready fasta file. ' +
                           'Is genotyping ongoing in another window?')
    seq_by_id = read_fasta(fasta_file)

    log.debug('Preparing info for run ' + str(run.id))
    info_by_project = dict()
    prs = sorted(run.projects.all(), key=lambda p_: p_.name)
    for i, p in enumerate(prs):
        info_by_project[p.name] = dict()
        info_by_project[p.name]['name'] = p.name
        info_by_project[p.name]['color'] = PROJ_COLORS[i % len(PROJ_COLORS)]
        info_by_project[p.name]['samples'] = dict()
        for s in p.samples.all():
            log.debug('Searching SNPs for sample ' + s.name + ' in ' + p.name)
            info_by_project[p.name]['samples'][s.name] = {
                'name': s.name,
                'id': s.id,
                'sex': s.sex,
                'seq': [nt for nt in seq_by_id[s.name + FASTA_ID_PROJECT_SEPARATOR + p.name]],
            }
    log.debug('Prepared info_by_sample_by_project. Counting all samples now')
    all_samples_count = sum(len(info_by_sample['samples']) for info_by_sample in info_by_project.values())
    log.debug('Total samples: ' + str(all_samples_count))
    locations = [dict(
            chrom=l.chrom.replace('chr', ''),
            pos=l.pos,
            rsid=l.rsid,
            gene=l.gene)
        for i, l in enumerate(run.locations)]

    tree_file = verify_file(run.tree_file_path())
    log.debug('Found the tree file: ' + tree_file)
    if not tree_file:
        raise RuntimeError('Run ' + project_names_line +
                           ' does not contain the tree file (probably failed building phylogeny)')

    return render_template(
        'tree.html',
        projects=[{
            'name': str(project_info['name']),
            'color': project_info['color'],
            'samples': [info['name'] for info in project_info['samples'].values()],
            'ids': [info['id'] for info in project_info['samples'].values()]
        } for i, project_info in enumerate(info_by_project.values())],
        title=', '.join(sorted(info_by_project.keys())),
        tree_newick=open(tree_file).read(),
        info_by_sample_by_project=json.dumps(info_by_project),
        samples_count=all_samples_count,
        locations=json.dumps(locations)
    )


