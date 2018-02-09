from __future__ import print_function
import json

import os
import re
import subprocess
import time

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
from ngs_utils.file_utils import safe_mkdir, file_transaction, can_reuse, verify_file
from ngs_utils.file_utils import can_reuse, safe_mkdir
from ngs_utils import logger as log

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
    ws = request.environ.get('wsgi.websocket', None)
    if not ws:
        raise RuntimeError('Environment lacks WSGI WebSocket support')

    run_log = join(DATA_DIR, str(project_names_line) + '.run_analysis.log')
    if isfile(run_log):
        _send_line(ws, f'Run for projects {project_names_line} already started. '
                       f'Please, wait until it funished. To restart, please remove {run_log} '
                       f'and reload the page.')
    else:
        log.debug(f'Recieved request to start analysis for {project_names_line}')

    manage_py = abspath(join(dirname(__file__), '..', 'manage.py'))
    cmdl = f'{sys.executable} {manage_py} analyse_projects {project_names_line}'
    log.debug(cmdl)
    _send_line(ws, f'\nStarting analysis:\n')
    _send_line(ws, f'{cmdl}\n')
    _send_line(ws, f'\nFollow the log at:\n')
    _send_line(ws, f'{run_log}\n')
    _send_line(ws, f'\nAnd reload the page when it\'s finished.')
    subprocess.Popen(cmdl.split(), stderr=subprocess.STDOUT, stdout=open(run_log, 'w'), env=os.environ,
                     close_fds=True)

    # run = Run.find_by_project_names_line(project_names_line)
    # if not run:
    #     _send_line(ws, 'Run ' + str(run.id) + ' for projects ' + project_names_line + ' cannot be found. Has genotyping been failed?', error=True)
    #
    # ws.send(json.dumps({'finished': True}))
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


def run_processing(project_names_line, redirect_to):
    pnames = project_names_line.split('--')
    return render_template(
        'processing.html',
        projects=pnames,
        title='Comparing projects ' + ', '.join(pnames),
        project_names_line=project_names_line,
        redirect_to=redirect_to
    )


def render_phylo_tree_page(project_names_line):
    run = Run.find_by_project_names_line(project_names_line)

    if not Run.is_ready(run) or run.rerun_on_usercall:
        return run_processing(project_names_line,
            redirect_to=url_for(
                'phylo_tree_page',
                project_names_line=project_names_line))

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


