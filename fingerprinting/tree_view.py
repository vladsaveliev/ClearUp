from __future__ import print_function
import json

import os
import re
import subprocess
import time
from sys import platform

import sys
from os.path import join, dirname
from collections import defaultdict
from os.path import abspath, join, dirname, splitext, basename

from Bio import SeqIO, Phylo
from flask import Flask, render_template, abort, request

from ngs_utils import logger as log
from ngs_utils.bed_utils import Region
from ngs_utils.file_utils import safe_mkdir, file_transaction, can_reuse, verify_file
from ngs_utils.file_utils import can_reuse, safe_mkdir

from fingerprinting.model import Project, db, Sample, Run, get_or_create_run
from fingerprinting.utils import read_fasta, FASTA_ID_PROJECT_SEPARATOR

suffix = 'lnx' if 'linux' in platform else 'osx'
prank_bin = join(dirname(__file__), 'prank', 'prank_' + suffix, 'bin', 'prank')


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


def run_analysis_socket_handler(run_id):
    log.debug('Recieved request to start analysis for ' + run_id)
    ws = request.environ.get('wsgi.websocket', None)
    if not ws:
        raise RuntimeError('Environment lacks WSGI WebSocket support')

    def _run_cmd(cmdl):
        log.debug(cmdl)
        proc = subprocess.Popen(cmdl.split(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, env=os.environ)
        # lines = []
        # prev_time = time.time()
        for stdout_line in iter(proc.stdout.readline, ''):
            # lines.append(stdout_line)
            # cur_time = time.time()
            # if cur_time - prev_time > 2:
            if '#(' not in stdout_line.strip():
                _send_line(ws, stdout_line)
            # lines = []
        log.debug('Exit from the subprocess')

    manage_py = abspath(join(dirname(__file__), '..', 'manage.py'))
    _run_cmd(sys.executable + ' ' + manage_py + ' analyse_projects ' + run_id)
    run = Run.query.get(run_id)
    if not run:
        _send_line(ws, 'Run ' + run_id + ' cannot be found. Has genotyping been failed?', error=True)

    fasta_file = verify_file(run.fasta_file_path())
    if not fasta_file:
        _send_line(ws, 'Run ' + run_id + ' does not contain ready fasta file. Is genotyping ongoing in another window?', error=True)

    prank_out = join(run.work_dir, splitext(basename(fasta_file))[0])
    _send_line(ws, '')
    _send_line(ws, 'Building phylogeny tree using prank...')
    _run_cmd(prank_bin + ' -d=' + fasta_file + ' -o=' + prank_out + ' -showtree')
    if not verify_file(prank_out + '.best.dnd'):
        _send_line(ws, 'Prank failed to run', error=True)
    
    os.rename(prank_out + '.best.dnd', run.tree_file_path())
    os.remove(prank_out + '.best.fas')
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


def render_phylo_tree_page(run_id):
    run = Run.query.filter_by(id=run_id).first()
    if not run or not can_reuse(verify_file(run.tree_file_path(), silent=True),
                                verify_file(run.fasta_file_path(), silent=True)):
        return render_template(
            'processing.html',
            projects=run_id.split(','),
            run_id=run_id,
            title='Processing ' + ', '.join(run_id.split(',')),
        )

    log.debug('Prank results found, rendering tree!')
    fasta_file = verify_file(run.fasta_file_path())
    if not fasta_file:
        raise RuntimeError('Run ' + run_id + ' does not contain ready fasta file. Is genotyping ongoing in another window?')
    seq_by_id = read_fasta(fasta_file)

    info_by_sample_by_project = dict()
    for i, p in enumerate(run.projects):
        info_by_sample_by_project[p.name] = dict()
        info_by_sample_by_project[p.name]['name'] = p.name
        info_by_sample_by_project[p.name]['color'] = PROJ_COLORS[i % len(PROJ_COLORS)]
        info_by_sample_by_project[p.name]['samples'] = dict()
        for s in p.samples:
            info_by_sample_by_project[p.name]['samples'][s.name] = {
                'name': s.name,
                'id': s.id,
                'sex': s.sex,
                'seq': [nt for nt in seq_by_id[s.name + FASTA_ID_PROJECT_SEPARATOR + p.name]],
                'snps': [snp.genotype for snp in s.snps_from_run(run)],
            }
    all_samples_count = sum(len(p.samples.all()) for p in run.projects)
    locations = [dict(
            chrom=l.chrom.replace('chr', ''),
            pos=l.pos,
            rsid=l.rsid,
            gene=l.gene,
            index=l.index)
        for l in run.locations]
    
    tree_file = verify_file(run.tree_file_path())
    if not tree_file:
        raise RuntimeError('Run ' + run_id + ' does not contain the tree file (probably failed building phylogeny)')
        
    return render_template(
        'tree.html',
        projects=[{
            'name': str(p.name),
            'color': PROJ_COLORS[i % len(PROJ_COLORS)],
            'samples': [str(sample.name) for sample in p.samples],
            'ids': [str(sample.id) for sample in p.samples]
        } for i, p in enumerate(run.projects)],
        title=', '.join(p.name for p in run.projects),
        tree_newick=open(tree_file).read(),
        info_by_sample_by_project=json.dumps(info_by_sample_by_project),
        samples_count=all_samples_count,
        locations=json.dumps(locations)
    )


