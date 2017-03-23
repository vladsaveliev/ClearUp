from __future__ import print_function
import json

import os
import re
import subprocess
import time
from sys import platform
from os.path import join, dirname
from collections import defaultdict
from os.path import abspath, join, dirname, splitext, basename

from Bio import SeqIO, Phylo
from flask import Flask, render_template, abort, request

from ngs_utils import logger as log
from ngs_utils.bed_utils import Region
from ngs_utils.file_utils import safe_mkdir, file_transaction, can_reuse, verify_file
from ngs_utils.file_utils import can_reuse, safe_mkdir

from fingerprinting import config
from fingerprinting.model import Project, db, Sample, PairedSample, Run, get_run
from fingerprinting.utils import read_fasta, get_sample_and_project_name, \
    FASTA_ID_PROJECT_SEPARATOR, calculate_distance_matrix
from fingerprinting.lookups import get_snp_record, get_sample_by_name

suffix = 'lnx' if 'linux' in platform else 'osx'
prank_bin = join(dirname(__file__), 'prank', 'prank_' + suffix, 'bin', 'prank')


PROJ_COLORS = ['#000000', '#90ed7d', '#f7a35c', '#8085e9', '#f15c80',
               '#e4d354', '#2b908f', '#f45b5b', '#91e8e1', '#7cb5ec']



def _get_run(run_id):
    run = get_run(run_id)
    if not run:
        abort(404)
    return run


def run_prank_socket_handler(run_id):
    log.debug('Recieved request to start prank for ' + run_id)

    ws = request.environ.get('wsgi.websocket', None)
    if not ws:
        raise RuntimeError('Environment lacks WSGI WebSocket support')

    run = _get_run(run_id)
    
    prank_out = join(run.work_dir, splitext(basename(run.fasta_file))[0])
    cmdl = prank_bin + ' -d=' + run.fasta_file + ' -o=' + prank_out + ' -showtree'
    _send_line(ws, '')
    _send_line(ws, 'Building phylogeny tree using prank...')
    log.debug(cmdl)
    proc = subprocess.Popen(cmdl.split(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    # lines = []
    # prev_time = time.time()
    for stdout_line in iter(proc.stdout.readline, ''):
        # lines.append(stdout_line)
        # cur_time = time.time()
        # if cur_time - prev_time > 2:
        if '#(' not in stdout_line.strip():
            _send_line(ws, stdout_line)
        # lines = []
    ws.send(json.dumps({
        'finished': True,
    }))
    if not verify_file(prank_out + '.best.dnd'):
        raise RuntimeError('Prank failed to run')
    os.rename(prank_out + '.best.dnd', run.tree_file)
    os.remove(prank_out + '.best.fas')
    return ''


def _send_line(ws, line):
    log.debug(line.rstrip())
    ws.send(json.dumps({
        'line': line.rstrip(),
    }))


def render_phylo_tree_page(run_id):
    run = _get_run(run_id)

    if not can_reuse(run.tree_file, run.fasta_file):
        log.debug('Tree for ' + run_id + ' was not found at ' + run.tree_file + ', renderring processing.html')
        return render_template(
            'processing.html',
            projects=[{'name': p.name} for p in run.projects],
            run_id=run_id,
            title='Processing ' + ', '.join(p.name for p in run.projects),
        )

    log.debug('Prank results found, rendering tree!')
    seq_by_id = read_fasta(run.fasta_file)
    # tree = next(Phylo.parse(run.tree_ffile, 'newick'))
    # tree_json = tree_to_json_for_d3(tree, seq_by_id, run.color_by_proj, run_id=run_id)
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
            }
    all_samples_count = sum(len(p.samples.all()) for p in run.projects)
    locations = [dict(chrom=l.chrom.replace('chr', ''), pos=l.pos, rsid=l.rsid, gene=l.gene) for l in run.locations]
    return render_template(
        'tree.html',
        projects=[{
            'name': str(p.name),
            'color': PROJ_COLORS[i % len(PROJ_COLORS)],
            'samples': [str(sample.name) for sample in p.samples],
            'ids': [str(sample.id) for sample in p.samples]
        } for i, p in enumerate(run.projects)],
        title=', '.join(p.name for p in run.projects),
        tree_newick=open(run.tree_file).read(),
        info_by_sample_by_project=json.dumps(info_by_sample_by_project),
        samples_count=all_samples_count,
        locations=json.dumps(locations)
    )


# def tree_to_json_for_d3(tree, seq_by_id, color_by_proj, run_id):
#     clade_dicts = []
#     distance_matrix = calculate_distance_matrix(tree)
#     _clade_to_json(tree.root, distance_matrix, cur_name='', clade_dicts=clade_dicts,
#                    seq_by_id=seq_by_id, color_by_proj=color_by_proj, run_name=run_id)
#     clade_dicts.sort(key=lambda c: c['id'])
#     json_str = json.dumps(clade_dicts)
#     return json_str
#
#
# def _clade_to_json(clade, distance_matrix, cur_name, clade_dicts, seq_by_id, color_by_proj, run_name):
#     if clade.name:  # leaf
#         sample_name, project = get_sample_and_project_name(clade.name, run_name)
#         sample = get_sample_by_name(sample_name, project)
#         clade_d = {
#             'id': cur_name + '.' + sample_name if cur_name else sample_name,
#             'sample_id': sample.id if sample else None,
#             'seq': [nt for nt in seq_by_id[clade.name]],
#             'sex': sample.sex if sample else None,
#             'project': project,
#             'color': color_by_proj.get(project, 'black'),
#         }
#         _add_paired_sample(run_name, clade, distance_matrix)
#     else:  # internal node
#         cur_name = cur_name + '.' + _clade_id(clade) if cur_name else _clade_id(clade)
#         clade_d = {
#             'id': cur_name,
#             'seq': None
#         }
#         for child in clade:
#             _clade_to_json(child, distance_matrix, cur_name, clade_dicts, seq_by_id, color_by_proj, run_name)
#     clade_dicts.append(clade_d)


def _clade_id(clade):
    if clade.name:
        return str(clade.name)
    else:
        return str(hash('__'.join(_clade_id(child) for child in clade)))


def _add_paired_sample(run_name, clade, distance_matrix):
    sample_name, project_name = get_sample_and_project_name(clade.name, run_name)
    project = Project.query.filter_by(name=project_name).first()
    if project:
        sample = project.samples.filter_by(name=sample_name).first()
        paired_clade = distance_matrix[clade][1]
        if paired_clade:
            paired_sample_name, paired_project_name = get_sample_and_project_name(paired_clade.name)
            if not paired_project_name:
                paired_project_name = project_name
            paired_sample_project = Project.query.filter_by(name=paired_project_name).first()
            paired_sample = paired_sample_project.samples.filter_by(name=paired_sample_name).first()
            ps = PairedSample(paired_sample.name, run_name, paired_sample.id, sample)
            db.session.add(ps)
            db.session.commit()


