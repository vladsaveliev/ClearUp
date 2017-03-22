from __future__ import print_function
import json
from genericpath import isfile

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
from ngs_utils.file_utils import safe_mkdir, file_transaction, can_reuse
from ngs_utils.file_utils import can_reuse, safe_mkdir

from fingerprinting import config
from fingerprinting.model import Project, db, Sample, PairedSample
from fingerprinting.utils import read_fasta, get_sample_and_project_name, \
    FASTA_ID_PROJECT_SEPARATOR, calculate_distance_matrix
from fingerprinting.lookups import get_snp_record, get_sample_by_name

suffix = 'lnx' if 'linux' in platform else 'osx'
prank_bin = join(dirname(__file__), 'prank', 'prank_' + suffix, 'bin', 'prank')


PROJ_COLORS = ['#000000', '#90ed7d', '#f7a35c', '#8085e9', '#f15c80',
               '#e4d354', '#2b908f', '#f45b5b', '#91e8e1', '#7cb5ec']


class TreeRun:
    def __init__(self, run_id):
        self.run_id = run_id
        project_names = run_id.split(',')
        self.projects = [Project.query.filter_by(name=pn).first() for pn in project_names]
        if not self.projects:
            log.err('Projects ' + ', '.join(project_names) + ' not found in database')
            abort(404)
        self.color_by_proj = {p.name: PROJ_COLORS[i % len(PROJ_COLORS)] for i, p in enumerate(self.projects)}
        self.work_dirpath = safe_mkdir(join(config.DATA_DIR, '__AND__'.join(project_names)))
        self.merged_fasta_fpath = merge_fasta(self.projects, self.work_dirpath)
        self.prank_out = os.path.join(self.work_dirpath,
             os.path.splitext(os.path.basename(self.merged_fasta_fpath))[0])
        self.tree_fpath = os.path.join(self.prank_out + '.best.dnd')


def run_prank_socket_handler(run_id):
    log.debug('Recieved request to start prank for ' + run_id)

    ws = request.environ.get('wsgi.websocket', None)
    if not ws:
        raise RuntimeError('Environment lacks WSGI WebSocket support')

    run = TreeRun(run_id)

    cmdl = prank_bin + ' -d=' + run.merged_fasta_fpath + ' -o=' + run.prank_out + ' -showtree'
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
    return ''


def _send_line(ws, line):
    log.debug(line.rstrip())
    ws.send(json.dumps({
        'line': line.rstrip(),
    }))


def render_phylo_tree_page(run_id):
    run = TreeRun(run_id)

    if not can_reuse(run.tree_fpath, run.merged_fasta_fpath):
        log.debug('Tree for ' + run_id + ' was not found at ' + run.tree_fpath + ', renderring processing.html')
        return render_template(
            'processing.html',
            projects=[{'name': p.name} for p in run.projects],
            run_id=run_id,
            title='Processing ' + ', '.join(p.name for p in run.projects),
        )

    log.debug('Prank results found, rendering tree!')
    seq_by_id = read_fasta(run.merged_fasta_fpath)
    # tree = next(Phylo.parse(run.tree_fpath, 'newick'))
    # tree_json = tree_to_json_for_d3(tree, seq_by_id, run.color_by_proj, run_id=run_id)
    info_by_sample_by_project = dict()
    for p in run.projects:
        info_by_sample_by_project[p.name] = dict()
        info_by_sample_by_project[p.name]['name'] = p.name
        info_by_sample_by_project[p.name]['color'] = run.color_by_proj.get(p.name, 'black')
        info_by_sample_by_project[p.name]['samples'] = dict()
        for s in p.samples:
            info_by_sample_by_project[p.name]['samples'][s.name] = {
                'name': s.name,
                'id': s.id,
                'sex': s.sex,
                'seq': [nt for nt in seq_by_id[s.name + FASTA_ID_PROJECT_SEPARATOR + p.name]],
            }
    all_samples_count = sum(len(p.samples.all()) for p in run.projects)
    return render_template(
        'tree.html',
        projects=[{
            'name': str(p.name),
            'color': str(run.color_by_proj[p.name]),
            'samples': [str(sample.name) for sample in p.samples],
            'ids': [str(sample.id) for sample in p.samples]
        } for i, p in enumerate(run.projects)],
        title=', '.join(p.name for p in run.projects),
        tree_newick=open(run.tree_fpath).read(),
        info_by_sample_by_project=json.dumps(info_by_sample_by_project),
        samples_count=all_samples_count
    )


def merge_fasta(projects, work_dirpath):
    fasta_fpaths = [p.fingerprints_fasta_fpath for p in projects]
    merged_fasta_fpath = join(work_dirpath, 'fingerprints.fasta')
    if not can_reuse(merged_fasta_fpath, fasta_fpaths):
        all_records = []
        for proj, fasta in zip(projects, fasta_fpaths):
            with open(fasta) as f:
                recs = SeqIO.parse(f, 'fasta')
                for rec in recs:
                    rec.id += FASTA_ID_PROJECT_SEPARATOR + proj.name
                    rec.name = rec.description = ''
                    all_records.append(rec)
        with file_transaction(None, merged_fasta_fpath) as tx:
            with open(tx, 'w') as out:
                SeqIO.write(all_records, out, 'fasta')
    return merged_fasta_fpath


def tree_to_json_for_d3(tree, seq_by_id, color_by_proj, run_id):
    clade_dicts = []
    distance_matrix = calculate_distance_matrix(tree)
    _clade_to_json(tree.root, distance_matrix, cur_name='', clade_dicts=clade_dicts,
                   seq_by_id=seq_by_id, color_by_proj=color_by_proj, run_name=run_id)
    clade_dicts.sort(key=lambda c: c['id'])
    json_str = json.dumps(clade_dicts)
    return json_str


def _clade_to_json(clade, distance_matrix, cur_name, clade_dicts, seq_by_id, color_by_proj, run_name):
    if clade.name:  # leaf
        sample_name, project = get_sample_and_project_name(clade.name, run_name)
        sample = get_sample_by_name(sample_name, project)
        clade_d = {
            'id': cur_name + '.' + sample_name if cur_name else sample_name,
            'sample_id': sample.id if sample else None,
            'seq': [nt for nt in seq_by_id[clade.name]],
            'sex': sample.sex if sample else None,
            'project': project,
            'color': color_by_proj.get(project, 'black'),
        }
        _add_paired_sample(run_name, clade, distance_matrix)
    else:  # internal node
        cur_name = cur_name + '.' + _clade_id(clade) if cur_name else _clade_id(clade)
        clade_d = {
            'id': cur_name,
            'seq': None
        }
        for child in clade:
            _clade_to_json(child, distance_matrix, cur_name, clade_dicts, seq_by_id, color_by_proj, run_name)
    clade_dicts.append(clade_d)


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


