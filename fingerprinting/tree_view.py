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

from Utils import logger as log
from Utils.file_utils import safe_mkdir, file_transaction, can_reuse
from Utils.file_utils import can_reuse, safe_mkdir

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
        self.work_dirpath = safe_mkdir(join(config.DATA_DIR, '_AND_'.join(project_names)))
        self.merged_fasta_fpath = merge_fasta(self.projects, self.work_dirpath)
        self.prank_out = os.path.join(self.work_dirpath,
             os.path.splitext(os.path.basename(self.merged_fasta_fpath))[0])
        self.tree_fpath = os.path.join(self.prank_out + '.best.dnd')


def run_prank_socket_handler(run_id):
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
        if not stdout_line.strip().startswith('#1#'):
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
        return render_template(
            'processing.html',
            projects=[{
                'name': p.name,
            } for i, p in enumerate(run.projects)],
            run_id=run_id,
            title='Processing ' + ', '.join(p.name for p in run.projects),
        )

    log.debug('Prank results found, rendering tree!')
    tree = next(Phylo.parse(run.tree_fpath, 'newick'))
    seq_by_id = read_fasta(run.merged_fasta_fpath)
    tree_json = tree_to_json_for_d3(tree, seq_by_id, run.color_by_proj, run_id=run_id)

    all_samples_count = sum(len(p.samples.all()) for p in run.projects)
    return render_template(
        'tree.html',
        projects=[{
            'name': p.name,
            'color': run.color_by_proj[p.name],
        } for i, p in enumerate(run.projects)],
        title=', '.join(p.name for p in run.projects),
        data=tree_json,
        tree_height=20 * all_samples_count,
        tree_width=5 * all_samples_count,
    )


def merge_fasta(projects, work_dirpath):
    fasta_fpaths = [p.fingerprints_fasta_fpath for p in projects]
    if len(projects) > 1:
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
    else:
        merged_fasta_fpath = fasta_fpaths[0]
    return merged_fasta_fpath



def tree_to_json_for_d3(tree, seq_by_id, color_by_proj, run_id):
    clade_dicts = []
    distance_matrix = calculate_distance_matrix(tree)
    _clade_to_json(tree.root, distance_matrix, cur_name='', clade_dicts=clade_dicts,
                   seq_by_id=seq_by_id, color_by_proj=color_by_proj, run_name=run_id)
    json_str = json.dumps(clade_dicts)
    return json_str


def _clade_to_json(clade, distance_matrix, cur_name, clade_dicts, seq_by_id, color_by_proj, run_name):
    if clade.name:  # leaf
        sample_name, project = get_sample_and_project_name(clade.name, run_name)
        sample = get_sample_by_name(sample_name, project)
        clade_d = {
            'id': cur_name + '.' + sample_name if cur_name else sample_name,
            'sample_id': sample.id if sample else None,
            'seq': _seq_to_json(seq_by_id[clade.name]),
            'project': project,
            'color': color_by_proj.get(project, 'black'),
        }
    else:  # internal node
        cur_name = cur_name + '.' + str(id(clade)) if cur_name else str(id(clade))
        clade_d = {
            'id': cur_name,
            'seq': None
        }
    clade_dicts.append(clade_d)
    if clade.clades:
        for child in clade:
            _clade_to_json(child, distance_matrix, cur_name, clade_dicts, seq_by_id, color_by_proj, run_name)
    else:
        _add_paired_sample(run_name, clade, distance_matrix)


def _seq_to_json(seq):
    return [(nuc, {
        'C': 'c_nuc',
        'A': 'a_nuc',
        'T': 't_nuc',
        'G': 'g_nuc',
        'N': 'n_nuc',
    }[nuc.upper()]) for nuc in seq]


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


