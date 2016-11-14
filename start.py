#!/usr/bin/env python2
import json

from Bio import SeqIO
from os.path import abspath, join, dirname
from flask import Flask, render_template, send_from_directory, abort, redirect

from Utils import logger as log
from Utils.file_utils import safe_mkdir, file_transaction, can_reuse

from fingerprinting import config
from fingerprinting.build_tree import build_tree
from fingerprinting.model import Project
from fingerprinting.utils import read_fasta

app = Flask(__name__)


@app.route('/favicon.ico/')
def send_favicon():
    return send_from_directory('static', 'dist/favicon.ico')


FASTA_ID_PROJECT_SEPARATOR = '____PROJECT_'


@app.route('/<project_names>/')
def project_page(project_names):
    project_names = project_names.split(',')
    projects = [Project.query.filter_by(name=pn).first() for pn in project_names]
    if not projects:
        log.err('Projects ' + project_names + ' not found in database')
        abort(404)

    fasta_fpaths = [p.fingerprints_fasta_fpath for p in projects]
    project_work_dirpath = safe_mkdir(join(config.DATA_DIR, '_AND_'.join(project_names)))
    safe_mkdir(project_work_dirpath)

    if len(projects) > 1:
        merged_fasta_fpath = join(project_work_dirpath, 'fingerprints.fasta')
        if not can_reuse(merged_fasta_fpath, fasta_fpaths):
            all_records = []
            for proj_name, fasta in zip(project_names, fasta_fpaths):
                with open(fasta) as f:
                    recs = SeqIO.parse(f, 'fasta')
                    for rec in recs:
                        rec.id += FASTA_ID_PROJECT_SEPARATOR + proj_name
                        rec.name = rec.description = ''
                        all_records.append(rec)
            with file_transaction(None, merged_fasta_fpath) as tx:
                with open(tx, 'w') as out:
                    SeqIO.write(all_records, out, 'fasta')
    else:
        merged_fasta_fpath = fasta_fpaths[0]

    tree = build_tree(merged_fasta_fpath, project_work_dirpath)
    seq_by_id = read_fasta(merged_fasta_fpath)

    project_colors = ['#000000', '#90ed7d', '#f7a35c', '#8085e9', '#f15c80',
                      '#e4d354', '#2b908f', '#f45b5b', '#91e8e1', '#7cb5ec']
    color_by_proj = {p.name: project_colors[i % len(project_colors)]
                     for i, p in enumerate(projects)}
    tree_json = _tree_to_json_for_d3(tree, seq_by_id, color_by_proj)

    all_samples_count = sum(len(p.samples.all()) for p in projects)
    t = render_template(
        'tree.html',
        projects=[{
            'name': p.name,
            'color': project_colors[i % len(project_colors)]
        } for i, p in enumerate(projects)],
        title=', '.join(project_names),
        data=tree_json,
        tree_height=20 * all_samples_count,
        tree_width=5 * all_samples_count,
    )
    return t


def _tree_to_json_for_d3(tree, seq_by_id, color_by_proj):
    clade_dicts = []
    _clade_to_json(tree.root, cur_name='', clade_dicts=clade_dicts,
                   seq_by_id=seq_by_id, color_by_proj=color_by_proj)
    json_str = json.dumps(clade_dicts)
    return json_str


def _clade_to_json(clade, cur_name, clade_dicts, seq_by_id, color_by_proj):
    if clade.name:  # leaf
        if FASTA_ID_PROJECT_SEPARATOR in clade.name:
            sample_name, project = clade.name.split(FASTA_ID_PROJECT_SEPARATOR)
        else:
            sample_name, project = clade.name, ''
        clade_d = {
            'id': cur_name + '.' + sample_name if cur_name else sample_name,
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
            _clade_to_json(child, cur_name, clade_dicts, seq_by_id, color_by_proj)


class_by_nuc = {
    'C': 'c_nuc',
    'A': 'a_nuc',
    'T': 't_nuc',
    'G': 'g_nuc',
    'N': 'n_nuc',
}
def _seq_to_json(seq):
    return [(nuc, class_by_nuc[nuc.upper()]) for nuc in seq]


@app.route('/')
def homepage():
    projects = Project.query.all()
    t = render_template(
        'index.html',
        projects=[{
            'name': p.name,
            'bcbio_final_path': p.bcbio_final_path,
            'genome': p.genome,
            'samples': [{
                'id': s.id,
                'name': s.name,
                'fingerprint': s.fingerprint,
            } for s in p.samples]
        } for p in projects],
    )
    return t


if __name__ == "__main__":
    # adds Flask command line options for setting host, port, etc.
    app.run(host=config.HOST_IP, debug=config.IS_DEBUG)
