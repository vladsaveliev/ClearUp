#!/usr/bin/env python2
import json
import os
import re

from Bio import SeqIO
from collections import defaultdict
from os.path import abspath, join, dirname
from flask import Flask, render_template, send_from_directory, abort, redirect
from flask import Response, request

from Utils import logger as log
from Utils.file_utils import safe_mkdir, file_transaction, can_reuse

from fingerprinting import config
from fingerprinting.build_tree import build_tree
from fingerprinting.model import Project, db, Sample, PairedSample
from fingerprinting.utils import read_fasta, get_sample_and_project_name, FASTA_ID_PROJECT_SEPARATOR, \
    calculate_distance_matrix
from fingerprinting.lookups import get_snp_record, get_sample_by_name

app = Flask(__name__)


@app.route('/favicon.ico/')
def send_favicon():
    return send_from_directory('static', 'dist/favicon.ico')


@app.route('/<project_names>/')
def project_page(project_names):
    fingerprint_project = project_names
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
    tree_json = _tree_to_json_for_d3(tree, seq_by_id, color_by_proj, fingerprint_project=fingerprint_project)

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


def _tree_to_json_for_d3(tree, seq_by_id, color_by_proj, fingerprint_project):
    clade_dicts = []
    distance_matrix = calculate_distance_matrix(tree)
    _clade_to_json(tree.root, distance_matrix, cur_name='', clade_dicts=clade_dicts,
                   seq_by_id=seq_by_id, color_by_proj=color_by_proj, fingerprint_project=fingerprint_project)
    json_str = json.dumps(clade_dicts)
    return json_str


def _clade_to_json(clade, distance_matrix, cur_name, clade_dicts, seq_by_id, color_by_proj, fingerprint_project):
    if clade.name:  # leaf
        sample_name, project = get_sample_and_project_name(clade.name, fingerprint_project)
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
            _clade_to_json(child, distance_matrix, cur_name, clade_dicts, seq_by_id, color_by_proj, fingerprint_project)
    else:
        add_paired_sample(fingerprint_project, clade, distance_matrix)


class_by_nuc = {
    'C': 'c_nuc',
    'A': 'a_nuc',
    'T': 't_nuc',
    'G': 'g_nuc',
    'N': 'n_nuc',
}
def _seq_to_json(seq):
    return [(nuc, class_by_nuc[nuc.upper()]) for nuc in seq]


def add_paired_sample(fingerprint_project, clade, distance_matrix):
    sample_name, project_name = get_sample_and_project_name(clade.name, fingerprint_project)
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
            ps = PairedSample(paired_sample.name, fingerprint_project, paired_sample.id, sample)
            db.session.add(ps)
            db.session.commit()


@app.route('/<project_name>/<sample_id>')
def sample_page(project_name, sample_id):
    sample = Sample.query.filter_by(id=sample_id).first()
    if not sample:
        log.err('Sample not found in project ' + project_name)
        abort(404)
    sample_name = sample.name
    paired_sample_id = sample.paired_samples.filter_by(project=project_name).first().sample_id
    if not paired_sample_id:
        log.err('No matched sample for ' + sample_name)
        abort(404)
    paired_sample = Sample.query.filter_by(id=paired_sample_id).first()
    snps_dict = defaultdict(int)
    snps_dict['total_score'] = 0
    snps_dict['confidence'] = 'Low'
    snp_tables = []
    snp_records = []
    for snp_a, snp_b in zip(sample.fingerprints, paired_sample.fingerprints):
        snp_records.append(get_snp_record(snps_dict, snp_a, snp_b))
        if snp_a.index % 50 == 0:
            snp_tables.append(snp_records)
            snp_records = []
    snp_tables.append(snp_records)

    bam_fpath_a = '/%s/bamfiles/%s' % (sample.project.name, sample_name + '.bam')
    bam_fpath_b = '/%s/bamfiles/%s' % (paired_sample.project.name, paired_sample.name + '.bam')
    sample_a = {
        'id': sample.id,
        'name': sample.name,
        'project': sample.project.name,
        'bam': bam_fpath_a
    }
    sample_b = {
        'id': paired_sample.id,
        'name': paired_sample.name,
        'project': paired_sample.project.name,
        'bam': bam_fpath_b
    }
    t = render_template(
        'sample.html',
        project_name=project_name,
        genome=sample.project.genome,
        sampleA=sample_a,
        sampleB=sample_b,
        snps_data=snps_dict,
        snp_tables=snp_tables
    )
    return t


@app.route('/<project_name>/<prev_sample_id>/<edit_sample_id>/<snp_index>/add_usercall', methods=['POST'])
def addUserCall(project_name, prev_sample_id, edit_sample_id, snp_index):
    sample = Sample.query.filter_by(id=edit_sample_id).first()
    if not sample:
        log.err('Sample not found')
        return redirect('/' + project_name + '/' + prev_sample_id)

    fingerprint = sample.fingerprints.filter_by(index=snp_index).first()
    fingerprint.usercall = request.form['usercall']
    db.session.commit()
    return redirect('/' + project_name + '/' + prev_sample_id)


@app.route('/<project_name>/bamfiles/<bam_fname>')
def bam_files(project_name, bam_fname):
    project_work_dirpath = safe_mkdir(join(config.DATA_DIR, project_name))
    bam_fpath = join(project_work_dirpath, bam_fname)
    full_path = abspath(bam_fpath)

    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get('Range', None)
    if not range_header:
        return send_from_directory(project_work_dirpath, bam_fname)

    m = re.search('(\d+)-(\d*)', range_header)
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        log.err(error_msg)
        return error_msg

    size = os.path.getsize(full_path)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset

    data = None
    with open(full_path, 'rb') as f:
        f.seek(offset)
        data = f.read(length)

    rv = Response(data, 206, mimetype="application/octet-stream", direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(offset, offset + length - 1, size))

    log.info("GET range request: %s-%s %s" % (m.group(1), m.group(2), full_path))
    return rv


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
                # 'fingerprint': s.fingerprint,
            } for s in p.samples]
        } for p in projects],
    )
    return t


if __name__ == "__main__":
    # adds Flask command line options for setting host, port, etc.
    app.run(host=config.HOST_IP, debug=config.IS_DEBUG)
