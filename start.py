#!/usr/bin/env python2
from Bio import SeqIO
from os.path import abspath, join, dirname
from flask import Flask, render_template, send_from_directory, abort, redirect

from Utils import logger as log
from Utils.file_utils import safe_mkdir, file_transaction

from fingerprinting import config
from fingerprinting.draw_tree import draw_tree
from fingerprinting.model import Project
from fingerprinting.utils import read_fasta

app = Flask(__name__)


@app.route('/favicon.ico/')
def send_favicon():
    return send_from_directory('static', 'dist/favicon.ico')


@app.route('/<project_names>/')
def project_page(project_names):
    project_names = project_names.split(',')
    projects = [Project.query.filter_by(name=pn).first() for pn in project_names]
    if not projects:
        log.err('Projects ' + project_names + ' not found in database')
        abort(404)

    fasta_fpaths = [p.fingerprints_fasta_fpath for p in projects]
    project_work_dirpath = safe_mkdir(join(config.DATA_DIR, '__'.join(project_names)))  # toto fix dir name
    if len(project_names) > 1:
        merged_fasta_fpath = join(project_work_dirpath, 'all_fingerprints.fasta')
        with file_transaction(None, merged_fasta_fpath) as tx:
            with open(tx, 'w') as out:
                for proj_name, fasta in zip(project_names, fasta_fpaths):
                    with open(fasta) as f:
                        recs = SeqIO.parse(f, 'fasta')
                        for rec in recs:
                            rec.id += '_' + proj_name
                        SeqIO.write(recs, out, 'fasta')
    else:
        merged_fasta_fpath = fasta_fpaths[0]

    json_tree = draw_tree(merged_fasta_fpath, project_work_dirpath)

    all_samples_count = sum(len(p.samples.all()) for p in projects)
    t = render_template(
        'tree.html',
        title=', '.join(project_names),
        data=json_tree,
        tree_height=20 * all_samples_count,
        tree_width=5 * all_samples_count,
    )
    return t


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
