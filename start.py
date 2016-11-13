#!/usr/bin/env python2
from os.path import abspath, join, dirname
from flask import Flask, render_template, send_from_directory, abort, redirect

from Utils import logger as log
from Utils.file_utils import safe_mkdir

from fingerprinting import config
from fingerprinting.draw_tree import draw_tree
from fingerprinting.model import Project


app = Flask(__name__)


@app.route('/favicon.ico/')
def send_favicon():
    return send_from_directory('static', 'dist/favicon.ico')


@app.route('/<project_name>/')
def project_page(project_name):
    project = Project.query.filter_by(name=project_name).first()
    if not project:
        log.err('Project ' + project_name + ' not found in database')
        abort(404)

    fasta_fpath = project.fingerprints_fasta_fpath
    project_work_dirpath = safe_mkdir(join(config.DATA_DIR, project_name))

    json_tree = draw_tree(fasta_fpath, project_work_dirpath)

    t = render_template(
        'tree.html',
        project_name=project_name,
        sample_names=[s.name for s in project.samples],
        data=json_tree,
        tree_height=20 * len(project.samples.all()),
        tree_width=5 * len(project.samples.all()),
    )
    return t


@app.route('/')
def homepage():
    return redirect('/Dev_0261_MiSeq_MCRC_PRCC/')


if __name__ == "__main__":
    # adds Flask command line options for setting host, port, etc.
    app.run(host=config.HOST_IP, debug=config.IS_DEBUG)
