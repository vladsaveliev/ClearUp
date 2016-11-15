#!/usr/bin/env python2
import json

import os
import subprocess
import time
import platform

from os.path import abspath, join, dirname, splitext, basename
from flask import Flask, render_template, send_from_directory, abort, redirect
from flask import Response, request

from geventwebsocket.handler import WebSocketHandler
from gevent.pywsgi import WSGIServer

from Utils.file_utils import safe_mkdir, file_transaction, can_reuse
from Utils import logger as log
log.is_debug = True

from fingerprinting import config
from fingerprinting.model import Project, db, Sample
from fingerprinting.sample_view import send_bam_files, render_sample_page
from fingerprinting.tree_view import run_prank_socket_handler, render_phylo_tree_page

app = Flask(__name__)


if log.is_local() and platform.system() == 'Darwin':
    PORT = 5002
else:
    PORT = 80


@app.route('/favicon.ico/')
def send_favicon():
    return send_from_directory('static', 'dist/favicon.ico')


@app.route("/<run_id>/run_prank/")
def run_prank(run_id):
    return run_prank_socket_handler(run_id)


@app.route('/<run_id>/')
def phylo_tree_page(run_id):
    return render_phylo_tree_page(run_id)


@app.route('/<run_id>/<sample_id>')
def sample_page(run_id, sample_id):
    return render_sample_page(run_id, sample_id)


@app.route('/<run_id>/<prev_sample_id>/<edit_sample_id>/<snp_index>/add_usercall', methods=['POST'])
def add_user_call(run_id, prev_sample_id, edit_sample_id, snp_index):
    sample = Sample.query.filter_by(id=edit_sample_id).first()
    if not sample:
        log.err('Sample not found')
        return redirect('/' + run_id + '/' + prev_sample_id)

    fingerprint = sample.fingerprints.filter_by(index=snp_index).first()
    fingerprint.usercall = request.form['usercall']
    db.session.commit()
    return redirect('/' + run_id + '/' + prev_sample_id)


@app.route('/<run_id>/bamfiles/<bam_fname>')
def bam_files_page(run_id, bam_fname):
    return send_bam_files(run_id, bam_fname)


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
    # app.run(host=config.HOST_IP, debug=config.IS_DEBUG)

    # if start_local_browser:
        # start server and web page pointing to it
        # url = "http://{HOST}:{PORT}".format(HOST=config.HOST_IP, PORT=PORT)
        # wb = webbrowser.get(None)  # instead of None, can be "firefox" etc
        # threading.Timer(1.25, lambda: wb.open(url)).start()

    http_server = WSGIServer((config.HOST_IP, PORT), app, handler_class=WebSocketHandler)
    log.info('Starting webserve at ' + config.HOST_IP + ':' + str(PORT))
    http_server.serve_forever()
