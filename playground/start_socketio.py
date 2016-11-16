#!/usr/bin/env python2
import json

import os
import subprocess
import time
import platform
import threading
import webbrowser
import random

from Bio import Phylo
from os.path import abspath, join, dirname, splitext, basename
from flask import Flask, render_template, send_from_directory, abort, redirect, jsonify
from flask import Response, request
from flask_socketio import SocketIO, send, emit

from geventwebsocket.handler import WebSocketHandler
from gevent.pywsgi import WSGIServer
import gevent

from Utils import logger as log

from fingerprinting.utils import read_fasta
log.is_debug = True
from Utils.file_utils import safe_mkdir, file_transaction, can_reuse

from fingerprinting import config
from fingerprinting.model import Project, db, Sample
from fingerprinting.sample_view import send_bam_files, render_sample_page
from fingerprinting.tree_view import tree_to_json_for_d3, \
    prank_bin, merge_fasta, PROJ_COLORS

app = Flask(__name__)
socketio = SocketIO(app)

start_local_browser = platform.system() == 'Darwin' and log.is_local()
if start_local_browser:
    PORT = 5000
else:
    PORT = 80


@app.route('/favicon.ico/')
def send_favicon():
    return send_from_directory('static', 'dist/favicon.ico')


#
# LONG_POLLING_VERSION
#
# @app.route('/processing')
# def processing():
#     def inner():
#         fasta_fpath = '/Users/vlad/vagrant/Fingerprinting/data/Dev_0261_MiSeq_MCRC_PRCC_AND_Dev_0261_MiSeq_COPY/fingerprints.fasta'
#         output_dirpath = '/Users/vlad/vagrant/Fingerprinting/data/Dev_0261_MiSeq_MCRC_PRCC_AND_Dev_0261_MiSeq_COPY'
#         output = join(output_dirpath, splitext(basename(fasta_fpath))[0])
#         cmdl = prank_bin + ' -d=' + fasta_fpath + ' -o=' + output + ' -showtree'
#         proc = subprocess.Popen(
#             cmdl,             #call something with a lot of output so we can see it
#             shell=True,
#             stderr=subprocess.STDOUT,
#             stdout=subprocess.PIPE,
#         )
#         for line in iter(proc.stdout.readline, ''):
#             log.info(line.rstrip())
#             line = '<span style="">' + line.rstrip() + '<br/>'
#             yield line + '\n'
#
#     return Response(inner(), mimetype='text/html')  # text/html is required for most browsers to show th$


# @app.route("/check_processing_log_updated/<run_id>")
# def check_processing_log_updated(run_id):
#     log.debug('Initializing websocket for ' + run_id)
#     ws = request.environ.get('wsgi.websocket', None)
#     if not ws:
#         raise RuntimeError('Environment lacks WSGI WebSocket support')
#     gevent.sleep(1)
#     log.debug('Initialized websocket for ' + run_id + ', sending ready')
#     ws.send('ready')


@app.route("/send_processing_log_updates/<run_id>")
def send_processing_log_updates(run_id):
    return send_processing_updates_ajax(run_id)


@socketio.on('connected')
def run_prank(run_id):
    project_names = run_id.split(',')
    projects = [Project.query.filter_by(name=pn).first() for pn in project_names]
    if not projects:
        log.err('Projects ' + ', '.join(project_names) + ' not found in database')
        abort(404)
    work_dirpath = safe_mkdir(join(config.DATA_DIR, '_AND_'.join(project_names)))
    safe_mkdir(work_dirpath)
    merged_fasta_fpath = merge_fasta(projects, work_dirpath)

    prank_out = os.path.join(work_dirpath, os.path.splitext(os.path.basename(merged_fasta_fpath))[0])
    cmdl = prank_bin + ' -d=' + merged_fasta_fpath + ' -o=' + prank_out + ' -showtree'
    log.debug('Starting prank ' + cmdl)
    proc = subprocess.Popen(cmdl.split(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    # lines = []
    # prev_time = time.time()
    for stdout_line in iter(proc.stdout.readline, ''):
        print stdout_line.rstrip()
        # lines.append(stdout_line)
        cur_time = time.time()
        # if cur_time - prev_time > 2:
        emit('running',
            json.dumps({
                'finished': False,
                'lines': [stdout_line.rstrip()],
            })
        )
        # lines = []
    emit('running',
        json.dumps({
            'finished': True,
            'lines': [],
        })
    )


@app.route('/<run_id>/')
def phylo_tree_page(run_id):
    project_names = run_id.split(',')
    projects = [Project.query.filter_by(name=pn).first() for pn in project_names]
    if not projects:
        log.err('Projects ' + ', '.join(project_names) + ' not found in database')
        abort(404)
    color_by_proj = {p.name: PROJ_COLORS[i % len(PROJ_COLORS)] for i, p in enumerate(projects)}
    work_dirpath = safe_mkdir(join(config.DATA_DIR, '_AND_'.join(project_names)))
    safe_mkdir(work_dirpath)
    merged_fasta_fpath = merge_fasta(projects, work_dirpath)

    prank_out = os.path.join(work_dirpath, os.path.splitext(os.path.basename(merged_fasta_fpath))[0])
    tree_fpath = os.path.join(prank_out + '.best.dnd')
    if not can_reuse(tree_fpath, merged_fasta_fpath):
        return render_template(
            'processing.html',
            projects=[{
                'name': p.name,
            } for i, p in enumerate(projects)],
            run_id=run_id,
            title='Processing ' + ', '.join(project_names),
        )

    log.debug('Prank results found, rendering tree!')
    tree = next(Phylo.parse(tree_fpath, 'newick'))
    seq_by_id = read_fasta(merged_fasta_fpath)
    tree_json = tree_to_json_for_d3(tree, seq_by_id, color_by_proj, run_id=run_id)

    all_samples_count = sum(len(p.samples.all()) for p in projects)
    return render_template(
        'tree.html',
        projects=[{
            'name': p.name,
            'color': color_by_proj[p.name],
        } for i, p in enumerate(projects)],
        title=', '.join(project_names),
        data=tree_json,
        tree_height=20 * all_samples_count,
        tree_width=5 * all_samples_count,
    )


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
    socketio.run(app, port=PORT, debug=True)

    # app.run(host=config.HOST_IP, debug=config.IS_DEBUG)

    # if start_local_browser:
        # start server and web page pointing to it
        # url = "http://{HOST}:{PORT}".format(HOST=config.HOST_IP, PORT=PORT)
        # wb = webbrowser.get(None)  # instead of None, can be "firefox" etc
        # threading.Timer(1.25, lambda: wb.open(url)).start()

    # http_server = WSGIServer(('', PORT), app, handler_class=WebSocketHandler)
    # http_server.serve_forever()
