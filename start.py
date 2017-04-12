#!/usr/bin/env python
from os.path import abspath, join, dirname, splitext, basename
from flask import Flask, render_template, send_from_directory, abort, redirect, url_for, send_file, request
from logging.handlers import RotatingFileHandler
import logging
from gevent.pywsgi import WSGIServer
from geventwebsocket.handler import WebSocketHandler

from ngs_utils import logger as log
from ngs_utils.file_utils import verify_file

from fingerprinting import app, DATA_DIR, HOST_IP, PORT
from fingerprinting.model import db, Sample, Project, Run, Location
from fingerprinting.sample_view import render_closest_comparison_page, send_file_for_igv
from fingerprinting.tree_view import run_analysis_socket_handler, render_phylo_tree_page


@app.route('/favicon.ico/')
def send_favicon():
    return send_from_directory('static', 'favicon.ico')


@app.route('/<project_names_line>/run_analysis/')
def run_analysis(project_names_line):
    return run_analysis_socket_handler(project_names_line)


@app.route('/<project_names_line>/tree/')
def phylo_tree_page(project_names_line):
    return render_phylo_tree_page(project_names_line)


@app.route('/<project_names_line>/tree/<int:sample_id>/')
def closest_comparison_page(project_names_line, sample_id):
    return render_closest_comparison_page(project_names_line, sample_id, request.args.get('snpIndex'))


@app.route('/<project_names_line>/tree/<int:sample_id>/add_usercall/', methods=['POST'])
def add_user_call(project_names_line, sample_id):
    log.info('Adding user call for ' + str(sample_id))
    edit_sample_id = request.form['editSampleId']
    sample = Sample.query.get(edit_sample_id)
    if not sample:
        log.error('Sample with ID=' + str(edit_sample_id) + ' not found')
        return redirect(url_for('closest_comparison_page', project_names_line=project_names_line, sample_id=sample_id))

    snp = sample.snps.join(Location).filter(Location.rsid==request.form['rsid']).first()
    snp.usercall = request.form['usercall']
    db.session.commit()
    return redirect(url_for('closest_comparison_page', project_names_line=project_names_line, sample_id=sample_id,
                            snpIndex=request.form['snpIndex']))


@app.route('/<run_id>/bamfiles/<bam_fname>/')
def bam_files_page(run_id, bam_fname):
    return send_file_for_igv(join(DATA_DIR, str(run_id), 'bams', bam_fname))


@app.route('/<project_names_line>/snps_bed/')
def locations_bed(project_names_line):
    run = Run.find_by_project_names_line(project_names_line)
    if not run:
        log.error('Run ' + project_names_line + ' not found')
        abort(404, {'message': 'Phylogenetic comparison for ' + project_names_line + ' is not found'})
    return send_file(run.snps_file)


@app.route('/<project_names_line>/')
def project_page(project_names_line):
    return redirect(url_for('phylo_tree_page', **locals()))


@app.route('/<project_names_line>/<sample_id>/')
def sample_page(project_names_line, sample_id):
    return redirect(url_for('closest_comparison_page', **locals()))


@app.route('/')
def homepage():
    projects = set()
    for run in Run.query.all():  # Finding projects with ready-to-view runs
        if verify_file(run.fasta_file_path(), silent=True) and run.projects.count() == 1:
            projects.add(run.projects[0])
    projects = sorted(projects, key=lambda p_: p_.name)
    t = render_template(
        'index.html',
        projects=[{
            'name': p.name,
            'data_dir': p.data_dir,
            'genome': p.genome,
            'samples': [{
                'id': s.id,
                'name': s.name,
            } for s in p.samples]
        } for p in projects],
    )
    return t


@app.errorhandler(404)
def page_not_found(error):
    return render_template('404.html', error=error.description), 404


if __name__ == "__main__":
    # app.run(host=config.HOST_IP, debug=config.IS_DEBUG)

    # if start_local_browser:
        # start server and web page pointing to it
        # url = "http://{HOST}:{PORT}".format(HOST=config.HOST_IP, PORT=PORT)
        # wb = webbrowser.get(None)  # instead of None, can be "firefox" etc
        # threading.Timer(1.25, lambda: wb.open(url)).start()
    
    # log_path = join(DATA_DIR, 'flask.log')
    # handler = RotatingFileHandler(log_path, maxBytes=10000, backupCount=10)
    # handler.setLevel(logging.INFO)
    # app.logger.addHandler(handler)
    
    http_server = WSGIServer((HOST_IP, PORT), app, handler_class=WebSocketHandler)
    log.info('Starting a webserver at ' + HOST_IP + ':' + str(PORT))
    http_server.serve_forever()
