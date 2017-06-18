#!/usr/bin/env python
import traceback
from collections import defaultdict
import os
from os.path import abspath, join, dirname, splitext, basename
import click
from flask import Flask, render_template, send_from_directory, abort, redirect, url_for, send_file, request
from logging.handlers import RotatingFileHandler
import logging
from gevent.pywsgi import WSGIServer
from geventwebsocket.handler import WebSocketHandler

from ngs_utils import logger as log
from ngs_utils.file_utils import verify_file, safe_mkdir
from ngs_utils.sambamba import index_bam
from ngs_utils.utils import is_local

from clearup import app, DATA_DIR, HOST_IP, PORT, get_version
from clearup.genotype import build_tree
from clearup.model import db, Sample, Project, Run, Location, SNP
from clearup.sample_view import render_closest_comparison_page, send_file_for_igv
from clearup.tree_view import run_analysis_socket_handler, render_phylo_tree_page


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--host',
              type=click.STRING,
              metavar='<IP address>',
              default=HOST_IP)
@click.option('--port',
              type=click.INT,
              metavar='<port>',
              default=PORT)
@click.version_option(version=get_version())
def main(host, port):
    os.environ['FLASK_DEBUG'] = '1'
    safe_mkdir(DATA_DIR)
    # log_path = join(DATA_DIR, 'flask.log')
    # handler = RotatingFileHandler(log_path, maxBytes=10000, backupCount=10)
    # handler.setLevel(logging.INFO)
    # app.logger.addHandler(handler)

    http_server = WSGIServer((host, port), app, handler_class=WebSocketHandler)
    log.info('Starting a webserver at ' + host + ':' + str(port))
    http_server.serve_forever()


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

    # snp = sample.snps.join(Location).filter(Location.rsid==request.form['rsid']).first()
    snp = sample.snps.filter(SNP.rsid==request.form['rsid']).first()
    snp.usercall = request.form['usercall']
    db.session.commit()

    # Forcing rebuilding the trees of affected runs
    for run in Run.query.all():
        if sample.project in run.projects:
            if any(l for l in run.locations if l.rsid == snp.rsid):
                # log.debug('Removing tree file ' + run.tree_file_path())
                # os.rename(run.fasta_file_path(), run.fasta_file_path() + '.bak')
                # os.rename(run.tree_file_path(), run.tree_file_path() + '.bak')
                run.rerun_on_usercall = True
                db.session.commit()

    return redirect(url_for('closest_comparison_page',
                            project_names_line=project_names_line,
                            sample_id=sample_id))


@app.route('/<run_id>/bamfiles/<fname>/')
def bam_files_page(run_id, fname):
    fpath = join(DATA_DIR, str(run_id), 'bams', fname)
    bam_fpath = fpath
    if fpath.endswith('.bai'):
        bam_fpath = fpath[:-4]
    index_bam(bam_fpath)
    return send_file_for_igv(fpath)


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
    proj_by_genome_build = defaultdict(list)
    for p in Project.query.all():
        proj_by_genome_build[p.genome].append(p)
    proj_by_genome_build = {
        g: sorted(ps, key=lambda p_: p_.name) for g, ps in proj_by_genome_build.items()
    }
    # for run in Run.query.all():  # Finding projects with ready-to-view runs
    #     if verify_file(run.fasta_file_path(), silent=True) and run.projects.count() == 1:
    #         projects.add(run.projects[0])
    t = render_template(
        'index.html',
        proj_by_genome_build=[(g, [{
                'name': p.name,
                'data_dir': p.data_dir,
                'genome': p.genome,
                'samples': [{
                    'id': s.id,
                    'name': s.name,
                } for s in p.samples]
            } for p in ps])
            for g, ps in sorted(proj_by_genome_build.items())
        ]
    )
    return t


@app.errorhandler(404)
def page_not_found(error):
    return render_template(
        'error.html',
        title='Page Not Found',
        lines=['Sorry, but the page you were trying to view does not exist.']), \
        404


@app.errorhandler(500)
def server_error(error):
    log.err('Error: ' + str(error))
    log.err(traceback.format_exc())

    lines = []
    for l in traceback.format_exc().split('\n'):
        if l.strip():
            lines.append(l.replace('    ', '&nbsp;'*4))

    return render_template(
        'error.html',
        title='Internal Server Error',
        error='Error: ' + str(error) + '',
        traceback=traceback.format_exc().split('\n')), \
        500


if __name__ == "__main__":
    main()
    # except Exception as e:
    #     log.err('Error ' + str(e) + ':\n' + traceback.print_exc())
    #     render_template('500.html', error=str(e)), 500
