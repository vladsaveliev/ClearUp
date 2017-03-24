#!/usr/bin/env python

from os.path import abspath, join, dirname, splitext, basename
from flask import Flask, render_template, send_from_directory, abort, redirect, url_for, send_file, request

from gevent.pywsgi import WSGIServer
from geventwebsocket.handler import WebSocketHandler

from fingerprinting import config
from fingerprinting.model import Project, db, Sample
from fingerprinting.sample_view import render_closest_comparison_page, send_file_for_igv
from fingerprinting.tree_view import run_analysis_socket_handler, render_phylo_tree_page

from ngs_utils import logger


app = Flask(__name__)
# import flask_assets
# assets = flask_assets.Environment()
# assets.init_app(app)


@app.route('/favicon.ico/')
def send_favicon():
    return send_from_directory('static', 'favicon.ico')


@app.route("/<run_id>/run_analysis/")
def run_analysis(run_id):
    return run_analysis_socket_handler(run_id)


@app.route('/<run_id>/tree/')
def phylo_tree_page(run_id):
    return render_phylo_tree_page(run_id)


@app.route('/<run_id>/tree/<int:sample_id>/')
def closest_comparison_page(run_id, sample_id):
    return render_closest_comparison_page(run_id, sample_id, request.args.get('snpIndex'))


@app.route('/<run_id>/tree/<int:sample_id>/add_usercall/', methods=['POST'])
def add_user_call(run_id, sample_id):
    logger.info('Adding user call for ' + str(sample_id))
    edit_sample_id = request.form['editSampleId']
    sample = Sample.query.filter_by(id=edit_sample_id).first()
    if not sample:
        logger.err('Sample with ID=' + str(edit_sample_id) + ' not found')
        return redirect(url_for('closest_comparison_page', run_id=run_id, sample_id=sample_id))

    snp = sample.snps.filter_by(index=request.form['snpIndex']).first()
    snp.usercall = request.form['usercall']
    db.session.commit()
    return redirect(url_for('closest_comparison_page', run_id=run_id, sample_id=sample_id,
                            snpIndex=snp.index))


@app.route('/<project_name>/bamfiles/<bam_fname>/')
def bam_files_page(project_name, bam_fname):
    return send_file_for_igv(join(config.DATA_DIR, project_name, 'bams', bam_fname))


@app.route('/snps/snps_bed/')
def locations_bed():
    return send_file(join('fingerprinting', 'snps', 'idt_snps.bed'))


@app.route('/<run_id>/')
def project_page(run_id):
    return redirect(url_for('phylo_tree_page', **locals()))


@app.route('/<run_id>/<sample_id>/')
def sample_page(run_id, sample_id):
    return redirect(url_for('closest_comparison_page', **locals()))


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
            } for s in p.samples]
        } for p in projects],
    )
    return t


@app.errorhandler(404)
def page_not_found(error):
    return render_template('404.html', error=error.description['message']), 404


if __name__ == "__main__":
    logger.init(True, join(dirname(__file__), 'data', 'log.txt'))
    # app.run(host=config.HOST_IP, debug=config.IS_DEBUG)

    # if start_local_browser:
        # start server and web page pointing to it
        # url = "http://{HOST}:{PORT}".format(HOST=config.HOST_IP, PORT=PORT)
        # wb = webbrowser.get(None)  # instead of None, can be "firefox" etc
        # threading.Timer(1.25, lambda: wb.open(url)).start()

    http_server = WSGIServer((config.HOST_IP, config.PORT), app, handler_class=WebSocketHandler)
    logger.info('Starting a webserver at ' + config.HOST_IP + ':' + str(config.PORT))
    http_server.serve_forever()
