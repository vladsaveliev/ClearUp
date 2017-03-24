#!/usr/bin/env python
import os
from collections import defaultdict
from os.path import join, abspath
from flask_script import Manager
from sqlalchemy.exc import OperationalError

from ngs_utils.utils import is_local
from start import app

import ngs_utils.logger as log
from ngs_reporting.bcbio.bcbio import BcbioProject

from fingerprinting.model import Project, Sample, db, SNP, get_or_create_run
from fingerprinting.utils import load_bam_file, get_snps_by_type

manager = Manager(app)


@manager.command
def load_project(bcbio_dir, project_name, genome, panel_type):
    log.init(is_debug_=True)
    
    get_snps_by_type(panel_type)  # Just to verify the type correctness

    log.info('-' * 70)
    log.info('Loading project ' + project_name + ' into the fingerprints database with SNPs ' + panel_type)
    log.info('-' * 70)
    log.info()

    log.info('Loading bcbio project from file system at ' + bcbio_dir)
    bcbio_proj = BcbioProject()
    bcbio_proj.load_from_bcbio_dir(bcbio_dir, project_name=project_name,
        proc_name='fingerprinting', need_coverage_interval=False)
    log.info('Loaded ' + bcbio_proj.final_dir)
    
    log.info('Loading fingerprints into the DB')
    fp_proj = Project(
        name=bcbio_proj.project_name,
        bcbio_final_path=bcbio_proj.final_dir,
        genome=genome,
        panel_type=panel_type)
    db.session.add(fp_proj)
    for s in bcbio_proj.samples:
        db_sample = Sample(s.name, fp_proj, s.bam)
        db.session.add(db_sample)
        
    log.info('Initializing run for single project')
    get_or_create_run(project_name)
    
    db.session.commit()
    log.info()
    log.info('Done.')


@manager.command
def analyse_projects(run_id):
    log.init(is_debug_=True)
    get_or_create_run(run_id)


@manager.command
def init_db():
    db.drop_all()
    db.create_all()


@manager.command
def reload_all_data():
    if is_local():
        db.drop_all()
        db.create_all()
        load_project(abspath('tests/Dev_0261_newstyle'), 'Dev_0261_newstyle', 'hg19', 'idt')
        load_project(abspath('tests/Dev_0261_newstyle_smallercopy'), 'Dev_0261_newstyle_smallercopy', 'hg19', 'idt')


if __name__ == "__main__":
    manager.run()

