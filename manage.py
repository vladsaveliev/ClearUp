#!/usr/bin/env python
import os
from collections import defaultdict
from os.path import join, abspath

from flask import Flask
from flask_script import Manager
from sqlalchemy.exc import OperationalError

import ngs_utils.logger as log
from ngs_utils.file_utils import safe_mkdir
from ngs_utils.utils import is_local
from ngs_utils.parallel import ParallelCfg, parallel_view
from ngs_reporting.bcbio.bcbio import BcbioProject

from fingerprinting.model import Project, Sample, db, SNP, get_or_create_run
from fingerprinting import app, DATA_DIR, parallel_cfg

manager = Manager(app)


@manager.command
def load_project(bcbio_dir, name=None):
    log.init(is_debug_=True)

    log.info('-' * 70)
    log.info('Loading project into the fingerprints database from ' + bcbio_dir)
    log.info('-' * 70)
    log.info()

    log.info('Loading bcbio project from file system at ' + bcbio_dir)
    bcbio_proj = BcbioProject()
    bcbio_proj.load_from_bcbio_dir(bcbio_dir, project_name=name,
        proc_name='fingerprinting', need_coverage_interval=False)
    log.info('Loaded ' + bcbio_proj.final_dir)
    log.info()

    work_dir = safe_mkdir(join(bcbio_proj.work_dir, 'fingerprinting'))
    with parallel_view(len(bcbio_proj.samples), parallel_cfg, work_dir) as parall_view:
        log.info('Genotyping sex')
        bcbio_summary_file = bcbio_proj.find_in_log('project-summary.yaml')
        sexes = parall_view.run(_determine_sex, [[s, work_dir, bcbio_summary_file]
                                                 for s in bcbio_proj.samples])
        log.info()
        log.info('Loading fingerprints into the DB')
        fp_proj = Project(
            name=bcbio_proj.project_name,
            bcbio_final_path=bcbio_proj.final_dir,
            genome=bcbio_proj.genome_build,
            panel=bcbio_proj.coverage_bed)
        db.session.add(fp_proj)
        for s, sex in zip(bcbio_proj.samples, sexes):
            db_sample = Sample(s.name, fp_proj, s.bam, sex=sex)
            db.session.add(db_sample)
            
        log.info('Initializing run for single project')
        get_or_create_run(bcbio_proj.project_name, parall_view=parall_view)
    
    db.session.commit()
    log.info()
    log.info('Done.')


def _determine_sex(s, work_dir, bcbio_summary_file=None):
    from ngs_reporting.coverage import get_avg_depth, determine_sex
    from os.path import join
    from ngs_utils.file_utils import safe_mkdir
    if bcbio_summary_file:
        avg_depth = get_avg_depth(bcbio_summary_file, s)
    else:
        avg_depth = 10
    sex = determine_sex(safe_mkdir(join(work_dir, s.name)), s.bam, avg_depth, s.genome_build,
                        target_bed=s.coverage_bed, min_male_size=1)
    return sex


@manager.command
def analyse_projects(run_id):
    log.init(is_debug_=True)
    get_or_create_run(run_id)


@manager.command
def init_db():
    safe_mkdir(DATA_DIR)
    db.init_app(app)
    db.drop_all()
    db.create_all()


@manager.command
def reload_all_data():
    if is_local():
        init_db()
        load_project(abspath('tests/Dev_0261_newstyle'), 'Dev_0261_newstyle')
        load_project(abspath('tests/Dev_0261_newstyle_smallercopy'), 'Dev_0261_newstyle_smallercopy')


if __name__ == "__main__":
    manager.run()

