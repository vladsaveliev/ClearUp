#!/usr/bin/env python2
from genericpath import isfile

from os.path import join, abspath
import sys
import glob
from flask_script import Manager

from fingerprinting import config
from fingerprinting.model import Project, Sample, db
from fingerprinting.utils import read_fasta, load_bam_file, get_fingerprints
from start import app

from Utils.file_utils import verify_file, verify_dir, safe_mkdir
import Utils.logger as log
IS_DEBUG = log.is_debug = True


manager = Manager(app)


@manager.command
def load_project(bcbio_final_path, project_name, genome):
    log.info('Loading project ' + project_name)

    if not verify_dir(bcbio_final_path):
        log.critical('Project path doesn\'t exist ' + bcbio_final_path)
    if 'final' not in bcbio_final_path.split('/')[-1]:
        log.critical('Directory must be a bcbio final dir ' + bcbio_final_path)

    fasta_glob = join(bcbio_final_path, '20??-??-??*', 'fingerprints', '*.fasta')
    fasta_fpath = next(iter(glob.glob(fasta_glob)), None)
    if not verify_file(fasta_fpath):
        log.critical('Fingerprints fasta file not found in ' + fasta_glob)

    project = Project(project_name, bcbio_final_path, fasta_fpath, genome)
    # project = Project.query.filter_by(name=project_name).first()
    seq_by_sample_id = read_fasta(fasta_fpath)
    db.session.add(project)
    project_work_dirpath = safe_mkdir(join(config.DATA_DIR, project_name))
    for sample_id, seq in seq_by_sample_id.items():
        bam_fpath = load_bam_file(bcbio_final_path, project_work_dirpath, sample_id)
        sample = Sample(sample_id, project, bam_fpath)
        db.session.add(sample)
        fingerprints = get_fingerprints(project, sample, seq._data.upper())
        for fingerprint in fingerprints:
            db.session.add(fingerprint)
    db.session.commit()

    log.info('Done')


@manager.command
def init_db():
    db.drop_all()
    db.create_all()


@manager.command
def reload_all_data():
    if log.is_local():
        db.drop_all()
        db.create_all()
        load_project(abspath('../test/analysis/dev/Dev_0261_MiSeq_MCRC_PRCC/bcbio/final'), 'Dev_0261_MiSeq_MCRC_PRCC', 'hg19')
        load_project(abspath('../test/analysis/dev/Dev_0261_MiSeq_COPY/bcbio/final'), 'Dev_0261_MiSeq_COPY', 'hg19')


if __name__ == "__main__":
    manager.run()

