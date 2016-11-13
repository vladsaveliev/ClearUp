#!/usr/bin/env python2
from genericpath import isfile

from os.path import join
import sys
import glob
from flask_script import Manager

from fingerprinting.model import Project, Sample, db
from fingerprinting.utils import read_fasta
from start import app

from Utils.file_utils import verify_file, verify_dir
import Utils.logger as log
IS_DEBUG = log.is_debug = True


manager = Manager(app)


@manager.command
def hello():
    print "hello"


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
    for sample_id, seq in seq_by_sample_id.items():
        sample = Sample(sample_id, project, seq)
        db.session.add(sample)
    db.session.commit()

    log.info('Done')



if __name__ == "__main__":
    manager.run()

