#!/usr/bin/env python2

from os.path import abspath, join, dirname
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from logging import info
from fingerprinting import config


app = Flask(__name__)
db_path = join(config.DATA_DIR, 'projects.db')
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:////' + db_path
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)


class Fingerprint(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    index = db.Column(db.Integer())
    pos = db.Column(db.Integer())
    chrom = db.Column(db.Integer())
    nucl = db.Column(db.String())
    usercall = db.Column(db.String())

    sample_id = db.Column(db.String, db.ForeignKey('sample.id'))
    sample = db.relationship('Sample', backref=db.backref('fingerprints', lazy='dynamic'))

    def __init__(self, sample, index, chrom, pos, nucl):
        self.sample = sample
        self.index = index
        self.chrom = chrom
        self.pos = pos
        self.nucl = nucl
        self.usercall = None

    def __repr__(self):
        return '<Fingerprint {}-{} for sample {}>'.format(str(self.pos), self.nucl, self.sample.name)


class Project(db.Model):
    name = db.Column(db.String(), primary_key=True)
    bcbio_final_path = db.Column(db.String())
    fingerprints_fasta_fpath = db.Column(db.String())
    genome = db.Column(db.String(20))

    def __init__(self, name, bcbio_final_path, fp_fpath, genome):
        self.name = name
        self.bcbio_final_path = bcbio_final_path
        self.fingerprints_fasta_fpath = fp_fpath
        self.genome = genome
        self.date = None
        self.region = None

    def __repr__(self):
        return '<Project {}>'.format(self.name)


class Sample(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String())
    bam_fpath = db.Column(db.String())
    paired_sample_id = db.Column(db.Integer())

    project_name = db.Column(db.String, db.ForeignKey('project.name'))
    project = db.relationship('Project', backref=db.backref('samples', lazy='dynamic'))

    def __init__(self, name, project, bam_fpath):
        self.name = name
        self.project = project
        self.bam_fpath = bam_fpath

    def __repr__(self):
        return '<Sample {} from project {}>'.format(self.name, self.project.name)


class PairedSample(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    project = db.Column(db.String())
    sample_id = db.Column(db.Integer())

    paired_sample_id = db.Column(db.String, db.ForeignKey('sample.id'))
    paired_sample = db.relationship('Sample', backref=db.backref('paired_samples', lazy='dynamic'))

    def __init__(self, name, project, sample_id, paired_sample):
        self.name = name
        self.project = project
        self.sample_id = sample_id
        self.paired_sample = paired_sample

    def __repr__(self):
        return '<Sample {} is matched with {}>'.format(self.name, self.sample.name)


if __name__ == '__main__':
    db.create_all()




