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
    fingerprint = db.Column(db.String())

    project_name = db.Column(db.String, db.ForeignKey('project.name'))
    project = db.relationship('Project', backref=db.backref('samples', lazy='dynamic'))

    def __init__(self, name, project, fingerprint=None):
        self.name = name
        self.project = project
        # self.fingerprint = fingerprint

    def __repr__(self):
        return '<Sample {} from project {}>'.format(self.name, self.project.name)


if __name__ == '__main__':
    db.create_all()




