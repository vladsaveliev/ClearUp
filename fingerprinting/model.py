#!/usr/bin/env python

from copy import copy
from os.path import abspath, join, dirname, splitext, basename

from cyvcf2 import VCF
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

from Bio import SeqIO
from pybedtools import BedTool
import vcf

from ngs_utils.Sample import BaseSample
from ngs_utils.file_utils import safe_mkdir, can_reuse, file_transaction
from ngs_utils.logger import info, debug, err
from ngs_utils.parallel import ParallelCfg, parallel_view

from fingerprinting.panel import build_snps_panel
from fingerprinting.genotype import genotype, vcfrec_to_seq, DEPTH_CUTOFF, write_fasta
from fingerprinting.utils import FASTA_ID_PROJECT_SEPARATOR, load_bam_file
from fingerprinting import app, db, DATA_DIR, parallel_cfg

import logging
logging.basicConfig()
logging.getLogger('sqlalchemy.engine').setLevel(logging.INFO)


run_to_project_assoc_table = db.Table(
    'run_to_project_association', db.Model.metadata,
    db.Column('run_id', db.Integer, db.ForeignKey('run.id')),
    db.Column('project_name', db.String, db.ForeignKey('project.name')))


class Location(db.Model):
    __tablename__ = 'location'
    id = db.Column(db.Integer, primary_key=True)
    rsid = db.Column(db.String)
    chrom = db.Column(db.String)
    pos = db.Column(db.Integer)
    gene = db.Column(db.String)
    # ref = db.Column(db.String)
    # alt = db.Column(db.String)
    
    run_id = db.Column(db.String, db.ForeignKey('run.id'))
    run = db.relationship('Run', backref=db.backref('locations', lazy='dynamic'))

    def __init__(self, rsid, chrom=None, pos=None, gene=None):
        self.rsid = rsid
        self.chrom = chrom
        self.pos = pos
        self.gene = gene
        # self.ref = None
        # self.alt = None

    def __repr__(self):
        return '<Location {}:{} {} at gene {}>'.format(self.chrom, str(self.pos), self.rsid, self.gene)


class SNP(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    rsid = db.Column(db.String)
    chrom = db.Column(db.String)
    pos = db.Column(db.Integer)
    gene = db.Column(db.String)

    genotype = db.Column(db.String)
    depth = db.Column(db.Integer)
    usercall = db.Column(db.String)

    location_id = db.Column(db.String, db.ForeignKey('location.id'))
    location = db.relationship('Location')

    sample_id = db.Column(db.String, db.ForeignKey('sample.id'))
    sample = db.relationship('Sample', backref=db.backref('snps', lazy='dynamic'))

    def __init__(self, loc):
        self.location = loc
        self.rsid = loc.rsid
        self.chrom = loc.chrom
        self.pos = loc.pos
        self.gene = loc.gene
        
    def __repr__(self):
        return '<SNP {}:{} {} {} for sample {}>'.format(
            str(self.location.chrom), str(self.location.pos), self.location.rsid,
            self.genotype, self.sample.name)


def _get_snps_not_called(snps_file, samples):
    # TODO: select and save only snps per sample that are not called in taht sample (special treatment for gender?)
    # lines_to_rerun = []
    # for i, interval in enumerate(BedTool(snps_file)):
    #     rsid, gene = interval.name.split('|')
    #     if all(s.snps.filter_by(snps=loc.id) for s in samples)
    #                 snp = .first()
    return snps_file


PROJ_COLORS = ['#000000', '#90ed7d', '#f7a35c', '#8085e9', '#f15c80',
               '#e4d354', '#2b908f', '#f45b5b', '#91e8e1', '#7cb5ec']

class Run(db.Model):
    __tablename__ = 'run'
    id = db.Column(db.Integer, primary_key=True)
    snps_file = db.Column(db.String)
    projects = db.relationship("Project", secondary=run_to_project_assoc_table,
                               backref=db.backref('runs', lazy='dynamic'), lazy='dynamic')
    
    def __init__(self):
        self.snps_file = None
     
    def work_dir_path(self):
        return safe_mkdir(join(DATA_DIR, str(self.id)))
     
    def fasta_file_path(self):
        return join(self.work_dir_path(), 'fingerprints.fasta')

    def tree_file_path(self):
        return join(self.work_dir_path(), 'fingerprints.newick')
    
    @staticmethod
    def create(projects, parall_view=None):
        run = Run()
        db.session.add(run)
        for p in projects:
            run.projects.append(p)
        db.session.commit()
        
        genome_builds = [p.genome for p in projects]
        assert len(set(genome_builds)) == 1, 'Error: different genome builds in projects'
        genome_build = genome_builds[0]
        
        snps_dir = safe_mkdir(join(run.work_dir_path(), 'snps'))
        run.snps_file = build_snps_panel(bed_files=[p.bed_fpath for p in projects if p.bed_fpath],
                                         output_dir=snps_dir, genome_build=genome_build)
        locations = extract_locations_from_file(run.snps_file)
        for loc in locations:
            db.session.add(loc)
        db.session.commit()
        location_by_rsid = {l.rsid: l for l in locations}
        
        info()
        info('Genotyping')
        samples = [s for p in projects for s in p.samples]
        if parall_view:
            fasta_file, vcf_by_sample = _genotype(run, samples, genome_build, parall_view)
        else:
            samples = [s for p in projects for s in p.samples]
            with parallel_view(len(samples), parallel_cfg, safe_mkdir(join(run.work_dir_path(), 'log'))) as parall_view:
                fasta_file, vcf_by_sample = _genotype(run, samples, genome_build, parall_view)

        info('Loading called SNPs into the DB')
        for s in samples:
            recs = [r for r in VCF(vcf_by_sample[s.long_name()])]
            for i, rec in enumerate(recs):
                loc = location_by_rsid[rec.ID]
                assert loc.pos == rec.POS
                # snp = s.snps.join(Location).filter(Location.rsid==loc.rsid).first()
                snp = s.snps.filter(SNP.rsid==rec.ID).first()
                if snp:
                    assert snp.depth == rec.INFO['DP']
                    assert snp.genotype == vcfrec_to_seq(rec, DEPTH_CUTOFF)
                else:
                    snp = SNP(loc)
                    snp.depth = rec.INFO['DP']
                    snp.genotype = vcfrec_to_seq(rec, DEPTH_CUTOFF)
                    s.snps.append(snp)
                    db.session.add(snp)

        info('Adding locations into the DB')
        for l in run.locations:
            l.run = None
        for snp in samples[0].snps:
            if snp.location.rsid in location_by_rsid:
                run.locations.append(snp.location)
        db.session.commit()
        info('Saved locations in the DB')
        return run
    
    @staticmethod
    def find_by_projects(projects):
        for r in Run.query.all():
            if sorted(tuple(p.name for p in projects)) == sorted(tuple(p.name for p in r.projects)):
                return r
        return None

    @staticmethod
    def find_by_project_names_line(project_names_line):
        pnames = project_names_line.split('--')
        projects = Project.query.filter(Project.name.in_(pnames))
        return Run.find_by_projects(projects)


def _genotype(run, samples, genome_build, parall_view):
    snps_left_to_call_file = _get_snps_not_called(run.snps_file, samples)

    vcf_dir = safe_mkdir(join(run.work_dir_path(), 'vcf'))
    work_dir = safe_mkdir(join(vcf_dir, 'work'))
    bs = [BaseSample(s.long_name(), bam=s.bam) for s in samples]
    vcf_by_sample = genotype(bs, snps_left_to_call_file, parall_view,
                             work_dir=work_dir, output_dir=vcf_dir, genome_build=genome_build)
    
    info()
    info('** Building fasta **')
    fasta_dir = safe_mkdir(join(run.work_dir_path(), 'fasta'))
    work_dir = safe_mkdir(join(fasta_dir, 'fasta'))
    fasta_file, vcf_by_sample = write_fasta(bs, vcf_by_sample, snps_left_to_call_file,
                                            parall_view, output_dir=fasta_dir, work_dir=work_dir, out_fasta=run.fasta_file_path())

    info('Loading BAMs sliced to fingerprints')
    parall_view.run(load_bam_file,
        [[s.bam, safe_mkdir(join(run.work_dir_path(), 'bams')), run.snps_file, s.long_name()]
         for s in samples])
    return fasta_file, vcf_by_sample
    

def get_or_create_run(projects, parall_view=None):
    run = Run.find_by_projects(projects)
    if not run:
        debug('Creating new run for projects ' + ', '.join(p.name for p in projects))
        run = Run.create(projects, parall_view)
        db.session.add(run)
        db.session.commit()
        debug('Done creating new run with ID ' + str(run.id))
    else:
        debug('Found run for ' + ', '.join([p.name for p in projects]) + ' with ID ' + str(run.id))
    return run


# def merge_fasta(projects, work_dirpath):
#     fasta_fpaths = [p.fingerprints_fasta_fpath for p in projects]
#     merged_fasta_fpath = join(work_dirpath, 'fingerprints.fasta')
#     if not can_reuse(merged_fasta_fpath, fasta_fpaths):
#         all_records = []
#         for proj, fasta in zip(projects, fasta_fpaths):
#             with open(fasta) as f:
#                 recs = SeqIO.parse(f, 'fasta')
#                 for rec in recs:
#                     rec.id += FASTA_ID_PROJECT_SEPARATOR + proj.name
#                     rec.name = rec.description = ''
#                     all_records.append(rec)
#         with file_transaction(None, merged_fasta_fpath) as tx:
#             with open(tx, 'w') as out:
#                 SeqIO.write(all_records, out, 'fasta')
#     return merged_fasta_fpath


class Project(db.Model):
    __tablename__ = 'project'
    name = db.Column(db.String(), primary_key=True)
    bcbio_final_path = db.Column(db.String)
    genome = db.Column(db.String(20))
    bed_fpath = db.Column(db.String)
    
    def __init__(self, name, bcbio_final_path, genome, bed_fpath):
        self.name = name
        self.bcbio_final_path = bcbio_final_path
        self.genome = genome
        self.bed_fpath = bed_fpath

    def __repr__(self):
        return '<Project {} {}>'.format(self.name, self.genome)


class Sample(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String)
    bam = db.Column(db.String)
    sex = db.Column(db.String)

    project_name = db.Column(db.String, db.ForeignKey('project.name'))
    project = db.relationship('Project', backref=db.backref('samples', lazy='dynamic'))

    def __init__(self, name, project, bam, sex=None):
        self.name = name
        self.project = project
        self.bam = bam
        self.sex = sex

    def long_name(self):
        return self.name + FASTA_ID_PROJECT_SEPARATOR + self.project.name

    def snps_from_run(self, run):
        locs_ids = set([l.rsid for l in run.locations.all()])
        info(locs_ids)
        # sq = run.locations.subquery()
        snps = self.snps.filter(SNP.rsid.in_(locs_ids))
        # snp = s.snps.join(Location).filter(Location.rsid==loc.rsid).first()
        # snps = [snp for snp in self.snps if snp.rsid in locs_ids]
        return snps
    
    def __repr__(self):
        return '<Sample {} from project {}>'.format(self.name, self.project.name)


def extract_locations_from_file(snps_file):
    locs = []
    for i, interval in enumerate(BedTool(snps_file)):
        pos = int(interval.start) + 1
        rsid, gene = interval.name.split('|')
        loc = Location(
            rsid=rsid,
            chrom=interval.chrom,
            pos=pos,
            gene=gene)
        locs.append(loc)
    return locs


if __name__ == '__main__':
    db.create_all()
