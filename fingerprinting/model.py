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
from fingerprinting.genotype import genotype, vcfrec_to_seq, DEPTH_CUTOFF, post_genotype
from fingerprinting.utils import FASTA_ID_PROJECT_SEPARATOR, load_bam_file
from fingerprinting import app, db, DATA_DIR, parallel_cfg

run_to_project_assoc_table = db.Table(
    'run_to_project_association', db.Model.metadata,
    db.Column('run_id', db.Integer, db.ForeignKey('run.id')),
    db.Column('project_name', db.String, db.ForeignKey('project.name')))


class Location(db.Model):
    __tablename__ = 'location'
    id = db.Column(db.Integer, primary_key=True)
    index = db.Column(db.Integer)
    rsid = db.Column(db.String)
    chrom = db.Column(db.String)
    pos = db.Column(db.Integer)
    gene = db.Column(db.String)
    # ref = db.Column(db.String)
    # alt = db.Column(db.String)
    
    run_id = db.Column(db.String, db.ForeignKey('run.id'))
    run = db.relationship('Run', backref=db.backref('locations', lazy='dynamic'))

    def __init__(self, rsid, index, chrom=None, pos=None, gene=None):
        self.rsid = rsid
        self.index = index
        self.chrom = chrom
        self.pos = pos
        self.gene = gene
        # self.ref = None
        # self.alt = None

    def __repr__(self):
        return '<Location {}:{} {} at gene {}>'.format(self.chrom, str(self.pos), self.rsid, self.gene)


class SNP(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    index = db.Column(db.Integer)
    genotype = db.Column(db.String)
    depth = db.Column(db.Integer)
    usercall = db.Column(db.String)

    location_id = db.Column(db.String, db.ForeignKey('location.id'))
    location = db.relationship('Location')

    sample_id = db.Column(db.String, db.ForeignKey('sample.id'))
    sample = db.relationship('Sample', backref=db.backref('snps', lazy='dynamic'))

    def __init__(self, index, location=None):
        self.index = index
        self.location = location
        self.genotype = None
        self.depth = None
        self.usercall = None

    def __repr__(self):
        return '<SNP {}:{} {} {} for sample {}>'.format(
            str(self.location.chrom), str(self.location.pos), self.location.rsid,
            self.genotype, self.sample.name)


def _get_snps_not_calls(snps_file, samples):
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
    id = db.Column(db.String, primary_key=True)
    snps_file = db.Column(db.String)
    work_dir = db.Column(db.String)
    projects = db.relationship("Project", secondary=run_to_project_assoc_table,
                               backref=db.backref('runs', lazy='dynamic'), lazy='dynamic')
    
    def __init__(self, project_names):
        self.id = ','.join(project_names)
        self.snps_file = None
        self.work_dir = safe_mkdir(join(DATA_DIR, '__AND__'.join(project_names)))
     
    def fasta_file_path(self):
        return join(self.work_dir, 'fingerprints.fasta')

    def tree_file_path(self):
        return join(self.work_dir, 'fingerprints.newick')
    
    @staticmethod
    def create(projects, parall_view=None):
        project_names = sorted([p.name for p in projects])
        run = Run(project_names)
        db.session.add(run)
        for p in projects:
            run.projects.append(p)
        db.session.commit()
        
        genome_builds = [p.genome for p in projects]
        assert len(set(genome_builds)) == 1, 'Error: different genome builds in projects'
        genome_build = genome_builds[0]
        
        snps_dir = safe_mkdir(join(run.work_dir, 'snps'))
        run.snps_file = build_snps_panel(bed_files=[p.bed_fpath for p in projects if p.bed_fpath],
                                         output_dir=snps_dir, genome_build=genome_build)
        locations = extract_locations_from_file(run.snps_file)
        for loc in locations:
            db.session.add(loc)
        db.session.commit()
        location_by_rsid = {l.rsid: l for l in locations}
        
        info()
        info('Genotyping')
        gt_work_dir = safe_mkdir(join(run.work_dir, 'genotyping'))
        samples = [s for p in projects for s in p.samples]
        if parall_view:
            fasta_file, vcf_by_sample = _genotype(run, samples, genome_build, parall_view, work_dir=gt_work_dir)
        else:
            samples = [s for p in projects for s in p.samples]
            with parallel_view(len(samples), parallel_cfg, gt_work_dir) as parall_view:
                fasta_file, vcf_by_sample = _genotype(run, samples, genome_build, parall_view, work_dir=gt_work_dir)

        info('Loading called SNPs into the DB')
        for s in samples:
            recs = [r for r in VCF(vcf_by_sample[s.long_name()])]
            for i, rec in enumerate(recs):
                loc = location_by_rsid[rec.INFO['rsid']]
                assert loc.pos == rec.POS
                snp = s.snps.join(Location).filter(Location.rsid==loc.rsid).first()
                if snp:
                    assert snp.depth == rec.INFO['DP']
                    assert snp.genotype == vcfrec_to_seq(rec, DEPTH_CUTOFF)
                else:
                    snp = SNP(index=i + 1, location=loc)
                    snp.depth = rec.INFO['DP']
                    snp.genotype = vcfrec_to_seq(rec, DEPTH_CUTOFF)
                    s.snps.append(snp)
                    db.session.add(snp)

        info('Adding locations into the DB')
        for snp in samples[0].snps:
            if snp.location.rsid in location_by_rsid:
                run.locations.append(snp.location)
        db.session.commit()
        return run


def _genotype(run, samples, genome_build, parall_view, work_dir=None):
    snps_left_to_call_file = _get_snps_not_calls(run.snps_file, samples)

    gt_work_dir = work_dir or safe_mkdir(join(run.work_dir, 'genotyping'))
    bs = [BaseSample(s.long_name(), bam=s.bam) for s in samples]
    vcf_by_sample = genotype(bs, snps_left_to_call_file, parall_view,
                             output_dir=gt_work_dir, genome_build=genome_build)
    
    info()
    info('** Post-genotyping **')
    fasta_file, vcf_by_sample = post_genotype(bs, vcf_by_sample, snps_left_to_call_file,
        parall_view, output_dir=run.work_dir, work_dir=gt_work_dir, out_fasta=run.fasta_file_path())

    info('Loading BAMs sliced to fingerprints')
    parall_view.run(load_bam_file,
        [[s.bam, safe_mkdir(join(run.work_dir, 'bams')), run.snps_file, s.long_name()]
         for s in samples])
    return fasta_file, vcf_by_sample
    

def get_or_create_run(run_id, parall_view=None):
    run = Run.query.filter_by(id=run_id).first()
    if not run:
        debug('Creating new run ' + run_id)
        project_names = run_id.split(',')
        projects = []
        not_found = []
        for pn in project_names:
            p = Project.query.filter_by(name=pn).first()
            if p:
                projects.append(p)
            else:
                not_found.append(pn)
        if not_found:
            err('Some projects not found in the database: ' + ', '.join(project_names))
            return None
        run = Run.create(projects, parall_view)
        db.session.add(run)
        db.session.commit()
        debug('Done creating new run ' + run_id)
    else:
        debug('Found run ' + run.id)
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
        snps = []
        for snp in self.snps:
            if snp.location.rsid in locs_ids:
                snps.append(snp)
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
            index=i + 1,
            chrom=interval.chrom,
            pos=pos,
            gene=gene)
        locs.append(loc)
    return locs


if __name__ == '__main__':
    db.create_all()
