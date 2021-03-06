#!/usr/bin/env python
import shutil
from os.path import abspath, join, dirname, splitext, basename, isfile, isdir
from collections import defaultdict
from cyvcf2 import VCF
from pybedtools import BedTool

from ngs_utils.Sample import BaseSample
from ngs_utils.file_utils import safe_mkdir, can_reuse, file_transaction, verify_file
from ngs_utils.parallel import ParallelCfg, parallel_view
from ngs_utils import logger as log

from clearup.panel import build_snps_panel
from clearup.genotype import genotype, build_tree, build_snp_from_records
from clearup.utils import FASTA_ID_PROJECT_SEPARATOR, load_bam_file
from clearup import app, db, DATA_DIR, parallel_cfg

# import logging
# logging.basicConfig()
# logging.getLogger('sqlalchemy.engine').setLevel(logging.INFO)
safe_mkdir(DATA_DIR)


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
    ref = db.Column(db.String)

    run_id = db.Column(db.String, db.ForeignKey('run.id'))
    run = db.relationship('Run', backref=db.backref('locations', lazy='dynamic', cascade='all, delete-orphan'))

    def __init__(self, rsid, chrom=None, pos=None, gene=None, ref=None):
        self.rsid = rsid
        self.chrom = chrom
        self.pos = pos
        self.gene = gene
        self.ref = ref

    def __repr__(self):
        return '<Location {}:{} {} at gene {}>'.format(self.chrom, str(self.pos), self.rsid, self.gene)


class SNP(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    rsid = db.Column(db.String)
    chrom = db.Column(db.String)
    pos = db.Column(db.Integer)
    gene = db.Column(db.String)

    allele1 = db.Column(db.String(1))
    allele2 = db.Column(db.String(1))
    depth = db.Column(db.Integer)
    allele1_depth = db.Column(db.Integer)
    allele2_depth = db.Column(db.Integer)
    usercall = db.Column(db.String)

    location_id = db.Column(db.String, db.ForeignKey('location.id'))
    location = db.relationship('Location')

    sample_id = db.Column(db.String, db.ForeignKey('sample.id'))
    sample = db.relationship('Sample', backref=db.backref('snps', lazy='dynamic', cascade='all, delete-orphan'))

    def __init__(self, loc):
        self.location = loc
        self.rsid = loc.rsid
        self.chrom = loc.chrom
        self.pos = loc.pos
        self.gene = loc.gene

    def __repr__(self):
        return '<SNP {}:{} {} {}|{} for sample {}>'.format(
            str(self.location.chrom), str(self.location.pos), self.location.rsid,
            self.allele1, self.allele2, self.sample.name)

    def get_gt(self):
        if self.usercall:
            return self.usercall
        assert self.allele1 and self.allele2, str(self)
        return ''.join(sorted([self.allele1, self.allele2]))


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
    projects = db.relationship("Project",
                               secondary=run_to_project_assoc_table,
                               backref=db.backref('runs', lazy='dynamic'),
                               lazy='dynamic')
    rerun_on_usercall = db.Column(db.Boolean, default=False)

    def __init__(self):
        self.snps_file = None
        self.rerun_on_usercall = False

    def work_dir_path(self):
        return safe_mkdir(join(DATA_DIR, str(self.id)))

    def fasta_file_path(self):
        return join(self.work_dir_path(), 'fingerprints.fasta')

    def tree_file_path(self):
        return join(self.work_dir_path(), 'fingerprints.newick')

    @staticmethod
    def is_ready(run):
        return run and \
               isfile(run.fasta_file_path()) and \
               isfile(run.tree_file_path())

    @staticmethod
    def create(projects, parall_view=None):
        run = Run()
        db.session.add(run)
        for p in projects:
            run.projects.append(p)
        db.session.commit()

        genomes = [p.genome for p in projects]
        if len(set(genomes)) > 1:
            log.critical('Error: multiple genomes in projects: ' + str(genomes))
        genome_build = genomes[0]

        snps_dir = safe_mkdir(join(run.work_dir_path(), 'snps'))
        run.snps_file = build_snps_panel(bed_files=[p.bed_fpath for p in projects if p.bed_fpath],
                                         output_dir=snps_dir, genome=genome_build)
        locations = extract_locations_from_file(run.snps_file)
        for loc in locations:
            db.session.add(loc)
        db.session.commit()

        log.info()
        log.info('Genotyping')
        samples = [s for p in projects for s in p.samples]
        snps_left_to_call_file = _get_snps_not_called(run.snps_file, samples)
        vcf_dir = safe_mkdir(join(run.work_dir_path(), 'vcf'))
        work_dir = safe_mkdir(join(vcf_dir, 'work'))
        bs = [BaseSample(s.long_name(), bam=s.bam) for s in samples]

        if parall_view:
            vcf_by_sample = genotype(bs, snps_left_to_call_file, parall_view,
                 work_dir=work_dir, output_dir=vcf_dir, genome_build=genome_build)
        else:
            n_threads = parallel_cfg.threads
            if len(samples) < n_threads:  # vardict is running in 1 thread
                parallel_cfg.threads = len(samples)
            with parallel_view(len(samples), parallel_cfg, safe_mkdir(join(run.work_dir_path(), 'log'))) as view:
                vcf_by_sample = genotype(bs, snps_left_to_call_file, view,
                     work_dir=work_dir, output_dir=vcf_dir, genome_build=genome_build)
            parallel_cfg.threads = n_threads


        # TODO: speed this up
        log.info('Loading called SNPs into the DB')
        for s in samples:
            recs = [r for r in VCF(vcf_by_sample[s.long_name()])]
            recs_by_rsid = defaultdict(list)
            for r in recs:
                recs_by_rsid[r.ID].append(r)
            for loc in locations:
                assert loc
                snp = s.snps.filter(SNP.rsid==loc.rsid).first()
                if not snp:
                    snp = SNP(loc)
                    build_snp_from_records(snp, recs_by_rsid[loc.rsid], s.project.min_depth)
                    s.snps.append(snp)
                    db.session.add(snp)

        log.info('Adding locations into the DB')
        run.locations.delete()
        for l in locations:
            run.locations.append(l)
        db.session.add(run)
        db.session.commit()
        log.info('Saved locations in the DB')

        log.info()
        log.info('Building tree')
        build_tree(run)

        log.info()
        log.info('Loading BAMs sliced to fingerprints')
        if parall_view:
            parall_view.run(load_bam_file,
                            [[s.bam, safe_mkdir(join(run.work_dir_path(), 'bams')), run.snps_file, s.long_name()]
                             for s in samples])
        else:
            with parallel_view(len(samples), parallel_cfg, safe_mkdir(join(run.work_dir_path(), 'log'))) as view:
                view.run(load_bam_file,
                         [[s.bam, safe_mkdir(join(run.work_dir_path(), 'bams')), run.snps_file, s.long_name()]
                          for s in samples])

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

    def delete(self):
        if isdir(self.work_dir_path()):
            shutil.rmtree(self.work_dir_path())
        db.session.delete(self)


def _genotype(run, samples, genome_build, parall_view):
    snps_left_to_call_file = _get_snps_not_called(run.snps_file, samples)

    vcf_dir = safe_mkdir(join(run.work_dir_path(), 'vcf'))
    work_dir = safe_mkdir(join(vcf_dir, 'work'))
    bs = [BaseSample(s.long_name(), bam=s.bam) for s in samples]
    vcf_by_sample = genotype(bs, snps_left_to_call_file, parall_view,
                             work_dir=work_dir, output_dir=vcf_dir, genome_build=genome_build)
    return vcf_by_sample


def get_or_create_run(projects, parall_view=None):
    genomes = set([p.genome for p in projects])
    if len(genomes) > 1:
        log.critical('Error: multiple genomes in projects: ' + str(genomes))
    run = Run.find_by_projects(projects)

    if run and run.rerun_on_usercall:
        log.info()
        log.info('Rebuilding tree on usercall')
        build_tree(run)
        run.rerun_on_usercall = False
        db.session.commit()
        return run

    if run and not Run.is_ready(run):
        log.debug('Tree files do not exist, recreating run for projects ' + ', '.join(p.name for p in projects))
        db.session.delete(run)
        db.session.commit()
        run = None

    if run:
        log.debug('Found run for ' + ', '.join([p.name for p in projects]) + ' with ID ' + str(run.id))
    else:
        log.debug('Creating new run for projects ' + ', '.join(p.name for p in projects))
        run = Run.create(projects, parall_view)
        log.debug('Done creating new run with ID ' + str(run.id))
    return run


class Project(db.Model):
    __tablename__ = 'project'
    name = db.Column(db.String(), primary_key=True)
    data_dir = db.Column(db.String)
    genome = db.Column(db.String(20))
    bed_fpath = db.Column(db.String)
    min_depth = db.Column(db.Integer)
    used_callable = db.Column(db.Boolean, default=False)

    def __init__(self, name, data_dir, genome, bed_fpath, min_depth, used_callable):
        self.name = name
        self.data_dir = data_dir
        self.genome = genome
        self.bed_fpath = bed_fpath
        self.min_depth = min_depth
        self.used_callable = used_callable

    def get_work_dir(self):
        return join(DATA_DIR, 'projects', self.name)

    def __repr__(self):
        return '<Project {} {}>'.format(self.name, self.genome)

    def delete(self, reuse_files=True):
        if isdir(self.get_work_dir()) and not reuse_files:
            shutil.rmtree(self.get_work_dir())
        for r in self.runs:
            r.delete()
        db.session.delete(self)


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
        """ Warning: returns unordered! Use run.locations to get order
        """
        locs_ids = set([l.rsid for l in run.locations.all()])
        # sq = run.locations.subquery()
        snps = self.snps.filter(SNP.rsid.in_(locs_ids))
        # snp = s.snps.join(Location).filter(Location.rsid==loc.rsid).first()
        # snps = [snp for snp in self.snps if snp.rsid in locs_ids]
        snps_by_rsid = {snp.rsid: snp for snp in snps}
        return snps_by_rsid

    def __repr__(self):
        return '<Sample {} from project {}>'.format(self.name, self.project.name)


def extract_locations_from_file(snps_file):
    locs = []
    for i, interval in enumerate(BedTool(snps_file)):
        pos = int(interval.start) + 1
        rsid, gene, ref, alts = interval.name.split('|')
        loc = Location(
            rsid=rsid,
            chrom=interval.chrom,
            pos=pos,
            gene=gene,
            ref=ref)
        locs.append(loc)
    return locs
