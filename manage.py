#!/usr/bin/env python
from collections import defaultdict
from os.path import join, abspath
import sys
import glob

import vcf
from flask_script import Manager
from ngs_utils.parallel import ParallelCfg

from fingerprinting import config
from fingerprinting.genotype import genotype_bcbio_dir, vcfrec_to_seq, DEPTH_CUTOFF, check_if_male
from fingerprinting.model import Project, Sample, db, Fingerprint
from fingerprinting.utils import read_fasta, load_bam_file, get_fingerprints, get_snps_file
from start import app

import az

from ngs_utils.file_utils import verify_file, verify_dir, safe_mkdir
import ngs_utils.logger as log
IS_DEBUG = log.is_debug = True

THREADS = 20


manager = Manager(app)


@manager.command
def load_project(bcbio_final_path, project_name, genome):
    log.info('-' * 70)
    log.info('Loading project ' + project_name)
    if not verify_dir(bcbio_final_path):
        log.critical('Project path doesn\'t exist ' + bcbio_final_path)
    if 'final' not in bcbio_final_path.split('/')[-1]:
        log.critical('Directory must be a bcbio final dir ' + bcbio_final_path)

    sys_cfg = az.init_sys_cfg()
    sys_cfg['threads'] = THREADS
    parallel_cfg = ParallelCfg(sys_cfg.get('scheduler'), sys_cfg.get('queue'),
                               sys_cfg.get('resources'), sys_cfg.get('threads'))
    parallel_cfg.set_tag('fingerprinting')

    project_work_dirpath = safe_mkdir(join(config.DATA_DIR, project_name))
    # TODO: use project_work_dirpath to write VCF and fasta
    all_fasta, vcf_by_sample = genotype_bcbio_dir(bcbio_final_path, get_snps_file(), sys_cfg, parallel_cfg)

    project = Project(project_name, bcbio_final_path, all_fasta, genome)
    db.session.add(project)
    seq_by_sample_name = read_fasta(all_fasta)
    for sname in seq_by_sample_name.keys():
        db.session.add(Sample(sname, project))

    log.info('')
    log.info('Reading SNPs...')
    fp_by_loc_by_sample = defaultdict(dict)
    fp_by_index_by_sample = defaultdict(dict)
    with open(get_snps_file()) as bed_f:
        for i, l in enumerate(l for l in bed_f if l[0] != '#'):
            chrom, pos0, pos1, ann = l.strip().split()
            for s in project.samples:
                fp = Fingerprint(index=i + 1, sample=s, chrom=chrom, pos=int(pos1), rsid=ann.split('-')[1])
                fp_by_loc_by_sample[s.name][(chrom, int(pos1))] = fp
                fp_by_index_by_sample[s.name][i + 1] = fp

    log.info('Populating genotypes and depth...')
    for s in project.samples:
        with open(vcf_by_sample[s.name]) as vcf_f:
            vcf_reader = vcf.Reader(vcf_f)
            recs = [r for r in vcf_reader]
            is_male = check_if_male(recs)
            for rec in recs:
                fp = fp_by_loc_by_sample[s.name][(rec.CHROM, rec.POS)]
                fp.depth = rec.samples[0]['DP']
                fp.genotype = vcfrec_to_seq(rec, is_male, DEPTH_CUTOFF)

    log.info('Pushing fingerprints into the database...')
    for s in project.samples:
        load_bam_file(bcbio_final_path, project_work_dirpath, s.name)
        fingerprints = sorted(fp_by_loc_by_sample[s.name].values(), key=lambda _fp: _fp.index)
        for fingerprint in fingerprints:
            db.session.add(fingerprint)

    db.session.commit()
    log.info()
    log.info('Done.')


@manager.command
def init_db():
    db.drop_all()
    db.create_all()


@manager.command
def reload_all_data():
    if log.is_local():
        db.drop_all()
        db.create_all()
        load_project(abspath('tests/Dev_0261_newstyle/final'), 'Dev_0261_newstyle', 'hg19')
        load_project(abspath('tests/Dev_0261_oldstyle_smallercopy/final'), 'Dev_0261_oldstyle_smallercopy', 'hg19')


if __name__ == "__main__":
    manager.run()

