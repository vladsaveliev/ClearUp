#!/usr/bin/env python

from collections import defaultdict
from os.path import join, abspath
import vcf
from flask_script import Manager
from sqlalchemy.exc import OperationalError
from start import app

import ngs_utils.logger as log
from ngs_utils.file_utils import verify_file, verify_dir, safe_mkdir
from ngs_utils.parallel import ParallelCfg
from ngs_reporting.bcbio.bcbio import BcbioProject
import az

from fingerprinting import config
from fingerprinting.genotype import genotype_bcbio_proj, vcfrec_to_seq, DEPTH_CUTOFF
from fingerprinting.model import Project, Sample, db, Fingerprint
from fingerprinting.utils import load_bam_file, get_snps_file


manager = Manager(app)


@manager.command
def load_project(bcbio_dir, project_name, genome):
    log.init(True)
    log.info('-' * 70)
    log.info('Loading project ' + project_name + ' into the fingerprints database')
    log.info('-' * 70)
    log.info()

    sys_cfg = az.init_sys_cfg()
    parallel_cfg = ParallelCfg(threads=sys_cfg.get('threads'))

    log.info('Loading bcbio project from file system at ' + bcbio_dir)
    bcbio_proj = BcbioProject()
    bcbio_proj.load_from_bcbio_dir(bcbio_dir, project_name=project_name,
                                   proc_name='fingerprinting', need_coverage_interval=False)
    log.info('Loaded ' + bcbio_proj.final_dir)
    
    log.info('Genotyping')
    fp_proj_dir = safe_mkdir(join(config.DATA_DIR, project_name))
    work_dir = safe_mkdir(join(fp_proj_dir, 'work'))
    all_fasta, vcf_by_sample, sex_by_sample = genotype_bcbio_proj(
        bcbio_proj, get_snps_file(), parallel_cfg,
        output_dir=fp_proj_dir,
        work_dir=work_dir)

    log.info('Loading BAMs sliced to fingerprints')
    for s in bcbio_proj.samples:
        load_bam_file(s.bam, safe_mkdir(join(fp_proj_dir, 'bams')), s.name)
    
    log.info('Loading fingerprints into the DB')
    fp_proj = Project(bcbio_proj.project_name, bcbio_proj.final_dir, all_fasta, genome)
    db.session.add(fp_proj)
    for s in bcbio_proj.samples:
        db_sample = Sample(s.name, fp_proj, sex=sex_by_sample[s.name])
        db.session.add(db_sample)
        with open(vcf_by_sample[s.name]) as vcf_f:
            vcf_reader = vcf.Reader(vcf_f)
            recs = [r for r in vcf_reader]
            for i, rec in enumerate(recs):
                fp = Fingerprint(index=i + 1,
                                 sample=db_sample,
                                 chrom=rec.CHROM,
                                 pos=rec.POS,
                                 rsid=rec.ID,
                                 gene=rec.INFO['GENE'])
                fp.depth = rec.samples[0]['DP']
                fp.genotype = vcfrec_to_seq(rec, DEPTH_CUTOFF)
                db.session.add(fp)
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

