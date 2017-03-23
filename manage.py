#!/usr/bin/env python

from collections import defaultdict
from os.path import join, abspath
from flask_script import Manager
from sqlalchemy.exc import OperationalError
from start import app

import ngs_utils.logger as log
from ngs_reporting.bcbio.bcbio import BcbioProject

from fingerprinting.model import Project, Sample, db, SNP, get_run
from fingerprinting.utils import load_bam_file, get_snps_by_type

manager = Manager(app)


@manager.command
def load_project(bcbio_dir, project_name, genome, panel_type):
    log.init(True)
    
    get_snps_by_type(panel_type)  # Just to verify the type correctness

    log.info('-' * 70)
    log.info('Loading project ' + project_name + ' into the fingerprints database with SNPs ' + panel_type)
    log.info('-' * 70)
    log.info()

    log.info('Loading bcbio project from file system at ' + bcbio_dir)
    bcbio_proj = BcbioProject()
    bcbio_proj.load_from_bcbio_dir(bcbio_dir, project_name=project_name,
        proc_name='fingerprinting', need_coverage_interval=False)
    log.info('Loaded ' + bcbio_proj.final_dir)
    
    log.info('Loading fingerprints into the DB')
    fp_proj = Project(
        name=bcbio_proj.project_name,
        bcbio_final_path=bcbio_proj.final_dir,
        genome=genome,
        panel_type=panel_type)
    db.session.add(fp_proj)
    for s in bcbio_proj.samples:
        db_sample = Sample(s.name, fp_proj, s.bam)
        db.session.add(db_sample)
        
    log.info('Initializing run for single project')
    get_run(project_name)
    
    db.session.commit()
    log.info()
    log.info('Done.')


# @manager.command
# def load_project(bcbio_dir, project_name, genome):
#     log.init(True)
#     log.info('-' * 70)
#     log.info('Loading project ' + project_name + ' into the fingerprints database')
#     log.info('-' * 70)
#     log.info()
#
#     sys_cfg = az.init_sys_cfg()
#     parallel_cfg = ParallelCfg(threads=sys_cfg.get('threads'))
#
#     log.info('Loading bcbio project from file system at ' + bcbio_dir)
#     bcbio_proj = BcbioProject()
#     bcbio_proj.load_from_bcbio_dir(bcbio_dir, project_name=project_name,
#                                    proc_name='fingerprinting', need_coverage_interval=False)
#     log.info('Loaded ' + bcbio_proj.final_dir)
#
#     log.info('Genotyping')
#     fp_proj_dir = safe_mkdir(join(config.DATA_DIR, project_name))
#     work_dir = safe_mkdir(join(fp_proj_dir, 'work'))
#     all_fasta, vcf_by_sample, sex_by_sample = genotype_bcbio_proj(
#         bcbio_proj, get_snps_file(), parallel_cfg,
#         output_dir=fp_proj_dir,
#         work_dir=work_dir)
#
#     log.info('Loading BAMs sliced to fingerprints')
#     for s in bcbio_proj.samples:
#         load_bam_file(s.bam, safe_mkdir(join(fp_proj_dir, 'bams')), s.name)
#
#     log.info('Loading fingerprints into the DB')
#     fp_proj = Project(bcbio_proj.project_name, bcbio_proj.final_dir, all_fasta, genome)
#     db.session.add(fp_proj)
#     for s in bcbio_proj.samples:
#         db_sample = Sample(s.name, fp_proj, sex=sex_by_sample[s.name])
#         db.session.add(db_sample)
#         with open(vcf_by_sample[s.name]) as vcf_f:
#             vcf_reader = vcf.Reader(vcf_f)
#             recs = [r for r in vcf_reader]
#             for i, rec in enumerate(recs):
#                 fp = SNP(index=i + 1,
#                          sample=db_sample,
#                          chrom=rec.CHROM,
#                          pos=rec.POS,
#                          rsid=rec.ID,
#                          gene=rec.INFO['GENE'])
#                 fp.depth = rec.samples[0]['DP']
#                 fp.genotype = vcfrec_to_seq(rec, DEPTH_CUTOFF)
#                 db.session.add(fp)
#     db.session.commit()
#     log.info()
#     log.info('Done.')


@manager.command
def init_db():
    db.drop_all()
    db.create_all()


@manager.command
def reload_all_data():
    if log.is_local():
        db.drop_all()
        db.create_all()
        load_project(abspath('tests/Dev_0261_newstyle'), 'Dev_0261_newstyle', 'hg19', 'idt')
        load_project(abspath('tests/Dev_0261_newstyle_smallercopy'), 'Dev_0261_newstyle_smallercopy', 'hg19', 'idt')


if __name__ == "__main__":
    manager.run()

