#!/usr/bin/env python
import os
from collections import defaultdict
from os.path import join, abspath

from cyvcf2 import VCF
from flask import Flask
from flask_script import Manager

import ngs_utils.logger as log
from ngs_utils import sambamba
from ngs_utils.bed_utils import get_total_bed_size
from ngs_utils.file_utils import safe_mkdir, file_transaction, intermediate_fname, can_reuse
from ngs_utils.utils import is_local, is_us
from ngs_utils.parallel import ParallelCfg, parallel_view
from ngs_reporting.bcbio.bcbio import BcbioProject

from fingerprinting.model import Project, Sample, db, SNP, get_or_create_run
from fingerprinting import app, DATA_DIR, parallel_cfg

manager = Manager(app)


@manager.command
def load_project(bcbio_dir, name=None):
    log.init(is_debug_=True)

    log.info('-' * 70)
    log.info('Loading project into the fingerprints database from ' + bcbio_dir)
    log.info('-' * 70)
    log.info()

    log.info('Loading bcbio project from file system at ' + bcbio_dir)
    bcbio_proj = BcbioProject()
    bcbio_proj.load_from_bcbio_dir(bcbio_dir, project_name=name,
        proc_name='fingerprinting', need_coverage_interval=False, need_vardict=False)
    work_dir = safe_mkdir(join(DATA_DIR, 'bcbio_projects', bcbio_proj.project_name))
    log.info('Loaded ' + bcbio_proj.final_dir)
    log.info()

    log.info()
    log.info('Loading fingerprints into the DB')
    fp_proj = Project(
        name=bcbio_proj.project_name,
        bcbio_final_path=bcbio_proj.final_dir,
        genome=bcbio_proj.genome_build,
        bed_fpath=bcbio_proj.coverage_bed)
    db.session.add(fp_proj)
    db_samples = []
    for s in bcbio_proj.samples:
        db_samples.append(Sample(s.name, fp_proj, s.bam))
        db.session.add(db_samples[-1])
    db.session.commit()
    
    log.info('Initializing run for single project')
    with parallel_view(len(bcbio_proj.samples), parallel_cfg, work_dir) as parall_view:
        get_or_create_run(bcbio_proj.project_name, parall_view=parall_view)
    
    log.info('Genotyping sex')
    sex_work_dir = safe_mkdir(join(work_dir, 'sex'))
    with parallel_view(len(bcbio_proj.samples), parallel_cfg, sex_work_dir) as parall_view:
        # sex_snp_bed = build_snps_panel(bcbio_projs=[bcbio_proj], output_dir=sex_work_dir,
        #                                genome_build=bcbio_proj.genome_build, take_autosomal=False, take_sex=True)
        # vcf_by_sample = genotype(bcbio_proj.samples, sex_snp_bed, parall_view, sex_work_dir, bcbio_proj.genome_build)
        # for s in db_samples:
        #     s.sex = _sex_from_x_snps(vcf_by_sample[s.name])
        bcbio_summary_file = bcbio_proj.find_in_log('project-summary.yaml')
        sexes = parall_view.run(_sex_from_bam, [[db_s, s, sex_work_dir, bcbio_summary_file]
                                                for db_s, s in zip(db_samples, bcbio_proj.samples)])
        for s, sex in zip(db_samples, sexes):
            s.sex = sex
    db.session.commit()
    
    log.info()
    log.info('Done.')


def _sex_from_x_snps(vcf_file):
    log.debug('Calling sex from ' + vcf_file)
    het_calls_num = 0
    hom_calls_num = 0
    for rec in VCF(vcf_file):
        if rec.CHROM == 'chrX':
            if rec.num_het > 0:
                het_calls_num += 1
            if rec.num_hom > 0:
                hom_calls_num += 1
    
    if het_calls_num + hom_calls_num > 10:
        if het_calls_num > 1.5 * hom_calls_num:
            return 'F'
        elif het_calls_num < 0.5 * hom_calls_num:
            return 'M'
        else:
            log.debug('het/hom ratio on chrX is ' + str(het_calls_num/hom_calls_num) +
                      ' - between 1.5 and 0.5, not confident enough to call sex.')
    else:
        log.debug('Total chrX calls number is ' + str(het_calls_num + hom_calls_num) +
                  ' - less than 10, not confident enough to call sex.')
    return None
    

def _sex_from_bam(db_sample, sample, work_dir, bcbio_summary_file=None):
    from ngs_reporting.coverage import get_avg_depth, determine_sex
    from os.path import join
    from ngs_utils.file_utils import safe_mkdir
    avg_depth = None
    if bcbio_summary_file:
        avg_depth = get_avg_depth(bcbio_summary_file, sample.name)
    if avg_depth is None:
        depths = [snp.depth for snp in db_sample.snps.all()]
        avg_depth = sum(depths) / len(depths)
    sex = determine_sex(safe_mkdir(join(work_dir, sample.name)), sample.bam, avg_depth,
                        sample.genome_build, target_bed=sample.coverage_bed, min_male_size=1)
    return sex


@manager.command
def analyse_projects(run_id):
    log.init(is_debug_=True)
    get_or_create_run(run_id)


@manager.command
def init_db():
    safe_mkdir(DATA_DIR)
    db.init_app(app)
    db.drop_all()
    db.create_all()


@manager.command
def reload_all_data():
    init_db()
    if is_local():
        load_project(abspath('tests/Dev_0261_newstyle'), 'Dev_0261_newstyle')
        load_project(abspath('tests/Dev_0261_newstyle_smallercopy'), 'Dev_0261_newstyle_smallercopy')
        load_project(abspath('/Users/vlad/vagrant/NGS_Reporting/tests/results/bcbio_postproc/dream_chr21/final'), 'dream_chr21')
    elif is_us():
        load_project(abspath(''), '')


if __name__ == "__main__":
    manager.run()

