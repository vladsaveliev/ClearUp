#!/usr/bin/env python
import glob
import os
from collections import defaultdict
from logging import critical
from os.path import join, abspath, basename, splitext

from cyvcf2 import VCF
from datetime import datetime
from flask import Flask
from flask_script import Manager

from ngs_utils.file_utils import safe_mkdir, file_transaction, intermediate_fname, can_reuse, verify_dir
from ngs_utils.utils import is_local, is_us
from ngs_utils import logger as log
from ngs_utils.parallel import ParallelCfg, parallel_view
from ngs_reporting.bcbio.bcbio import BcbioProject

from fingerprinting.model import Project, Sample, db, SNP, get_or_create_run
from fingerprinting import app, DATA_DIR, parallel_cfg

manager = Manager(app)


@manager.command
def _add_project(bam_by_sample, bed_file, project_name, data_dir='', genome_build='hg19', bcbio_summary_file=None):
    work_dir = safe_mkdir(join(DATA_DIR, 'projects', project_name))

    log.info('Loading fingerprints into the DB')
    fp_proj = Project(
        name=project_name,
        data_dir=data_dir,
        genome=genome_build,
        bed_fpath=bed_file)
    db.session.add(fp_proj)
    db_samples = []
    for sname, bam_file in bam_by_sample.items():
        db_samples.append(Sample(sname, fp_proj, bam_file))
        db.session.add(db_samples[-1])
    db.session.commit()
    
    log.info('Initializing run for single project')
    with parallel_view(len(bam_by_sample), parallel_cfg, work_dir) as parall_view:
        get_or_create_run([fp_proj], parall_view=parall_view)
    
    log.info('Genotyping sex')
    sex_work_dir = safe_mkdir(join(work_dir, 'sex'))
    with parallel_view(len(bam_by_sample), parallel_cfg, sex_work_dir) as parall_view:
        sexes = parall_view.run(_sex_from_bam, [
            [db_s, bam_file, bed_file, sex_work_dir, genome_build, bcbio_summary_file]
            for db_s, bam_file in zip(db_samples, bam_by_sample.values())])
        for s, sex in zip(db_samples, sexes):
            s.sex = sex
    db.session.commit()
    
    log.info()
    log.info('Done.')


@manager.command
def load_data(data_dir, project_name):
    data_dir = verify_dir(data_dir, is_critical=True)
    bam_files = glob.glob(join(data_dir, '*.bam'))
    assert bam_files, 'No BAM files in ' + data_dir
    bed_files = glob.glob(join(data_dir, '*.bed'))
    assert len(bed_files) == 1, 'Multuple BED files in ' + data_dir + ': ' + str(bed_files)
    
    _add_project(
        bam_by_sample={splitext(basename(bf))[0]: bf for bf in bam_files},
        bed_file=bed_files[0],
        project_name=project_name,
        data_dir=data_dir)
    

@manager.command
def load_bcbio_project(bcbio_dir, name=None):
    log.info('-' * 70)
    log.info('Loading project into the fingerprints database from ' + bcbio_dir)
    log.info('-' * 70)
    log.info()

    bcbio_proj = BcbioProject()
    bcbio_proj.load_from_bcbio_dir(bcbio_dir, project_name=name,
        proc_name='fingerprinting', need_coverage_interval=False, need_vardict=False)

    _add_project(
        bam_by_sample={s.name: s.bam for s in bcbio_proj.samples},
        bed_file=bcbio_proj.coverage_bed,
        project_name=name,
        data_dir=bcbio_proj.final_dir,
        genome_build=bcbio_proj.genome_build,
        bcbio_summary_file=bcbio_proj.find_in_log('project-summary.yaml'))


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
    

def _sex_from_bam(db_sample, bam_file, bed_file, work_dir, genome_build, bcbio_summary_file=None):
    from os.path import join
    from ngs_utils.file_utils import safe_mkdir
    from ngs_reporting.coverage import get_avg_depth, determine_sex
    avg_depth = None
    if bcbio_summary_file:
        avg_depth = get_avg_depth(bcbio_summary_file, db_sample.name)
    if avg_depth is None:
        depths = [snp.depth for snp in db_sample.snps.all()]
        if not depths:
            critical('Error: no SNPs in sample ' + db_sample.long_name())
        avg_depth = sum(depths) / len(depths)
    sex = determine_sex(safe_mkdir(join(work_dir, db_sample.name)), bam_file, avg_depth,
                        genome_build, target_bed=bed_file)
    return sex


@manager.command
def analyse_projects(project_names_line):
    log.init(is_debug_=True)
    project_names = project_names_line.split('--')
    projects = Project.query.filter(Project.name.in_(project_names))
    if projects.count() < len(project_names):
        raise RuntimeError('Some projects in ' + str(project_names) + ' are not found in the database: ' +
                           str(set(project_names) - set(p.name for p in projects)))
    get_or_create_run(projects)


@manager.command
def init_db():
    safe_mkdir(DATA_DIR)
    db.init_app(app)
    db.drop_all()
    db.create_all()


@manager.command
def reload_all_data():
    os.rename(DATA_DIR, DATA_DIR + '.bak' + datetime.now().strftime("%Y_%m_%d_%H_%M_%S"))
    safe_mkdir(DATA_DIR)
    init_db()
    if is_local():
        load_bcbio_project(abspath('tests/Dev_0261_newstyle'), 'Dev_0261_newstyle')
        load_bcbio_project(abspath('tests/Dev_0261_newstyle_smallercopy'), 'Dev_0261_newstyle_smallercopy')
        load_bcbio_project(abspath('/Users/vlad/vagrant/NGS_Reporting/tests/results/bcbio_postproc/dream_chr21/final'), 'dream_chr21')
    elif is_us():
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/bcbio/final'), 'EXT_070_Plasma_Seq_Pilot_Resolution')
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/Foundation/bcbio/final'), 'EXT_070_Plasma_Seq_Pilot_Foundation')
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/Foundation/plasma/bcbio/final'), 'EXT_070_Plasma_Seq_Pilot_Foundation_plasma')
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/Foundation/tissue/OurType/bcbio_complete'), 'EXT_070_Plasma_Seq_Pilot_Foundation_tissue')
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/PGDx/bcbio/final'), 'EXT_070_Plasma_Seq_Pilot_PGDx')

        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0288_HiSeq4000_PlasmaSeqExome/bcbio/final'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0287_HiSeq4000_PlasmaSeqAZ100/bcbio_umi/final_vardict_1.4.8'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0300_HiSeq4000_PlasmaSeq_Tissue_AZ100/bcbio_umi/final'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0306_HiSeq4000_PlasmaSeq_Tissue_Exome/bcbio_umi/final'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0308_HiSeq4000_MerckFFPE_AZ100/bcbio_umi_102/final'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0309_HiSeq4000_PlasmaSeq_Tissue_Exome/bcbio_umi/final'))
        # load_project(abspath('/ngs/oncology/Analysis/dev/Dev_0310_HiSeq4000_PlasmaSeq_Tissue_AZ100/bcbio_umi_102/final'))


if __name__ == "__main__":
    manager.run()
