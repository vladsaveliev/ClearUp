#!/usr/bin/env python
import glob
import traceback
from os.path import join, abspath, basename, splitext, isdir

import sys
from cyvcf2 import VCF
from flask_script import Manager

from ngs_utils.file_utils import safe_mkdir, file_transaction, intermediate_fname, can_reuse, verify_dir
from ngs_utils.utils import is_local, is_us, is_uk
from ngs_utils import logger as log, call_process
from ngs_utils.parallel import ParallelCfg, parallel_view

from clearup.panel import get_dbsnp
from clearup.callable import batch_callable_bed
from clearup.model import Project, Sample, db, SNP, get_or_create_run, Run
from clearup import app, DATA_DIR, parallel_cfg, DEPTH_CUTOFF
from clearup.utils import bam_samplename, get_ref_fasta

manager = Manager(app)
log.init(True)


def _add_project(bam_by_sample, project_name, bed_file=None, use_callable=False,
                 data_dir='', genome='hg19', min_depth=DEPTH_CUTOFF, depth_by_sample=None):
    work_dir = safe_mkdir(join(DATA_DIR, 'projects', project_name))

    with parallel_view(len(bam_by_sample), parallel_cfg, work_dir) as p_view:
        log.info('Initializing run for single project')
        if use_callable:
            log.info('No BED file specified for project ' + project_name + ', calculating callable regions.')
            bed_file = join(work_dir, 'callable_regions.bed')

            genome_fasta_file = get_ref_fasta(genome)
            batch_callable_bed(bam_by_sample.values(), bed_file, work_dir, genome_fasta_file, min_depth,
                               parall_view=p_view)

        fp_proj = Project.query.filter(Project.name==project_name).first()
        if fp_proj:
            db_samples = fp_proj.samples.all()
        else:
            log.info('Loading fingerprints into the DB')
            fp_proj = Project(
                name=project_name,
                data_dir=data_dir,
                genome=genome,
                bed_fpath=bed_file,
                min_depth=min_depth,
                used_callable=use_callable,
            )
            db.session.add(fp_proj)
            db_samples = []
            for sname, bam_file in bam_by_sample.items():
                db_samples.append(Sample(sname, fp_proj, bam_file))
            db.session.add_all(db_samples)
            db.session.commit()

        # get_or_create_run([fp_proj], parall_view=p_view)

        # _add_to_ngb(work_dir, project_name, bam_by_sample, genome, bed_file, p_view)

        log.info('Genotyping sex')
        sex_work_dir = safe_mkdir(join(work_dir, 'sex'))

        # sexes = p_view.run(_sex_from_bam, [
        #     [db_s.name, bam_by_sample[db_s.name], bed_file, sex_work_dir, genome,
        #      depth_by_sample.get(db_s.name) if depth_by_sample else None,
        #      [snp.depth for snp in db_s.snps.all()]]
        #     for db_s in db_samples])
        # for s, sex in zip(db_samples, sexes):
        #     s.sex = sex
    db.session.commit()

    log.info()
    log.info('Done.')


def _add_to_ngb(work_dir, project_name, bam_by_sample, genome_build, bed_file, p_view):
    if is_us() or is_uk():
        try:
            from az.ngb import add_bcbio_project_to_ngb, add_data_to_ngb, add_file_to_ngb
        except ImportError:
            log.warn('If you want to, install NGS Reporting with `conda install -v vladsaveliev ngs_reporting`')
        else:
            log.info('Exposing project to NGB...')
            try:
                dataset = project_name + '_Fingerprints'
                add_data_to_ngb(work_dir, p_view, bam_by_sample, dict(), dataset,
                                bed_file=bed_file, genome=genome_build)
                add_file_to_ngb(work_dir, get_dbsnp(genome_build), genome_build, dataset, dataset,
                                skip_if_added=True)
            except Exception:
                traceback.print_exc()
                log.err('Error: cannot export to NGB')
            log.info('*' * 70)


@manager.command
def load_data(data_dir, name, genome):
    data_dir = verify_dir(data_dir, is_critical=True)
    bam_files = glob.glob(join(data_dir, '*.bam'))
    assert bam_files, 'No BAM files in ' + data_dir

    if name.startswith('--name'):
        name = name.split('--name')[1]

    bed_file = None
    bed_files = glob.glob(join(data_dir, '*.bed'))
    if bed_files:
        assert len(bed_files) == 1, 'Multuple BED files in ' + data_dir + ': ' + str(bed_files)
        bed_file = bed_files[0]

    sample_by_bam = dict()
    for bam_file in bam_files:
        # sample_by_bam[bam_file] = check_output('goleft samplename ' + bam_file)
        sample_by_bam[bam_file] = bam_samplename(bam_file)

    _add_project(
        bam_by_sample={sample_by_bam[bf]: bf for bf in bam_files},
        project_name=name,
        bed_file=bed_file,
        use_callable=not bed_file,
        data_dir=data_dir,
        genome=genome,
        min_depth=DEPTH_CUTOFF)


@manager.command
def load_bcbio_project(bcbio_dir, name=None, use_callable=False):
    try:
        from ngs_reporting.bcbio.bcbio import BcbioProject
    except ImportError:
        log.critical('Error: cannot import ngs_reporting, needed to load a bcbio project. '
                     'Please, install it with `conda install -v vladsaveliev ngs_reporting')
    else:
        log.info('-' * 70)
        log.info('Loading project into the fingerprints database from ' + bcbio_dir)
        log.info('-' * 70)
        log.info()

        bcbio_proj = BcbioProject()
        bcbio_proj.load_from_bcbio_dir(bcbio_dir, project_name=name,
            proc_name='clearup', need_coverage_interval=False, need_vardict=False)

        _add_project(
            bam_by_sample={s.name: s.bam for s in bcbio_proj.samples},
            project_name=name or bcbio_proj.project_name,
            bed_file=bcbio_proj.coverage_bed,
            use_callable=use_callable,
            data_dir=bcbio_proj.final_dir,
            genome=bcbio_proj.genome_build,
            min_depth=DEPTH_CUTOFF,
            depth_by_sample={s.name: s.get_avg_depth() for s in bcbio_proj.samples},
        )


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


# def get_gender(genome, bam_fpath, bed_fpath, sample, avg_depth):
#     gender = None
#     chrom_lengths = ref.get_chrom_lengths(genome)
#     chrom_names = [chrom for chrom, length in chrom_lengths]
#     if 'Y' in chrom_names or 'chrY' in chrom_names:
#         gender = determine_sex(sample.work_dir, bam_fpath, avg_depth, genome, bed_fpath)
#         if gender:
#             with open(join(safe_mkdir(sample.dirpath), 'gender.txt'), 'w') as f:
#                 f.write(gender[0].upper())
#     return gender


        # with open(self.find_in_log('project-summary.yaml')) as f:
        #     data = yaml.load(f)


def _sex_from_bam(sname, bam_file, bed_file, work_dir, genome_build, avg_depth=None, snp_depths=None):
    from os.path import join
    from ngs_utils.file_utils import safe_mkdir
    from ngs_utils.sex import determine_sex
    if avg_depth is None:
        if not snp_depths:
            log.critical('Error: avg_depth is NOT provided and no SNPs in sample ' + sname)
        avg_depth = sum(snp_depths) / len(snp_depths)
    sex = determine_sex(safe_mkdir(join(work_dir, sname)), bam_file, avg_depth,
                        genome_build, target_bed=bed_file)
    return sex


@manager.command
def analyse_projects(project_names_line):
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
    # if verify_dir(DATA_DIR):
    #     os.rename(DATA_DIR, DATA_DIR + '.bak' + datetime.now().strftime("%Y_%m_%d_%H_%M_%S"))
    safe_mkdir(DATA_DIR)
    init_db()
    if is_local():
        load_data(abspath('tests/test_project'), 'SA-1826796__SA-30853', genome='hg38')
        load_bcbio_project(abspath('tests/Dev_0261_newstyle'), 'Dev_0261_newstyle')
        load_bcbio_project(abspath('tests/Dev_0261_newstyle_smallercopy'), 'Dev_0261_newstyle_smallercopy')
        load_bcbio_project(abspath('/Users/vlad/vagrant/NGS_Reporting/tests/results/bcbio_postproc/dream_chr21/final'), 'dream_chr21')
    elif is_us():
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/Resolution/bcbio/final'), 'EXT_070_Plasma_Seq_Pilot_Resolution', use_callable=True)
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/Foundation/bcbio/final'), 'EXT_070_Plasma_Seq_Pilot_Foundation')
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/Foundation/plasma/bcbio/final'), 'EXT_070_Plasma_Seq_Pilot_Foundation_plasma', use_callable=True)
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/Foundation/tissue/OurType/bcbio_complete'), 'EXT_070_Plasma_Seq_Pilot_Foundation_tissue')
        load_bcbio_project(abspath('/ngs/oncology/analysis/external/EXT_070_Plasma_Seq_Pilot/PGDx/bcbio/final'), 'EXT_070_Plasma_Seq_Pilot_PGDx')
        load_bcbio_project(abspath('/ngs/oncology/analysis/dev/Dev_0327_MiSeq_SNP251/bcbio_preprint/final'), 'Dev_0327_MiSeq_SNP251_initial_preprint')
        load_bcbio_project(abspath('/ngs/oncology/analysis/dev/Dev_0320_HiSeq4000_PARPiResistant_Exome/bcbio/final'), 'Dev_0320_HiSeq4000_PARPiResistant_Exome')

        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0288_HiSeq4000_PlasmaSeqExome/bcbio/final'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0287_HiSeq4000_PlasmaSeqAZ100/bcbio_umi/final_vardict_1.4.8'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0300_HiSeq4000_PlasmaSeq_Tissue_AZ100/bcbio_umi/final'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0306_HiSeq4000_PlasmaSeq_Tissue_Exome/bcbio_umi/final'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0308_HiSeq4000_MerckFFPE_AZ100/bcbio_umi_102/final'))
        # load_project(abspath('/ngs/oncology/analysis/dev/Dev_0309_HiSeq4000_PlasmaSeq_Tissue_Exome/bcbio_umi/final'))
        # load_project(abspath('/ngs/oncology/Analysis/dev/Dev_0310_HiSeq4000_PlasmaSeq_Tissue_AZ100/bcbio_umi_102/final'))


@manager.command
def dump_projects(output_file):
    projects = set()
    for run in Run.query.all():  # Finding projects with ready-to-view runs
        for p in run.projects:
            projects.add(p)
    for p in sorted(list(projects), key=lambda p: p.name):
        if 'final' in p.data_dir:
            cmdl = './manage.py load_bcbio_project ' + p.data_dir + ' --name=' + p.name
        else:
            cmdl = './manage.py load_data ' + p.data_dir + ' ' + p.name + ' --genome=' + p.genome
        if 'used_callable' in p.__dict__ and p.used_callable:
            cmdl += ' --used_callable'
        with open(output_file, 'w') as f:
            f.write(cmdl + '\n')


if __name__ == '__main__':
    manager.run()
