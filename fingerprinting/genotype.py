from collections import OrderedDict

import os
from os.path import join, dirname, splitext
import sys
import vcf
from subprocess import check_output
from variant_filtering.filtering import combine_vcfs
from vcf.parser import field_counts, _Info

from ngs_utils.bed_utils import bgzip_and_tabix
from ngs_utils.call_process import run
from ngs_utils.utils import is_local, is_us
from ngs_utils.parallel import ParallelCfg, parallel_view
from pybedtools import BedTool

from ngs_utils.file_utils import file_transaction, safe_mkdir, chdir, which, adjust_path, can_reuse, add_suffix, \
    verify_file, intermediate_fname, splitext_plus
from ngs_utils.logger import info, err, critical, debug
from ngs_utils.sambamba import index_bam
from ngs_reporting.coverage import get_gender, determine_sex
import az

from fingerprinting.utils import is_sex_chrom


DEPTH_CUTOFF = 5


def genotype(samples, snp_bed, parallel_cfg, output_dir, work_dir, genome_build,
             depth_cutoff=DEPTH_CUTOFF):
    autosomal_bed, sex_bed = _split_bed(snp_bed, work_dir)
    
    genome_cfg = az.get_refdata(genome_build)
    info('** Running VarDict ** ')
    with parallel_view(len(samples), parallel_cfg, work_dir) as view:
        vcfs = view.run(_vardict_pileup_sample,
            [[s, safe_mkdir(join(output_dir, 'vcf')), genome_cfg, view.cores_per_job, autosomal_bed]
             for s in samples])
    vcf_by_sample = OrderedDict(zip([s.name for s in samples], vcfs))
    info('** Finished running VarDict **')
    
    info('** Annotate variants with gene names and rsIDs **')
    for s in samples:
        vcf_by_sample[s.name] = bgzip_and_tabix(_annotate_vcf(vcf_by_sample[s.name], snp_bed))

    info('** Writing fasta **')
    fasta_by_sample = OrderedDict((s.name, join(work_dir, s.name + '.fasta')) for s in samples)
    for s in samples:
        info('Writing fasta for sample ' + s.name)
        vcf_to_fasta(s, vcf_by_sample[s.name], fasta_by_sample[s.name], depth_cutoff)

    info('** Merging fasta **')
    all_fasta = join(output_dir, 'fingerprints.fasta')
    if not can_reuse(all_fasta, fasta_by_sample.values()):
        with open(all_fasta, 'w') as out_f:
            for s in samples:
                with open(fasta_by_sample[s.name]) as f:
                    out_f.write(f.read())
    info('All fasta saved to ' + all_fasta)

    info('** Determining sexes **')
    sex_by_sample = dict()
    for s in samples:
        avg_depth = calc_avg_depth(vcf_by_sample[s.name])
        sex = determine_sex(safe_mkdir(join(work_dir, s.name)), s.bam, avg_depth, genome_build,
                            target_bed=sex_bed, min_male_size=1)
        sex_by_sample[s.name] = sex
        
    info('Combining VCFs')
    combined_vcf = combine_vcfs(vcf_by_sample.values(), join(output_dir, 'fingerprints.vcf'))
    assert combined_vcf

    # ped_file = splitext_plus(combined_vcf)[0]
    # cmdl = 'vcftools --plink --gzvcf {combined_vcf}'.format(**locals())
    # cmdl += ' --out ' + splitext(ped_file)[0]
    # run(cmd=cmdl)
    info('Converting VCFs to PED')
    ped_file = vcf_to_ped(vcf_by_sample, join(output_dir, 'fingerprints.ped'), sex_by_sample, depth_cutoff)
    
    return all_fasta, vcf_by_sample, sex_by_sample


def vcf_to_ped(vcf_by_sample, ped_file, sex_by_sample, depth_cutoff):
    if can_reuse(ped_file, vcf_by_sample.values()):
        return ped_file

    with file_transaction(None, ped_file) as tx:
        with open(tx, 'w') as out_f:
            for sname, vcf_file in vcf_by_sample.items():
                with open(vcf_file) as vcf_f:
                    vcf_reader = vcf.Reader(vcf_f)
                    recs = [r for r in vcf_reader]
                    sex = {'M': '1', 'F': '2'}.get(sex_by_sample[sname], '0')
                    out_fields = [sname, sname, '0', '0', sex, '0']
                    for rec in recs:
                        gt = vcfrec_to_seq(rec, depth_cutoff).replace('N', '0')
                        out_fields.extend([gt[0], gt[1]])
                    out_f.write('\t'.join(out_fields) + '\n')

    info('PED saved to ' + ped_file)
    return ped_file


def genotype_bcbio_proj(proj, snp_bed, parallel_cfg, depth_cutoff=DEPTH_CUTOFF,
                        output_dir=None, work_dir=None):
    output_dir = output_dir or safe_mkdir(join(proj.date_dir, 'fingerprints'))
    work_dir = work_dir or safe_mkdir(proj.work_dir)
    return genotype(proj.samples, snp_bed, parallel_cfg, output_dir, work_dir,
                    proj.genome_build, depth_cutoff=depth_cutoff)


def _split_bed(bed_file, work_dir):
    """ Splits into autosomal and sex chromosomes
    """
    autosomal_bed = intermediate_fname(work_dir, bed_file, 'autosomal')
    sex_bed = intermediate_fname(work_dir, bed_file, 'sex')
    if not can_reuse(autosomal_bed, bed_file) or not can_reuse(sex_bed, bed_file):
        with open(bed_file) as f, open(autosomal_bed, 'w') as a_f, open(sex_bed, 'w') as s_f:
            for l in f:
                chrom = l.split()[0]
                if is_sex_chrom(chrom):
                    s_f.write(l)
                else:
                    a_f.write(l)
    return autosomal_bed, sex_bed


def _vardict_pileup_sample(sample, output_dir, genome_cfg, threads, snp_file):
    vardict_snp_vars = join(output_dir, sample.name + '_vars.txt')
    vardict_snp_vars_vcf = join(output_dir, sample.name + '.vcf')
    
    if can_reuse(vardict_snp_vars, [sample.bam, snp_file]) and can_reuse(vardict_snp_vars_vcf, vardict_snp_vars):
        return vardict_snp_vars_vcf
    
    if is_local():
        vardict_dir = '/Users/vlad/vagrant/VarDict/'
    elif is_us():
        vardict_dir = '/group/cancer_informatics/tools_resources/NGS/bin/'
    else:
        vardict_pl = which('vardict.pl')
        if not vardict_pl:
            critical('Error: vardict.pl is not in PATH')
        vardict_dir = dirname(vardict_pl)

    # Run VarDict
    index_bam(sample.bam)
    vardict = join(vardict_dir, 'vardict.pl')
    ref_file = adjust_path(genome_cfg['seq'])
    cmdl = '{vardict} -G {ref_file} -N {sample.name} -b {sample.bam} -p -D {snp_file}'.format(**locals())
    run(cmdl, output_fpath=vardict_snp_vars)
    
    # Convert to VCF
    test_strandbias = join(vardict_dir, 'teststrandbias.R')
    var2vcf_valid = join(vardict_dir, 'var2vcf_valid.pl')
    cmdl = ('cut -f-34 {vardict_snp_vars}'
            ' | {test_strandbias}'
            ' | {var2vcf_valid}'
            ' | grep "^#\|TYPE=SNV\|TYPE=REF" '
            ).format(**locals())
    run(cmdl, output_fpath=vardict_snp_vars_vcf)
    _fix_vcf(vardict_snp_vars_vcf, ref_file)
    
    return vardict_snp_vars_vcf


def _fix_vcf(vardict_snp_vars_vcf, ref_file):
    """ Fixes VCF generated by VarDict in puleup debug mode:
        - Fix non-call records with empty REF and LAT, and "NA" values assigned to INFO's SN and HICOV
    :param vardict_snp_vars_vcf: VarDict's VCF in pileup debug mode
    """
    vardict_snp_vars_fixed_vcf = add_suffix(vardict_snp_vars_vcf, 'fixed')
    info('Fixing VCF, writing to ' + vardict_snp_vars_fixed_vcf)
    with open(vardict_snp_vars_vcf) as inp, open(vardict_snp_vars_fixed_vcf, 'w') as out_f:
        for l in inp:
            if not l.startswith('#'):
                fs = l.split('\t')
                chrom, start, ref, alt = fs[0], fs[1], fs[3], fs[4]
                if ref in ['.', '']:
                    samtools = which('samtools')
                    if not samtools: sys.exit('Error: samtools not in PATH')
                    cmdl = '{samtools} faidx {ref_file} {chrom}:{start}-{start}'.format(**locals())
                    faidx_out = check_output(cmdl, shell=True)
                    fasta_ref = faidx_out.split('\n')[1].strip().upper()
                    assert faidx_out
                    fs[3] = fasta_ref  # '.'  # REF
                    fs[4] = fasta_ref  # '.'  # ALT
                    l = '\t'.join(fs)
                    l = l.replace('=NA;', '=.;')
                    l = l.replace('=;', '=.;')
            out_f.write(l)
            
    assert verify_file(vardict_snp_vars_fixed_vcf) and \
           len(open(vardict_snp_vars_vcf).readlines()) == len(open(vardict_snp_vars_fixed_vcf).readlines()), \
        vardict_snp_vars_fixed_vcf
    os.rename(vardict_snp_vars_fixed_vcf, vardict_snp_vars_vcf)
    return vardict_snp_vars_vcf


def _annotate_vcf(vcf_file, snp_bed):
    gene_by_snp = dict()
    rsid_by_snp = dict()
    for interval in BedTool(snp_bed):
        rsid, gene = interval.name.split('|')
        pos = int(interval.start) + 1
        gene_by_snp[(interval.chrom, pos)] = gene
        rsid_by_snp[(interval.chrom, pos)] = rsid
    
    annotated_vcf = add_suffix(vcf_file, 'ann')
    with open(vcf_file) as f, open(annotated_vcf, 'w') as out:
        vcf_reader = vcf.Reader(f)
        vcf_writer = vcf.Writer(out, vcf_reader)
        vcf_reader.infos['GENE'] = _Info('GENE', '1', 'String', 'Gene name', '', '')
        for rec in vcf_reader:
            rec.INFO['GENE'] = gene_by_snp[(rec.CHROM, rec.POS)]
            rec.ID = rsid_by_snp[(rec.CHROM, rec.POS)]
            vcf_writer.write_record(rec)
    assert verify_file(annotated_vcf) and \
           len(open(vcf_file).readlines()) == len(open(annotated_vcf).readlines()), annotated_vcf
    annotated_vcf = bgzip_and_tabix(annotated_vcf)

    annotated_header_vcf = add_suffix(annotated_vcf, 'hdr')
    header = '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">'
    cmdl = 'bcftools annotate -h <(echo \'{header}\') {annotated_vcf}'.format(**locals())
    run(cmdl, output_fpath=annotated_header_vcf)
    assert verify_file(annotated_header_vcf), annotated_header_vcf
    os.rename(annotated_header_vcf, vcf_file)

    return vcf_file


def __p(rec):
    s = rec.samples[0]
    gt_type = {
        0: 'hom_ref',
        1: 'het',
        2: 'hom_alt',
        None: 'uncalled'
    }.get(s.gt_type)
    gt_bases = s.gt_bases
    return str(rec) + ' FILTER=' + str(rec.FILTER) + ' gt_bases=' + str(gt_bases) + ' gt_type=' + gt_type + ' GT=' + str(s['GT']) + ' Depth=' + str(s['VD']) + '/' + str(s['DP'])


def vcfrec_to_seq(rec, depth_cutoff):
    var = rec.samples[0]

    depth_failed = rec.INFO['DP'] < depth_cutoff
    filter_failed = any(v in ['MSI12', 'InGap'] for v in rec.FILTER)
    if depth_failed or filter_failed:
        var.called = False

    if is_sex_chrom(rec.CHROM):  # We cannot confidentelly determine sex, and thus determine X heterozygocity,
        gt_bases = ''            # so we can't promise constant fingerprint length across all samples
    elif var.called:
        gt_bases = ''.join(sorted(var.gt_bases.split('/')))
    else:
        gt_bases = 'NN'
    
    return gt_bases


def calc_avg_depth(vcf_file):
    with open(vcf_file) as f:
        vcf_reader = vcf.Reader(f)
        recs = [r for r in vcf_reader]
    depths = [r.INFO['DP'] for r in recs]
    return float(sum(depths)) / len(depths)
    

# def check_if_male(recs):
#     y_total_depth = 0
#     for rec in recs:
#         depth = rec.INFO['DP']
#         if 'Y' in rec.CHROM:
#             y_total_depth += depth
#     return y_total_depth >= 5


def vcf_to_fasta(sample, vcf_file, fasta_file, depth_cutoff):
    if can_reuse(fasta_file, vcf_file):
        return fasta_file

    info('Parsing VCF ' + vcf_file)
    with open(vcf_file) as f:
        vcf_reader = vcf.Reader(f)
        recs = [r for r in vcf_reader]

    with open(fasta_file, 'w') as fhw:
        fhw.write('>' + sample.name + '\n')
        fhw.write(''.join(vcfrec_to_seq(rec, depth_cutoff) for rec in recs) + '\n')

    info('Fasta saved to ' + fasta_file)
