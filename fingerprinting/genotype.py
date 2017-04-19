import shutil
from collections import OrderedDict, defaultdict
import os
from os.path import join, dirname, splitext, basename
import sys
import vcf
from subprocess import check_output

from ngs_utils.bed_utils import bgzip_and_tabix
from ngs_utils import call_process
from ngs_utils.utils import is_local, is_us
from ngs_utils.parallel import ParallelCfg, parallel_view
from ngs_utils.file_utils import file_transaction, safe_mkdir, chdir, which, adjust_path, can_reuse, add_suffix, \
    verify_file, intermediate_fname, splitext_plus
from ngs_utils.logger import info, err, critical, debug, warn
from ngs_utils.sambamba import index_bam

import az

from fingerprinting.utils import is_sex_chrom
from fingerprinting import DEPTH_CUTOFF, AF_CUTOFF


class Allele:
    def __init__(self, rec=None, nt=None, af=None, depth=None):
        self.nt = nt
        self.af = af
        self.depth = depth
        self.passing = True
        if rec:
            self.af = rec.INFO['AF']
            self.depth = rec.INFO['VD']
            if rec.INFO['TYPE'] == 'REF':
                self.af = 1.0
                self.depth = rec.INFO['DP']
                self.nt = rec.ALT[0]
            else:
                is_complex = len(rec.REF) > 1 or any(len(a) > 1 for a in rec.ALT)
                self.passing = rec.num_called > 0 and not rec.FILTER and not is_complex
                if self.passing:
                    self.nt = rec.ALT[0]
                else:
                    self.nt = 'N'


def build_snp_from_records(snp, records, min_depth):
    # TODO:
    # keep alt in the dbSNP file, and select only the variant here with matching alt. Then report AF of that
    if not records:
        snp.depth = snp.allele1_depth = snp.allele2_depth = 0
    else:
        snp.depth = records[0].INFO['DP']

    alleles = [Allele(rec) for rec in records]
    sum_af = sum(c.af for c in alleles)
    if sum_af < 1:
        ref_af = 1.0 - sum_af
        alleles.append(Allele(nt=snp.location.ref, af=ref_af,
                              depth=ref_af * snp.depth))  # adding the reference allele
    alleles = [a for a in alleles if a.af >= AF_CUTOFF]
    alleles.sort(key=lambda a_: a_.af)
    
    if len(alleles) == 1:  # only 1 allele with AF>AF_CUTOFF - homozygous
        snp.allele1 = snp.allele2 = alleles[-1].nt
        snp.allele1_depth = snp.allele2_depth = alleles[-1].depth
    else:  # multiple alleles with AF>AF_CUTOFF - heterozygous
        snp.allele1 = alleles[-1].nt
        snp.allele2 = alleles[-2].nt
        snp.allele1_depth = alleles[-1].depth
        snp.allele2_depth = alleles[-2].depth

    if snp.depth < min_depth:  # Not enough depth on location to call variation
        snp.allele1, snp.allele2 = 'N', 'N'
        
        # for i, rec in enumerate(high_af_calls):
        #     called = rec.num_called > 0
        #     filter_failed = rec.FILTER
        #     is_complex = len(rec.REF) > 1 or any(len(a) > 1 for a in rec.ALT)
        #     called = called and not filter_failed and not is_complex
        #     if called:
        #         alleles[i] = rec.ALT[0]
        #     else:
        #         alleles[i] = 'N'
    return snp


def build_tree(run):
    info('Writing fasta to ' + run.fasta_file_path())
    samples = [s for p in run.projects for s in p.samples]
    with open(run.fasta_file_path(), 'w') as fhw:
        for s in samples:
            snps_by_rsid = s.snps_from_run(run)
            fhw.write('>' + s.long_name() + '\n')
            fhw.write(''.join(snps_by_rsid[loc.rsid].get_gt() for loc in run.locations.all()) + '\n')
    info('All fasta saved to ' + run.fasta_file_path())

    info()
    info('Building phylogeny tree using prank...')
    prank_out = join(run.work_dir_path(), splitext(basename(run.fasta_file_path()))[0])
    suffix = 'lnx' if 'linux' in sys.platform else 'osx'
    prank_bin = join(dirname(__file__), 'prank', 'prank_' + suffix, 'bin', 'prank')
    call_process.run(prank_bin + ' -d=' + run.fasta_file_path() + ' -o=' + prank_out + ' -showtree')
    if not verify_file(prank_out + '.best.dnd'):
        critical('Prank failed to run')
    os.rename(prank_out + '.best.dnd', run.tree_file_path())
    os.remove(prank_out + '.best.fas')
    
    return run.fasta_file_path()


def genotype(samples, snp_bed, parall_view, work_dir, output_dir, genome_build):
    genome_cfg = az.get_refdata(genome_build)
    info('** Running VarDict ** ')
    vcfs = parall_view.run(_vardict_pileup_sample,
        [[s, work_dir, output_dir, genome_cfg, parall_view.cores_per_job, snp_bed]
         for s in samples])
    vcf_by_sample = OrderedDict(zip([s.name for s in samples], vcfs))
    info('** Finished running VarDict **')
    return vcf_by_sample
    

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


def _vardict_pileup_sample(sample, work_dir, output_dir, genome_cfg, threads, snp_file):
    vardict_snp_vars = join(work_dir, sample.name + '_vars.txt')
    vcf_file = join(output_dir, sample.name + '.vcf')
    if can_reuse(vardict_snp_vars, [sample.bam, snp_file]) and can_reuse(vcf_file, vardict_snp_vars):
        return vcf_file
    
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
    call_process.run(cmdl, output_fpath=vardict_snp_vars)

    # Complex variants might have a shifted start positions with respect to rsid so we are
    # associating starts with rsid for futher snp identification
    ann_by_var = defaultdict(list)
    with open(vardict_snp_vars) as f:
        for l in f:
            fs = l.split('\t')
            ann, chrom, start = fs[1], fs[2], fs[3]
            ann_by_var[(chrom, start)] = ann
    
    info()
    info('Converting to VCF')
    work_vcf_file = join(work_dir, sample.name + '_vars.vcf')
    cmdl = ('cut -f-34 ' + vardict_snp_vars +
            ' | awk -F"\\t" -v OFS="\\t" \'{for (i=1;i<=NF;i++) { if ($i=="") $i="0" } print $0 }\''
            ' | ' + join(vardict_dir, 'teststrandbias.R') +
            ' | ' + join(vardict_dir, 'var2vcf_valid.pl') + ' -A -f 0.2' +
            '')
    call_process.run(cmdl, output_fpath=work_vcf_file)
    
    # Fix non-call records with empty REF and LAT, and "NA" values assigned to INFO's SN and HICOV
    fixed_vcf_file = add_suffix(work_vcf_file, 'fixed')
    info('Fixing VCF for parsing, writing to ' + fixed_vcf_file)
    with open(work_vcf_file) as inp, open(fixed_vcf_file, 'w') as out_f:
        for l in inp:
            if l.startswith('#'):
                out_f.write(l)
            else:
                fs = l.split('\t')
                chrom, pos, _, ref, alt = fs[0], int(fs[1]), fs[2], fs[3], fs[4]
                if alt in ['.', '']:
                    fs[4] = fs[3] = _get_fasta_ref(ref_file, chrom, pos)  # Reading the reference allele from fasta
                l = '\t'.join(fs)
                l = l.replace('=NA;', '=.;')
                l = l.replace('=;', '=.;')
                l = l.replace('TYPE=0', 'TYPE=REF')
                out_f.write(l)
    assert verify_file(fixed_vcf_file)

    info('Annotating VCF with gene names and rsIDs')
    ann_vcf_file = add_suffix(fixed_vcf_file, 'ann')
    with open(fixed_vcf_file) as f, open(ann_vcf_file, 'w') as out:
        vcf_reader = vcf.Reader(f)
        vcf_writer = vcf.Writer(out, vcf_reader)
        for rec in vcf_reader:
            ann = ann_by_var[(rec.CHROM, str(rec.POS))]
            rsid, gene, ref = ann.split('|')
            rec.INFO['GENE'] = gene
            rec.ID = rsid
            vcf_writer.write_record(rec)
    assert verify_file(ann_vcf_file), ann_vcf_file

    ann_hdr_vcf_file = add_suffix(ann_vcf_file, 'hdr')
    cmdl = 'bcftools annotate -h <(echo ' \
           '\'##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">\') ' + \
           bgzip_and_tabix(ann_vcf_file)
    call_process.run(cmdl, output_fpath=ann_hdr_vcf_file)
    
    debug('Renaming ' + ann_hdr_vcf_file + ' -> ' + vcf_file)
    os.rename(ann_hdr_vcf_file, vcf_file)
    return vcf_file


def _get_fasta_ref(ref_file, chrom, pos):
    samtools = which('samtools')
    if not samtools: sys.exit('Error: samtools not in PATH')
    cmdl = '{samtools} faidx {ref_file} {chrom}:{pos}-{pos}'.format(**locals())
    faidx_out = check_output(cmdl, shell=True)
    fasta_ref = faidx_out.split('\n')[1].strip().upper()
    assert faidx_out
    return fasta_ref
    

# def vcf_to_ped(vcf_by_sample, ped_file, sex_by_sample, depth_cutoff):
#     if can_reuse(ped_file, vcf_by_sample.values()):
#         return ped_file
#
#     from cyvcf2 import VCF
#     with file_transaction(None, ped_file) as tx:
#         with open(tx, 'w') as out_f:
#             for sname, vcf_file in vcf_by_sample.items():
#                 recs = [v for v in VCF(vcf_file)]
#                 sex = {'M': '1', 'F': '2'}.get(sex_by_sample[sname], '0')
#                 out_fields = [sname, sname, '0', '0', sex, '0']
#                 for rec in recs:
#                     gt = vcfrec_to_seq(rec, depth_cutoff).replace('N', '0')
#                     out_fields.extend([gt[0], gt[1]])
#                 out_f.write('\t'.join(out_fields) + '\n')
#
#     info('PED saved to ' + ped_file)
#     return ped_file


# def genotype_bcbio_proj(proj, snp_bed, parallel_cfg, depth_cutoff=DEPTH_CUTOFF,
#                         output_dir=None, work_dir=None):
#     output_dir = output_dir or safe_mkdir(join(proj.date_dir, 'fingerprints'))
#     work_dir = work_dir or safe_mkdir(proj.work_dir)
#     with parallel_view(len(proj.samples), parallel_cfg, work_dir) as parall_view:
#         vcf_by_sample = genotype(proj.samples, snp_bed, parall_view, work_dir, output_dir, proj.genome_build)
#         fasta_file, vcf_by_sample = write_fasta(proj.samples, vcf_by_sample,
#                                                 snp_bed, parall_view, output_dir, work_dir,
#                                                 out_fasta=join(output_dir, 'fingerprints.fasta'), depth_cutoff=depth_cutoff)
#     return vcf_by_sample


# def calc_avg_depth(vcf_file):
#     with open(vcf_file) as f:
#         vcf_reader = vcf.Reader(f)
#         recs = [r for r in vcf_reader]
#     depths = [r.INFO['DP'] for r in recs]
#     return float(sum(depths)) / len(depths)


# def check_if_male(recs):
#     y_total_depth = 0
#     for rec in recs:
#         depth = rec.INFO['DP']
#         if 'Y' in rec.CHROM:
#             y_total_depth += depth
#     return y_total_depth >= 5
