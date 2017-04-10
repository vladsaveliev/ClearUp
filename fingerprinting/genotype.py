import shutil
from collections import OrderedDict, defaultdict
import os
from os.path import join, dirname, splitext, basename
import sys
import vcf
from subprocess import check_output

from ngs_utils.bed_utils import bgzip_and_tabix
from ngs_utils.call_process import run
from ngs_utils.utils import is_local, is_us
from ngs_utils.parallel import ParallelCfg, parallel_view
from ngs_utils.file_utils import file_transaction, safe_mkdir, chdir, which, adjust_path, can_reuse, add_suffix, \
    verify_file, intermediate_fname, splitext_plus
from ngs_utils.logger import info, err, critical, debug, warn
from ngs_utils.sambamba import index_bam

import az

from fingerprinting import DEPTH_CUTOFF
from fingerprinting.utils import is_sex_chrom


def genotype(samples, snp_bed, parall_view, work_dir, output_dir, genome_build):
    genome_cfg = az.get_refdata(genome_build)
    info('** Running VarDict ** ')
    vcfs = parall_view.run(_vardict_pileup_sample,
        [[s, work_dir, output_dir, genome_cfg, parall_view.cores_per_job, snp_bed]
         for s in samples])
    vcf_by_sample = OrderedDict(zip([s.name for s in samples], vcfs))
    info('** Finished running VarDict **')
    return vcf_by_sample
    
    
def write_fasta(samples, vcf_by_sample, snp_bed, parall_view, output_dir, work_dir, out_fasta, depth_cutoff=DEPTH_CUTOFF):
    info('Removing sex chromosomes')
    autosomal_vcf_by_sample = OrderedDict()
    for sn, vf in vcf_by_sample.items():
        autosomal_vcf_by_sample[sn] = add_suffix(vf, 'autosomal')
        run('grep -v chrX ' + vf + ' | grep -v chrY', output_fpath=autosomal_vcf_by_sample[sn])
    
    info('Writing fasta')
    fasta_by_sample = OrderedDict((s.name, join(work_dir, s.name + '.fasta')) for s in samples)
    parall_view.run(_vcf_to_fasta, [[s, autosomal_vcf_by_sample[s.name], fasta_by_sample[s.name], depth_cutoff]
                                    for s in samples])
    info('Merging fasta')
    if not can_reuse(out_fasta, fasta_by_sample.values()):
        with open(out_fasta, 'w') as out_f:
            for s in samples:
                with open(fasta_by_sample[s.name]) as f:
                    out_f.write(f.read())
    info('All fasta saved to ' + out_fasta)

    info('Saving VCFs into the output directory')
    for sn, vf in autosomal_vcf_by_sample.items():
        shutil.copy(vf, join(output_dir, basename(autosomal_vcf_by_sample[sn]) + '.gz'))
        shutil.copy(vf, join(output_dir, basename(autosomal_vcf_by_sample[sn]) + '.gz.tbi'))

    # info('Combining VCFs')
    # combined_vcf = combine_vcfs(vcf_by_sample.values(), join(output_dir, 'fingerprints.vcf'))
    # assert combined_vcf

    # ped_file = splitext_plus(combined_vcf)[0]
    # cmdl = 'vcftools --plink --gzvcf {combined_vcf}'.format(**locals())
    # cmdl += ' --out ' + splitext(ped_file)[0]
    # run(cmd=cmdl)
    # info('Converting VCFs to PED')
    # ped_file = vcf_to_ped(vcf_by_sample, join(output_dir, 'fingerprints.ped'), sex_by_sample, depth_cutoff)
    
    return out_fasta, autosomal_vcf_by_sample


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


def genotype_bcbio_proj(proj, snp_bed, parallel_cfg, depth_cutoff=DEPTH_CUTOFF,
                        output_dir=None, work_dir=None):
    output_dir = output_dir or safe_mkdir(join(proj.date_dir, 'fingerprints'))
    work_dir = work_dir or safe_mkdir(proj.work_dir)
    with parallel_view(len(proj.samples), parallel_cfg, work_dir) as parall_view:
        vcf_by_sample = genotype(proj.samples, snp_bed, parall_view, work_dir, output_dir, proj.genome_build)
        fasta_file, vcf_by_sample = write_fasta(proj.samples, vcf_by_sample,
                                                snp_bed, parall_view, output_dir, work_dir,
                                                out_fasta=join(output_dir, 'fingerprints.fasta'), depth_cutoff=depth_cutoff)
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
    run(cmdl, output_fpath=vardict_snp_vars)

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
            ' | ' + join(vardict_dir, 'var2vcf_valid.pl') +
            # ' | grep "^#\|TYPE=SNV\|TYPE=REF" ' +
            '')
    run(cmdl, output_fpath=work_vcf_file)
    
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
                out_f.write(l)
    assert verify_file(fixed_vcf_file)

    info('Annotating VCF with gene names and rsIDs')
    ann_vcf_file = add_suffix(fixed_vcf_file, 'ann')
    with open(fixed_vcf_file) as f, open(ann_vcf_file, 'w') as out:
        vcf_reader = vcf.Reader(f)
        vcf_writer = vcf.Writer(out, vcf_reader)
        for rec in vcf_reader:
            # if any(len(a) > 0 for a in rec.ALT) or len(ref) > 0:
            #     rec.REF = rec.ALT = get_fasta_ref(ref_file, chrom, pos)  # Reading the reference allele from fasta
            #     rsid = rsid_by_var[(chrom, pos, ref, alt)]
                # fs[1] = snp_by_rsid[rsid]  # Resetting the POS for complex variants to the ogirinal SNP location
                # fs[6] = 'complex' if fs[6]  # Filter
                # debug('Complex variant ' + chrom + ':' + str(pos) + ' ' + ref + '>' + alt + '; resetting POS to ' + fs[1])

            # if len(ref) > 1:   # Deletion or complex variant
            #     if loc.pos >= pos:  # SNP is within the variant
            #         fs[3] = ref[loc.pos - pos]
            #     else:               # SNP is outside the variant
            #         fs[3] = fasta_ref
            # if len(alt) > 1:   # Insertion or complex variant
            #     if loc.pos >= pos:
            #         fs[4] = alt[loc.pos - pos]
            #     else:
            #         fs[4] = fasta_ref
            # if len(ref) > 1 or len(alt) > 1:
            #     fs[3] = fs[4] = ''
            ann = ann_by_var[(rec.CHROM, str(rec.POS))]
            rsid, gene = ann.split('|')
            rec.INFO['GENE'] = gene
            rec.ID = rsid
            vcf_writer.write_record(rec)
    assert verify_file(ann_vcf_file), ann_vcf_file

    info('Selecting unique records per rsid (prioritizing SNV)')
    unq_vcf_file = add_suffix(ann_vcf_file, 'unq')
    with open(ann_vcf_file) as f, open(unq_vcf_file, 'w') as out:
        vcf_reader = vcf.Reader(f)
        vcf_writer = vcf.Writer(out, vcf_reader)
        recs = [r for r in vcf_reader]
        recs_by_rsid = defaultdict(list)
        for r in recs:
            recs_by_rsid[r.ID].append(r)
        for r in recs:
            rs = recs_by_rsid[r.ID]
            if not rs:
                continue
            del recs_by_rsid[r.ID]
            if len(rs) > 1:
                debug('Multiple records found for id ' + r.ID)
            if len(rs) == 1:
                vcf_writer.write_record(rs[0])
            else:
                snp_rs = [r for r in rs if r.INFO['TYPE'] == 'SNV']
                if len(snp_rs) > 0:
                    vcf_writer.write_record(snp_rs[0])
                else:
                    vcf_writer.write_record(rs[0])
    assert verify_file(unq_vcf_file), unq_vcf_file

    ann_hdr_vcf_file = add_suffix(unq_vcf_file, 'hdr')
    cmdl = 'bcftools annotate -h <(echo \'##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">\') ' + bgzip_and_tabix(unq_vcf_file)
    run(cmdl, output_fpath=ann_hdr_vcf_file)
    
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
    

# def _annotate_vcf(vcf_file, snp_bed):
#     gene_by_snp = dict()
#     rsid_by_snp = dict()
#     for interval in BedTool(snp_bed):
#         rsid, gene = interval.name.split('|')
#         pos = int(interval.start) + 1
#         gene_by_snp[(interval.chrom, pos)] = gene
#         rsid_by_snp[(interval.chrom, pos)] = rsid
#
#     ann_vcf_file = add_suffix(vcf_file, 'ann')
#     with open(vcf_file) as f, open(ann_vcf_file, 'w') as out:
#         vcf_reader = vcf.Reader(f)
#         vcf_writer = vcf.Writer(out, vcf_reader)
#         for rec in vcf_reader:
#             if (rec.CHROM, rec.POS) in gene_by_snp:
#                 rec.INFO['GENE'] = gene_by_snp[(rec.CHROM, rec.POS)]
#                 rec.ID = rsid_by_snp[(rec.CHROM, rec.POS)]
#                 vcf_writer.write_record(rec)
#     ann_vcf_file = bgzip_and_tabix(ann_vcf_file)
#     vcf_file += '.gz'
#     # debug('Tabixed, annotating into ' + ann_vcf_file)
#     # from cyvcf2 import VCF, Writer
#     # vcf = VCF(vcf_file)
#     # debug('Created reader')
#     # vcf.add_info_to_header({'ID': 'GENE', 'Description': 'Overlapping gene', 'Type': 'String', 'Number': '1'})
#     # vcf.add_info_to_header({'ID': 'rsid', 'Description': 'dbSNP rsID', 'Type': 'String', 'Number': '1'})
#     # debug('Added header to reader')
#     # w = Writer(ann_vcf_file, vcf)
#     # debug('Creted writer')
#     # for rec in vcf:
#     #     debug('Annotating rec ' + str(rec))
#     #     if (rec.CHROM, rec.POS) in gene_by_snp:
#     #         debug('Rec is in gene_by_snp')
#     #         rec.INFO['GENE'] = gene_by_snp[(rec.CHROM, rec.POS)]
#     #         rec.INFO['rsid'] = rsid_by_snp[(rec.CHROM, rec.POS)]
#     #         w.write_record(rec)
#     #         debug('Written rec')
#     # debug('Annotated, saved into ' + ann_vcf_file)
#     # w.close()
#     # debug('Closed ' + ann_vcf_file)
#     assert verify_file(ann_vcf_file), ann_vcf_file
#
#     ann_hdr_vcf_file = add_suffix(ann_vcf_file, 'hdr')
#     header = '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">'
#     cmdl = 'bcftools annotate -h <(echo \'{header}\') {ann_vcf_file}'.format(**locals())
#     run(cmdl, output_fpath=ann_hdr_vcf_file)
#
#     debug('Renaming ' + ann_hdr_vcf_file + ' -> ' + vcf_file)
#     os.rename(ann_hdr_vcf_file, vcf_file)
#     return vcf_file


# def __p(rec):
#     s = rec.samples[0]
#     gt_type = {
#         0: 'hom_ref',
#         1: 'het',
#         2: 'hom_alt',
#         None: 'uncalled'
#     }.get(s.gt_type)
#     gt_bases = s.gt_bases
#     return str(rec) + ' FILTER=' + str(rec.FILTER) + ' gt_bases=' + str(gt_bases) + ' gt_type=' + gt_type + ' GT=' + str(s['GT']) + ' Depth=' + str(s['VD']) + '/' + str(s['DP'])


def vcfrec_to_seq(rec, depth_cutoff):
    called = rec.num_called > 0
    depth_failed = rec.INFO['DP'] < depth_cutoff
    filter_failed = rec.FILTER and any(v in ['MSI12', 'InGap'] for v in rec.FILTER)
    non_snp = len(rec.REF) > 1 or any(len(a) > 1 for a in rec.ALT)
    if depth_failed or filter_failed or non_snp:
        called = False

    if is_sex_chrom(rec.CHROM):  # We cannot confidentelly determine sex, and thus determine X heterozygocity,
        gt_bases = 'NN'          # so we can't promise constant fingerprint length across all samples
    elif called:
        gt_bases = ''.join(sorted(rec.gt_bases[0].split('/')))
    else:
        gt_bases = 'NN'
    
    return gt_bases


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


def _vcf_to_fasta(sample, vcf_file, fasta_file, depth_cutoff):
    if can_reuse(fasta_file, vcf_file):
        return fasta_file

    from cyvcf2 import VCF
    info('Parsing VCF ' + vcf_file)
    recs = [v for v in VCF(vcf_file)]
    with open(fasta_file, 'w') as fhw:
        fhw.write('>' + sample.name + '\n')
        fhw.write(''.join(vcfrec_to_seq(rec, depth_cutoff) for rec in recs) + '\n')

    info('Fasta saved to ' + fasta_file)
