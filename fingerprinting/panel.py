#!/usr/bin/env python
import random
import shutil
from collections import OrderedDict, defaultdict
import click
from os.path import join, dirname, isfile, isdir, basename

import math
from pybedtools import BedTool

from ngs_utils import logger as log
from ngs_utils import call_process
from ngs_utils.file_utils import safe_mkdir, can_reuse, verify_file, file_transaction, verify_dir, add_suffix, splitext_plus
from ngs_utils.logger import debug
from ngs_utils.reference_data import get_chrom_order
from ngs_reporting.bcbio.bcbio import BcbioProject

from fingerprinting.utils import is_sex_chrom
from fingerprinting import get_version, DEPTH_CUTOFF


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('paths',
                type=click.Path(exists=True),
                metavar="<bcbio directories or bed files>",
                nargs=-1)
@click.option('-o', '--output-dir',
              type=click.Path(dir_okay=True),
              metavar='<output bed path>',
              required=True)
@click.option('-g', '--genome',
              type=click.STRING,
              metavar='<genome>',
              required=True)
@click.option('-D', '--depth',
              type=int,
              help='Minimum coverage depth for calls.',
              default=DEPTH_CUTOFF,
              )
@click.version_option(version=get_version())
def main(paths, output_dir, genome, depth):
    bed_files = [verify_file(f, is_critical=True) for f in paths if isfile(f)]
    bcbio_projs = []
    for d in [verify_dir(f, is_critical=True) for f in paths if isdir(f)]:
        proj = BcbioProject()
        proj.load_from_bcbio_dir(d, proc_name='fingerprinting', need_coverage_interval=False)
        bcbio_projs.append(proj)
    
    log.init(True)
    build_snps_panel(bcbio_projs, bed_files, safe_mkdir(output_dir), genome)


def build_snps_panel(bcbio_projs=None, bed_files=None, output_dir=None, genome_build=None):
    selected_snps_file = join(output_dir, 'snps.bed')
    if can_reuse(selected_snps_file, bed_files):
        return selected_snps_file
    
    work_dir = safe_mkdir(join(output_dir, 'work'))

    all_bed_files = set()
    for proj in bcbio_projs or []:
        if proj.coverage_bed:
            log.info(proj.project_name + ': selecting ' + proj.coverage_bed)
            all_bed_files.add(proj.coverage_bed)
    all_bed_files |= set(bed_files or [])
    
    overlapped_bed = join(work_dir, 'merged_bed_files.bed')
    _overlap_bed_files(all_bed_files, overlapped_bed)
    
    # Selecting SNPs from dbSNP
    dbsnp_file = get_dbsnp(genome_build)
    dbsnp_snps_file = join(work_dir, 'snps_in_merged_bed_files.bed')
    if not can_reuse(dbsnp_snps_file, [dbsnp_file, overlapped_bed]):
        cmdl = 'bedtools intersect -header -a ' + dbsnp_file + ' -b ' + overlapped_bed
        call_process.run(cmdl, dbsnp_snps_file)

    subset_bed_file = add_suffix(dbsnp_snps_file, 'subset')
    _make_snp_file(dbsnp_snps_file, genome_build, subset_bed_file)
    
    shutil.copyfile(subset_bed_file, selected_snps_file)
    return selected_snps_file


def _make_snp_file(dbsnp_snps_file, genome_build, output_file,
                   autosomal_locations_limit=175, min_snp_amount=30):
    if can_reuse(output_file, dbsnp_snps_file):
        return output_file

    locs_by_gene = defaultdict(list)
    total_locs = 0
    for i, interval in enumerate(BedTool(dbsnp_snps_file)):
        if is_sex_chrom(interval.chrom):
            continue
        pos = int(interval.start) + 1
        annots = interval.name.split('|')
        if len(annots) == 2:
            rsid, gene = interval.name.split('|')
            ref = interval[4]
        else:
            rsid, gene, ref = interval.name.split('|')
        loc = (interval.chrom, pos, rsid, gene, ref)
        locs_by_gene[gene].append(loc)
        total_locs += 1
    
    random.seed(1234)  # seeding random for reproducability

    # Selecting random genes
    gnames = random.sample(locs_by_gene.keys(), min(len(locs_by_gene), autosomal_locations_limit))
    locs_by_gene = {g: locs_by_gene[g] for g in gnames}
    # Selecting random SNPs in each gene
    # min_locs_per_gene = min(len(locs) for locs in locs_by_gene.values())
    # if pick_unclustered:
    #     locs_per_gene = min(autosomal_locations_limit / len(gnames), min_locs_per_gene)
    #     while locs_per_gene * len(gnames) < min_snp_amount:
    #         locs_per_gene = math.ceil(float(min_snp_amount) / len(gnames))
    #     selected_locs_by_gene = {g: random.sample(locs_by_gene[g], locs_per_gene) for g in gnames}
    #     selected_locs = [l for gene_locs in selected_locs_by_gene.values() for l in gene_locs]
    # else:
    all_locs = [l for gene_locs in locs_by_gene.values() for l in gene_locs]

    # Selecting unclustered SNPs within genes
    non_clustered_locs = []
    prev_pos = 0
    for (chrom, pos, rsid, gene, ref) in all_locs:
        if 0 < pos - prev_pos < 500:
            continue
        else:
            prev_pos = pos
            non_clustered_locs.append((chrom, pos, rsid, gene, ref))

    # Selecting random SNPs within the limit
    selected_locs = random.sample(non_clustered_locs, min(len(non_clustered_locs), autosomal_locations_limit))

    # Sorting final locations
    chrom_order = get_chrom_order(genome_build)
    selected_locs.sort(key=lambda a: (chrom_order.get(a[0], -1), a[1:]))

    log.debug('Selected the following autosomal SNPs:')
    for (chrom, pos, rsid, gene, ref) in selected_locs:
        log.debug('  ' + chrom + ':' + str(pos) + '\t' + rsid + '\t' + gene)
    
    with file_transaction(None, output_file) as tx:
        with open(tx, 'w') as out:
            for (chrom, pos, rsid, gene, ref) in selected_locs:
                out.write('\t'.join([chrom, str(pos-1), str(pos), rsid + '|' + gene + '|' + ref]) + '\n')
    return output_file


def _overlap_bed_files(bed_files, output_bed_file):
    if can_reuse(output_bed_file, bed_files):
        return output_bed_file
    if len(bed_files) == 1:
        shutil.copy(bed_files.pop(), output_bed_file)
        return output_bed_file
    cmdl = 'bedops --intersect' + ''.join([' <(sort-bed ' + bf + ')' for bf in bed_files])
    call_process.run(cmdl, output_bed_file)
    return output_bed_file


# def _filter_snps(snp_vcf):
#     # filt_snp_vcf = add_suffix(snp_vcf, 'maf10pct')
#     snps = []
#     for v in VCF(snp_vcf):
#         caf = v.INFO.get('CAF')
#         if caf:
#             cafs = [float(a) if a != '.' else 0 for a in caf.split(',')[1:]]
#             alts = v.ALT
#             assert len(cafs) == len(alts)
#             caf_by_alt = dict({a: c for a, c in zip(alts, cafs) if c > 0.1})
#             assert len(caf_by_alt) <= 1, str(v)
#             if len(caf_by_alt) == 1:
#                 alt, caf = caf_by_alt.items()[0]
#                 log.info(str(alt) + ' ' + str(caf))
#                 snps.append(v)
#
#     log.info('Found SNPs with MAF>10%: ' + str(len(snps)))
#     return snps
#


def get_snps_file(fname):
    return verify_file(join(dirname(__file__), 'snps', fname), is_critical=True)


def get_dbsnp(genome):
    # return get_snps_file('dbsnp.autosomal.bed.gz')
    return get_snps_file('dbsnp_maf10pct.no_selfchain_gc25-30_65-70_lowcomp50.autosomal.bed.gz')


if __name__ == '__main__':
    main()
