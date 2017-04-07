#!/usr/bin/env python
import random
import shutil
from collections import OrderedDict, defaultdict

import click
from os.path import join, dirname, isfile, isdir, basename

from fingerprinting.utils import is_sex_chrom
from pybedtools import BedTool

from ngs_utils import logger as log
from ngs_utils.call_process import run
from ngs_utils.file_utils import safe_mkdir, can_reuse, verify_file, file_transaction, verify_dir, add_suffix, splitext_plus
from ngs_utils.logger import debug
from ngs_utils.reference_data import get_chrom_order

from fingerprinting import get_version
from fingerprinting.genotype import DEPTH_CUTOFF

from ngs_reporting.bcbio.bcbio import BcbioProject



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
    # TODO:
    #   if a BED is found:
    #       if only exomes or wgs: select get_snps_by_type('exome')
    #       if panel present (check if is_small_bed()): select from the dbSNP_v149.mafs10pct.vcf.gz.tx. If too many, select less clustered

    all_bed_files = set()
    for proj in bcbio_projs or []:
        if proj.coverage_bed:
            log.info(proj.project_name + ': selecting ' + proj.coverage_bed)
            all_bed_files.add(proj.coverage_bed)
    all_bed_files |= set(bed_files or [])
    
    if not all_bed_files:  # Empty list? Using exome
        all_bed_files.add(get_snps_by_type('exome'))
    
    overlapped_bed = join(output_dir, 'merged_panel.bed')
    _overlap_bed_files(all_bed_files, overlapped_bed)
    
    # Selecting SNPs from dbSNP
    dbsnp_file = get_dbsnp(genome_build)
    dbsnp_snps_in_bed = join(output_dir, 'snps.bed')
    if not can_reuse(dbsnp_snps_in_bed, [dbsnp_file, overlapped_bed]):
        cmdl = 'bedtools intersect -header -a ' + dbsnp_file + ' -b ' + overlapped_bed
        run(cmdl, dbsnp_snps_in_bed)
    
    return _reduce_number_of_locations(dbsnp_snps_in_bed, genome_build)


def _reduce_number_of_locations(dbsnp_snps_in_bed, genome_build, autosomal_locations_limit=200):
    out_file = add_suffix(dbsnp_snps_in_bed, str(autosomal_locations_limit))
    if can_reuse(out_file, dbsnp_snps_in_bed):
        return out_file

    # TODO: split at <max_number> same-size clusters with at least 1 snp in each, and select 1 snp from each
    locs_by_gene = defaultdict(list)
    total_locs = 0
    for i, interval in enumerate(BedTool(dbsnp_snps_in_bed)):
        if is_sex_chrom(interval.chrom):
            continue
        pos = int(interval.start) + 1
        rsid, gene = interval.name.split('|')
        loc = (interval.chrom, pos, rsid, gene)
        locs_by_gene[gene].append(loc)
        total_locs += 1
    
    gnames = random.sample(locs_by_gene.keys(), min(len(locs_by_gene), autosomal_locations_limit))
    min_locs_per_gene = min(len(locs) for locs in locs_by_gene.values())
    locs_per_gene = min(autosomal_locations_limit / len(gnames), min_locs_per_gene)
    selected_locs_by_gene = {g: random.sample(locs, locs_per_gene) for g, locs in locs_by_gene.items()}
    selected_locs = [l for locs in selected_locs_by_gene.values() for l in locs]
    selected_locs.sort(key=lambda a: (get_chrom_order(genome_build).get(a[0], -1), a[1:]))

    # non_clustered_locs = []
    # prev_pos = 0
    # for (chrom, pos, rsid, gene) in selected_locs:
    #     if 0 < pos - prev_pos < 200:
    #         continue
    #     else:
    #         prev_pos = pos
    #         non_clustered_locs.append((chrom, pos, rsid, gene))
    
    log.debug('Autosomal locations:')
    for l in selected_locs:
        log.debug('  ' + str(l))
    
    with file_transaction(None, out_file) as tx:
        with open(tx, 'w') as out:
            for (chrom, pos, rsid, gene) in selected_locs:
                out.write('\t'.join([chrom, str(pos-1), str(pos), rsid + '|' + gene]) + '\n')
    return out_file


def _overlap_bed_files(bed_files, output_bed_file):
    if can_reuse(output_bed_file, bed_files):
        return output_bed_file
    if len(bed_files) == 1:
        shutil.copy(bed_files.pop(), output_bed_file)
        return output_bed_file
    cmdl = 'bedops --intersect' + ''.join([' <(sort-bed ' + bf + ')' for bf in bed_files])
    run(cmdl, output_bed_file)
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

BED_BY_TYPE = {
    'idt': 'idt_snps.bed',
    'exome': 'exome_snps.bed',
}
def get_snps_by_type(panel_type):
    snps_fname = BED_BY_TYPE.get(panel_type)
    if not panel_type:
        log.critical('SNPs for panel type ' + panel_type + ' is not defined')
    return get_snps_file(snps_fname)

def get_snps_file(fname):
    return verify_file(join(dirname(__file__), 'snps', fname), is_critical=True)

def get_dbsnp(genome, take_autosomal=True, take_sex=False):
    if take_autosomal:
        return get_snps_file('dbsnp.autosomal.bed.gz')
    elif take_sex:
        return get_snps_file('dbsnp.chrX.bed.gz')


if __name__ == '__main__':
    main()
