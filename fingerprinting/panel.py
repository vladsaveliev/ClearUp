#!/usr/bin/env python

import shutil
import click
from os.path import join, dirname, isfile, isdir, basename

from cyvcf2 import VCF
from pybedtools import BedTool

import az

from ngs_utils import logger as log
from ngs_utils.call_process import run
from ngs_utils.file_utils import safe_mkdir, can_reuse, verify_file, file_transaction, verify_dir, add_suffix, splitext_plus
from ngs_utils.logger import debug

from fingerprinting import get_version
from fingerprinting.genotype import DEPTH_CUTOFF

from ngs_reporting.bcbio.bcbio import BcbioProject

from variant_filtering import reference_data as filt_ref_data, get_anno_config
from variant_filtering.utils import parse_gene_blacklists, check_gene_in_a_blacklist, all_blacklisted_genes


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
    
    # Empty list? Using exome
    if not all_bed_files:
        all_bed_files.add(get_snps_by_type('exome'))
    
    overlapped_bed = join(output_dir, 'merged_panel.bed')
    _overlap_bed_files(all_bed_files, overlapped_bed)
    
    # Selecting SNPs from dbSNP
    dbsnp_file = get_dbsnp(genome_build)
    dbsnp_snps_in_bed = join(output_dir, 'snps.bed')
    cmdl = 'bedtools intersect -header -a ' + dbsnp_file + ' -b ' + overlapped_bed
    run(cmdl, dbsnp_snps_in_bed)
    
    dbsnp_snps_in_bed = _reduce_number_of_locations(dbsnp_snps_in_bed, 200)
    return dbsnp_snps_in_bed


def _reduce_number_of_locations(dbsnp_snps_in_bed, max_number):
    # TODO: split at <max_number> same-size clusters with at least 1 snp in each, and select 1 snp from each
    locs = []
    for i, interval in enumerate(BedTool(dbsnp_snps_in_bed)):
        pos = int(interval.start) + 1
        rsid, gene = interval.name.split('|')
        loc = (interval.chrom, pos, rsid, gene)
        locs.append(loc)
    locs = locs[:max_number]
    out_file = add_suffix(dbsnp_snps_in_bed, '200')
    with file_transaction(None, out_file) as tx:
        with open(tx, 'w') as out:
            for (chrom, pos, rsid,gene) in locs:
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


def _subset_dbsnp(genome_build):
    az.init_sys_cfg()
    genome_cfg = az.get_refdata(genome_build)
    dbsnp_file = verify_file(genome_cfg['dbsnp'], is_critical=True)

    incidentalome_dir = verify_dir(filt_ref_data.incidentalome_dir(), 'incidentalome')
    anno_cfg = get_anno_config()
    blacklisted_genes = all_blacklisted_genes(anno_cfg['blacklist']['genes'], incidentalome_dir)

    subset_dbsnp = join(dirname(__name__), 'snps', 'dbsnp', splitext_plus(dbsnp_file)[0] + 'mafs10pct.bed')
    if verify_file(subset_dbsnp):
        return subset_dbsnp
    total_snps_written = 0
    current_chrom = None
    with file_transaction(None, subset_dbsnp) as tx:
        with open(tx, 'w') as out:
            for v in VCF(dbsnp_file):
                caf = v.INFO.get('CAF')
                if caf and len(v.REF) == 1:
                    gene_val = v.INFO.get('GENEINFO')
                    if gene_val:
                        cafs = [float(a) if a != '.' else 0 for a in caf.split(',')[1:]]
                        alts = v.ALT
                        assert len(cafs) == len(alts)
                        caf_by_snp_alt = dict({a: c for a, c in zip(alts, cafs) if c > 0.1 and len(a) == 1})
                        if len(caf_by_snp_alt) > 0:
                            gene = gene_val.split(':')[0]
                            if gene not in blacklisted_genes:
                                rsid = v.ID
                                fields = [v.CHROM, str(v.POS - 1), str(v.POS), rsid + '|' + gene, v.REF, ','.join(v.ALT)]
                                out.write('\t'.join(fields) + '\n')
                                if current_chrom != v.CHROM:
                                    current_chrom = v.CHROM
                                    debug(v.CHROM + ', written ' + str(total_snps_written))
                                total_snps_written += 1

    # exclude those from LCR and tricky regions:
    ''' bedtools intersect -header -a dbsnp_maf10pct.bed -b ~/bcbio/genomes/Hsapiens/hg19/coverage/problem_regions/GA4GH/*.bed.gz ~/bcbio/genomes/Hsapiens/hg19/coverage/problem_regions/ENCODE/wgEncodeDacMapabilityConsensusExcludable.bed.gz ~/bcbio/genomes/Hsapiens/hg19/coverage/problem_regions/repeats/LCR.bed.gz ~/bcbio/genomes/Hsapiens/hg19/coverage/problem_regions/repeats/sv_repeat_telomere_centromere.bed -v > dbsnp_maf10pct.good.bed
    '''
    # exclude duplicates - remove each second from this output:
    ''' uniq -d dbsnp_maf10pct.good.bed
    '''
    # TODO: check clustered in HapMap?
    
    return subset_dbsnp


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

def get_dbsnp(genome):
    return get_snps_file('dbsnp.bed.gz')


if __name__ == '__main__':
    main()
