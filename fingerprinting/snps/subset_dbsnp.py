#!/usr/bin/env python
import sys

from ngs_utils.utils import is_chihua
from os.path import join, dirname, isfile, isdir, basename
from cyvcf2 import VCF

from ngs_utils.logger import critical, info, warn
from ngs_utils.call_process import run
from ngs_utils.file_utils import safe_mkdir, can_reuse, verify_file, file_transaction, verify_dir, add_suffix, splitext_plus

from variant_filtering import reference_data as filt_ref_data, get_anno_config
from variant_filtering.utils import parse_gene_blacklists, check_gene_in_a_blacklist, all_blacklisted_genes


MIN_MAF = 0.1  # 10%


if is_chihua():
    dbsnp_loc_by_genome = {
        'hg19': '/home/vsaveliev/bcbio/genomes/Hsapiens/hg19/variation/dbSNP_v150.vcf.gz',
        'hg38': '/home/vsaveliev/bcbio/genomes/Hsapiens/hg38/variation/dbSNP_v150.vcf.gz',
    }
else:
    dbsnp_loc_by_genome = {
        'hg19': '/Users/vlad/googledrive/az/reference_data/genomes/Hsapiens/hg19/variation/dbSNP_v149.vcf.gz',
        'hg38': '/Users/vlad/googledrive/az/reference_data/genomes/Hsapiens/hg38/variation/dbSNP_v150.vcf.gz',
    }


def subset_dbsnp(genome):
    dbsnp_location = dbsnp_loc_by_genome.get(genome)
    if not dbsnp_location:
        critical('Genome ' + genome + ' is not supported')
    dbsnp_file = verify_file(dbsnp_location, is_critical=True)

    incidentalome_dir = verify_dir(filt_ref_data.incidentalome_dir(), 'incidentalome')
    anno_cfg = get_anno_config()
    blacklisted_genes = all_blacklisted_genes(anno_cfg['blacklist']['genes'], incidentalome_dir)

    total_snps_written = 0
    current_chrom = None
    output_path = join(dirname(__file__), 'dbsnp.maf10pct.' + genome + '.1.bed')
    info('Writing to ' + output_path)
    prev_pos = 0
    with file_transaction(None, output_path) as tx:
        with open(tx, 'w') as out:
            for v in VCF(dbsnp_file):
                caf = v.INFO.get('CAF')
                if not caf or not len(v.REF) == 1:
                    continue
                gene_val = v.INFO.get('GENEINFO')
                if not gene_val:
                    continue
                if 0 < v.POS - prev_pos < 1:  # Checking if clustered (1 means basically turn off, i.e. removing only exact duplicates)
                    continue
                prev_pos = v.POS
                cafs = [float(a) if a != '.' else 0 for a in caf.split(',')[1:]]
                alts = v.ALT
                assert len(cafs) == len(alts)
                caf_by_snp_alt = dict({a: c for a, c in zip(alts, cafs) if c > MIN_MAF and len(a) == 1})
                if len(caf_by_snp_alt) == 0:
                    continue
                gene = gene_val.split(':')[0]
                if gene in blacklisted_genes:
                    continue
                chrom = v.CHROM if v.CHROM.startswith('chr') else 'chr' + v.CHROM.replace('MT', 'M')
                # for alt, caf in caf_by_snp_alt.items():
                fields = [chrom, str(v.POS - 1), str(v.POS), v.ID + '|' + gene + '|' + v.REF + '|' + ','.join(caf_by_snp_alt.keys())]
                out.write('\t'.join(fields) + '\n')
                if current_chrom != chrom:
                    current_chrom = chrom
                    info(chrom + ', written ' + str(total_snps_written))
                total_snps_written += 1
    info('Done, saved to ' + output_path)

    # exclude those from LCR and tricky regions:
    ''' scp dbsnp.maf10pct.hg38.bed chi:~/tmp
        ssh chi:~/tmp
        bedtools intersect -v -header -a dbsnp.maf10pct.hg38.bed -b \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/bad_promoter.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc0to15.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc70to75.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc75to80.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc80to85.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc85to100.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_51to200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_gt200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/repeats.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/heng_universal_mask.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/LCR.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/sv_repeat_telomere_centromere.bed.gz \
        > dbsnp_maf10pct.hg38.no_selfchain_gc25-30_65-70_lowcomp50.bed

        cat dbsnp_maf10pct.hg38.no_selfchain_gc25-30_65-70_lowcomp50.bed | grep -v chrY | grep -v chrX \
        > dbsnp_maf10pct.hg38.no_selfchain_gc25-30_65-70_lowcomp50.autosomal.bed
        bgzip dbsnp_maf10pct.hg38.no_selfchain_gc25-30_65-70_lowcomp50.autosomal.bed
        tabix dbsnp_maf10pct.hg38.no_selfchain_gc25-30_65-70_lowcomp50.autosomal.bed.gz
    '''
    # exclude duplicates - remove each second from this output:
    ''' uniq -d dbsnp_maf10pct.good.bed
    '''
    # split into autosomal (for fingerprinting) and chrX (for sex check)

    # TODO: check clustered in HapMap?

    # TODO: include with some of more moderate GC numbers, and self chain?
    '''
  - bad_promoter
  - gc0to15
  - gc15to20
  - gc20to25
  - gc25to30                    - keep?
  - gc65to70                    - keep?
  - gc70to75
  - gc75to80
  - gc80to85
  - gc85to100
  - low_complexity_lt51bp       - keep?
  - low_complexity_51to200bp
  - low_complexity_gt200bp
  - repeats
  - self_chain                  - keep?
  - heng_universal_mask         - keep?

TODO: annotated, and prioritize when small number left

bedtools intersect -v -header -a dbsnp.maf10pct.bed -b \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/bad_promoter.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc0to15.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc70to75.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc75to80.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc80to85.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc85to100.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_51to200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_gt200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/repeats.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/heng_universal_mask.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/LCR.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/sv_repeat_telomere_centromere.bed.gz \
        > dbsnp_maf10pct.no_selfchain_gc25-30_65-70_lowcomp50.bed

bedtools intersect -v -header -a dbsnp.maf10pct.bed -b \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/self_chain.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/bad_promoter.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc0to15.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc25to30.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc70to75.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc75to80.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc80to85.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc85to100.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_51to200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_gt200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/repeats.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/heng_universal_mask.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/LCR.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/sv_repeat_telomere_centromere.bed.gz \
        > dbsnp_maf10pct.no_lowcomplexity1to50.bed

bedtools intersect -v -header -a dbsnp.maf10pct.bed -b \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/self_chain.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/bad_promoter.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc0to15.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc25to30.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc75to80.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc80to85.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc85to100.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_lt51bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_51to200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_gt200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/repeats.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/heng_universal_mask.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/LCR.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/sv_repeat_telomere_centromere.bed.gz \
        > dbsnp_maf10pct.no_gc70tp75.bed

bedtools intersect -v -header -a dbsnp.maf10pct.bed -b \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/self_chain.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/bad_promoter.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc0to15.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc70to75.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc75to80.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc80to85.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc85to100.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_lt51bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_51to200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_gt200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/repeats.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/heng_universal_mask.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/LCR.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/sv_repeat_telomere_centromere.bed.gz \
        > dbsnp_maf10pct.no_gc25tp30.bed

bedtools intersect -v -header -a dbsnp.maf10pct.bed -b \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/self_chain.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/bad_promoter.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc0to15.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc75to80.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc80to85.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc85to100.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_lt51bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_51to200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_gt200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/repeats.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/heng_universal_mask.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/LCR.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/sv_repeat_telomere_centromere.bed.gz \
        > dbsnp_maf10pct.no_gc.bed

bedtools intersect -v -header -a dbsnp.maf10pct.bed -b \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/bad_promoter.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc0to15.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc75to80.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc80to85.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc85to100.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_51to200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/low_complexity_gt200bp.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/repeats.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/heng_universal_mask.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/LCR.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/sv_repeat_telomere_centromere.bed.gz \
        > dbsnp_maf10pct.no_selfchain_lowcomplexity1to50_gc.bed


        1772719   dbsnp_maf10pct.good.bed
        2045141   dbsnp_maf10pct.no_gc25tp30.bed
        1817857   dbsnp_maf10pct.no_gc70tp75.bed
        2059411   dbsnp_maf10pct.no_gc.bed
        1819416   dbsnp_maf10pct.no_lowcomplexity1to50.bed
        1853549   dbsnp_maf10pct.noselfchain.bed
        2134550   dbsnp_maf10pct.no_selfchain_lowcomplexity1to50_gc.bed
        2578547   dbsnp.maf10pct.bed
    '''

    return subset_dbsnp


if __name__ == '__main__':
    subset_dbsnp(sys.argv[1])
