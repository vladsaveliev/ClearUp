#!/usr/bin/env python
from os.path import join, dirname, isfile, isdir, basename
from cyvcf2 import VCF

from ngs_utils import logger as log
from ngs_utils.call_process import run
from ngs_utils.file_utils import safe_mkdir, can_reuse, verify_file, file_transaction, verify_dir, add_suffix, splitext_plus

from variant_filtering import reference_data as filt_ref_data, get_anno_config
from variant_filtering.utils import parse_gene_blacklists, check_gene_in_a_blacklist, all_blacklisted_genes


dbsnp_location = '/Users/vlad/googledrive/az/reference_data/genomes/Hsapiens/hg19/variation/dbSNP_v149.vcf.gz'


def subset_dbsnp():
    dbsnp_file = verify_file(dbsnp_location, is_critical=True)

    incidentalome_dir = verify_dir(filt_ref_data.incidentalome_dir(), 'incidentalome')
    anno_cfg = get_anno_config()
    blacklisted_genes = all_blacklisted_genes(anno_cfg['blacklist']['genes'], incidentalome_dir)

    total_snps_written = 0
    current_chrom = None
    output_path = join(dirname(__name__), 'dbsnp.mafs10pct.bed')
    met_locations = set()
    with file_transaction(None, output_path) as tx:
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
                                locid = v.CHROM, v.POS
                                if locid not in met_locations:
                                    met_locations.add(locid)
                                    fields = [v.CHROM, str(v.POS - 1), str(v.POS), v.ID + '|' + gene, v.REF, ','.join(v.ALT)]
                                    out.write('\t'.join(fields) + '\n')
                                    if current_chrom != v.CHROM:
                                        current_chrom = v.CHROM
                                        print(v.CHROM + ', written ' + str(total_snps_written))
                                    total_snps_written += 1
                                
    # exclude those from LCR and tricky regions:
    ''' scp dbsnp.mafs10pct.bed chi:~/tmp
        ssh chi:~/tmp
        bedtools intersect -v -header -a dbsnp.mafs10pct.bed -b \
        ~/bcbio/genomes/Hsapiens/hg19/coverage/problem_regions/GA4GH/*.bed.gz \
        ~/bcbio/genomes/Hsapiens/hg19/coverage/problem_regions/ENCODE/wgEncodeDacMapabilityConsensusExcludable.bed.gz \
        ~/bcbio/genomes/Hsapiens/hg19/coverage/problem_regions/repeats/LCR.bed.gz \
        ~/bcbio/genomes/Hsapiens/hg19/coverage/problem_regions/repeats/sv_repeat_telomere_centromere.bed \
        > dbsnp_maf10pct.good.bed
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

bedtools intersect -v -header -a dbsnp.mafs10pct.bed -b \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/bad_promoter.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc0to15.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc25to30.bed.gz \
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
        > dbsnp_maf10pct.noselfchain.bed
        
bedtools intersect -v -header -a dbsnp.mafs10pct.bed -b \
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

bedtools intersect -v -header -a dbsnp.mafs10pct.bed -b \
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

bedtools intersect -v -header -a dbsnp.mafs10pct.bed -b \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/self_chain.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/bad_promoter.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc0to15.bed.gz \
        /home/vsaveliev/NGS_Reporting/az/reference_data/tricky_regions/hg19/gc20to25.bed.gz \
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

bedtools intersect -v -header -a dbsnp.mafs10pct.bed -b \
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

bedtools intersect -v -header -a dbsnp.mafs10pct.bed -b \
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
        2578547   dbsnp.mafs10pct.bed
    '''
    
    return subset_dbsnp


if __name__ == '__main__':
    subset_dbsnp()
