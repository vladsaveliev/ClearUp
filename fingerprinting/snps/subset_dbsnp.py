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
  - gc25to30 - keep?
  - gc65to70 - keep?
  - gc70to75
  - gc75to80
  - gc80to85
  - gc85to100
  - low_complexity_lt51bp
  - low_complexity_51to200bp
  - low_complexity_gt200bp
  - repeats
  - self_chain - keep?
  - heng_universal_mask - keep?
    '''
    
    return subset_dbsnp


if __name__ == '__main__':
    subset_dbsnp()
