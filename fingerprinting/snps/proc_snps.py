# import sys

# with open(sys.argv[1]) as bed_rsid_f:   # chr1  123123  123124  ?-rs12312-G
#     rsid_lines = [l for l in bed_rsid_f if l.strip() and not l.startswith('#')]
# with open(sys.argv[2]) as bed_gene_f:   # chr1  123123  123124  GENE
#     gene_lines = [l for l in bed_gene_f if l.strip() and not l.startswith('#')]

# i = 0
# for l_rsid, l_gene in zip(rsid_lines, gene_lines):
#     chrom, pos0, pos1, ann = l_rsid.split()[:4]
#     _, rsid, ref = ann.split('-')

#     chrom2, pos02, pos12, gene = l_gene.split()[:4]

#     assert (chrom, pos0, pos1) == (chrom2, pos02, pos12), str(i) + ' ' + rsid

#     i += 1

#     new_ann = rsid + '|'
#     if gene != '.':
#         new_ann += gene

#     sys.stdout.write('\t'.join([chrom, pos0, pos1, new_ann]) + '\n')

import sys

with open('exome_snps.txt') as f:
    for l in f:
        chrom, pos1, alt, ref, rsid = l.strip().split('\t')
        print '\t'.join([chrom, str(int(pos1)-1), pos1, rsid])
