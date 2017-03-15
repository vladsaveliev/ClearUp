# from pyfaidx import Fasta
# records = Fasta('/Users/vlad/googledrive/az/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa')

import sys

with open(sys.argv[1]) as bed_rsid_f:
    rsid_lines = [l for l in bed_rsid_f if l.strip() and not l.startswith('#')]
with open(sys.argv[2]) as bed_gene_f:
    gene_lines = [l for l in bed_gene_f if l.strip() and not l.startswith('#')]

i = 0
for l_rsid, l_gene in zip(rsid_lines, gene_lines):
    chrom, pos0, pos1, ann = l_rsid.split()[:4]
    _, rsid, ref = ann.split('-')

    chrom2, pos02, pos12, gene = l_gene.split()[:4]

    assert (chrom, pos0, pos1) == (chrom2, pos02, pos12), str(i) + ' ' + rsid

    i += 1

    new_ann = rsid + '|'
    if gene != '.':
        new_ann += gene

    sys.stdout.write('\t'.join([chrom, pos0, pos1, new_ann]) + '\n')
