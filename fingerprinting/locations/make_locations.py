# from pyfaidx import Fasta
# records = Fasta('/Users/vlad/googledrive/az/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa')

with open('locations.bed') as f, \
        open('locations.igv.bed', 'w') as out_bed:
    i = 0
    for l in f:
        if not l.strip() or l[0] == '#':
            pass
        else:
            fs = l.strip().split()
            chrom, pos0, pos1, ann = l.strip().split()[:4]
            _, rsid, ref = ann.split('-')
            i += 1
            out_bed.write('\t'.join([chrom, pos0, pos1, str(i) + '-' + rsid]) + '\n')
                #, str(records[chrom][int(pos0)].seq).upper()])
