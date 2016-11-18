# from pyfaidx import Fasta
# records = Fasta('/Users/vlad/googledrive/az/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa')

met = set()
with open('/Users/vlad/vagrant/Fingerprinting/test/2016-10-11_bcbio/fingerprints/sambambadepth.txt') as f, \
        open('locations.bed', 'w') as out_bed:
    i = 0
    for l in f:
        if l[0] == '#':
            pass
        else:
            fs = l.strip().split()
            chrom, pos0, pos1, rsid = l.strip().split()[:4]
            if (chrom, pos1) in met:
                pass
            else:
                met.add((chrom, pos1))
                i += 1
                out_bed.write('\t'.join([chrom, pos0, pos1, str(i) + '-' + rsid]) + '\n')  #, str(records[chrom][int(pos0)].seq).upper()])
