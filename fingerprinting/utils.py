from Bio import SeqIO


def read_fasta(fasta_fpath):
    seq_by_id = dict()
    with open(fasta_fpath) as f:
        reference_records = SeqIO.parse(f, 'fasta')
        for record in reference_records:
            seq_by_id[record.id] = record.seq
    return seq_by_id


