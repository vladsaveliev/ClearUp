import glob
import shutil

import itertools
from collections import defaultdict, OrderedDict
import os
from os.path import isfile

from Bio import SeqIO
from ngs_utils.call_process import run
from os.path import join, dirname

from ngs_utils.file_utils import verify_file, safe_mkdir, file_transaction, which
import ngs_utils.logger as log

from fingerprinting.model import Fingerprint


FASTA_ID_PROJECT_SEPARATOR = '____PROJECT_'


def read_fasta(fasta_fpath):
    seq_by_id = OrderedDict()
    with open(fasta_fpath) as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            seq_by_id[record.id] = record.seq
    return seq_by_id


def load_bam_file(bcbio_final_path, project_work_dirpath, sample_id):
    """
    Assuming the BAM files are sliced to fingerprints locations
    """
    bam_glob = join(bcbio_final_path, sample_id, '*-ready.bam')
    bam_fpath = next(iter(glob.glob(bam_glob)), None)
    if not verify_file(bam_fpath):
        log.critical('BAM file not found in ' + bam_glob)
    bam_index_fpath = bam_fpath + '.bai'
    if not verify_file(bam_index_fpath):
        log.critical('BAM index file not found in ' + bam_index_fpath)
    bam_copy_fpath = join(safe_mkdir(join(project_work_dirpath, 'bams')), sample_id + '.bam')
    bam_index_copy_fpath = bam_copy_fpath + '.bai'
    shutil.copy(bam_fpath, bam_copy_fpath)
    shutil.copy(bam_index_fpath, bam_index_copy_fpath)
    return bam_copy_fpath


def get_snps_file():
    return verify_file(join(dirname(__file__), 'snps', 'idt_snps.bed'), is_critical=True)


def get_fingerprints(project, seq_by_sample_name):
    fp_by_loc_by_sample = defaultdict(dict)
    fp_by_index_by_sample = defaultdict(dict)
    with open(get_snps_file()) as f:
        for i, l in enumerate(l for l in f if l[0] != '#'):
            chrom, pos0, pos1, ann = l.strip().split()
            for s in project.samples:
                fp = Fingerprint(index=i + 1, sample=s, chrom=chrom, pos=int(pos1), rsid=ann.split('-')[1])
                fp_by_loc_by_sample[s.name][(chrom, int(pos1))] = fp
                fp_by_index_by_sample[s.name][i + 1] = fp

    for (chrom, _1, pos, rsid, depth, _2, sname) in parse_sambambadepth(project.bcbio_final_path):
        if sname in fp_by_loc_by_sample:
            assert (chrom, int(pos)) in fp_by_loc_by_sample[sname], (chrom + ':' + str(pos) + ' not found in ' + sname)
            fp_by_loc_by_sample[sname][(chrom, int(pos))].depth = depth

    for sname, seq in seq_by_sample_name.items():
        for i in range(0, len(seq), 2):
            index0 = i / 2
            genotype = str(seq[i: i + 2]).upper()
            fp_by_index_by_sample[sname][index0 + 1].genotype = genotype

    return fp_by_loc_by_sample


#     alt_by_loc = dict()
#     for (sname, chrom, pos, ids, ref, alt) in parse_vardict_txt(bcbio_final_path):
#         if sname == sample.name:
#             fp = fp_by_loc.get((chrom, int(pos)))
#             if fp and fp.genotype is not None:
#                 assert fp.ref == ref, fp.ref + ' ' + ref
#                 alt_by_loc[(chrom, int(pos))] = alt
#
#     fingerprints = sorted(fp_by_loc.values(), key=lambda _fp: _fp.index)
#     vcf_fpath = join(safe_mkdir(join(project_work_dirpath, 'vcfs')), sample.name + '.vcf')
#     with file_transaction(None, vcf_fpath) as tx:
#         with open(tx, 'w') as out:
#             out.write('##fileformat=VCFv4.2\n')
#             out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
#             for fp in fingerprints:
#                 fs = [fp.chrom, str(fp.pos), fp.rsid, 'N', alt_by_loc.get((fp.chrom, fp.pos), fp.ref), '.', '.', '.']
#                 out.write('\t'.join(fs) + '\n')
#     igvtools_index(vcf_fpath)
#     return vcf_fpath
#
#
# def igvtools_index(vcf_fpath):
#     if not which('igvtools'):
#         log.critical('igvtools is not installed!')
#
#     cmdline = 'igvtools index {vcf_fpath}'.format(**locals())
#     run(cmdline)
#     if isfile('igv.log'):
#         try:
#             os.remove('igv.log')
#         except OSError:
#             pass
#     return vcf_fpath + '.idx'


def parse_sambambadepth(bcbio_final_path):
    # sambamba_glob = join(bcbio_final_path, '20??-??-??*', 'fingerprints', 'sambambadepth.txt')
    # sambamba_fpath = next(iter(glob.glob(sambamba_glob)), None)
    # if not verify_file(sambamba_fpath):
    #     log.critical('Fingerprints sambambadepth file not found in ' + sambamba_glob)
    with open(sambamba_fpath) as in_f:
        for line in in_f:
            if not line or line.startswith('#'):
                continue
            chrom, _1, pos, id, depth, _2, sample = line.strip().split('\t')
            yield chrom, _1, int(pos), id, int(depth), _2, sample


def parse_vardict_txt(bcbio_final_path):
    vardict_glob = join(bcbio_final_path, '20??-??-??*', 'var', 'vardict.txt')
    vardict_fpath = next(iter(glob.glob(vardict_glob)), None)
    if not verify_file(vardict_fpath):
        log.critical('Fingerprints vardict.txt file not found in ' + vardict_glob)
    with open(vardict_fpath) as in_f:
        for line in in_f:
            if not line or line.startswith('Sample'):
                continue
            sname, chrom, pos, rsid, ref, alt = line.strip('\n').split('\t')[:6]
            yield sname, chrom, int(pos), rsid, ref, alt


def get_sample_and_project_name(name, fingerprint_project=None):
    project_names = fingerprint_project.split(',') if fingerprint_project else []
    if FASTA_ID_PROJECT_SEPARATOR in name:
        sample_name, project_name = name.split(FASTA_ID_PROJECT_SEPARATOR)
    else:
        sample_name, project_name = name, ''
    if len(project_names) != 1:
        return sample_name, project_name
    else:
        return sample_name, project_names[0]


def calculate_distance_matrix(tree):
    clades = tree.get_terminals()
    distance_matrix = defaultdict(lambda: (1, None))
    for clade, clade2 in list(itertools.combinations(clades, 2)):
        distance = tree.distance(clade, clade2)
        if distance < distance_matrix[clade][0]:
            distance_matrix[clade] = (distance, clade2)
            distance_matrix[clade2] = (distance, clade)
    return distance_matrix
