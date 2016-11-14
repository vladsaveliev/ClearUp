import glob
import shutil

import itertools
from collections import defaultdict

from Bio import SeqIO
from os.path import join

from Utils.file_utils import verify_file
import Utils.logger as log
from fingerprinting.model import Fingerprint

FASTA_ID_PROJECT_SEPARATOR = '____PROJECT_'

def read_fasta(fasta_fpath):
    seq_by_id = dict()
    with open(fasta_fpath) as f:
        reference_records = SeqIO.parse(f, 'fasta')
        for record in reference_records:
            seq_by_id[record.id] = record.seq
    return seq_by_id


def load_bam_file(bcbio_final_path, project_work_dirpath, sample_id):
    bam_glob = join(bcbio_final_path, sample_id, '*-ready.bam')
    bam_fpath = next(iter(glob.glob(bam_glob)), None)
    if not verify_file(bam_fpath):
        log.critical('BAM file not found in ' + bam_glob)
    bam_copy_fpath = join(project_work_dirpath, sample_id + '.bam')
    bam_index_fpath = bam_fpath + '.bai'
    bam_index_copy_fpath = bam_copy_fpath + '.bai'
    shutil.copy(bam_fpath, bam_copy_fpath)
    shutil.copy(bam_index_fpath, bam_index_copy_fpath)
    return bam_copy_fpath


def get_fingerprints(project, sample, seq):
    fingerprints = []
    fingerprints_positions = get_fingerprints_positions(project.bcbio_final_path)
    for i in range(0, len(seq), 2):
        index0 = i / 2
        snp_index = index0 + 1
        snp = seq[i: i+2]
        chrom, pos = fingerprints_positions[index0]
        fingerprint = Fingerprint(sample, snp_index, chrom, pos, snp)
        fingerprints.append(fingerprint)
    return fingerprints


def get_fingerprints_positions(bcbio_final_path):
    sambamba_glob = join(bcbio_final_path, '20??-??-??*', 'fingerprints', 'sambambadepth.txt')
    sambamba_fpath = next(iter(glob.glob(sambamba_glob)), None)
    if not verify_file(sambamba_fpath):
        log.critical('Fingerprints fasta file not found in ' + sambamba_glob)
    fingerprints_positions = []
    with open(sambamba_fpath) as in_f:
        prev_pos = 0
        for line in in_f:
            if not line or line.startswith('#'):
                continue
            fs = line.split()
            chrom, pos = fs[0], fs[1]
            if pos != prev_pos:
                fingerprints_positions.append((chrom, pos))
                prev_pos = pos
    return fingerprints_positions


def get_sample_and_project_name(name, fingerprint_project=None):
    project_names = fingerprint_project.split(',') if fingerprint_project else []
    if FASTA_ID_PROJECT_SEPARATOR in name:
        sample_name, project = name.split(FASTA_ID_PROJECT_SEPARATOR)
    else:
        sample_name, project = name, ''
    if len(project_names) != 1:
        return sample_name, project
    else:
        return sample_name, project_names[0]


def calculate_distance_matrix(tree):
    clades = tree.get_terminals()
    distance_matrix = defaultdict(lambda: (1, None))
    for (clade, clade2) in list(itertools.combinations(clades, 2)):
        distance = tree.distance(clade, clade2)
        if distance < distance_matrix[clade][0]:
            distance_matrix[clade] = (distance, clade2)
            distance_matrix[clade2] = (distance, clade)
    return distance_matrix