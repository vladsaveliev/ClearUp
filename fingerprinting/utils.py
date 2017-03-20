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


FASTA_ID_PROJECT_SEPARATOR = '____PROJECT_'


def read_fasta(fasta_fpath):
    seq_by_id = OrderedDict()
    with open(fasta_fpath) as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            seq_by_id[record.id] = record.seq
    return seq_by_id


def load_bam_file(bam_fpath, bams_dir, sample_name):
    """ Assuming the BAM files are sliced to fingerprint locations
    """
    # TODO: slice BAM here
    bam_index_fpath = bam_fpath + '.bai'
    if not verify_file(bam_index_fpath):
        log.critical('BAM index file not found in ' + bam_index_fpath)
    bam_copy_fpath = join(bams_dir, sample_name + '.bam')
    bam_index_copy_fpath = bam_copy_fpath + '.bai'
    shutil.copy(bam_fpath, bam_copy_fpath)
    shutil.copy(bam_index_fpath, bam_index_copy_fpath)
    return bam_copy_fpath


def get_snps_file():
    return verify_file(join(dirname(__file__), 'snps', 'idt_snps.bed'), is_critical=True)


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


def is_sex_chrom(chrom):
    # type: (object) -> object
    return chrom in ['X', 'Y', 'chrX', 'chrY']
