import glob
import shutil

import itertools
from collections import defaultdict, OrderedDict
import os
from os.path import isfile

from Bio import SeqIO
from pybedtools import BedTool

from ngs_utils.call_process import run
from os.path import join, dirname

from ngs_utils.file_utils import verify_file, safe_mkdir, file_transaction, which, can_reuse
import ngs_utils.logger as log
from ngs_utils.sambamba import index_bam

FASTA_ID_PROJECT_SEPARATOR = '____PROJECT_'


def read_fasta(fasta_fpath):
    seq_by_id = OrderedDict()
    with open(fasta_fpath) as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            seq_by_id[record.id] = record.seq
    return seq_by_id


def load_bam_file(bam_file, bams_dir, snp_bed, sample_name):
    """ Slicing to fingerprints locations
    """
    bam_index_file = bam_file + '.bai'
    if not verify_file(bam_index_file):
        log.critical('BAM index file not found in ' + bam_index_file)
    sliced_bam_file = join(bams_dir, sample_name + '.bam')
    if not can_reuse(sliced_bam_file, [bam_file, snp_bed]):
        cmdl = 'sambamba view {bam_file} -L {snp_bed} -f bam -o {sliced_bam_file}'.format(**locals())
        run(cmdl, output_fpath=sliced_bam_file, stdout_to_outputfile=False)
        index_bam(sliced_bam_file)
    return sliced_bam_file


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


def is_sex_chrom(chrom):
    return chrom in ['X', 'Y', 'chrX', 'chrY']

