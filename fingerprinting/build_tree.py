#!/usr/bin/env python2

import sys
from sys import platform
import os
from os.path import join, dirname

from Bio import Phylo
import json

from Utils.call_process import run
from Utils.file_utils import can_reuse, safe_mkdir
from Utils import logger

from fingerprinting.utils import read_fasta


suffix = 'lnx' if 'linux' in platform else 'osx'
prank_bin = join(dirname(__file__), 'prank', 'prank_' + suffix, 'bin', 'prank')


def build_tree(fasta_fpath, output_dirpath):
    output = os.path.join(output_dirpath, os.path.splitext(os.path.basename(fasta_fpath))[0])
    tree_fpath = os.path.join(output + '.best.dnd')
    if not can_reuse(tree_fpath, cmp_f=fasta_fpath):
        run(prank_bin + ' -d=' + fasta_fpath + ' -o=' + output + ' -showtree')
    tree = next(Phylo.parse(tree_fpath, 'newick'))
    # print(tree)
    # tree.ladderize()
    return tree


