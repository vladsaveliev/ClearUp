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


def draw_tree(fasta_fpath, output_dirpath):
    safe_mkdir(output_dirpath)
    tree = _build_tree(fasta_fpath, output_dirpath)
    seq_by_id = read_fasta(fasta_fpath)
    return _make_csv_for_d3(tree, seq_by_id)


def _make_csv_for_d3(tree, seq_by_id):
    clade_dicts = []
    _clade_to_json(tree.root, cur_name='', clade_dicts=clade_dicts, seq_by_id=seq_by_id)
    json_str = json.dumps(clade_dicts)
    print json_str
    return json_str


def _build_tree(fasta_fpath, output_dirpath):
    output = os.path.join(output_dirpath, os.path.splitext(os.path.basename(fasta_fpath))[0])
    tree_fpath = os.path.join(output + '.best.dnd')
    if not can_reuse(tree_fpath, cmp_f=fasta_fpath):
        run(prank_bin + ' -d=' + fasta_fpath + ' -o=' + output + ' -showtree')
    tree = next(Phylo.parse(tree_fpath, 'newick'))
    # print(tree)
    # tree.ladderize()
    return tree


def _draw_clade(clade, cur_name, out, seq_by_id):
    separator = '\t'
    if clade.name:  # leaf
        leaf_name = clade.name
        line = cur_name + separator + leaf_name if cur_name else leaf_name
        seq = str(seq_by_id[clade.name])
        line += ',' + seq
    else:  # internal node
        cur_name = cur_name + separator + str(id(clade)) if cur_name else str(id(clade))
        line = cur_name + ','

    out.write(line + '\n')
    if clade.clades:
        for child in clade:
            _draw_clade(child, cur_name, out, seq_by_id)


def _clade_to_json(clade, cur_name, clade_dicts, seq_by_id):
    separator = '.'
    if clade.name:  # leaf
        leaf_name = clade.name
        clade_d = {
            'id': cur_name + separator + leaf_name if cur_name else leaf_name,
            'seq': _seq_to_json(seq_by_id[clade.name]),
        }
    else:  # internal node
        cur_name = cur_name + separator + str(id(clade)) if cur_name else str(id(clade))
        clade_d = {
            'id': cur_name,
            'seq': None
        }
    clade_dicts.append(clade_d)
    if clade.clades:
        for child in clade:
            _clade_to_json(child, cur_name, clade_dicts, seq_by_id)


class_by_nuc = {
    'C': 'c_nuc',
    'A': 'a_nuc',
    'T': 't_nuc',
    'G': 'g_nuc',
    'N': 'n_nuc',
}
def _seq_to_json(seq):
    return [(nuc, class_by_nuc[nuc.upper()]) for nuc in seq]


if __name__ == '__main__':
    draw_tree(sys.argv[1], sys.argv[2])


