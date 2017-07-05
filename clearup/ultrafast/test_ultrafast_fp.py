#!/usr/bin/env python

import sys
import os
from os import listdir
from os.path import join, basename, dirname
from glob import glob
from collections import defaultdict
from enum import Enum
import hashlib

import itertools as it
import numpy as np
import scipy.stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from cyvcf2 import VCF
from ngs_utils.file_utils import splitext_plus, can_reuse, safe_mkdir, \
    intermediate_fname, adjust_path, str_to_filename
from ngs_utils.call_process import run
from ngs_utils import logger as log
from clearup.panel import overlap_bed_files

# os.environ['REUSE'] = '0'


class Params:
    L = 20
    NORMALIZE_VAR = True
    NORMALIZE_DIST = True
    SLICE = True
    SKIP_DAMAGE = False
    SKIP_REJECT = False
    SKIP_NOCALL = False
    MIN_AF = 0.25


basedir = '/Users/vladsaveliev/vagrant/ClearUp/clearup/ultrafast'
# input_subdirs = 'platinum_hg19 platinum_hg38'.split()
input_subdirs = 'umi_onseq az600_probes az600_targets platinum_hg19 platinum_hg38'.split()
# input_subdirs = 'az600_probes az600_targets platinum_hg19 platinum_hg38'.split()
# input_subdirs = 'umi_onseq'.split()


def main():
    title = ', '.join(input_subdirs) + ', L=' + str(Params.L) + ''
    if Params.SLICE: title += ', sliced to BED files'
    else:            title += ', not sliced to BED files'

    if not Params.NORMALIZE_VAR:  title += ', not normalized by variants'
    if not Params.NORMALIZE_DIST: title += ', not normalized by distance'
    if Params.SKIP_DAMAGE: title += ', skipped damage'
    if Params.SKIP_REJECT: title += ', skipped REJECT'
    if Params.SKIP_NOCALL: title += ', skipped num called = 0'
    if Params.MIN_AF:      title += ', min AF=' + str(Params.MIN_AF)

    run_id = '__'.join(input_subdirs)
    if not Params.SLICE: run_id += '_not_sliced'

    run_dir = safe_mkdir(join(basedir, 'runs', run_id))
    work_dir = safe_mkdir(join(run_dir, 'work'))
    png_file = join(run_dir, str_to_filename(title) + '.png')

    log.init(True, join(run_dir, 'log.txt'))

    vcf_by_label = dict()
    all_bed_files = []
    for sd in input_subdirs:
        dp = join(basedir, sd)
        all_bed_files.extend(glob(join(dp, '*.bed')))
        for vcf_fpath in glob(join(dp, '*.vcf.gz')):
            label = join(sd, basename(splitext_plus(vcf_fpath)[0]))
            assert label not in vcf_by_label
            vcf_by_label[label] = vcf_fpath

    if Params.SLICE and all_bed_files:
        log.info('Slicing VCF files to regions in ' + str(all_bed_files))
        fnames = [basename(splitext_plus(fp)[0]) for fp in all_bed_files]
        target_id = '__'.join(fnames)
        overlapped_bed = join(work_dir, target_id + '.bed')
        if not can_reuse(overlapped_bed, all_bed_files):
            overlap_bed_files(all_bed_files, overlapped_bed)
        for label, vcf_fpath in list(vcf_by_label.items()):
            sliced_vcf_fpath = intermediate_fname(work_dir, vcf_fpath, 'sliced')
            if not can_reuse(sliced_vcf_fpath, vcf_fpath):
                run(f'bcftools view {vcf_fpath} --targets-file {overlapped_bed} -o {sliced_vcf_fpath}')
            vcf_by_label[label] = sliced_vcf_fpath

    log.info('Calculating fingerprints for individual samples')
    all_prints = [(label, make_fingerprint(vcf_fpath, work_dir, label, fp_size=Params.L))
                  for label, vcf_fpath in sorted(vcf_by_label.items())]

    log.info('Comparing fingerprints pairwise')
    pairwise = defaultdict(dict)
    for ((label1, print1), (label2, print2)) in it.combinations_with_replacement(all_prints, 2):
        dist, pvalue = compare(print1, print2)
        pairwise[label1][label2] = dist
        pairwise[label2][label1] = dist
        log.info(f'   {label1} VS {label2}: {dist:.2f}, Pvalue={pvalue:.2f}')

    log.info('Plotting comparison')
    df = pd.DataFrame(data=pairwise)
    log.info(df)

    # Generate a mask for the upper triangle + main diagonale
    mask = np.zeros_like(df, dtype=np.bool)
    mask[np.triu_indices_from(mask, k=1)] = True

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(10, 7))
    if title:
        ax.set_title(title)

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    g = sns.heatmap(df, vmin=0.5, vmax=1, center=0.75, mask=mask,
                    cmap=cmap, annot=True, fmt='.2f', ax=ax,
                    annot_kws={'size': 8})
    g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=7)
    g.set_xticklabels(g.get_xticklabels(), rotation=90, fontsize=7)
    sns.set(font_scale=2)
    sns.plt.subplots_adjust(left=0.2, right=1, top=0.95, bottom=0.26)

    sns.plt.savefig(png_file)
    log.info('')
    log.info('Saved heatmap into ' + adjust_path(png_file))


def compare(fp1, fp2):
    try:
        res = scipy.stats.spearmanr(fp1.flatten(), fp2.flatten())
    except ValueError as e:
        log.err(e)
        log.err('Error calculating correlation between fingerpirnts, '
                'likely too small numbrer of mutations. Try encreasing target '
                'size or filtering criteria, or decrease L.')
        sys.exit(1)
    else:
        return res.correlation, res.pvalue


possible_keys = [''.join(v[0]+v[1]) for v in it.product(it.permutations('ACGT', 2), repeat=2)]
index_by_key = {sys.intern(k): possible_keys.index(k) for k in possible_keys}


def make_fingerprint(vcf_file, work_dir=None, label=None, fp_size=20):
    log.info('Starting processing file ' + vcf_file)

    if label: print_name = label.replace('/', '_')
    else:     print_name = splitext_plus(basename(vcf_file))[0]
    print_name += '.print' + str(fp_size)
    work_dir = work_dir or dirname(vcf_file)

    raw_print_file = join(work_dir, print_name)
    if can_reuse(raw_print_file, vcf_file):
        with open(raw_print_file) as f:
            raw = np.fromfile(f).reshape((len(index_by_key), fp_size))
    else:
        raw = raw_fingerprint(vcf_file, fp_size=fp_size)
        with open(raw_print_file, 'w') as f:
            raw.tofile(f)
        log.info('Saved raw fingerprints into ' + raw_print_file)

    norm_print_file = raw_print_file
    if Params.NORMALIZE_DIST: norm_print_file += '_normdist'
    if Params.NORMALIZE_VAR:  norm_print_file += '_normvar'

    if can_reuse(norm_print_file, raw_print_file):
        with open(norm_print_file) as f:
            norm = np.fromfile(f).reshape((len(index_by_key), fp_size))
    else:
        norm = normalize_fingerprint(raw)
        with open(norm_print_file, 'w') as f:
            norm.tofile(f)
        log.info('Saved normalised fingerprints into ' + norm_print_file)

    return norm


def raw_fingerprint(vcf_file, fp_size=20):
    raw = np.zeros((len(index_by_key), fp_size))
    prev_rec = None

    total_vars = 0
    rejected_af = 0
    rejected_other = 0
    passed = 0
    pairs = 0

    for rec in VCF(vcf_file):
        total_vars += 1
        if Params.SKIP_REJECT and rec.FILTER:
            continue
        if Params.SKIP_NOCALL and rec.num_called == 0:
            continue
        if Params.SKIP_DAMAGE and rec.INFO.get('DKFZBias'):
            continue
        if 'X' in rec.CHROM or 'Y' in rec.CHROM:
            rejected_other += 1
            continue  # remove sex chromosome
        if len(rec.ALT) > 1:
            rejected_other += 1
            continue  # many ALTs
        if len(rec.REF) > 1 or any(len(a) > 1 for a in rec.ALT):
            rejected_other += 1
            continue  # indel or complex variant
        if rec.INFO.get('AF') and rec.INFO.get('AF') < Params.MIN_AF:
            rejected_af += 1
            continue
        if 'N' in rec.REF or 'N' in rec.ALT[0]:
            rejected_other += 1
            continue  # likely REF or ALT is N
        passed += 1
        if prev_rec is None or prev_rec.CHROM != rec.CHROM:
            prev_rec = rec
            continue  # looping into a new chromosome
        pair_dist = rec.POS - prev_rec.POS
        if pair_dist < 5:  # minimal allowed distance between variants (to avoid complex variants)
            continue  # too closely located variants
        pair_key = sys.intern(prev_rec.REF + prev_rec.ALT[0] + rec.REF + rec.ALT[0])
        raw[index_by_key[pair_key], pair_dist % fp_size] += 1
        prev_rec = rec
        pairs += 1

    log.info('total_vars: ' + str(total_vars) +
             ', rejected_af: ' + str(rejected_af) +
             ', rejected_other: ' + str(rejected_other) +
             ', passed: ' + str(passed) +
             ', pairs: ' + str(pairs))

    return raw


def normalize_fingerprint(arr):
    if Params.NORMALIZE_DIST:
        arr = arr - np.mean(arr, axis=0)
        arr = arr / np.std(arr, axis=0)

    if Params.NORMALIZE_VAR:
        arr = arr - np.mean(arr, axis=1)[:, np.newaxis]
        arr = arr / np.std(arr, axis=1)[:, np.newaxis]

    return arr


def _hash(s, l=10):
    return str(int(hashlib.sha256(s.encode('utf-8')).hexdigest(), 16) % 10 ** l)


main()
