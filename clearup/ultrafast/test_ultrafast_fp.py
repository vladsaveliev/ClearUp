#!/usr/bin/env python

import sys
import os
from os.path import join, basename, dirname, isfile
from glob import glob
from collections import defaultdict
import hashlib

import click
import itertools as it
import numpy as np
import scipy.stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt; plt.switch_backend('agg')
from pybedtools import BedTool
from cyvcf2 import VCF
import matplotlib

from ngs_utils.file_utils import splitext_plus, can_reuse, safe_mkdir, adjust_path, str_to_filename
from ngs_utils.call_process import run
from ngs_utils import logger as log
from ngs_utils.parallel import ParallelCfg, parallel_view
from ngs_utils.bcbio import BcbioProject
from ngs_utils.bed_utils import bgzip_and_tabix, sort_bed


class Params:
    L = 20
    NORMALIZE_VAR = False
    NORMALIZE_DIST = False
    SKIP_DAMAGE = False
    SKIP_REJECT = False
    SKIP_NOCALL = False
    MIN_AF = 0.25
    MIN_DIST = 20
    INTERREGION_PAIRS = True  # when regional sequencing, allow snp pairs between regions


code_dir = dirname(__file__)
# input_subdirs = 'platinum_hg19 platinum_hg38'.split()
# input_subdirs = 'umi_onseq az600_probes az600_targets platinum_hg19 platinum_hg38'.split()
# input_subdirs = 'az600_probes az600_targets platinum_hg19 platinum_hg38'.split()
# input_subdirs = 'umi_onseq'.split()


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('subdirs',
                metavar='<Either directories with VCF files (and optional BED files), or bcbio projects locations. '
                        'If providing BED files, add the genome build name to the end of the directory name '
                        '(e.g. platinum_hg19:genome>',
                nargs=-1)
@click.option('-o', '--output_dir',
              type=click.Path(dir_okay=True),
              metavar='<output directory>')
@click.option('-t', '--threads',
              type=int,
              help='Number of threads. Default is in corresponding system_info_*.yaml or 1. ' +
                   'If set to 1, skip starting cluster even if scheduler is specified.',
              default=1)
@click.option('-L', '--fp_size',
              type=int,
              help='Fingerprint size. For small targets, use values below 20.',
              default=Params.L)
@click.option('-d', '--isdebug',
              is_flag=True,
              default=False)
@click.pass_context
def main(ctx, subdirs, output_dir, threads=1, fp_size=Params.L, isdebug=False):
    """ Generates a PNG image with a relatedness heatmap.
    """
    if not subdirs:
        ctx.fail('Provide at least on input directory.')

    datasets = _load_datasets(subdirs)

    title = ', '.join(d.name for d in datasets) + '\nL=' + str(fp_size) + ''
    if not Params.NORMALIZE_VAR:  title += ', not norm by var'
    if not Params.NORMALIZE_DIST: title += ', not norm by dist'
    if Params.SKIP_DAMAGE:        title += ', skipped damage'
    if Params.SKIP_REJECT:        title += ', skipped REJECT'
    if Params.SKIP_NOCALL:        title += ', skipped num called = 0'
    if Params.MIN_AF:             title += ', min AF=' + str(Params.MIN_AF)
    if Params.MIN_DIST:           title += ', min dist=' + str(Params.MIN_DIST)
    if Params.INTERREGION_PAIRS:  title += ', used SNP pairs between regions'
    else:                         title += ', skiped SNP pairs between regions'

    run_id = '__'.join(d.name for d in datasets)

    run_dir = safe_mkdir(join((output_dir or join(code_dir, 'runs')), run_id))
    log.init(isdebug, join(run_dir, 'log.txt'), save_previous=True)
    work_dir = safe_mkdir(join(run_dir, 'work'))

    all_vcf_by_label = dict()
    bed_files_by_genome = defaultdict(set)
    for d in datasets:
        all_vcf_by_label.update(d.vcf_by_label)
        if d.bed_file:
            bed_files_by_genome[d.genome].add(d.bed_file)  # d.bed_file=None for WGS

    genome_by_label = dict()
    for d in datasets:
        for label in d.vcf_by_label:
            genome_by_label[label] = d.genome

    parallel_cfg = ParallelCfg(threads=threads)
    log.info(f'Starting using {parallel_cfg.threads} threads')

    with parallel_view(len(all_vcf_by_label), parallel_cfg, work_dir) as parall_view:
        overlap_bed_file_by_genome = dict()
        if bed_files_by_genome:
            overlap_bed_file_by_genome = _prep_bed(work_dir, bed_files_by_genome, overlap_bed_file_by_genome)
            log.info('Slicing VCFs to regions in BED files')
            out = parall_view.run(_slice_vcf_fn,
               [[work_dir, label, vcf, overlap_bed_file_by_genome.get(genome_by_label[label])]
                for label, vcf in all_vcf_by_label.items()])
            all_vcf_by_label = dict(out)
            log.info()

        log.info('Calculating fingerprints for individual samples')
        out = parall_view.run(make_fingerprint,
           [[vcf, work_dir, label, fp_size, overlap_bed_file_by_genome.get(genome_by_label[label])]
            for label, vcf in all_vcf_by_label.items()])
        print_label_pairs = dict(out)
        log.info()

    log.info('Comparing fingerprints pairwise')
    pairwise_dict = defaultdict(dict)
    for ((label1, print1), (label2, print2)) in it.combinations_with_replacement(print_label_pairs.items(), 2):
        dist, pvalue = compare(print1, print2)
        if dist:
            log.info(f'   {label1} VS {label2}: {dist:.2f}, Pvalue={pvalue:.2f}')
        else:
            log.info(f'   {label1} VS {label2}: failed to calculate')
            dist = float('NaN')
        pairwise_dict[label1][label2] = dist
        pairwise_dict[label2][label1] = dist

    log.info('Plotting comparison heatmap')
    plot_heatmap(pairwise_dict, run_dir, title)


def _prep_bed(work_dir, bed_files_by_genome, overlap_bed_file_by_genome=None):
    overlap_bed_file_by_genome = dict()
    log.info(f'Found BED files: {bed_files_by_genome}')
    for genome, bed_files in bed_files_by_genome.items():
        bed_files = [b for b in bed_files if b]
        log.info(f'Overlapping BED files for genome {genome}')
        overlap_bed_file_by_genome[genome] = _overlap_bed_files(bed_files, work_dir, genome) \
            if bed_files else None

    primary_genome = sorted(bed_files_by_genome.items(), key=lambda kv: len(kv[1]))[-1][0]
    lifted_bed_files = []
    for genome, overlap_bed_file in overlap_bed_file_by_genome.items():
        if overlap_bed_file and genome != primary_genome:
            from clearup.panel import lift_over
            lifted_bed_file = lift_over(overlap_bed_file, genome, primary_genome)
            lifted_bed_files.append(lifted_bed_file)
    if lifted_bed_files:
        primary_bed_files = [b for b in lifted_bed_files + [overlap_bed_file_by_genome[primary_genome]] if b]
        overlap_bed_file_by_genome[primary_genome] = _overlap_bed_files(
            primary_bed_files, work_dir, primary_genome)

    log.info('Lifting BED files back')
    for genome in overlap_bed_file_by_genome:
        if genome != primary_genome:
            from clearup.panel import lift_over
            overlap_bed_file_by_genome[genome] = lift_over(
                overlap_bed_file_by_genome[primary_genome], primary_genome, genome)
    log.info()

    log.info('Sorting, bgzipping and tabixing BED files')
    for g, bed in overlap_bed_file_by_genome.items():
        overlap_bed_file_by_genome[g] = bgzip_and_tabix(sort_bed(bed, genome=g))
    log.info()
    return overlap_bed_file_by_genome


def compare(fp1, fp2):
    try:
        res = scipy.stats.spearmanr(fp1.flatten(), fp2.flatten())
    except ValueError as e:
        log.err(e)
        log.err('Error calculating correlation between fingerpirnts, '
                'likely too small numbrer of mutations. Try encreasing target '
                'size or filtering criteria, or decrease L.')
        return None, None
    else:
        return res.correlation, res.pvalue


def _overlap_bed_files(bed_files, work_dir, genome):
    from clearup.panel import overlap_bed_files

    fnames = [basename(splitext_plus(fp)[0]) for fp in bed_files]
    overlapped_file = join(work_dir, f'{"__".join(fnames)}.{genome}.bed')
    if not can_reuse(overlapped_file, bed_files):
        overlap_bed_files(bed_files, overlapped_file)
    return overlapped_file


def _slice_vcf_fn(work_dir, label, vcf_file, overlapped_bed):
    sliced_vcf_file = join(work_dir, label + '.sliced.vcf')
    if not can_reuse(sliced_vcf_file, [vcf_file]):
        run(f'bcftools view {vcf_file} --targets-file {overlapped_bed} -o {sliced_vcf_file}')

    # ann_vcf_file = join(work_dir, label + '.sliced.ann.vcf')
    # if not can_reuse(ann_vcf_file, [sliced_vcf_file]):
    #     vcf_header = join(work_dir, label + '.vcf_header')
    #     with open(vcf_header, 'w') as f:
    #         f.write('##INFO=<ID=CHROM,Number=1,Type=String,Description="Region chromosome">\n')
    #         f.write('##INFO=<ID=FROM,Number=1,Type=String,Description="Region start">\n')
    #         f.write('##INFO=<ID=TO,Number=1,Type=String,Description="Region end">\n')
    #     run(f'bcftools annotate -c CHROM,FROM,TO -a {overlapped_bed} {sliced_vcf_file} '
    #         f'-h {vcf_header} -o {ann_vcf_file}')

    return label, sliced_vcf_file


class Dataset:
    def __init__(self):
        self.name = None
        self.genome = None
        self.bed_file = None
        self.vcf_by_label = dict()


def _load_datasets(subdirs):
    vcf_by_project_by_genome = defaultdict(dict)
    # vcf_by_label = dict()
    # all_bed_files = []
    # project_names = []
    datasets = []

    for subdir in subdirs:
        dataset = Dataset()

        if ':' in subdir:
            subdir, dataset.genome = subdir.split(':')
        else:
            dataset.genome = 'hg19'

        dir_path = subdir
        if glob(join(dir_path, '*.vcf.gz')):
            log.info(f'Found .vcf.gz files in directory {dir_path}')
            # Simple directory with VCF files and an optional BED file?
            dataset.name = subdir.replace('/', '__')
            if glob(join(dir_path, '*.bed')):
                dataset.bed_file = glob(join(dir_path, '*.bed'))[0]
            for vcf_fpath in glob(join(dir_path, '*.vcf.gz')):
                label = join(subdir, basename(splitext_plus(vcf_fpath)[0])).replace('/', '__')
                dataset.vcf_by_label[label] = vcf_fpath
        else:
            log.info(f'Not found any .vcf.gz files in directory {dir_path}. Checking if that\'s a bcbio folder.')
            # Bcbio directory?
            bcbio_proj = BcbioProject()
            bcbio_proj.load_from_bcbio_dir(subdir, proc_name='clearup')
            dataset.name = bcbio_proj.project_name
            dataset.genome = bcbio_proj.genome_build
            for s in bcbio_proj.samples:
                vcf_file = s.find_raw_vcf()
                if vcf_file:
                    dataset.vcf_by_label[bcbio_proj.project_name + '__' + s.name] = vcf_file
            if bcbio_proj.coverage_bed:
                dataset.bed_file = bcbio_proj.coverage_bed

        datasets.append(dataset)
    return datasets


possible_keys = [''.join(v[0]+v[1]) for v in it.product(it.permutations('ACGT', 2), repeat=2)]
index_by_key = {sys.intern(k): possible_keys.index(k) for k in possible_keys}


def make_fingerprint(vcf_file, work_dir=None, label=None, fp_size=20, bed_file=None):
    log.info('Starting processing file ' + vcf_file)
    work_dir = work_dir or dirname(vcf_file)

    if label: print_name = label
    else:     print_name = splitext_plus(basename(vcf_file))[0]
    print_name += '.print' + str(fp_size)
    print_name += '_dist' + str(Params.MIN_DIST)
    print_name += '_af' + str(Params.MIN_AF)
    if not Params.INTERREGION_PAIRS:
        print_name += '_skip_interregion_pairs'

    raw_print_file = join(work_dir, print_name)
    if can_reuse(raw_print_file, vcf_file):
        with open(raw_print_file) as f:
            raw = np.fromfile(f).reshape((len(index_by_key), fp_size))
    else:
        raw = _raw_fingerprint(vcf_file, fp_size=fp_size, bed_file=bed_file)
        with open(raw_print_file, 'w') as f:
            raw.tofile(f)
        log.info(f'Saved raw fingerprints into {raw_print_file}')

    norm_print_name = print_name
    if Params.NORMALIZE_DIST: norm_print_name += '_normdist'
    if Params.NORMALIZE_VAR:  norm_print_name += '_normvar'

    norm_print_file = join(work_dir, norm_print_name)
    if can_reuse(norm_print_file, raw_print_file):
        with open(norm_print_file) as f:
            norm = np.fromfile(f).reshape((len(index_by_key), fp_size))
    else:
        norm = _normalize_fingerprint(raw)
        with open(norm_print_file, 'w') as f:
            norm.tofile(f)
        log.info(f'Saved normalised fingerprints into {norm_print_file}')

    return label, norm


def _raw_fingerprint(vcf_file, fp_size=20, bed_file=None):
    raw = np.zeros((len(index_by_key), fp_size))
    prev_rec = None

    total_vars = 0
    rejected_af = 0
    rejected_other = 0
    passed = 0
    pairs = 0

    regions = BedTool(bed_file).features() if bed_file else None
    if regions:
        region = next(regions)
    def in_region(_rec, _reg):
        return _rec.CHROM == _reg.chrom and _reg.start < _rec.POS <= _reg.stop

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
            continue  # multiple ALTs
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

        if prev_rec is None:
            prev_rec = rec
            continue
        # looping into a new chromosome?
        if prev_rec.CHROM != rec.CHROM:
            prev_rec = rec
            continue  # looping into a new chromosome
        # looping into a new region?
        if regions and Params.INTERREGION_PAIRS:
            assert region
            while not in_region(rec, region):
                try:
                    region = next(regions)  # finding the region for current variant
                except:
                    log.critical(f'Rec {rec.CHROM}:{rec.POS} outside of regions')
            if not in_region(prev_rec, region):  # if the previos is not at the same region, continue
                prev_rec = rec
                continue  # looping into a new region

        pair_dist = rec.POS - prev_rec.POS
        if pair_dist < Params.MIN_DIST:  # minimal allowed distance between variants
            continue  # too closely located variants
        log.debug(f'{prev_rec.CHROM}:{prev_rec.POS} {prev_rec.REF}>{prev_rec.ALT[0]} '
                  f'- {rec.CHROM}:{rec.POS} {rec.REF}>{rec.ALT[0]} - {pair_dist}')
        pair_key = sys.intern(prev_rec.REF + prev_rec.ALT[0] + rec.REF + rec.ALT[0])
        raw[index_by_key[pair_key], pair_dist % fp_size] += 1
        prev_rec = rec
        pairs += 1

    log.info(f'total_vars: {total_vars}'
             f', rejected_af: {rejected_af}'
             f', rejected_other: {rejected_other}'
             f', passed: {passed}'
             f', pairs: {pairs}')
    log.info()

    return raw


def _normalize_fingerprint(arr):
    if Params.NORMALIZE_DIST:
        arr = arr - np.mean(arr, axis=0)
        arr = arr / np.std(arr, axis=0)

    if Params.NORMALIZE_VAR:
        arr = arr - np.mean(arr, axis=1)[:, np.newaxis]
        arr = arr / np.std(arr, axis=1)[:, np.newaxis]

    return arr


def plot_heatmap(pairwise, run_dir, title):
    df = pd.DataFrame(data=pairwise)
    log.info(df)

    # Generate a mask for the upper triangle + main diagonale
    mask = np.zeros_like(df, dtype=np.bool)
    mask[np.triu_indices_from(mask, k=1)] = True

    # Set up the matplotlib figure
    n = len(pairwise)
    figsize = (n / 2, n * 7 / 20)
    log.info(f'Saving figure of size {figsize}')
    f, ax = plt.subplots(figsize=figsize)  # For 20 samples, take 10x7
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
    matplotlib.pyplot.subplots_adjust(left=0.2, right=1, top=0.93, bottom=0.29)

    png_file = join(run_dir, str_to_filename(title) + '.png')
    if isfile(png_file):
        os.remove(png_file)
    matplotlib.pyplot.savefig(png_file)
    if isfile(png_file):
        log.info('')
        log.info('Saved heatmap into ' + adjust_path(png_file))
        try:
            from az.webserver.exposing import convert_gpfs_path_to_url
        except ImportError:
            pass
        else:
            url = convert_gpfs_path_to_url(png_file)
            if url:
                log.info('    url: ' + url)
                return url
        return png_file


def _hash(s, l=10):
    return str(int(hashlib.sha256(s.encode('utf-8')).hexdigest(), 16) % 10 ** l)


if __name__ == '__main__':
    main()
