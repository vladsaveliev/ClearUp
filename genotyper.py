#!/usr/bin/env python

import os
import click

from ngs_utils import logger as log
from ngs_utils.file_utils import safe_mkdir, can_reuse, verify_file
from ngs_utils.parallel import ParallelCfg

from ngs_reporting.bcbio.bcbio import BcbioProject

from fingerprinting.config import get_version
from fingerprinting.genotype import genotype_bcbio_proj, DEPTH_CUTOFF

import az


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('bcbio_dir',
                type=click.Path(exists=True, dir_okay=True),
                default=os.getcwd(),
                metavar="<bcbio directory>"
                )
@click.option('-b', '--bed',
              type=click.Path(exists=True),
              help='Sorted BED file for all the snps, in the following format.',
              required=True,
              )
@click.option('-D', '--depth',
              type=int,
              help='Minimum coverage depth for calls.',
              default=DEPTH_CUTOFF,
              )
@click.option('-t', '--threads',
              type=int,
              help='Number of threads. Default is in corresponding system_info_*.yaml or 1. '
                   'If set to 1, skip starting cluster even if scheduler is specified.',
              )
@click.option('-d', '--isdebug',
              is_flag=True
              )
@click.version_option(version=get_version())
def main(bcbio_dir, bed, depth, threads=None, isdebug=True):
    snp_file = verify_file(bed)
    depth_cutoff = depth

    log.init(isdebug)

    sys_cfg = az.init_sys_cfg()
    parallel_cfg = ParallelCfg(
        scheduler=sys_cfg.get('scheduler'),
        queue=sys_cfg.get('queue'),
        resources=sys_cfg.get('resources'),
        threads=threads or sys_cfg.get('threads'),
        tag='fingerprinting')

    log.info('Loading bcbio project from ' + bcbio_dir)
    log.info('-' * 70)
    proj = BcbioProject()
    proj.load_from_bcbio_dir(bcbio_dir, proc_name='fingerprinting', need_coverage_interval=False)
    log.info('Loaded ' + proj.final_dir)
    genotype_bcbio_proj(proj, snp_file, parallel_cfg, depth_cutoff)


if __name__ == '__main__':
    main()
