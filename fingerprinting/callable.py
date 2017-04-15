import random
import tempfile
import os
import shutil
from os.path import basename, splitext, join, dirname
import pybedtools
from contextlib import contextmanager

from ngs_utils.call_process import run
from ngs_utils.logger import info, debug, err, warn, critical
from ngs_utils.file_utils import tx_tmpdir, can_reuse, file_transaction, safe_mkdir, adjust_path, chdir, splitext_plus

from fingerprinting.utils import bam_samplename


def sample_callable_bed(bam_file, output_bed_file, work_dir, genome_cfg, min_depth):
    """Retrieve callable regions for a sample subset by defined analysis regions.
    """
    with bedtools_tmpdir(work_dir):
        callable_bed = _calculate(bam_file, work_dir, genome_cfg, min_depth)
        if not can_reuse(output_bed_file, callable_bed):
            with file_transaction(work_dir, output_bed_file) as tx_out_file:
                callable_regions = pybedtools.BedTool(callable_bed).filter(lambda x: x.name == 'CALLABLE')
                callable_regions.saveas(tx_out_file)
    return output_bed_file


def batch_callable_bed(bam_files, output_bed_file, work_dir, genome_cfg, min_depth, parall_view):
    """ Picking random 3 samples and getting a callable for them.
        Trade off between looping through all samples in a huge batch,
        and hitting an sample with outstanding coverage.
    """
    if can_reuse(output_bed_file, bam_files):
        return output_bed_file
        
    random.seed(1234)  # seeding random for reproducability
    bam_files = random.sample(bam_files, min(len(bam_files), 3))

    callable_beds = parall_view.run(_calculate, [[bam_file, work_dir, genome_cfg, min_depth]
         for bam_file in bam_files])

    with bedtools_tmpdir(work_dir):
        with file_transaction(work_dir, output_bed_file) as tx:
            pybedtools.BedTool(callable_beds[0])\
                .cat(*callable_beds[1:])\
                .filter(lambda x: x.name == 'CALLABLE')\
                .merge()\
                .saveas(tx)
    return output_bed_file
    

@contextmanager
def bedtools_tmpdir(work_dir):
    with tx_tmpdir(work_dir) as tmpdir:
        orig_tmpdir = tempfile.gettempdir()
        pybedtools.set_tempdir(tmpdir)
        yield
        if orig_tmpdir and os.path.exists(orig_tmpdir):
            pybedtools.set_tempdir(orig_tmpdir)
        else:
            tempfile.tempdir = None


def _calculate(bam_file, work_dir, genome_cfg, min_depth):
    """Calculate coverage in parallel using samtools depth through goleft.

    samtools depth removes duplicates and secondary reads from the counts:
    if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
    """
    params = {'window_size': 5000, 'parallel_window_size': 1e5, 'high_multiplier': 20}
    prefix = os.path.join(work_dir, bam_samplename(bam_file) + '-coverage')
    depth_file = prefix + '.depth.bed'
    callable_file = prefix + '.callable.bed'
    if can_reuse(callable_file, bam_file):
        return callable_file
        
    ref_file = adjust_path(genome_cfg['seq'])
    cmdl = 'goleft depth --q 1 --mincov {min_depth} --reference {ref_file} --ordered'
    with file_transaction(work_dir, depth_file) as tx_depth_file:
        with chdir(dirname(tx_depth_file)):
            tx_callable_file = tx_depth_file.replace('.depth.bed', '.callable.bed')
            prefix = tx_depth_file.replace('.depth.bed', '')
            cmdl += '--prefix {prefix} {bam_file}'
            info('Calculating coverage at ' + bam_file)
            run(cmdl.format(**locals()))
            shutil.move(tx_callable_file, callable_file)

    return callable_file


