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
from ngs_utils.parallel import parallel_view, ParallelCfg

from clearup.utils import bam_samplename


def sample_callable_bed(bam_file, output_bed_file, work_dir, genome_fasta_file, min_depth):
    """Retrieve callable regions for a sample subset by defined analysis regions.
    """
    callable_bed = _calculate(bam_file, work_dir, genome_fasta_file, min_depth)
    if not can_reuse(output_bed_file, callable_bed):
        with file_transaction(work_dir, output_bed_file) as tx_out_file:
            callable_regions = pybedtools.BedTool(callable_bed).filter(lambda x: x.name == 'CALLABLE')
            callable_regions.saveas(tx_out_file)
    return output_bed_file


def batch_callable_bed(bam_files, output_bed_file, work_dir, genome_fasta_file, min_depth,
                       parall_view=None):
    """ Picking random 3 samples and getting a callable for them.
        Trade off between looping through all samples in a huge batch,
        and hitting an sample with outstanding coverage.
    """
    if can_reuse(output_bed_file, bam_files):
        return output_bed_file

    work_dir = safe_mkdir(join(work_dir, 'callable_work'))
    # random.seed(1234)  # seeding random for reproducability
    # bam_files = random.sample(bam_files, min(len(bam_files), 3))

    if parall_view:
        callable_beds = parall_view.run(_calculate, [
            [bf, work_dir, genome_fasta_file, min_depth]
            for bf in bam_files])
    else:
        with parallel_view(len(bam_files), ParallelCfg(threads=len(bam_files)), work_dir) as parall_view:
            callable_beds = parall_view.run(_calculate, [
                [bf, work_dir, genome_fasta_file, min_depth]
                for bf in bam_files])

    good_overlap_sample_fraction = 0.8  # we want to pick those regions that have coverage at 80% of samples
    good_overlap_count = min(1, good_overlap_sample_fraction * len(callable_beds))
    info(f'Intersecting callable regions and picking good overlaps with >={good_overlap_count} ' +
         f'samples ({100 * good_overlap_sample_fraction}% of {len(callable_beds)})')
    with file_transaction(work_dir, output_bed_file) as tx:
        pybedtools.set_tempdir(safe_mkdir(join(work_dir, 'pybedtools_tmp')))
        intersection = pybedtools.BedTool() \
            .multi_intersect(i=callable_beds) \
            .filter(lambda r: len(r[4].split(',')) >= good_overlap_count)
        intersection.saveas(tx)
    info(f'Saved to {output_bed_file}')
    return output_bed_file


def _calculate(bam_file, work_dir, genome_fasta_file, min_depth):
    """Calculate coverage in parallel using samtools depth through goleft.

    samtools depth removes duplicates and secondary reads from the counts:
    if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
    """
    output_prefix = os.path.join(work_dir, bam_samplename(bam_file))
    callability_annotation_file = output_prefix + '.callable.bed'
    callable_file = output_prefix + '.callable.CALLABLE.bed'
    if can_reuse(callable_file, bam_file):
        return callable_file

    info(f'Calculating coverage at {bam_file}')
    run(f'goleft depth --q 1 --mincov {min_depth} --reference {genome_fasta_file} --ordered'
        f' --prefix {output_prefix} {bam_file}')

    with file_transaction(None, callable_file) as tx:
        pybedtools.BedTool(callability_annotation_file)\
            .filter(lambda x: x.name == 'CALLABLE')\
            .saveas(tx)

    return callable_file


