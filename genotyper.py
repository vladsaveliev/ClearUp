#!/usr/bin/env python

import os
import sys
import glob
import click
from os.path import join, dirname, basename
from variant_filtering.vcf import bgzip_and_tabix

from ngs_utils import logger
from ngs_utils.file_utils import safe_mkdir, can_reuse, verify_file
from ngs_utils.logger import err, info, debug, critical
from ngs_utils.call_process import run
from ngs_utils.parallel import ParallelCfg

from ngs_reporting.bcbio.bcbio import BcbioProject

from fingerprinting.config import get_version
from fingerprinting.genotype import fingerprint, fingerwriter, paircompare, run_vardict

import az


""" Rapidly genotype the samples after BCBio runs.
    Make sure you have bcftools loaded"
"""

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
# @click.option('-v', '--vcf',
#               type=click.Path(exists=True),
#               help='Sorted VCF file only for sites for the snps.',
#               required=True,
#               )
@click.option('-D', '--depth',
              type=int,
              help='Minimum coverage depth for calls.',
              default=5,
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
    """ Step one gather up the needed files and set up the data structures."""
    snp_file = verify_file(bed)
    depth_cutoff = depth

    logger.init(isdebug)

    info('Loading bcbio project from ' + bcbio_dir)
    info('-' * 70)
    proj = BcbioProject()
    proj.load_from_bcbio_dir(bcbio_dir, proc_name='fingerprinting', need_coverage_interval=False)
    safe_mkdir(proj.work_dir)
    info('Loaded ' + proj.final_dir)
    bam_files = [s.bam for s in proj.samples]
    debug('Found BAM files: ' + str(bam_files))
    
    bigwig_file = set()
    genelist = []

    sys_cfg = az.init_sys_cfg()
    if threads:
        sys_cfg['threads'] = threads
    parallel_cfg = ParallelCfg(sys_cfg.get('scheduler'), sys_cfg.get('queue'),
                               sys_cfg.get('resources'), sys_cfg.get('threads'))
    parallel_cfg.set_tag('fingerprinting')
    genome_cfg = az.get_refdata(proj.genome_build)

    info('** Running VarDict ** ')
    run_vardict(proj, proj.work_dir, snp_file, parallel_cfg, genome_cfg)
    info('** Finished running VarDict **')
    
    sys.exit(1)
    
    # NEW fingerprintvcf creation #################
    # Step 2) quickly run vardict and create a new vcf for fingerprinting, then use bcftools to query => filter out the DP and AF (needed for tumor normal samples)
    # Written to each dir/bcbio/final/subdir/ as idtfingervar.txt, idtfingervar.vcf.gz, and finally to fingerprintvcf.txt
    for sample in proj.samples:
        if not sample.bam:
            critical('BAM file not found for sample ' + sample.nam)
        fingerprintvcf_txt = join(safe_mkdir(join(proj.work_dir, sample.name)), 'fingerprintvcf.txt')
        if not can_reuse(fingerprintvcf_txt, cmp_f=sample.bam):
            info('Use vardict.pl to quickly scan the ready.bam for files')
            vardict_out = join(proj.work_dir, sample.name, 'idtfingervar.txt')
            run('vardict.pl -b {sample.bam} -N {sample.name} -f 0.005 {snp_bed_file} > {vardict_out}'.format(**locals()))
            vardict_vcf = join(proj.work_dir, sample.name, 'idtfingervar.vcf')
            run('var2vcf_valid.pl {vardict_out} > {vardict_vcf}'.format(**locals()))
            vardict_vcf_gz = bgzip_and_tabix(vardict_vcf)
            run('bcftools query -H -R {vcf_file} -f "%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE AF=%%AF][\\t%%SAMPLE DP=%%DP]\\t%%ID\\t%%LINE"'
                ' {vardict_vcf_gz} > {fingerprintvcf_txt}'.format(**locals()))

    info('Running sambamba to get depth at the positions for all samples simultaneously output to dir/bcbio/final/sambambadepth.txt')
    fingerpints_dir = safe_mkdir(join(proj.date_dir, 'fingerprints'))
    depth_txt = join(fingerpints_dir, 'sambambadepth.txt')
    if not can_reuse(depth_txt, cmp_f=bam_files):
        bam_lines = ' '.join(bam_files)
        run('sambamba depth region -L {snp_bed_file} -o {depth_txt} {bam_lines}'.format(**locals()))
    
    # Step 4) parse this info into dictionaries so taht you can check if it has a variatn called and what the depth is.
    with open(depth_txt) as fh:
        fh.readline()  # gather infromatino and indexing from here if need be later
        bwdict = {}  # key this on sample names
        for line in fh:
            line2 = line.strip().split()
            chrom = line2[0]
            start = int(line2[1])
            end = int(line2[2])
            sample_name = line2[-1]
            readcount = int(line2[-3])
            avgcount = float(line2[-2])
            if sample_name not in bigwig_file:  # add it and create a dictionary of dictionaries to contain the write outs at the end for depth file
                bigwig_file.add(sample_name)
                bwdict[sample_name] = {}  # key this on chrid
            if max(int(start), int(end)) - min(int(start), int(end)) > 1:
                err(line + " WAS LONGER IN DEPTH")  # just identify the longer than a SNP positions (make this an actual error later or save for gender calling)
            bwdict[sample_name][chrom + "_" + str(start)] = readcount
    
    # Now parse the dictionary according to the snp file (later delete this step and just output the sambamba directly)
    for sample_name in bwdict:
        fingerprint_depth_file = join(proj.work_dir, sample_name, 'fingerprintdepth.txt')
        if not can_reuse(fingerprint_depth_file, snp_file):
            with open(snp_file) as fh, open(fingerprint_depth_file, 'w') as fhw:
                ytotal = 0
                for line in fh:
                    line2 = line.split('\t')
                    chrom = line2[0]
                    start = int(line2[1])
                    fhw.write(line.strip()+'\t1\t'+str(bwdict[sample_name][chrom + "_" + str(start)]) + '\n')
                    # this has a 1 in it because bedtools puts the coverage fraction in the depth output so I'm matching in... going forward can remove and decrease the position later in teh code
                    if chrom == "chrY":
                        ytotal += bwdict[sample_name][chrom + "_" + str(start)]
    
    # Ok now parse through the checklist for each base directory that was analyzed and kick out the final files. (parallelize these (defs for each of the above) pool and set threads to 6(?) Do I need this here? maybe for the checking? single pipeline that? It's here now so I can print out the order and make sure things are going through correctly for the output.
    refdict = {}  # to process the checklist
    with open(snp_file) as fh:
        for line in fh:
            line2 = line.split('\t')
            chrid = line2[0] + "_" + line2[2]
            name = line2[-1].strip()
            refdict[chrid] = name
            genelist.append(name)
    
    info('Printing any files missing either the bigwig or vardict VCF')
    if bigwig_file ^ vcf_file:
        for bdir in bigwig_file ^ vcf_file:
            err("lacking depth or vardict vcf in: %s" % bdir)
    
    info(str(bigwig_file), str(vcf_file))
    finaldatefolder = glob.glob(join(bcbio_dir, 'final', '*_bcbio'))
    with open(join(finaldatefolder[0], "Paired_FingerPrint_Comparison.txt"), 'w') as fhw2:
        info('Running through things normally and check for inclusion of both bigwig and VCF (tumor normal will have all in one VCF')
        for eachdir in bigwig_file & vcf_file:
            info('Sending to the def and let it return a dictionary with keys')
            info(eachdir + " eachdir")
            data_by_sample = fingerprint(bcbio_dir, eachdir)  # this is a dictionary of dictionaries keyed on the sample name and then the chrid
            
            info('Writing files for samples')
            fingerwriter(vcf_file, bcbio_dir, data_by_sample, depth_cutoff, refdict)  # this will write out the files
        
            # if paired then we want to make the immediate comparison. Else we want to run comparison later so def the comparison stuff into this.
            if len(data_by_sample) > 1:  # we have a paired analysis.
                if len(data_by_sample) > 2:
                    info("Have more than pairs need to figure things out on pairing")
                    # although if the first is the main sample the rest can be iterated through (Go back and take out of genotypecompareV2 to do this)
                info("Paired sample")
                if fhw2:
                    paircompare(data_by_sample, fhw2, genelist)
                # do something to tag the two files in the folder or
            
    # new 7/11/16
    # now generate the comprehensive blast dictionary
    # takes the -Fingerprint.txt output and creates a fasta
    fingerprint_texts = glob.glob(join(proj.final_dir, '*', '*-Fingerprint.txt'))
    longprints_fasta = join(fingerpints_dir, 'longprints.fasta')
    if not can_reuse(longprints_fasta):
        with open(longprints_fasta, 'w') as blastfasta:
            for fpfile in fingerprint_texts:
                sample_name = fpfile.split('/')[1]
                blastfasta.write(">" + sample_name + "\n")
                with open(fpfile) as fh:
                    seqline = fh.readline()
                    blastfasta.write(seqline.replace("\t", ''))
                with open('final/%s/%s-longprint.fasta' % (sample_name, sample_name), 'w') as fhw:
                    fhw.write(">%s\n" % sample_name)
                    fhw.write(seqline.replace("\t",''))


if __name__ == '__main__':
    main()
