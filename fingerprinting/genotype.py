from itertools import izip

from collections import defaultdict

from ngs_utils.sambamba import index_bam

from ngs_utils.call_process import run
from ngs_utils.utils import is_local, is_us
from os.path import join

from ngs_utils.parallel import ParallelCfg, parallel_view
from ngs_utils.file_utils import file_transaction, safe_mkdir, chdir, which, adjust_path
from ngs_utils.logger import info, err, critical


def run_vardict(proj, work_dir, snp_file, parallel_cfg, genome_cfg):
    with parallel_view(len(proj.samples), parallel_cfg, safe_mkdir(work_dir)) as view:
        view.run(_vardict_one_sample, [[s, work_dir, genome_cfg, view.cores_per_job, snp_file] for s in proj.samples])


def _vardict_one_sample(sample, work_dir, genome_cfg, threads, snp_file):
    if is_local():
        vardict_pl = '/Users/vlad/vagrant/VarDict/vardict.pl'
    elif is_us():
        vardict_pl = '/group/cancer_informatics/tools_resources/NGS/bin/vardict.pl'
    else:
        vardict_pl = which('vardict.pl')
        if not vardict_pl:
            critical('Error: vardict.pl is not in PATH')

    index_bam(sample.bam)

    ref_file = adjust_path(genome_cfg['seq'])
    sample_work_dir = safe_mkdir(join(work_dir, sample.name))
    output_file = join(sample_work_dir, 'vardict_snp_vars.txt')
    cmdl = '{vardict_pl} -G {ref_file} -c 1 -S 2 -E 2 -N {sample.name} -b {sample.bam} -p -D {snp_file}'.format(**locals())
    run(cmdl, output_fpath=output_file)
    return output_file


def fingerprint(bcbio_dir, sample_name):  # the purpose of this is to rip through the vcf and create dictionaries
                                    # keyed on sample ID and containing chrid keyed ref,alt,AF,DP information
    data_by_sample = defaultdict(dict())
    with open('%s/work/%s/fingerprintvcf.txt' % (bcbio_dir, sample_name)) as fh:
        for line in fh:
            line2 = line.split('\t')
            if "SNV" in line:
                chrid = line2[0] + "_" + line2[1]
                ref = line2[2]
                alt = line2[3]
                for item in line2:
                    if " AF=" in item:  # figure out the sample names here
                        samplename = item.split(" AF")[0]
                        info(samplename + " sample name from AF")
                        af = item.split('=')[1]
                        if chrid not in data_by_sample[samplename]:
                            data_by_sample[samplename][chrid] = [ref, alt, af, 0]
                        else:
                            info("DOUBLE IDED CHRID IN SYSTEM", line, data_by_sample[samplename][chrid])
                    if " DP=" in item:
                        samplename = item.split(" DP")[0]
                        dp = item.split('=')[1]
                        data_by_sample[samplename][chrid][3] = dp
    return data_by_sample


def fingerwriter(reffile, bcbio_dir, vcfdict, depthcut, refdict):  # this takes in a dictionary keyed on sample names each keyed on chrid for each identified snp
    for sample_name in vcfdict:  # for each sample go through the snps IDed #basedir is the key for samples
        info(sample_name)
        depthdict = {}
        ismale = False
        ytotal = 0
        with open('%s/work/%s/fingerprintdepth.txt' % (bcbio_dir, sample_name)) as fh:  # should open the depth and we shoudl have these because of bigwigs for each T/N pair
            for line in fh:
                line2 = line.split('\t')
                chrid = line2[0] + "_" + line2[2]
                if not depthdict.has_key(chrid):
                    depthdict[chrid] = (line2[-3], line2[-1].strip())  # identifier and depth
                else:
                    err("FAIL DEPTH!!!!!! " + line + ' ' + str(depthdict[chrid]))
                if line2[0] == "chrY":
                    ytotal += int(line2[-1].strip())
        if ytotal > depthcut:
            ismale = True
        
        with open(reffile) as fh, \
             open("%s/final/%s/%s-Fingerprint.txt" % (bcbio_dir, sample_name, sample_name), "w") as fhw, \
             open("%s/work/%s/%s-printcoverage.log" % (bcbio_dir, sample_name, sample_name), "w") as fhw2, \
             open("%s/work/%s/%s-error.log" % (bcbio_dir, sample_name, sample_name), "w") as fherror:
            outstr = ""
            info(sample_name + ": processing...")
            rlc = 0  # refline count
            woc = 0  # write out count
            for line in fh:
                rlc += 1
                line2 = line.split('\t')
                chrid = line2[0] + "_" + line2[1]
                snpid = line2[3]
                if chrid in depthdict:
                    if int(depthdict[chrid][1]) > depthcut:  # check if it's below the set depth.
                        if chrid in vcfdict[sample_name]:  # see if we have a variant called
                            if line2[0] == "chrX" and ismale and 0.25 <= float(vcfdict[sample_name][chrid][2]):
                                woc += 1
                                outstr += vcfdict[sample_name][chrid][1] + "\t"
                            elif line2[0] == 'chrY' and 0.25 <= float(vcfdict[sample_name][chrid][2]):
                                woc += 1
                                outstr += vcfdict[sample_name][chrid][1] + "\t"
                            else:  # is not a sex chromosome or is chrX and needs to be two alleles
                                if 0.25 <= float(vcfdict[sample_name][chrid][2]) < 0.75:  # heterozygous
                                    woc += 1
                                    outstr += vcfdict[sample_name][chrid][0] + vcfdict[sample_name][chrid][1] + "\t"
                                elif 0.75 <= float(vcfdict[sample_name][chrid][2]):  # homozygous ALT
                                    woc += 1
                                    outstr += vcfdict[sample_name][chrid][1] + vcfdict[sample_name][chrid][1] + "\t"
                                else:  # lower than 25% so we'll call this homozygous reference
                                    fhw2.write("Variant Called but below 25%%.\t%s\t%s\n" % (
                                    vcfdict[sample_name][chrid], refdict[chrid]))
                                    woc += 1
                                    outstr += vcfdict[sample_name][chrid][0] + vcfdict[sample_name][chrid][0] + "\t"
                        else:  # no variant called so we just report reference if we have sufficient depth. Homozygous reference
                            if line2[0] == 'chrY':
                                woc += 1
                                outstr += line2[3].strip() + "\t"
                            else:  # not a chrY
                                woc += 1
                                outstr += line2[3].strip() + line2[3].strip() + "\t"
                    else:  # lower depth for now print out NN because of low depth
                        if vcfdict[sample_name].has_key(chrid):  # see if we have a variant called
                            fhw2.write("variant called but LOW depth (<%s): %s\n" % (depthcut, line))
                        fhw2.write("low depth\t%s\t%s\n" % ('\t'.join(depthdict[chrid]), line))
                        if line2[0] == "chrY":
                            woc += 1
                            outstr += "N" + "\t"
                        else:  # Not a Y chr
                            if line2[0] == 'chrX' and ismale:  # is there a Y somewhere?
                                woc += 1
                                outstr += "N" + "\t"
                            else:
                                woc += 1
                                outstr += "NN" + "\t"
                else:  # check for variant call but without depth this will require some deep dive
                    if vcfdict[sample_name].has_key(chrid):
                        fhw2.write("variant called but without depth: %s\n" % (line))
                    else:  # no variant was called adn we lack depth
                        fhw2.write("NO depth detected and no variant (likely chrY): %s\n" % (line))
                    if line2[0] == "chrY":
                        woc += 1
                        outstr += "N" + "\t"
                    else:  # X or non-sex chromosome
                        if line2[0] == 'chrX' and ismale:
                            woc += 1
                            outstr += "N" + "\t"
                        else:
                            woc += 1
                            outstr += "NN" + "\t"
                if not woc == rlc:
                    # the odd case
                    if woc > rlc:
                        fherror.write("woc is higher than the ref line, %s\n" % (line))
                        while rlc < woc:
                            rlc += 1
                    elif woc < rlc:  # more common
                        fherror.write("refline didn't get analyzed, %s\n" % line)
                        while woc < rlc:
                            woc += 1
            fhw.write(outstr + '\n')


def paircompare(samples, fhw2, genelist):
    # This is the tool to compare a T/N vcf and provide the statistics back for the pair identified takes
    # in the dictionary of dictoionaries keyed on sample name (basedirs)
    # initially the length should only be two. so I'm writing this for that (although I'm going to iterate so that
    # the first is the main tumor sample and the rest are normals in this case)
    normals = {}
    
    tumor = join('final', samples.keys()[0], samples.keys()[0] + '-Fingerprint.txt')
    tumorbasedir = samples.keys()[0]
    for samplekey in samples.keys()[1:]:
        # commented out assuming no duplication here?
        # if not normals.has_key(samplekey):
        # 	normals[samplekey]=[]
        normals[samplekey] = join('final', samplekey, samplekey + '-Fingerprint.txt')
    
    for normalsample in normals:
        fh = open(tumor)
        a = fh.readline()
        fh.close()
        info(str(normals))
        info(str(normalsample))
        with open(normals[normalsample]) as fh:
            b = fh.readline()
        
        if a == b:
            info("Genotype matches 100 Percent between: %s and %s" % (tumorbasedir, normalsample))
            info("%s\n%s" % (tumorbasedir, '\t'.join(a.strip().split())))
            info("%s\n%s" % (normalsample, '\t'.join(b.strip().split())))
            fhw2.write("Genotype matches 100 Percent between: %s and %s\n" % (tumorbasedir, normalsample))
            fhw2.write("%s\n%s\n" % (tumorbasedir, '\t'.join(a.strip().split())))
            fhw2.write("%s\n%s\n" % (normalsample, '\t'.join(b.strip().split())))
        
        else:  # if not an exact match print out metrics of how close we are
            a2 = a.strip().split()
            b2 = b.strip().split()
            count = 0  # init matches
            totalcount = 0
            for c1, c2 in izip(a2, b2):
                if c1 == c2:
                    count += 1
                    totalcount += 1
                else:
                    totalcount += 1
            info("NO MATCH between: %s and %s" % (tumorbasedir, normalsample))
            info('\t'.join(genelist))
            info("%s\n%s" % (tumorbasedir, '\t'.join(a.strip().split())))
            info("%s\n%s" % (normalsample, '\t'.join(b.strip().split())))
            info("However we matched %s Percent of the reads (%s out of %s total)" %
                 ((float(count) / float(totalcount)), count, totalcount))
            info("Genes where we mismatch: " + ",".join(
                [genelist[i] for i, (c1, c2) in enumerate(izip(a2, b2)) if c1 != c2]))
            fhw2.write("NO MATCH between: %s and %s\n" % (tumorbasedir, normalsample))
            fhw2.write('\t'.join(genelist) + '\n')
            fhw2.write("%s\n%s\n" % (tumorbasedir, '\t'.join(a.strip().split())))
            fhw2.write("%s\n%s\n" % (normalsample, '\t'.join(b.strip().split())))
            fhw2.write("However we matched %s Percent of the reads (%s out of %s total)\n" % (
            (float(count) / float(totalcount)), count, totalcount))
            fhw2.write("Genes where we mismatch: " + ",".join(
                [genelist[i] for i, (c1, c2) in enumerate(izip(a2, b2)) if c1 != c2]) + '\n')
            info("=============================================================")
            fhw2.write("=============================================================\n\n")
