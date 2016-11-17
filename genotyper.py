#! python

#Author Tristan Lubinski
#Date 3-15-16
#purpose to rapidly genotype the samples after BCBio runs
#date 6/21/16
#modifications to speed up using bcftools instead of bedtools and now integrating the information from tumor normal/multiple samples when availble.

import sys, glob, subprocess, argparse, os
from bx.bbi.bigwig_file import BigWigFile
from itertools import izip


print "make sure you have bcftools loaded" 
print "usage: python genotyperswift.py sortedsnp.bed desiredSNPs.vcf ##fordepthcutoff<defaultto0please>"
print "now with argparse please use python genotyperswiftV12longprint.py -h for help"

parser = argparse.ArgumentParser(prog='genotyperswift', usage='%(prog)s --bed path/to/file.bed --vcf path/to/file.vcf --depth mindepth --denovo True/False')
parser.add_argument('-b', '--bed', help='sorted bedfile for all the snps', required=True)
parser.add_argument('-r', '--ref', help='sorted reference file for all the mutations', required=True)
parser.add_argument('-d', '--depth', default=5, type=int, help='minimum coverage depth for calls. Default = 5', required=True)
#long term project will scan the annotated vcfs for snps and create a new fingerprint across the experiement
#parser.add_argument('-D', '--denovo', action='store_true', help='Use if creating a new fingerprint (not the production print)')

def fingerprint(homedir, basedir): #the purpose of this is to rip through the vcf and create dictionaries keyed on sample ID and containing chrid keyed ref,alt,AF,DP information
	samplesdict = {}
	fh = open('%s/work/%s/fingerprintvcf.txt'%(homedir, basedir))
	print "these aren't annotated so they won't have the usual markers CODE IS MODIFIED"
	for line in fh:
		line2 = line.split('\t')
		#if "rs" in line2[6] and "TYPE=SNV" in line:
		if "SNV" in line:
			chrid = line2[0]+"_"+line2[1]
			ref = line2[2]
			alt = line2[3]
			vid = line2[6]
			for item in line2:
				if " AF=" in item: #figure out the sample names here
					samplename = item.split(" AF")[0]
					print samplename, "sample name from AF"
					if not samplesdict.has_key(samplename):
						samplesdict[samplename]={}
					AF = item.split('=')[1]
					if not samplesdict[samplename].has_key(chrid):
						samplesdict[samplename][chrid]=[ref,alt,AF,0]
					else:
						print "DOUBLE IDED CHRID IN SYSTEM", line, samplesdict[samplename][chrid]
				if " DP=" in item:
					samplename = item.split(" DP")[0]
					DP = item.split('=')[1]
					samplesdict[samplename][chrid][3]=DP
	fh.close()
	return samplesdict
		
def fingerwriter(homedir, vcfdict): #this takes in a dictionary keyed on sample names each keyed on chrid for each identified snp
	for basedir in vcfdict: #for each sample go through the snps IDed #basedir is the key for samples
		print basedir
		depthdict = {}
		ismale = False
		ytotal = 0
		fh = open('%s/work/%s/fingerprintdepth.txt'%(homedir, basedir)) #should open the depth and we shoudl have these because of bigwigs for each T/N pair
		for line in fh:
			line2 = line.split('\t')
			chrid = line2[0]+"_"+line2[2]
			#print line
			#print chrid, line2[4],line2[6]
			if not depthdict.has_key(chrid):
				#print line2, "TRISTAN NEW FIND"
				depthdict[chrid]=(line2[-3],line2[-1].strip()) #identifier and depth
			else:
				print "FAIL DEPTH!!!!!!", line, depthdict[chrid]
			if line2[0]=="chrY":
				ytotal += int(line2[-1].strip())
		fh.close()
		if ytotal > depthcut: 
			ismale = True
		
		fh =open(reffile)
		fhw=open("%s/final/%s/%s-Fingerprint.txt"%(homedir, basedir, basedir),"w")
		fhw2=open("%s/work/%s/%s-printcoverage.log"%(homedir, basedir, basedir),"w")
		fherror=open("$s/work/%s/%s-error.log"%(homedir, basedir, basedir),"w")
		outstr = ""
		print basedir+" Processing..."
		rlc = 0 #refline count
		woc = 0 #write out count
		for line in fh:
			rlc+=1
			line2 = line.split('\t')
			chrid = line2[0]+"_"+line2[1]
			snpid = line2[3]
			#print chrid, refdict[chrid]
			#print snpid
			#outstr +=refdict[chrid]
			if depthdict.has_key(chrid):
				if int(depthdict[chrid][1])>depthcut: #check if it's below the set depth.
					if vcfdict[basedir].has_key(chrid): #see if we have a variant called
						if line2[0] == "chrX" and ismale and 0.25<=float(vcfdict[basedir][chrid][2]):
							woc +=1
							outstr += vcfdict[basedir][chrid][1]+"\t"
						elif line2[0] == 'chrY' and 0.25<=float(vcfdict[basedir][chrid][2]):
							woc +=1
							outstr += vcfdict[basedir][chrid][1]+"\t"
						
						else: #is not a sex chromosome or is chrX and needs to be two alleles
							if 0.25<=float(vcfdict[basedir][chrid][2])<0.75: #heterozygous
								woc +=1
								outstr += vcfdict[basedir][chrid][0]+vcfdict[basedir][chrid][1]+"\t"
							elif 0.75<=float(vcfdict[basedir][chrid][2]): #homozygous ALT
								woc +=1
								outstr += vcfdict[basedir][chrid][1]+vcfdict[basedir][chrid][1]+"\t"
							else: #lower than 25% so we'll call this homozygous reference
								print "variant Called but below 25%", vcfdict[basedir][chrid], refdict[chrid]
								fhw2.write("Variant Called but below 25%%.\t%s\t%s\n"% (vcfdict[basedir][chrid], refdict[chrid]))
								woc +=1
								outstr += vcfdict[basedir][chrid][0]+vcfdict[basedir][chrid][0]+"\t"
					else: #no variant called so we just report reference if we have sufficient depth. Homozygous reference
						if line2[0] == 'chrY':
							woc +=1
							outstr += line2[3].strip()+"\t"
						else: #not a chrY 
							woc +=1
							outstr+= line2[3].strip()+line2[3].strip()+"\t"
				else: #lower depth for now print out NN because of low depth
					print "low depth", line, depthdict[chrid]
					if vcfdict[basedir].has_key(chrid): #see if we have a variant called
						fhw2.write("variant called but LOW depth (<%s): %s\n"%(depthcut, line))
						print vcfdict[basedir][chrid][3], "with variant low depth"
					fhw2.write("low depth\t%s\t%s\n"%('\t'.join(depthdict[chrid]), line))
					if line2[0] == "chrY":
						woc +=1
						outstr +="N"+"\t"
					else:#Not a Y chr 
						if line2[0]=='chrX' and ismale: #is there a Y somewhere?
							woc +=1
							outstr += "N"+"\t"
						else:
							woc +=1
							outstr +="NN"+"\t"
			else: #check for variant call but without depth this will require some deep dive
				if vcfdict[basedir].has_key(chrid):
					print "Dive Deep on ", line
					fhw2.write("variant called but without depth: %s\n"%(line))
				else: #no variant was called adn we lack depth
					print "NO depth", line
					fhw2.write("NO depth detected and no variant (likely chrY): %s\n"%(line))
				if line2[0] == "chrY":
					woc +=1
					outstr +="N"+"\t"
				else:#X or non-sex chromosome
					if line2[0]=='chrX' and ismale:
						woc +=1
						outstr += "N"+"\t"
					else:
						print basedir, line
						woc +=1
						outstr +="NN"+"\t"
			if not woc == rlc:
				#the odd case
				if woc >rlc:
					fherror.write("woc is higher than the ref line, %s\n"%(line))
					print woc, rlc
					while rlc < woc:
						rlc += 1
				elif woc < rls: #more common
					fherror.write("refline didn't get analyzed, %s\n"%line)
					print woc, rlc
					while woc <rlc:
						woc +=1
		fh.close()
		fhw.write(outstr+'\n')
		fhw.close()
		fhw2.close()
		fherror.close()

def paircompare(samples): #this is the tool to compare a T/N vcf and provide the statistics back for the pair identified takes in the dictionary of dictoionaries keyed on sample name (basedirs)
	#initially the length should only be two. so I'm writing this for that (although I'm going to iterate so that the first is the main tumor sample and the rest are normals in this case)
	normals={}

	tumor = 'final/%s/%s-Fingerprint.txt'%(samples.keys()[0], samples.keys()[0])
	tumorbasedir = samples.keys()[0]
	for samplekey in samples.keys()[1:]:
		#commented out assuming no duplication here?
		# if not normals.has_key(samplekey):
		# 	normals[samplekey]=[]
		normals[samplekey]='final/%s/%s-Fingerprint.txt'%(samplekey, samplekey)

	for normalsample in normals:
		fh = open(tumor)
		a = fh.readline()
		fh.close()
		print normals, normalsample
		fh = open(normals[normalsample])
		b = fh.readline()
		fh.close()
		
		if a == b:
			print "Genotype matches 100 Percent between: %s and %s"%(tumorbasedir, normalsample)
			print "%s\n%s"%(tumorbasedir, '\t'.join(a.strip().split()))
			print "%s\n%s"%(normalsample, '\t'.join(b.strip().split()))
			fhw2.write("Genotype matches 100 Percent between: %s and %s\n"%(tumorbasedir, normalsample))
			fhw2.write("%s\n%s\n"%(tumorbasedir, '\t'.join(a.strip().split())))
			fhw2.write("%s\n%s\n"%(normalsample, '\t'.join(b.strip().split())))

		else: #if not an exact match print out metrics of how close we are
			a2 = a.strip().split()
			b2 = b.strip().split()
			count = 0 #init matches
			totalcount = 0
			for c1,c2  in izip(a2,b2):
				if c1==c2:
					count +=1
					totalcount+=1
				else:
					totalcount+=1
			print "NO MATCH between: %s and %s"%(tumorbasedir, normalsample)
			print '\t'.join(genelist)
			print "%s\n%s"%(tumorbasedir, '\t'.join(a.strip().split()))
			print "%s\n%s"%(normalsample, '\t'.join(b.strip().split()))
			print "However we matched %s Percent of the reads (%s out of %s total)" % ((float(count)/float(totalcount)),count, totalcount)
			print "Genes where we mismatch: "+",".join([genelist[i] for i,(c1,c2)  in enumerate(izip(a2,b2)) if c1!=c2])
			fhw2.write("NO MATCH between: %s and %s\n"%(tumorbasedir, normalsample))
			fhw2.write('\t'.join(genelist)+'\n')
			fhw2.write("%s\n%s\n"%(tumorbasedir, '\t'.join(a.strip().split())))
			fhw2.write("%s\n%s\n"%(normalsample, '\t'.join(b.strip().split())))
			fhw2.write("However we matched %s Percent of the reads (%s out of %s total)\n" % ((float(count)/float(totalcount)),count, totalcount))
			fhw2.write("Genes where we mismatch: "+",".join([genelist[i] for i,(c1,c2)  in enumerate(izip(a2,b2)) if c1!=c2])+'\n')
			print "============================================================="
			fhw2.write("=============================================================\n\n")


##main#####
#step one gather up the needed files (some old stuff grabbed here bigwigs and vcf files (now using sambamba and vardict to quickly custom call the sites))
args=parser.parse_args()
basedirs =[]
if args.bed:
	snpbedfile = args.bed()
if args.ref:
	reffile = args.ref()
if args.depth:
	depthcut = args.depth()
homedir = os.getcwd()
#now making individual vcf files and not relying on the annotated output don't use bigwigs either. 
# vcffiles = glob.glob('final/*/*-vardict*.anno.filt.vcf.gz')
# if len(vcffiles) ==0:
# 	vcffiles = glob.glob('final/*/varFilter/*-vardict*.anno.filt.vcf.gz')
# bigwigfiles = glob.glob('final/*/*-ready.bigwig')
bamfiles = glob.glob('final/*/*-ready.bam')
bdvcf=set()
bdbw=set()
genelist=[]

#next check if folder structure for output is there. if not make it
if not os.path.exists(homedir+"/work"):
	os.mkdir(homedir+"/work")

###NEW fingerprintvcf creation#################
##Step 2 ) quickly run vardict and create a new vcf for fingerprinting, then use bcftools to query => filter out the DP and AF (needed for tumor normal samples)
#written to each dir/bcbio/final/subdir/ as idtfingervar.txt, idtfingervar.vcf.gz, and finally to fingerprintvcf.txt
print "this code is running its own vardict to create a vcf for the IDT fingerprinting panel and will not use the annotated vcf file. "
print bamfiles
for bamfile in bamfiles:
	basedir = bamfile.split('/')[1]
	if not os.path.exists(homedir+"/work/"+basedir):#check if there is a basedir in the work for temp files
		os.mkdir(homedir+"/work/"+basedir)
	if not glob.glob(homedir+'/final/%s/fingerprintvcf.txt'%(basedir)): #added to skip completed runs (should just add reuse option)
		#use vardict.pl to quickly scan the ready.bam for files
		cline = "vardict.pl -b %s -N %s -f 0.005 %s > %s/work/%s/idtfingervar.txt"%(bamfile, basedir, snpbedfile, homedir, basedir)
		print cline
		pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
		cline = 'var2vcf_valid.pl final/%s/idtfingervar.txt > %s/work/%s/idtfingervar.vcf'%(basedir, homedir, basedir)
		print cline
		pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
		cline = 'bgzip -f %s/work/%s/idtfingervar.vcf'%(homedir, basedir)
		print cline
		pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
		cline = 'tabix %s/work/%s/idtfingervar.vcf.gz'%(homedir, basedir)
		print cline
		pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
		cline = "bcftools query -H -R %s -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE AF=%%AF][\\t%%SAMPLE DP=%%DP]\\t%%ID\\t%%LINE' %s/work/%s/idtfingervar.vcf.gz >%s/work/%s/fingerprintvcf.txt"%(reffile, homedir, basedir, homedir, basedir)
		#cline = "bcftools query -H -R %s -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE AF=%%AF][\\t%%SAMPLE DP=%%DP][\\t%%SAMPLE GT=%%GT]\\t%%LINE' %s >final/%s/fingerprintvcf.txt"%(reffile, vcffile, basedir)
		print cline
		pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
	bdvcf.add(basedir)


##Step 3 ) run sambamba to get depth at the positions for all samples simultaneously output to dir/bcbio/final/sambambadepth.txt
cline = 'sambamba depth region -L %s -o %s/final/sambambadepth.txt final/*/*-ready.bam'%(snpbedfile, homedir)
print cline
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()

##step 4) parse this info into dictionaries so taht you can check if it has a variatn called and what the depth is.
fh =open(homedir+'/final/sambambadepth.txt')
fh.readline() #gather infromatino and indexing from here if need be later
bwdict={} #key this on sample names
for line in fh:
	line2 = line.strip().split()
	chrom = line2[0]
	start = int(line2[1])
	end = int(line2[2])
	basedir = line2[-1]
	readcount = int(line2[-3])
	avgcount = float(line2[-2])
	if not basedir in bdbw: #add it and create a dictionary of dictionaries to contain the write outs at the end for depth file
		bdbw.add(basedir)
		bwdict[basedir]={} #key this on chrid
	if max(int(start),int(end)) - min(int(start),int(end))>1:
		print line, "WAS LONGER IN DEPTH" #just identify the longer than a SNP positions (make this an actual error later or save for gender calling)
	bwdict[basedir][chrom+"_"+str(start)]=readcount
fh.close()
#now parse the dictionary according to the snp file (later delete this step and just output the sambamba directly )
for basedir in bwdict:
	#if not glob.glob('final/%s/fingerprintdepth.txt'%(basedir)):
	fh = open(snpbedfile) #open the bed file and parse it 
	fhw = open('%s/work/%s/fingerprintdepth.txt'%(homedir, basedir),'w')
	ytotal = 0
	for line in fh:
		line2 = line.split('\t')
		chrom = line2[0]
		start = int(line2[1])
		fhw.write(line.strip()+'\t1\t'+str(bwdict[basedir][chrom+"_"+str(start)])+'\n') #this has a 1 in it because bedtools puts the coverage fraction in the depth output so I'm matching in... going forward can remove and decrease the position later in teh code
		if chrom == "chrY":
			ytotal += bwdict[basedir][chrom+"_"+str(start)]		
	fhw.close()
	fh.close()



			

#ok now parse through the checklist for each base directory that was analyzed and kick out the final files. (parallelize these (defs for each of the above) pool and set threads to 6(?) Do I need this here? maybe for the checking? single pipeline that? It's here now so I can print out the order and make sure things are going through correctly for the output.
refdict = {} # to process the checklist
fh=open(snpbedfile)
for line in fh:
	line2 = line.split('\t')
	chrid = line2[0]+"_"+line2[2]
	name = line2[-1].strip()
	refdict[chrid]= name
	genelist.append(name)
fh.close()
	
#print any files missing either the bigwig or vardict vcf
if bdbw^bdvcf:
	for bdir in bdbw^bdvcf:
		print "lacking depth or vardict vcf in: %s" %(bdir)

#old print loop
# for item,y in enumerate(genelist): #iterate through the gene list and print out the genes in the analysis
# 	print item,y 

print bdbw, bdvcf
printflag = True
#run through things normally and check for inclusion of both bigwig and vcf (tumor normal will have all in one vcf)
for eachdir in bdbw&bdvcf:
	#send to the def and let it return a dictionary with keys 
	print eachdir, "eachdir"
	samples= fingerprint(homedir, eachdir) #this is a dictionary of dictionaries keyed on the sample name and then the chrid
	
	#print samples
	fingerwriter(homedir, samples) #this will write out the files
	

	#if paired then we want to make the immediate comparison. Else we want to run comparison later so def the comparison stuff into this. 
	if len(samples)>1: # we have a paired analysis. 
		if printflag:
			finaldatefolder = glob.glob(homedir+'/final/*_bcbio')
			fhw2 =open(finaldatefolder[0]+"/Paired_FingerPrint_Comparison.txt",'w')
			printflag = False
		if len(samples)>2:
			print "have more than pairs need to figure things out on pairing"
			#although if the first is the main sample the rest can be iterated through (Go back and take out of genotypecompareV2 to do this)
		print "paired sample"
		paircompare(samples)
		#do something to tag the two files in the folder or
if not printflag:
	fhw2.close()

###new 7/11/16
#now generate the comprehensive blast dictionary	
#takes the -Fingerprint.txt output and creates a fasta 
blastfasta = open('final/longprints.fasta','w')
fingerprinttexts = glob.glob('final/*/*-Fingerprint.txt')
for fpfile in fingerprinttexts:
	basedir = fpfile.split('/')[1]
	blastfasta.write(">%s\n"%(basedir))
	fh=open(fpfile)
	seqline=fh.readline()
	blastfasta.write(seqline.replace("\t",''))
	fh.close()
	fhw=open('final/%s/%s-longprint.fasta'%(basedir, basedir),'w')
	fhw.write(">%s\n"%(basedir))
	fhw.write(seqline.replace("\t",''))
	fhw.close()
blastfasta.close()