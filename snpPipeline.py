#!/usr/bin/env python3
import argparse, sys, os, time, threading, subprocess as sp
from collections import deque

try:
    import Queue
except ImportError:
    import queue as Queue


"""
   _____  _____  _          ___         
  / __/ |/ / _ \(_)__  ___ / (_)__  ___ 
 _\ \/    / ___/ / _ \/ -_) / / _ \/ -_)
/___/_/|_/_/  /_/ .__/\__/_/_/_//_/\__/ 
               /_/                      
#Python by Ryan Culligan
"""

def trimmomatic():
	forwardPairedReads = []
	reversePairedReads = []
	forwardUnpairedReads = []
	reverseUnpairedReads = []

	for forward, reverse in zip(forSamples, revSamples):
		forBase = forward.split(".")[0]
		revBase = reverse.split(".")[0]
		cmd =  "java -jar trimmomatic.jar {readType} -{phredType} ".format(**globals())
		
		for baseName in (forBase, revBase):
			paired   = "{baseName}.paired.fq.gz".format(**locals())
			unpaired = "{baseName}.unpaired.fq.gz ".format(**locals())
			
			if baseName is forBase:
				forwardPairedReads.append(paired)
				forwardUnpairedReads.append(unpaired)
			else:
				forwardUnpairedReads.append(paired)
				reverseUnpairedReads.append(unpaired)

			cmd += "{paired} {unpaired} ".format(**locals())

		sp.call(cmd, shell=True)

	return (forwardPairedReads, reversePairedReads)


#BEGIN STEP 0
def bwaIndex():
	sp.call("bwa index {referenceFileName}; ".format(**globals()), shell=True)
	return True

def samIndexFa():
	cmd = "samtools faidx {referenceFileName}; ".format(**globals())
	sp.call(cmd, shell=True)
	return True

#Begin Step 1
#SCATTER GATHER!!!
def bwaAlign(f, r):
	cmd =  "bwa mem {referenceFileName} ".format(**globals())
	cmd += "{f} {r} > ".format(**locals())
	cmd += "{bwaOutputFileName}; ".format(**globals())
	return cmd

def samToBam():
	cmd =  'samtools view {bwaOutputFileName} '.format(**globals())
	cmd += '-b -t {referenceFileName} '.format(**globals())
	cmd += '-o {bamOutputFileName}; '.format(**globals())
	return cmd

def picardSort():
	cmd =  "java -jar bin/SortSam.jar "
	cmd += "INPUT={bamOutputFileName} ".format(**globals())
	cmd += "OUTPUT={bamFileWorking} ".format(**globals())
	cmd += "SORT_ORDER=coordinate; "
	return cmd

def picardMark():
	cmd =  "java -jar bin/MarkDuplicates.jar "
	cmd += "INPUT={bamFileWorking} ".format(**globals())
	cmd += "OUTPUT={markedBamFile} ".format(**globals())
	cmd += "REMOVE_DUPLICATES=true METRICS_FILE=dup.txt ASSUME_SORTED=true; "
	return cmd

def addReadGroups():
	cmd =  "java -jar bin/AddOrReplaceReadGroups.jar "
	cmd += "INPUT={markedBamFile} ".format(**globals())
	cmd += "OUTPUT={taggedBamFile} ".format(**globals())
	cmd += "SO=coordinate "
	cmd += "RGID={taggedBamFile} ".format(**globals())
	cmd += "RGLB=1 "
	cmd += "RGPL=illumina "
	cmd += "RGPU=1 "
	cmd += "RGSM={taggedBamFile}; ".format(**globals())
	return cmd

def picardDict():
	cmd =  "java -jar bin/CreateSequenceDictionary.jar "
	cmd += "REFERENCE= {referenceFileName} ".format(**globals())
	cmd += "OUTPUT= {referenceDictName}; ".format(**globals())
	return cmd

def samIndexBam():
	cmd = "samtools index {taggedBamFile}; ".format(**globals())
	return cmd

def createTargets():
	cmd =  "java -jar bin/GenomeAnalysisTK.jar "
	cmd += "-T RealignerTargetCreator "
	cmd += "-nt 10 "
	cmd += "-R {referenceFileName} ".format(**globals()) 
	cmd += "-I {taggedBamFile} ".format(**globals())
	cmd += "-o {intervalOut}; ".format(**globals())
	return cmd

def gatkIndel():
	cmd =  "java -jar bin/GenomeAnalysisTK.jar "
	cmd += "-T IndelRealigner "
	cmd += "--filter_bases_not_stored "
	cmd += "-R {referenceFileName} ".format(**globals()) 
	cmd += "-I {taggedBamFile} ".format(**globals())
	cmd += "-targetIntervals {intervalOut} ".format(**globals())
	cmd += "-o {indelOut}; ".format(**globals())
	return cmd

def gatkHaplo():
	"""
	Not used anymore.  Replaced with gatkDiagnose Targets
	"""
	cmd =  "java -jar bin/GenomeAnalysisTK.jar "
	cmd += "-T HaplotypeCaller "
	cmd += "-R {referenceFileName} ".format(**globals())
	cmd += "-PF out.log "
	cmd += "-I {indelOut} ".format(**globals())
	cmd += "--emitRefConfidence GVCF "
	cmd += "-o {gVcfOut}; ".format(**globals())
	return cmd

#BEGIN STEP 2
def gatkGVCF():
	"""
	Not used anymore.  Replaced with gatkDiagnose Targets
	"""
	cmd =  "java -jar bin/GenomeAnalysisTK.jar "
	cmd += "-T GenotypeGVCFs "
	cmd += "-R {referenceFileName} ".format(**globals())
	cmd += "".join(["--variant {vcf} ".format(**locals()) for vcf in gVcfList]) + " "
	cmd += "-o {outVCF}; ".format(**globals())
	sp.call(cmd, shell=True)
	return True

def gatkDiagnoseTargets():
	"""
	"""
	cmd =  "java -jar bin/GenomeAnalysisTK.jar "
	cmd += "-T DiagnoseTargets "
	cmd += "-R {referenceFileName} ".format(**globals())
	cmd += "".join(["-I {bam} ".format(**locals()) for bam in bamList]) + " "
	cmd += "-L {intervalOut} ".format(**globals())
	cmd += "--minimum_coverage {coverageDepth} ".format(**globals())
	cmd += "-o {outVCF}; ".format(**globals())
	sp.call(cmd, shell=True)
	exit()
	return True

def gatkSelectVariants():
	cmd =  "java -jar bin/GenomeAnalysisTK.jar "
	cmd += "-T SelectVariants "
	cmd += "-R {referenceFileName} ".format(**globals())
	cmd += "--variant {outVCF} ".format(**globals()) #TreeshrewsGT100815.vcf
	cmd += "-o  {variantOut} ".format(**globals()) #TreeshrewsSV101215.vcf 
	cmd += "-selectType SNP "
	cmd += "-selectType MNP "
	cmd += "-env ;"
	sp.call(cmd, shell=True)
	return True

def gatkVariantFilt():
	cmd =  "java -jar bin/GenomeAnalysisTK.jar "
	cmd += "-T VariantFiltration "
	cmd += "-R {referenceFileName} ".format(**globals())
	cmd += "-V {variantOut} ".format(**globals())
	cmd += '--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" '
	cmd += '--filterName "my_snp_filter" '
	cmd += "-o {outFilteredVCF}; ".format(**globals()) #CHANGE
	sp.call(cmd, shell=True)
	return True

#STEP 3
def vcftools(): 
	cmd =  "vcftools --vcf {outFilteredVCF} ".format(**globals())
	cmd += "--remove-filtered-all "
	cmd += "--thin {basePairSep} ".format(**globals())
	cmd += "--max-missing {maxMissing} ".format(**globals())
	cmd += "--recode "
	cmd += "--out {refBaseName};".format(**globals()) #CHANGE
	sp.call(cmd, shell=True)
	return True

def gatkCovered():
	cmd =  "java -jar bin/gatkOld.jar "
	cmd += "-T CoveredByNSamplesSites "
	cmd += "-R {referenceFileName} ".format(**globals())
	cmd += "-V {outFilteredVCF} ".format(**globals()) #TreeshrewsVF101215.vcf #CHANGE
	cmd += "-minCov {coverageDepth} ".format(**globals())
	cmd += "-nt 10 "
	cmd += "-percentage 0.001 "
	cmd += "-out {coveredVCF}; ".format(**globals())
	sp.call(cmd, shell=True)
	return True

def vcfOutputTable():
	cmd =  "vcftools "
	cmd += "--vcf  {coveredVCF} ".format(**globals()) #Treeshrews_2kb_SNP_spacing_rmFail.vcf
	cmd += "--freq "
	cmd += "--out {vcfTable}; ".format(**globals()) #Treeshrew_allelefreqs_spaced 
	sp.call(cmd, shell=True)
	return True

def gatkVarToTable():
	cmd = "java -jar bin/GenomeAnalysisTK.jar "
	cmd += "-R {referenceFileName} ".format(**globals())
	cmd += "-T VariantsToTable "
	cmd += "-V {vcfTable}".format(**globals()) #CHANGE Treeshrews_2kb_SNP_spacing_rmFail.vcf.recode.vcf 
	cmd += "-F CHROM -F POS -F ID -F QUAL "
	cmd += "-o {vcfTable2}; ".format(**globals()) # Treeshrews2275SNPs.table
	sp.call(cmd, shell=True)
	return True

#start pipeline
def runPipeline_step0():
	"""This is untested but should be a viable first step"""
	if untrimmed is True:
		forSamples, revSamples = trimmomatic()	
	bwaIndex()
	samIndexFa()

def runPipeline_step1(forSamples, revSamples):
	commands = []
	commands.append(bwaAlign(forSamples, revSamples))
	commands.append(samToBam())
	commands.append(picardSort())
	commands.append(picardMark())
	commands.append(addReadGroups())
	commands.append(picardDict())
	commands.append(samIndexBam())
	commands.append(createTargets())
	commands.append(gatkIndel())
	#commands.append(gatkHaplo()) #Removed due to deprication of CoveredByNSampleSites
	return commands

def runPipeline_step2():
	#gatkGVCF() #Removed due to deprication of CoveredByNSampleSites
	gatkDiagnoseTargets()
	gatkSelectVariants()
	gatkVariantFilt()

def runPipeline_step3():
	vcftools()
	gatkCovered()
	vcfOutputTable()
	gatkVarToTable()
	return True

if __name__ == "__main__":
	'''
	"""SETUP VARIABLES"""

	parser = argparse.ArgumentParser()

	parser.add_argument('--referenceFileName', action='store', type=str, required=True, help="A reference file to align against")
	parser.add_argument('--forward', nargs="+", action="store", required=True, help="An in order list of forward fastq samples")
	parser.add_argument('--reverse', nargs="+", action="store", required=True, help="An in order list of reverse fastq samples")
	
	parser.add_argument('--numCPU', action='store', type=int, default=1, help="Number of CPUS to use on step 2")
	parser.add_argument('--untrimmed', action='store', type=bool, default=False, help="Run trimmomatic")
	parser.add_argument('--readType', action='store', type= str, default="PE", help="PE or SE")
	parser.add_argument('--phredType', action='store', type=str, default="phred33", help="Type of phred scoring to use.")
	parser.add_argument('--basePairSepLow', action='store', type=int, default=350, help="Low end filter to SNPs that are separated by at least n base pairs")
	parser.add_argument('--basePairSepHigh', action='store', type=int, default=350, help="High end filter to SNPs that are separated by at least n base pairs")
	parser.add_argument('--maxMissingLow', action='store', type=float, default=0.7, help="")
	parser.add_argument('--maxMissingHigh', action='store', type=float, default=0.7, help="")
	parser.add_argument('--maxMissingStep', action='store', type=float, default=0.05, help="")
	parser.add_argument('--coverageDepthLow', action='store', type=int, default=5, help="Low end select only SNPs with a coverage of at least n per individual")
	parser.add_argument('--coverageDepthHigh', action='store', type=int, default=5, help="High end select  only SNPs with a coverage of at least 5 per individual")
	parser.add_argument('--coverageDepthStep', action='store', type=int, default=1, help="High end select  only SNPs with a coverage of at least 5 per individual")
	parser.add_argument('--step0', action='store_true', help="Start pipeline at step 0")
	parser.add_argument('--step1', action='store_true', help="Start pipeline at step 1")
	parser.add_argument('--step2', action='store_true', help="Start pipeline at step 2")
	parser.add_argument('--step3', action='store_true', help="Start pipeline at step 3")
	args = parser.parse_args()
	
	"""
	if args.referenceFileName:
		referenceFileName = args.referenceFileName
	if args.forward:
		forSamples = args.forward
	if args.reverse:
		revSamples = args.reverse
	"""
	if args.numCPU:
		numCPU = args.numCPU
	if args.untrimmed:
		untrimmed = args.untrimmed
	if args.readType:
		readType = args.readType
	if args.phredType:
		phredType = args.phredType
	if args.basePairSepLow:
		basePairSepLow = args.basePairSepLow
	if args.basePairSepHigh:
		basePairSepHigh = args.basePairSepHigh
	if args.maxMissingLow:
		maxMissingLow = args.maxMissingLow
	if args.maxMissingHigh:
		maxMissingHigh = args.maxMissingHigh
	if args.maxMissingStep:
		maxMissingStep = args.maxMissingStep
	if args.coverageDepthLow:
		coverageDepthLow = args.coverageDepthLow
	if args.coverageDepthHigh:
		coverageDepthHigh = args.coverageDepthHigh
	if args.coverageDepthStep:
		coverageDepthStep = args.coverageDepthStep


	if args.step0 is not None:
		step0 = args.step0
	if args.step1 is not None:
		step1 = args.step1
	if args.step2 is not None:
		step2 = args.step2
	if args.step3 is not None:
		step3 = args.step3

	step2 =True

	#Run the pipeline starting at a step
	if step0 == True:
		step1 = True
		step2 = True
		step3 = True
	elif step1 == True:
		step2 = True
		step3 = True
	elif step2 == True:
		step3 = True
	else:
		step0 = True
		step1 = True
		step2 = True
		step3 = True

	step3 = False
	'''

	untrimmed = False
	readType = "PE" #PE OR SE
	phredType = "phred33" #-phred33|-phred64
	basePairSep = 350
	maxMissing = 0.7
	coverageDepth = 5 #CHANGE
	numCPU = 2 #CHANGE TO NUMBER FORSAMPLES

	referenceFileName = "kian8.fa" #CHANGE

	refBaseName = referenceFileName.split(".")[0]
	referenceDictName = refBaseName + ".dict"
	outVCF = refBaseName + ".vcf"
	variantOut = refBaseName + ".variants.vcf"
	outFilteredVCF = refBaseName + ".filtered.vcf"
	recodedVCF = refBaseName + ".recode.vcf"
	coveredVCF = refBaseName + ".covered.vcf"
	vcfTable = refBaseName + ".table.vcf"
	vcfTable2 = refBaseName + ".table2.vcf"

	forSamples = [ "<FILEPATH>/<FILE>.fastq"] #CHANGE
	revSamples = [ "<FILEPATH>/<FILE>.fastq"] #CHANGE 

	gVcfList = []
	bamList  = []
	commands = []

	if step0 == True:
		"""
		STEP 0:
			Run sequence indexing commands, nothing too interesting.
		"""
		runPipeline_step0()

	
	"""
	STEP 1:
		create bamfiles, filter sort, and run GATK for GVCF
	"""
	for f,r in zip(forSamples, revSamples):
		"""GET FILENAMES"""
		baseName = f.split(".")[2].strip("/")
		bwaOutputFileName = baseName + ".sam"
		bamOutputFileName = baseName + ".bam"
		bamFileWorking = baseName + ".sorted.sam"
		markedBamFile = baseName + ".marked.bam"
		taggedBamFile = baseName + ".marked.tag.bam"
		intervalOut = baseName + ".intervals"
		indelOut = baseName + ".realign.bam"
		gVcfOut = baseName + ".g.vcf"

		bamList.append(indelOut)
		gVcfList.append(gVcfOut)
		commands.append(runPipeline_step1(f, r))

	if step1 == True:
		"""
		THREADING STUFF
		Basically there is no need to run prelimanary steps in order 
		this means that we can run the commands concurrently.  
		In other words, if you have more less samples than cpu cores,
		your results will take the same amount of time as if you had 
		1 sample.
		"""
		threadList=[]
		for i in range(numCPU): #NUMBER CPU
			threadList.append("Thread-{}".format(i+1))

		queueLock = threading.Lock()
		workQueue = Queue.Queue(len(commands))
		threads = []
		threadID = 1
		exitFlag = 0

		class myThread (threading.Thread):
			def __init__(self, threadID, name, q):
				threading.Thread.__init__(self)
				self.threadID = threadID
				self.name = name
				self.q = q
			def run(self):
				print("Starting " + self.name)
				process_data(self.name, self.q)#,self.alpha)
				print("Exiting " + self.name)

		def process_data(threadName, q):
			while not exitFlag:
				queueLock.acquire()
				if not workQueue.empty():
					data = q.get()
					queueLock.release()
					for item in data:
						sp.call('{}'.format(item), shell=True)
				else:
				    queueLock.release()
				time.sleep(1)

		# Create new threads
		for tName in threadList:
			thread = myThread(threadID, tName, workQueue)
			thread.start()
			threads.append(thread)
			threadID += 1

		# Fill the queue
		queueLock.acquire()
		for command in commands:
			workQueue.put(command)
		queueLock.release()


		# Wait for queue to empty
		while not workQueue.empty():
			pass
		# Notify threads it's time to exit
		exitFlag = 1

		# Wait for all threads to complete
		for thread in threads:
			thread.join()
		print("Exiting Main Thread")
		exitFlag = 0

	if step2 == True:
		"""
		STEP 2:
			Now that there are results for each read, we can 
			merge all of the gVcf files together to determine 
			snps for all the individual samples.  This cant 
			be done concurrently...
		"""
		#after the loop we can merge the groups together
		runPipeline_step2()
	
	if step3 == True:
		for basePairSep in range(basePairSepLow, basePairSepHigh, basePairSepStep):
			coverageDepth = 1
			maxMissing = maxMissingLow
			while maxMissing <= maxMissingHigh:
				recodedVCF = "ranges/{refBaseName}.{maxMissing}.{basePairSep}.{maxMissing}.recode.vcf".format(**globals())
				coveredVCF = "ranges/{refBaseName}.{maxMissing}.{basePairSep}.{maxMissing}.covered.vcf".format(**globals())
				vcfTable = "ranges/{refBaseName}.{maxMissing}.{basePairSep}.{maxMissing}.table.vcf".format(**globals())
				vcfTable2 = "ranges/{refBaseName}.{maxMissing}.{basePairSep}.{maxMissing}.table2.vcf".format(**globals())
				runPipeline_step3()
				maxMissing += maxMissingStep





	
