__author__ = 'akoziol'

# Import the necessary modules
# OS is used for file/folder manipulations
import os
# Subprocess->call is used for making system calls
import subprocess
# Errno is used in the file creation command  - I think it's similar to the $! variable in Perl
import errno
# Glob finds all the path names matching a specified pattern according to the rules used by the Unix shell
import glob
# Shutil is useful for file moving/copying
import shutil
# prints variables in a form which can be used as input to the interpreter - similar to Data::Dumper?
#import pprint
# Regex
import re
# System tools
import sys
# Can import date, and calculate length of run, etc.
import time
# Multiprocessing module
from multiprocessing import Pool
# Numerical python - used in the parsing of vcf files
import numpy
# Math module - used in the parsing of vcf files
import math
# Argument parser for user-inputted values, and a nifty help menu
from argparse import ArgumentParser

#Parser for arguments
parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
parser.add_argument('-v', '--version', action='version', version='%(prog)s v1.0')
parser.add_argument('-i', '--input', required=True, help='Specify input directory')
parser.add_argument('-l', '--readLength', required=True, help='Specify list of read lengths to be used e.g. 18, 19, 20, 21, 22')
parser.add_argument('-f', '--foldCoverage', required=True, help='Specify list of fold coverage values to be used e.g. 1, 1.5, 2, 2.5, 5')
parser.add_argument('-k', '--kmerLength', required=True, help='Specify list of kmer lengths to be used e.g. 5, 7, 11')

# Get the arguments into a list
args = vars(parser.parse_args())

# Define variables from the arguments - there may be a more streamlined way to do this
path = args['input']

# Since the following three variables need to be lists, the string entered is split on the ,
readLength = args['readLength'].split(",")
foldCoverage = args['foldCoverage'].split(",")
kmer = args['kmerLength'].split(",")

# maxRL will be used in the fold coverage calculations during simulated file generation
maxRL = readLength[-1]

# Initialize the required dictionaries
vcfData = {}

# Target files are in scriptPath/targets
scriptPath = os.getcwd()
os.chdir("%s/targets" % scriptPath)
# Must have a .fa extension
targetFiles = glob.glob("*.fa")

# Add targets to appropriate list
targets = ["%s/targets/" % scriptPath + fastaFile for fastaFile in targetFiles]

# Genome reference files need to be in "path/reference"
os.chdir("%s/reference" % path)

# File extension needs to be .fa .fas .fasta, but not .ffn .fsa, etc.
referenceFile = glob.glob("*.fa*")

# Add all the reference to the references list
references = ["%s/reference/" % path + fastaFile for fastaFile in referenceFile]

# Define the path for the outputs
outPath = "%s/outputs" % path


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def createSimulatedFilesProcesses(reference):
    """Creates a pool of processes, and maps data in a parallel fashion to createSimulatedFiles"""
    print "Creating simulated files"
    # Initialise the args list
    simulatedArgs = []
    # Every Python module has its __name__ defined and if this is '__main__',
    # it implies that the module is being run standalone by the user and we can do corresponding appropriate actions.
    # http://ibiblio.org/g2swap/byteofpython/read/module-name.html
    if __name__ == '__main__':
        # Initialise the pool of processes - it defaults to the number of processors
        simulatedFilepool = Pool()
        # Create a tuple of the appropriate read lengths and fold coverages
        # eg. (30, 1), (30, 2), ... (30, 100), (35, 1), (35, 2), ... (150, 100)
        for rLength in readLength:
            for fCov in foldCoverage:
                simulatedArgs.append((rLength, fCov, reference))
                # Use the map function and the tuple created above to process the data rapidly
        simulatedFilepool.map(createSimulatedFiles, simulatedArgs)


def createSimulatedFiles((rLength, fCov, reference)):
    """Iterates over the readLength and foldCoverage lists to create folders (if necessary)\
     and perform analyses"""
    os.chdir(path)
    # Create a new folder(if necessary) at the appropriate location
    newPath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
    newFile = "%s/%s_%s" % (newPath, rLength, fCov)
    # The adjusted coverage keeps the number of reads constant for each readLength value supplied.
    # a modelling experiment with a readLength of 20 will have a adjCov that is 40% the value of
    # one with a readLength of 50
    adjCov = float(fCov) * float(rLength)/float(maxRL)
    # If using Sakai as the reference, then multiplying the foldCoverage by the constant below
    # will allow for the use of precise read lengths - using a foldCoverage value of 5 will yield almost
    # exactly 500 000 reads
    # adjCov = float(fCov) * 0.90935049 * float(rLength)/float(maxRL)

    # Call art_illumina to simulate the reads into the appropriate folders - general format of system call:
    # art_illumina -i /path-to-file/Escherichia_coli_O157_H7_str_Sakai.fas -l "readLength" -f "foldCoverage" \
    # -m 225 -s 60 -o /path-to-folder/Appropriate_name
    artIlluminaCall = "art_illumina -i %s -l %s -f %s -o %s" % (reference, rLength, adjCov, newFile)
    # If not using an adjusted coverage value, then uncomment the line below
    # artIlluminaCall = "art_illumina -i %s -l %s -f %s -o %s" % (reference, rLength, float(fCov), newFile)

    make_path(newPath)

    if not os.path.isfile("%s.fq" % newFile):
        sys.stdout.write('.')
        # Subprocess.call requires that the command be finished before the loop can continue
        # this ensures that processes will not be started, and continue running, while the
        # script believes that it is "safe" to start more processes, eventually leading to problems
        subprocess.call(artIlluminaCall, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    else:
        print sys.stdout.write('.')


def faidxTargetsProcesses():
    """Faidx multiprocessing helper function"""
    print '\nProcessing targets with faidx'
    # Initialise the args list
    if __name__ == '__main__':
        # Initialise the pool of processes - it defaults to the number of processors
        faidxPool = Pool()
        faidxPool.map(faidxTargets, targets)


def faidxTargets(file):
    """Creates .fai index files of the targets, which are necessary for the conversion
    of sorted BAM files to fastq files."""
    # print json.dumps(file, sort_keys=True, indent=4, separators=(',', ': '))
    fileName = os.path.split(file)[1]
    faidxFile = "%s.fai" % file
    # print faidxFile
    faidxPath = "%s/targets/faidxFiles" % path
    make_path(faidxPath)
    if not os.path.isfile("%s/%s.fai" % (faidxPath, fileName)):
        faidxCommand = "samtools faidx %s" % file
        subprocess.call(faidxCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        # Move the file and faidx-processed file to the appropriate folder for further processing
        shutil.move(faidxFile, faidxPath)
        shutil.copy(file, faidxPath)
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def indexTargets():
    """Performs smalt index on the targets using the range of k-mers stored in the variable kmer"""
    print '\nIndexing targets'
    for target in targets:
        for size in kmer:
            # Format the target names properly
            filename = os.path.split(target)[1]
            fileNoExt = filename.split(".")[0]
            # Create a new path to be created (if necessary) for the generation of the range of k-mers
            indexPath = "%s/targets/%s/%s_%s" % (path, fileNoExt, fileNoExt, size)
            # Call the make_path function to make folders as necessary
            make_path(indexPath)
            indexFileSMI = "%s.smi" % fileNoExt
            indexFileSMA = "%s.sma" % fileNoExt
            # Index the appropriate files
            if not os.path.isfile("%s/%s" % (indexPath, indexFileSMI)):
                indexCommand = "smalt index -k %s -s 1 %s %s/targets/faidxFiles/%s" % (size, fileNoExt, path, filename)
                subprocess.call(indexCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
                shutil.move(indexFileSMI, indexPath)
                shutil.move(indexFileSMA, indexPath)
                sys.stdout.write('.')
            else:
                sys.stdout.write('.')


def mappingProcesses():
    """Mapping threads!"""
    os.chdir(path)
    print '\nPerforming reference mapping'
    mappingProcessesArgs = []
    if __name__ == '__main__':
        mappingProcessesPool = Pool()
        # uses kmer, targets, readLength, foldCoverage
        for rLength in readLength:
            for fCov in foldCoverage:
                for target in targets:
                    for size in kmer:
                        mappingProcessesArgs.append((rLength, fCov, target, size))
        mappingProcessesPool.map(mapping, mappingProcessesArgs)


def mapping((rLength, fCov, target, size)):
    """Performs the mapping of the simulated reads to the targets"""
    filename = os.path.split(target)[1]
    fileNoExt = filename.split(".")[0]
    megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, fileNoExt, size)
    filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
    newPath = "%s/%s" % (filePath, megaName)
    make_path(newPath)
    targetPath = "%s/targets/%s/%s_%s" % (path, fileNoExt, fileNoExt, size)
    # Map the files to the reference
    if not os.path.isfile("%s/%s.bam" % (newPath, megaName)):
        smaltMap = "smalt map -o %s/%s.bam -f bam -x %s/%s %s/%s_%s.fq" \
                   % (newPath, megaName, targetPath, fileNoExt, filePath, rLength, fCov)
        subprocess.call(smaltMap, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def sortingProcesses():
    """Multiprocessing for sorting bam files"""
    print "\nSorting bam files"
    sortingProcessesArgs = []
    if __name__ == '__main__':
        sortingProcessesPool = Pool()
        # uses kmer, targets, readLength, foldCoverage
        for rLength in readLength:
            for fCov in foldCoverage:
                for target in targets:
                    for size in kmer:
                        sortingProcessesArgs.append((rLength, fCov, target, size))
        sortingProcessesPool.map(sorting, sortingProcessesArgs)


def sorting((rLength, fCov, target, size)):
    """Performs samtools sort to return a sorted bam file"""
    filename = os.path.split(target)[1]
    fileNoExt = filename.split(".")[0]
    megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, fileNoExt, size)
    sorted = megaName + "_sorted"
    sortedMegaName = megaName + "_sorted.bam"
    filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
    newPath = "%s/%s" % (filePath, megaName)
    #Sort the BAM file
    if not os.path.isfile("%s/%s" % (newPath, sortedMegaName)):
        bamSort = "samtools sort %s/%s.bam %s/%s" % (newPath, megaName, newPath, sorted)
        subprocess.call(bamSort, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def bamIndexingProcesses():
    print '\nIndexing bam files'
    bamIndexingArgs = []
    if __name__ == '__main__':
        bamIndexingPool = Pool()
        # uses kmer, targets, readLength, foldCoverage
        for rLength in readLength:
            for fCov in foldCoverage:
                for target in targets:
                    for size in kmer:
                        bamIndexingArgs.append((rLength, fCov, target, size))
        bamIndexingPool.map(bamIndexing, bamIndexingArgs)


def bamIndexing((rLength, fCov, target, size)):
    """Indexes the sorted bam files in order to visualize the assemblies with tablet - note this is OPTIONAL"""
    filename = os.path.split(target)[1]
    fileNoExt = filename.split(".")[0]
    megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, fileNoExt, size)
    sortedMegaName = megaName + "_sorted.bam"
    filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
    newPath = "%s/%s" % (filePath, megaName)
    indexedName = megaName + "_sorted.bai"
    if not os.path.isfile("%s/%s" % (newPath, indexedName)):
        bamIndex = "samtools index %s/%s %s/%s" % (newPath, sortedMegaName, newPath, indexedName)
        subprocess.call(bamIndex, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def createVCFProcesses():
    print '\nCreating vcf files'
    createVCFArgs = []
    if __name__ == '__main__':
        createVCFPool = Pool()
        # uses kmer, targets, readLength, foldCoverage
        for rLength in readLength:
            for fCov in foldCoverage:
                for target in targets:
                    for size in kmer:
                        createVCFArgs.append((rLength, fCov, target, size))
        createVCFPool.map(createVCF, createVCFArgs)


def createVCF((rLength, fCov, target, size)):
    """Creates the variant calling format files from which all relevant data can be pulled"""
    filename = os.path.split(target)[1]
    fileNoExt = filename.split(".")[0]
    megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, fileNoExt, size)
    sortedMegaName = megaName + "_sorted.bam"
    filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
    vcfFile = megaName + "_sorted.vcf"
    newPath = "%s/%s" % (filePath, megaName)
    faidxTarget = "%s/targets/faidxFiles/%s" % (path, filename)
    # Read this to understand why certain flags were used
    # http://samtools.sourceforge.net/mpileup.shtml
    if not os.path.isfile("%s/%s" % (newPath, vcfFile)):
        vcfPipe = "samtools mpileup -A -BQ0 -d 1000000 -uf %s %s/%s | bcftools view -cg - > %s/%s" \
                  % (faidxTarget, newPath, sortedMegaName, newPath, vcfFile)
        subprocess.call(vcfPipe, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def createOutputFiles():
    """Parses the vcf files created above to create a handy summary table of mapping stats"""
    print "\nCreating outputs"
    make_path(outPath)
    os.chdir(outPath)
    outFile = open("SipprModelling_%s.csv" % start, "wb")
    outFile.write("readLength\tfoldCoverage\ttarget\tkmerLength\tMedianQualityScore\t"
                  "QualityScoreSD\tMedianFoldCoverage\tFoldCoverageSD\tMedianPercentID\tqualityMetric\n")
    for rLength in readLength:
            for fCov in foldCoverage:
                for target in targets:
                    for size in kmer:
                        total1 = 0
                        sys.stdout.write('.')
                        filename = os.path.split(target)[1]
                        fileNoExt = filename.split(".")[0]
                        megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, fileNoExt, size)
                        filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
                        vcfFile = megaName + "_sorted.vcf"
                        newPath = "%s/%s" % (filePath, megaName)
                        outputFile = "%s/%s" % (newPath, vcfFile)
                        # Initialise the counter, which will be used to track lines in the vcf file - if positions in the
                        # target are not mapped, then the position field will jump ahead of the counter
                        count = 1
                        # Initialise the arrays, which will keep track of the appropriate values for each dataset
                        arrQual = []
                        arrCov = []
                        arrSum = []
                        output = open(outputFile, "r")
                        for line in output:
                            # vcf files have 36 commented out lines at the top of each file - these are not necessary
                            if re.search('#', line):
                                pass
                            else:
                                total1 += 1
                                # Format of file
                                # CHROM	    POS	ID	REF	ALT	QUAL FILTER	INFO	                                   FORMAT
                                # adk-12	8	.	G	.	32.7	.	DP=1;AF1=0;AC1=0;DP4=0,1,0,0;MQ=29;FQ=-30.3	PL	0
                                # data[0] [1]  [2] [3]  [4] [5]    [6]  [7]
                                data = line.split("\t")
                                #target = data[0]
                                pos = data[1]
                                refSeq = data[3]
                                mapSeq = data[4]
                                qual = data[5]
                                # Depth of coverage is reported prior to the first ";"
                                dpLine = data[7].split(";")[0]
                                # For now, I'm skipping lines that indicated the presence of a possible indel
                                # - I may return to this later
                                if re.search("INDEL", dpLine):
                                    pass
                                else:
                                    # If the called base (mapSeq) is identical to the reference base (refSeq)
                                    # - denoted by a ".", then set seq to equal refSeq, otherwise, pull the
                                    # value of mapSeq for seq
                                    avgQual = sum(arrQual)/total1
                                    if mapSeq == ".":
                                        seq = refSeq
                                        match = 1
                                    # This section corrects for the fact that during the conversion of bam files to vcf
                                    # files, SNP calls and ambiguous calls look identical, except for the fact that for
                                    # SNPs, the qualityScore (qual) tends to be higher than the surrounding bases,
                                    # while ambiguous calls have a lower qualityScore - this loop screens for quality
                                    # scores that are at least 10 lower than the score of the previous base
                                    else:
                                        if float(arrQual[-1] - 10) >= 0:
                                            prevValue = float(arrQual[-1] - 10)
                                        else:
                                            prevValue = 0
                                        if float(qual) <= prevValue:
                                            seq = refSeq
                                            match = 1
                                        else:
                                            # This attempts to catch if there are two ambiguous bases in a row;
                                            # they will hopefully have the same value
                                            if float(qual) == prevValue:
                                                seq = refSeq
                                                match = 1
                                            else:
                                                # "True" SNPs seem to have increased qualityScore compared to the
                                                # surrounding values, this will catch that
                                                if float(qual) > prevValue:
                                                    seq = mapSeq
                                                    match = 0
                                    # Strip the "DP=" from dpLine
                                    DP = dpLine.split("=")[1]
                                    #vcfData[pos] = (fileName, target, refSeq, mapSeq, DP)
                                    # If pos > count, then there is a gap in the mapping (or a deletion, but ignoring
                                    # this possibility for now). For my purposes, I want to have data reported for
                                    # every position, whether it is present in the vcf file or not, so I will use count
                                    # as the position, "-" as the seq, and 0 as the quality and depth of coverage
                                    if int(pos) > count:
                                        #print int(pos) - count, pos, count, range(count, int(pos))
                                        # the number of skipped positions is equal to the value for pos - count
                                        # For each skipped position (i), set appropriate variables to appropriate values
                                        for i in range(count, int(pos)):
                                            posAdj = count
                                            seqAdj = "-"
                                            matchAdj = 0
                                            qualAdj = 0
                                            DPAdj = 0
                                            #vcfData[fileName][rL][fC][target][size][int(posAdj)][seqAdj][matchAdj][qualAdj] = DP
                                            arrQual.append(float(qualAdj))
                                            arrCov.append(float(DPAdj))
                                            arrSum.append(float(matchAdj))
                                            count += 1
                                            if int(pos) == count:
                                                #vcfData[fileName][rL][fC][target][size][int(pos)][seq][match][qual] = DP
                                                arrQual.append(float(qual))
                                                arrCov.append(float(DP))
                                                arrSum.append(float(match))
                                                count += 1
                                    else:
                                        #vcfData[fileName][rL][fC][target][size][int(pos)][seq][match][qual] = DP
                                        arrQual.append(float(qual))
                                        arrCov.append(float(DP))
                                        arrSum.append(float(match))
                                        count += 1
                        # In the case of no data being present in a file,
                        total = count - 1
                        if total == 0:
                            avgQual = 0
                            stdQual = 0
                            avgCov = 0
                            stdCov = 0
                            avgID = 0
                            qualMet = 0
                        else:
                            avgQual = sum(arrQual)/total
                            stdQual = numpy.std(arrQual)
                            avgCov = sum(arrCov)/total
                            stdCov = numpy.std(arrCov)
                            avgID = sum(arrSum)/total * 100
                            qualMet = avgQual * avgCov

                        outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                    % (rLength, fCov, fileNoExt, size, avgQual, stdQual, avgCov, stdCov, avgID, qualMet))

                        output.close()
    outFile.close()


def pipeline():
    """Calls all the functions in a way that they can be multi-processed"""
    for reference in references:
        createSimulatedFilesProcesses(reference)
        faidxTargetsProcesses()
        indexTargets()
        #Start the mapping operations
        mappingProcesses()
        sortingProcesses()
        bamIndexingProcesses()
        createVCFProcesses()
        createOutputFiles()

start = time.time()
pipeline()
print "\nElapsed Time: %s seconds" % (time.time() - start)