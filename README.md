# geneSippr README

This repository contains two major scripts and one accessory script, as well as all necessary targets for GeneSipping

###modellingMultiprocessing.py 

Python script used to simulate fastq files with controlled read length and fold coverage values from a reference genome. FASTQ files are mapped 
to gene targets with an adjustable kmer value. Creates a summary report of mapping results 


###rawReadSippingUpdate.pl

Perl script used to convert bcl basecall files to FASTQ files, and subsequently map FASTQ files to user-defined targets. Creates a summary
report of mapping results. Optionally can call gapGraphing.R to create histogram summary of gap size and frequency 


This software was designed and run on Linux. 

## Requirements

The following software is necessary in order for GeneSippr to function

- Python, Perl, R
- [BCL2FASTQ Conversion Software](http://support.illumina.com/downloads/bcl2fastq_conversion_software.html "bcl2fastq")
- [smalt](https://www.sanger.ac.uk/resources/software/smalt/)
- [samtools](http://samtools.sourceforge.net/)
- [bcftools](http://samtools.github.io/bcftools/bcftools.html)
- [art_illumina](http://wiki.techfak.uni-bielefeld.de/cmg/Art_Illumina) - Only required for the modelling script

## Use

###Installation

Install the repo using the method of your choice

For instance, in your terminal, navigate to the local directory of your choice and enter
`git clone https://github.com/OLC-LOC-Bioinformatics/geneSippr.git`

Make sure the scripts have execute permissions

###Inputs

####modellingMultiprocessing.py 

Requires four arguments:

* -i, --input: Input directory. Sets the path variable in the script
* -l, --readLength: List of read lengths to be used e.g. 18, 19, 20, 21, 22
* -f, --foldCoverage: List of fold coverage values to be used e.g. 1, 1.5, 2, 2.5, 5
* -k, --kmerLength: List of kmer lengths to be used e.g. 5, 7, 11

Additionally, the scrip requires sequence files in two locations:

* a reference genome in the "path/reference" folder. This genome must have a file extension of .fa* (.fa, .fasta, .fas, etc.)
* target files in the supplied "script path/targets" folder. The files in this folder must have a .fa extension. The folder supplied
in the repository contains 19 gene targets specific to *Escherichia coli*, strain Sakai

####rawReadSipping.pl

Ideally, the MiseqOutput folder on the MiSeq should be mounted as a network drive. Either that, or copy the entire "MiSeqOutput/[run Name]" folder to the location of your choice.

As with modellingMultiprocessing.py, a folder of target sequences is required. This folder is supplied in the repository as "script path/target/allTargets" 

Required arguments:

* -m, --miseqPath: The path of the folder that contains the run data
* -f, --folderPath: The path of the folder in which to place the output folder
* -o, --outPath: The name of the output folder
* -t, --targetPath: The path of the folder that contains the target sequences

Optional arguments:

* -r, --readLength: Length of forward reads to use
* -p, --project: Name of the project

###Outputs

####modellingMultiprocessing.py 

* "path/tmp" folder with simulated FASTQ, mapped, sorted, and index BAM files, as well as VCF (variant calling format) files - these files can be deleted
  once you are sure you won't need them anymore
* "path/targets" folder with the targets
* "path/outputs" folder with a .csv report of the target mapping results

####rawReadSipping.pl

* "folderPath/outPath/project/query" containing the FASTQ files generated for the analyses
* "folderPath/outPath/project/results" containing the temporary files generated in the mapping portion of the script - these files can be deleted, though if
an analysis need to be re-performed, these files can be reused
* "folderPath/outPath/project/reports" containing a .csv report of the target mapping results


