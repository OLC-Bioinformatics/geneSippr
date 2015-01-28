#!usr/bin/perl

use warnings;
#use strict;
use Cwd;
use Time::Piece;
use threads;
use File::Path qw(make_path remove_tree);
use File::Find;
use Bio::SeqIO;
use File::Copy;
use Data::Dumper qw(Dumper);
use Statistics::Basic qw(:all);
use List::MoreUtils qw(uniq);
use Getopt::ArgParse;

# This start time will be used in calculating the total time of the run
my $start_time = time;

# Initialize variables
my ($sequenceName, $forwardReverse, $geneName, @folders, $fastaTitle, $reads, $investigator, $flowCell, $project, @stats, @presence, @sequence, @categories, @miSeqDirectories);
my (%results, %sampleSheet, %sippr, %targets, %targetLength, %targetCategories, %targetPresence);

my $lane = 1;

# Determine the number of threads present in the system
my @cpus = `awk '/^processor/ { N++} END { print N }' /proc/cpuinfo`;
chomp @cpus;

# Adds an argument parser to allow the user to easily specify certain variables
$ap = Getopt::ArgParse->new_parser(
        prog        => $0,
        description => 'Creates per-read FASTQ files using per-cycle BCL basecall files with the bcl2fastq package from Illumina. Reference maps FASTQ reads to target files in order to
determine presence/absence of pathogenicity and quality assurance genes',
    epilog      => 'Requires bcl2fastq, smalt, samtools, as well as bcftools',
 );

# , default=> 21
$ap->add_arg('--miseqPath', '-m', required => 1, help => 'The path of the folder that contains the run data');
$ap->add_arg('--folderPath', '-f', required => 1, help => 'The path of the folder in which to place the output folder');
$ap->add_arg('--outPath', '-o', required => 1, help => 'The name of the output folder');
$ap->add_arg('--targetPath', '-t', required => 1, help => 'The path of the folder that contains the target sequences');
$ap->add_arg('--readLength', '-r', required => 0, help => 'Optional. Specify the the length of forward reads to use');
$ap->add_arg('--project', '-p', required => 0, help => 'Optional. Specify the name of the project');


$ns = $ap->parse_args();

my $miSeqPath = $ns->miseqPath;
my $folderPath = $ns->folderPath;
my $outPath = $ns->outPath;
my $filesPath = $ns->targetPath;

# Parse the optional arguments
if ($ns->readLength) {
	$reads = $ns->readLength;
}

if ($ns->project) {
	$project = $ns->project;
}

print "$folderPath $miSeqPath $reads\n";

# This section finds the most recent run folder on the MiSeq and determines which cycle the run is currently on
# NB: This script must be run after the run has been initialised, or the wrong folder will be identified as the current folder, and this will not work!

chdir($miSeqPath);

# Grab all the files in $miSeqPath
my @miSeqFolders = glob("*");

# Ensure that you are in fact working with a directory
foreach (@miSeqFolders) {
	if (-d $_) {
		push(@miSeqDirectories, $_);
	}
}

# Since the folders are all named starting with the date, they can be sorted (from highest to lowest) and the first folder will be the current run
@miSeqFolders = sort{$b cmp $a}(@miSeqDirectories);

my $folder = $miSeqFolders[0];

print "$folder\n";
# Get the flowcell ID from the end of the folder name
$folder =~ /.+_.+_.+_(\S+)/;
$flowCell = $1;

# Make the appropriate folder in the required location - append the folder name to $folderPath
$folderPath .= "/" . $folder;
make_path($folderPath);

# Open the modified sample sheet and write the headers
open(OUTPUT, ">", "$folderPath/SampleSheet_modified.csv") or die $!;
print OUTPUT "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n";

# Change to the appropriate directory
chdir($miSeqPath . "/" . $folder);

# Open the parameters file to determine the run length - this is necessary, as we will be pulling data once the first index has been processed (just over halfway through the run)
open(SHEET, "<", "SampleSheet.csv");
while (<SHEET>) {
	# Grab the investigator name - perhaps because the input file is a CSV, the chop function (x2) seems to be required to format the cells properly
	if (/Investigator Name/) {
		($investigator = $_) =~ s/Investigator Name,//;
		$investigator =~ s/,//g;
		chomp $investigator;
		chop $investigator;
	}
	# The reads counts are set up as follows:
	# [Reads]\n
	# reads1\n
	# reads2\n
	# To get only reads1 - find [Reads]
	if (/\[Reads\]/) {
		# Then start going through SHEET again
		while (<SHEET>) {
			# reads1 will be the first cell
			if (not $ns->readLength) {
				$reads = $_;
				# Get rid of \n
				chomp $reads;
				$reads =~ s/,//g;
			}
			# Exit the loop, as we have the value for reads1
			last;
		}
	}
	# As above - find [Data], then proceed reading through SHEET
	if (/Sample_ID/) {
		while (<SHEET>) {
			my @line = split(/,/);
			# $line[5], $line[7]: index 1 and 2, respectively. $line[9]: description
			###
			chop $line[9]; chop $line[9];
			if (not $ns->project) {
				$project = $line[8];	
			}
			
			my $index = $line[5] . "-" . $line[7];
			print OUTPUT "$flowCell,", 	# flowcell
						 "$lane,",			# lane 1
						 "$line[0],",  	# sample ID
						 "no_ref,",    	# no_ref sampleRef
						 "$index,",    	# index
						 "$line[9],",  	# description
						 "N,",      		# N control
						 "NA,",		 	# NA recipe
						 "$investigator,", # operator
						 "$project\n";		# sample_project
		}
	}
}

# Close the open files
close SHEET;
close OUTPUT;

# Find out how many cycles have been completed
chdir($miSeqPath . "/" . $folder . "/Thumbnail_Images/L001");
my @cycleNum = glob("*C*");
my $cycles = scalar @cycleNum;

# As the number of cycles required is the number of forward reads + the index(8) + the second index(8) + 1 (just to be safe)
my $readsNeeded = $reads + 17;

# A while loop that waits until the required number of cycles has been achieved
while ($cycles < $readsNeeded) {
	running_time("Currently at cycle $cycles. Need to wait until the MiSeq reaches cycle $readsNeeded.");
	sleep(300);
	chdir($miSeqPath . "/" . $folder . "/Data/Intensities/BaseCalls/L001");
	@cycleNum = glob("*C*");
	$cycles = scalar @cycleNum;
}

# Adding 0 to $reads converts it from a string to a number
my $numReads = $reads + 0;

# Sets the base-mask string, which is important for determining which bases to use from each run
my $baseMask = "Y" . $numReads . "n*,I8,I8,n*";

# Call configureBclToFastq.pl - in order to prevent the compression of the fastq files, I had to manually edit the Config.mk file in /usr/local/share/bcl2fastq-1.8.3/makefiles to not include the compression and compression suffix (commented out lines 174 and 175)
unless(-e("$folderPath/$outPath")) {
	print "Calling script\n";
	system("configureBclToFastq.pl --input-dir $miSeqPath/$folder/Data/Intensities/BaseCalls/ --output-dir $folderPath/$outPath --force --sample-sheet $folderPath/SampleSheet_modified.csv --mismatches 1 --no-eamss --fastq-cluster-count 0 --compression none --use-bases-mask $baseMask");
	chdir("$folderPath/$outPath");
	# If you decide to change $baseMask to include reads from both directions, remove the r1 in the command below
	system("nohup make -j 16 r1");
}

# Chdir to the working directory
chdir("$folderPath/$outPath/Project_$project") or die $!;

my $path = getcwd;

# This uncompresses the fastq files - I'm looking at editing the configureBclToFastq.pl script or one of its modules to try and avoid compression in the first place
opendir(DIR, $path) or die "can't find $path: $!";
while (defined(my $file = readdir(DIR))) {
	# Ignore special files
	next if $file =~ /^\.\.?$/;
	# Ignore folders without sequence data
	if ($file =~ m/Sample_/) {
		if (-d $path ."/". $file){
			chdir $path ."/". $file;
			(my $folder = $file) =~ s/Sample_//g;
			my @files = glob("*.gz");
			foreach my $fileZ (@files){
				(my $folder = $fileZ) =~ s/_.*//g;
				(my $unzip = $fileZ) =~ s/.gz//g;
				(my $filename = $fileZ) =~ s/_L001|.gz//g;
				system("gzip -d $fileZ");
				if (-f $unzip){unlink "$fileZ"};
			}
			my @fastqFiles = glob("*001.fastq");			
			make_path("$path/query");
			(my $filename1 = $fastqFiles[0]) =~ s/_L001//g;
			(my $filename2 = $fastqFiles[0]) =~ s/_L001_R1_001.fastq//g;
			unless (-e ("$path/query/$filename2.fastq")) {
				copy("$path/$file/$fastqFiles[0]", "$path/query/$filename2.fastq");
			}
		}
	}
}
print "Processing fastq files\n";
# Get the target files from the folder

#my $filesPath = "/home/blais/git/geneSippr/target";
chdir ("$filesPath");

# Get the names of the two folder with target sequences
opendir(DIR, $filesPath) or die "can't find $path: $!";
while (defined(my $file = readdir(DIR))) {
	# Ignore special files
	next if $file =~ /^\.\.?$/;
	# Ignore folders without sequence data
	if (-d $filesPath . "/" . $file) {
		push(@categories, $file);
		chdir $filesPath . "/" . $file;
		my @targetGenes = glob("*.fa");
		$targetCategories{$file} = \@targetGenes;
	}
}

# Perform necessary indexing and manipulations on target files
while (my ($category, $targetGenes) = each (%targetCategories)) {
	chdir $filesPath . "/" . $category;
	foreach my $gene (@$targetGenes) {
		
		# Perform smalt indexing of the target genes (if necessary) - it might eventually be advisable to play around with the k-mer and step sizes in this command
		(my $geneName = $gene) =~ s/.fa//;
		unless (-e ("$geneName.smi")) {
			system ("smalt index -k 5 -s 1 $geneName $gene");
		}
		# Create the index .fai file for use in mapping
		unless (-e ("$geneName.fai")) {
			system ("samtools faidx $gene");
		}
		# Get the length of the full-length target gene
		open(INDEX, "<", "$gene" . ".fai") or die $!;
		while (<INDEX>) {
			my ($a, $length, $b) = split("\t", $_);
			$targetLength{$geneName} = $length;
		}
		close INDEX;
	}
}

# The files must be in the "query" subfolder. This subfolder must only have sequences that you wish to examine, or the program won't be able to find them
chdir ("$path/query");

#Grab the fastq files for manipulations
my @fastq = glob("*.fastq");

# Use a multithreading approach
while (scalar(@threads) < @cpus) {
	# Loop through each file to be subtyped
	foreach my $file (@fastq) {
		while (my ($category, $targetGenes) = each (%targetCategories)) {
			foreach my $gene (@$targetGenes) {
				# Make a new thread using the blast subroutine - pass appropriate values to subroutine
				my $r = threads->new(\&rawMapping, $path, $filesPath, $file, $category, $gene);
				# Output data to @threads. Right now, all this does is increase the size for the scalar(@threads) portion of the while loop
				push(@threads,$r);
			}
		}
	}
}

# This loop ensures that each thread is complete before terminating
foreach (@threads) {
	my $num = $_->join;
}

# Parses the VCF file and returns arrays of the mapping stats, presence/absence of gene targets, and the sequence of each mapped gene target
foreach my $file (@fastq) {
	foreach my $category (sort keys %targetCategories) {
		$targetGenes = $targetCategories{$category};
		foreach my $gene (@$targetGenes) {	
			my ($stats, $presence, $seq) = parseVCF($file, $category, $gene, \%targetLength);
			push (@stats, $stats);
			push (@presence, $presence);
			push (@sequence, $seq);
		}
	}
}


# Print results to report - make appropriate path, and open file
my $reportPath = $path . "/reports";
make_path($reportPath);
chdir("$reportPath");

open(REPORT, ">", "GeneSipprReport" . $start_time . ".csv") or die $!;


# Custom hash ordering courtesy of:
# http://stackoverflow.com/questions/8171528/in-perl-how-can-i-sort-hash-keys-using-a-custom-ordering
my @order = qw (O26 O45 O103 O111 O121 O145 O157 eae VT1 VT2 hlyA adk dinB fumC gyrB icd icdA mdh pabB polB purA putP recA trpA trpB uidA);
my %order_map = map { $order[$_] => $_ } 0 .. $#order;
my $pat = join '|', @order;


# Use @presence to properly format %targetPresence for the creation of a report
# $targetPresence{$fileName}{$category}{$geneName} = "+";
foreach my $component (sort @presence) {
	my %presence = %$component;
	foreach my $fName (sort keys %presence) {
		foreach my $cat (sort keys % {$presence{$fName}}) {
			while (my ($gName, $pres) = each ( % {$presence{$fName}{$cat}})) {
				$targetPresence{$fName}{$cat}{$gName} = $pres;
			}
		}
	}
}

# Similar to custom ordering above, but the quality and pathogen genes are separated into separate arrays
my @qualityArray = qw (adk dinB fumC gyrB icd icdA mdh pabB polB purA putP recA trpA trpB uidA);
my @pathoArray = qw (O26 O45 O103 O111 O121 O145 O157 eae VT1 VT2 hlyA);
my @pathoPres;
my $countQuality = 0;
print REPORT "Strain\t";
print REPORT join("\t", @order);
print REPORT "\tPathotype\t";
print REPORT "Quality Pass/Total\t";

foreach my $fName (sort keys %targetPresence) {
	print REPORT "\n$fName\t";
	foreach my $cat (sort keys % {$targetPresence{$fName}}) {
		$countQuality = 0;
		foreach my $gName (sort { my ($x, $y) = map /^($pat)/, $a, $b; $order_map{$x} <=> $order_map{$y}} keys %{ $targetPresence{$fName}{$cat} }) { # see URL above for some idea what's going on
			my $pres = $targetPresence{$fName}{$cat}{$gName};
			print REPORT "$pres\t";
			if ($pres eq "+" and $gName ~~ @qualityArray) {
				$countQuality++;
			} elsif ($pres eq "+" and $gName ~~ @pathoArray) {
				push(@pathoPres, $gName);
			}
		}
		if (scalar @pathoPres) {
			print REPORT join(";", @pathoPres);
		} else {
			print REPORT "NA"
		}
		print REPORT "\t$countQuality" . "/" . scalar @qualityArray . "\t";
		undef @pathoPres;
	}
}

close REPORT;

# Return the run time
my $end_time = time;
my $total_time = $end_time - $start_time;
my $legible_time = sprintf("%.1d", $total_time);

print "-------------------------------------------\n";
print "The total run time was $legible_time seconds.\n\n";

exit;

##########################################################
sub rawMapping {
	my ($path, $filesPath, $rawFile, $category, $gene) = @_;
	(my $geneName = $gene) =~ s/.fa//;
	(my $file = $rawFile) =~ s/_R1_001.fastq//;
	(my $rawFile1 = $rawFile) =~ s/R1_001/R2_001/;
	my $name = "$file" . "_$geneName";
	
	make_path("$path/results/tmp");
	chdir ("$path/results/tmp");
	
	# Perform smalt mapping
	unless (-e ("$path/results/tmp/$name.bam")) {
		system("smalt map -f bam -n 24 -o $path/results/tmp/$name.bam $filesPath/$category/$geneName $path/query/$rawFile");
	}

	# Sort the bam file
	unless (-e("$name" . "_sorted.bam")) {
		system ("samtools sort $name.bam $name" . "_sorted");
	}
	
	#Index the BAM file
	unless (-e("$name" . "_sorted.bai")) {
		system("samtools index $name" . "_sorted.bam $name" . "_sorted.bai");		
	}
	
	# Creates a vcf file from which all relevant sequence data/metadata can be extracted
	unless (-e("$name" . "_sorted.vcf")) {
		print "Creating vcf files\n";
		system("samtools mpileup -A -BQ0 -d 1000000 -uf $filesPath/$category/$gene $name" . "_sorted.bam | bcftools view -cg - > $name" . "_sorted.vcf")
	}
}

##########################################################
sub parseVCF {
	my ($rawFile, $category, $gene, $tLength) = @_;
	my %targetLength = %$tLength;
	# Initialise variables
	my ($count, $counts, $matchBoolean) = 0;
	my ($targetLength, $sequence, $match, $depth, $prevValue, $avgQual, $stdQual, $avgCov, $stdCov, $avgID, $lengthRatio);
	my (@quality, @coverage, @identity, @depth);
	my (%rawStats, %targetPresence, %sequence, %matchCount, %matchTally);
	my $curPath = getcwd;
	
	# Manipulate variable data for easy formatting
	(my $geneName = $gene) =~ s/.fa//;
	(my $file = $rawFile) =~ s/_R1_001.fastq//;
	(my $fileName = $file) =~ s/_.+//g;
	my $name = "$file" . "_$geneName";
	my $useSeq;
	# Parse the vcf file
	# Remove the index sequence and gene target which follow the strain name. These were added to the file names in the mapping process (e.g. 2014-SEQ-121_AAGAGGCA_adk)
	(my $noIndex = $name) =~ s/_.+//g;
	open(VCF, "<", "$name" . "_sorted.vcf");
	while (<VCF>) {
		if (/\#/) {
			next;
		} else {
			# Format of file
			### 35 lines of information on file format, samtools version, reference file used, contig ID, and abbreviation information
            # CHROM	    POS	ID	REF		ALT		QUAL 	FILTER	INFO	                                   FORMAT
            # adk-12	8	.	G		.		32.7	.		DP=1;AF1=0;AC1=0;DP4=0,1,0,0;MQ=29;FQ=-30.3	PL	0
            # data[0]  [1] 	[2] [3]  	[4] 	[5]    	[6]  	[7]
            # $a	   $pos	$b 	$refSeq $mapSeq $qual 	$c 		$info
			$count++;
			# Assign variables to the values split on tab
			(my $a, $pos, $a, $refSeq, $mapSeq, $qual, $a, $info) = split("\t");
			# For now, the only important data from the INFO category is located prior to the first semi-colon
			my @splitInfo = split(";", $info);
			# For now, I'm skipping lines that indicated the presence of a possible indel
            # - I may return to this later
			if ($splitInfo[0] =~ /INDEL/) {
				next;
			} else {
				# Depth of coverage is reported prior to the first ";"
				($depth = $splitInfo[0]) =~ s/DP=//;
				push(@depth, $depth);
				# If $pos > $count, then there is a gap in the mapping (or a deletion, but ignoring
	            # this possibility for now). For my purposes, I want to have data reported for
	            # every position, whether it is present in the vcf file or not, so I will use count
	            # as the position, "-" as the seq, and 0 as the quality and depth of coverage
				if ($pos > $count) {
					# the number of skipped positions is equal to the value for pos - count
	                # For each skipped position (i), set appropriate variables to appropriate values
	                for (my $range = $count; $range < $pos; $range++) {
	                	my $posAdj = $count;
	                    my $seqAdj = "-";
	                    my $matchAdj = 0;
	                    my $qualAdj = 0;
	                    my $DPAdj = 0;
	                    
	                    push(@quality, $qualAdj);
	                    push(@coverage, $DPAdj);
	                    push(@identity, $matchAdj);
	                    
	                    $useSeq = $seqAdj;
	                    $sequence .= $seqAdj;
	                    $count++;
	                    if ($pos == $count) {
	                    	# This match shouldn't technically be 1, but I don't want to recreate the entire 30 line loop right now
	                    	# I'll eventually make it into a subroutine, and chances are, 99.9% of the time, $mapSeq eq "."
	                    	$match = 1;
							$useSeq = $refSeq;
	                    	$sequence .= $refSeq;
	                    	push(@quality, $qual);
	                    	push(@coverage, $depth);
	                    	push(@identity, $match);
	                    }
	                }
				} else {
					# If the called base ($mapSeq) is identical to the reference base ($refSeq)
            		# - denoted by a ".", then append $refSeq to $sequence, and set $match to 1
					if ($mapSeq eq ".") {
						$useSeq = $refSeq;
						$sequence .= $refSeq;
						$match = 1;
					} else {
					# This section corrects for the fact that during the conversion of bam files to vcf
	                # files, SNP calls and ambiguous calls look identical, except for the fact that for
	                # SNPs, the qualityScore ($qual) tends to be higher than the surrounding bases,
	                # while ambiguous calls have a lower qualityScore - this loop screens for quality
	                # scores that are at least 10 lower than the score of the previous base
	                	if (($quality[-1] - 10) >= 0) {
	                		$prevValue = $quality[-1] - 10;
	                	} else {
	                		$prevValue = 0;
	                	}
	                	if ($qual <= $prevValue) {
	                		$useSeq = $refSeq;
	                		$sequence .= $refSeq;
	                		$match = 1;
	                	} else {
	                		# This attempts to catch if there are two ambiguous bases in a row;
	                    	# they will hopefully have the same value
	                    	if ($qual == $prevValue) {
	                    		$useSeq = $refSeq;
	                    		$sequence .= $refSeq;
	                    		$match = 1;
	                    	} else {
	                    		# "True" SNPs seem to have increased qualityScore compared to the
	                            # surrounding values, this will catch that
	                            if ($qual > $prevValue) {
	                            	if (length($mapSeq) > 1) {
	                					my $firstBase = split(",", $mapSeq);
	                					$useSeq = $firstBase;
	                					$sequence .= $firstBase;
	                					$match = 0;
	                				} else {
	                					$useSeq = $mapSeq;
	                            		$sequence .= $mapSeq;
	                            		$match = 0;
	                				}
	                            }
	                    	}
	                	}
					}
					# Populate the arrays with the appropriate variables and increment the count
					push(@quality, $qual);
	                push(@coverage, $depth);
	                push(@identity, $match);
				}
			}
		}
	}
	close VCF;
	
	# This loop checks to see if $pos exists - ensures that some mapping occurred - and subsequently if the mapping terminated early
	# It then adds an appropriate number of "-"s to $sequence and 0s to the arrays, so the arrays will have the same length as the target sequence
	if ($pos) {
		if ($pos < $targetLength{$geneName}) {
			for ($pos; $pos < $targetLength{$geneName}; $pos++) {
				$sequence .= "-";
	            push(@quality, 0);
	            push(@coverage, 0);
	            push(@identity, 0);
			}
		}
	} else {
		$pos = 0;
	}
	
	# In the case of no data being present in a file, set all the values to 0
	if ($count == 0) {
	   	$avgQual = 0;
        $stdQual = 0;
        $avgCov = 0;
        $stdCov = 0;
	    $avgID = 0;
	    $lengthRatio = 0;
	} else {
		# Calculate the average to two decimals (sprintf) as follows:
		# add all values from an appropriate array together (eval) and divide this number by the total length of the target sequence ($targetLength{$geneName})
		$avgQual = sprintf("%.2f", (eval join "+", @quality) / $targetLength{$geneName});
		# Calculate the standard deviation using the Statistics::Basic module
	    $stdQual = sprintf("%.2f", (stddev(@quality)));
	    $avgCov = sprintf("%.2f", (eval join "+", @coverage) / $targetLength{$geneName});
	    $stdCov = sprintf("%.2f", (stddev(@coverage)));
	    $avgID = sprintf("%.2f", ((eval join "+", @identity) / $targetLength{$geneName}) * 100);
	    # Calculate the ratio of the length of the mapped sequence to the length of the target
	    $lengthRatio = sprintf("%.2f", (length($sequence) / $targetLength{$geneName}) * 100);
	}
	
	# Populate the data hash
	$rawStats{$fileName}{$category}{$geneName}{$avgID}{$lengthRatio}{$avgCov} = $stdCov;
	
	# This array search for unique depth values - this shows that two or more unique reads have been mapped to 
	# the sequence. This will be used below in calculating the presence/absence of gene targets 
	my @uniqueDepth = uniq @depth;

	# Calculates the gap size
	foreach my $position (@coverage) {
		if ($position >= 1) {
			$matchBoolean = 0;
		} elsif ($position =~ /0/ and $matchBoolean == 0) {
			$counts++;
			$matchBoolean = 1;
			$matchCount{$fileName}{$counts}{$position}++;
		} else {
			$matchCount{$fileName}{$counts}{$position}++;
			}
		}
	# Calculates gap frequency
	foreach my $fName (sort keys %matchCount) {
		foreach my $counted (sort keys %{ $matchCount{$fName} }) {
			while (my ($pos, $size) = each (%{ $matchCount{$fName}{$counted} })) {
				$matchTally{$fName}{$size}++;
			}
		}
	}
	
	# Print the calculated gap sizes and gap frequencies to file
	open (MATCH, ">", $fileName . "_gaps.csv");
	print MATCH "$fileName,GapSize,GapFrequency\n"; 
	foreach my $fName (sort keys %matchTally) {
		while (my ($size, $frequency) = each ( %{ $matchTally{$fName} })) {
			print MATCH "$fName,$size,$frequency\n";
		}
	}
	
	close MATCH;

	# Creates a figure of the coverage of each target -uncomment if desired
	#system("/home/blais/git/RawReadSipping/gapGraphing.R $curPath $fileName");

	# Populate the sequence hash
	# Ensure that $sequence exists, otherwise, just indicate the lack of sequence as "" in the hash
	if ($sequence) {
		$sequence{$fileName}{$category}{$geneName} = $sequence;
	} else {
		$sequence{$fileName}{$category}{$geneName} = "";
	}
	
	# Populate the presence/absence hash
	# Using double cut-off now: either greater than 55% or there are two unique reads present in the vcf file
	if (($avgID < 55 and scalar @uniqueDepth > 1) or ($avgID > 55)) {
		$targetPresence{$fileName}{$category}{$geneName} = "+";
	} else {
		$targetPresence{$fileName}{$category}{$geneName} = "-";
	}
	
	# Return the hashes
	return(\%rawStats, \%targetPresence, \%sequence);
}

##########################################################
# This subroutine allows the printing of the local time each time it is invoked - it allows for the user to see how much time has passed
sub running_time {
	my $time =  localtime;
	my $hms = "[" . $time->hms . "] @_\n";
	print $hms;
}
