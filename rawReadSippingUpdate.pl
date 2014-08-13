#usr/bin/perl

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

# This start time will be used in calculating the total time of the run
my $start_time = time;

# Initialize variables
my ($sequenceName, $forwardReverse, $geneName, @folders, $fastaTitle, $reads, $investigator, $flowCell, $project, @stats, @presence, @sequence, @categories);
my (%results, %sampleSheet, %sippr, %targets, %targetLength, %targetCategories);

my $lane = 1;

# Determine the number of threads present in the system
my @cpus = `awk '/^processor/ { N++} END { print N }' /proc/cpuinfo`;
chomp @cpus;

# This section finds the most recent run folder on the MiSeq and determines which cycle the run is currently on
# NB: This script must be run after the run has been initialised, or the wrong folder will be identified as the current folder, and this will not work!
###
#my $miSeqPath = "/media/miseq/MiSeqOutput";
#my $miSeqPath = "/media/nas/akoziol/Raw_read_sipping/runData";
my $miSeqPath = "/media/nas";
chdir($miSeqPath);

# Grab all the folders
my @miSeqFolders = glob("*");

# Since the folders are all named starting with the date, they can be sorted (from highest to lowest) and the first folder will be the current run
@miSeqFolders = sort{$b cmp $a}(@miSeqFolders);
#my $folder = $miSeqFolders[0];
my $folder = "140325_M02466_0010_000000000-A8MVB";

# Get the flowcell ID from the end of the folder name
$folder =~ /.+_.+_.+_(\S+)/;
$flowCell = $1;

# As this run will not be performed connected to the NAS, the path will have to be modified
#my $folderPath = "/media/nas/akoziol/WGS/RawReadSipping/$folder";
#my $folderPath = "/media/nas/akoziol/Raw_read_sipping/Sandbox/$folder";
my $folderPath = "/home/blais/git/RawReadSipping/Sandbox/$folder";


# Make the appropriate folder in the required location - I'm not sure whether the permissions section at the end is necessary - it worked when I was trying to make it work, so, for now, it stays
make_path($folderPath, {owner => "blais", group => "blais", mode => 0777});

open(OUTPUT, ">", "$folderPath/SampleSheet_modified.csv") or die $!;
#system("chown blais:blais $folderPath/SampleSheet_modified.csv");
#system ("chmod 777 $folderPath/SampleSheet_modified.csv");
print OUTPUT "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n";

chdir($miSeqPath . "/" . $folder);

# Open the parameters file to determine the run length - this is necessary, as we will be pulling data once the first index has been processed (just over halfway through the run)
open(SHEET, "<", "SampleSheet.csv");
while (<SHEET>) {
	# Grab the investigator name - perhaps because the input file is a CSV, the chop function (x2) seems to be required to format the cells properly
	if (/Investigator Name/) {
		($investigator = $_) =~ s/Investigator Name,//;
		$investigator =~ s/,//g;
		chomp $investigator;
		#chop $investigator;
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
			$reads = $_;
			# Get rid of \n
			chomp $reads;
			$reads =~ s/,//g;
			# Exit the loop, as we have the value for reads1
			last;
		}
	}
	# As above - find [Data], then proceed reading through SHEET
	if (/Sample_ID/) {
		while (<SHEET>) {
			#print "$_";
			my @line = split(/,/);
			# $line[5], $line[7]: index 1 and 2, respectively. $line[9]: description
			###
			#chop $line[5]; chop $line[7];
			chop $line[9]; chop $line[9];
			$project = $line[8];
			my $index = $line[5];
			#my $index = $line[5] . "-" . $line[7];
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

close SHEET;
close OUTPUT;

# Find out how many cycles have been completed
chdir($miSeqPath . "/" . $folder . "/Thumbnail_Images/L001");

my @cycleNum = glob("*C*");

my $cycles = scalar @cycleNum;

# As the number of cycles required is the number of forward reads + the index(8)
my $readsNeeded = $reads + 8;

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
#my $baseMask = "Y" . $numReads . ",I8,n*,n*";
my $baseMask = "Y" . $numReads . "n*,I8,n*,n*";

# Call configureBclToFastq.pl - in order to prevent the compression of the fastq files, I had to manually edit the Config.mk file in /usr/local/share/bcl2fastq-1.8.3/makefiles to not include the compression and compression suffix (commented out lines 174 and 175)
unless(-e("$folderPath/Unaligned")) {
	print "Calling script\n";
	system("configureBclToFastq.pl --input-dir $miSeqPath/$folder/Data/Intensities/BaseCalls/ --output-dir $folderPath/Unaligned --force --sample-sheet $folderPath/SampleSheet_modified.csv --mismatches 1 --no-eamss --fastq-cluster-count 0 --compression none --use-bases-mask $baseMask");
	chdir("$folderPath/Unaligned");
	system("nohup make -j 16 r1");
}


# Chdir to the working directory
#$project = "000000000-AAERG";
chdir("$folderPath/Unaligned/Project_$project") or die $!;

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
			#(my $unzip = $files[0]) =~ s/.gz//g;
			#
			
			my @files = glob("*.gz");
#			if (scalar @files == 0) {
#				next;
#			}
			foreach my $fileZ (@files){
				(my $folder = $fileZ) =~ s/_.*//g;
				(my $unzip = $fileZ) =~ s/.gz//g;
				(my $filename = $fileZ) =~ s/_L001|.gz//g;
				#running_time("Now extracting $filename");
				#mkdir $folder;
				#move($file,"$folder/$file");
				system("gzip -d $fileZ");
				if (-f $unzip){unlink "$fileZ"};
				#move("$folder/$unzip","$folder/$filename");
			}
			my @fastqFiles = glob("*.fastq");			
			make_path("$path/query");
			(my $filename1 = $fastqFiles[0]) =~ s/_L001//g;
			copy("$path/$file/$fastqFiles[0]", "$path/query/$filename1");
			# Right now, the GeneSippr doesn't work with paired-end data. I left this in because maybe someday it will. Just uncomment the lines with $filename2
			
#			if (scalar @fastqFiles == 2) {
#				(my $filename1 = $fastqFiles[0]) =~ s/_L001//g;
#				#(my $filename2 = $files[1]) =~ s/_L001//g;
#				copy("$path/$file/$filename1", "$path/query/$filename1");
#				#copy("$path/$file/$filename2", "$path/query/$filename2");
#				chdir("$path");
#			}
#			elsif (scalar @fastqFiles == 1){
#				(my $filename = $fastqFiles[0]) =~ s/_L001//g;
#				copy("$path/$file/$filename", "$path/query/$filename");
#				chdir("$path");
#			}
		}
	}
}

print "Processing fastq files\n";
# Get the target files from the folder
#my $filesPath = "/media/nas/akoziol/WGS/RawReadSipping/target";
my $filesPath = "/home/blais/git/RawReadSipping/target";
chdir ("$filesPath");

# Get the names of the two folder with target sequences
opendir(DIR, $filesPath) or die "can't find $path: $!";
while (defined(my $file = readdir(DIR))) {
	# Ignore special files
	next if $file =~ /^\.\.?$/;
	# Ignore folders without sequence data
	next if $file =~ /holding/;
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
			system ("smalt index -k 20 -s 1 $geneName $gene");
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
	while (my ($category, $targetGenes) = each (%targetCategories)) {
			foreach my $gene (@$targetGenes) {	
			my ($stats, $presence, $seq) = parseVCF($file, $category, $gene, \%targetLength);
			push (@stats, $stats);
			push (@presence, $presence);
			push (@sequence, $seq);
		}
	}
}

print Dumper(\@sequence);

# Return the run time
my $end_time = time;
my $total_time = $end_time - $start_time;
my $legible_time = sprintf("%.1d", $total_time);

print "-------------------------------------------\n";
print "The total run time was $legible_time seconds.\n\n";

exit;

###############################################################################
sub rawMapping {
	my ($path, $filesPath, $rawFile, $category, $gene) = @_;
	(my $geneName = $gene) =~ s/.fa//;
	(my $file = $rawFile) =~ s/_R1_001.fastq//;
	make_path("$path/results/tmp");
	chdir ("$path/results/tmp");
	my $name = "$file" . "_$geneName";
	
	# Perform smalt mapping
	unless (-e ("$path/results/tmp/$name.bam")) {
		system("smalt map -f bam -o $path/results/tmp/$name.bam -n 24 $filesPath/$category/$geneName $path/query/$rawFile");
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
		system("samtools mpileup -A -BQ0 -d 1000000 -uf $filesPath/$gene $name" . "_sorted.bam | bcftools view -cg - > $name" . "_sorted.vcf")
	}
}

##########################################################
sub parseVCF {
	my ($rawFile, $category, $gene, $tLength) = @_;
	my %targetLength = %$tLength;
	# Initialise variables
	my $count = 0;
	my ($targetLength, $sequence, $match, $depth, $prevValue, $avgQual, $stdQual, $avgCov, $stdCov, $avgID, $lengthRatio);
	my (@quality, @coverage, @identity);
	my (%rawStats, %targetPresence, %sequence);
	
	# Manipulate variable data for easy formatting
	(my $geneName = $gene) =~ s/.fa//;
	(my $file = $rawFile) =~ s/_R1_001.fastq//;
	(my $fileName = $file) =~ s/_.+//g;
	my $name = "$file" . "_$geneName";
	
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
	                    $sequence .= $seqAdj;
	                    $count++;
	                    if ($pos == $count) {
	                    	# This match shouldn't technically be 1, but I don't want to recreate the entire 30 line loop right now
	                    	# I'll eventually make it into a subroutine, and chances are, 99% of the time, $mapSeq eq "."
	                    	$match = 1;
	                    	push(@quality, $qual);
	                    	push(@coverage, $depth);
	                    	push(@identity, $match);
	                    	$count++;
	                    }
	                }
				} else {
					# If the called base ($mapSeq) is identical to the reference base ($refSeq)
            		# - denoted by a ".", then append $refSeq to $sequence, and set $match to 1
					if ($mapSeq eq ".") {
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
	                		$sequence .= $refSeq;
	                		$match = 1;
	                	} else {
	                		# This attempts to catch if there are two ambiguous bases in a row;
	                    	# they will hopefully have the same value
	                    	if ($qual == $prevValue) {
	                    		$sequence .= $refSeq;
	                    		$match = 1;
	                    	} else {
	                    		# "True" SNPs seem to have increased qualityScore compared to the
	                            # surrounding values, this will catch that
	                            if ($qual > $prevValue) {
	                            	if (length($mapSeq) > 1) {
	                					my $firstBase = split(",", $mapSeq);
	                					$sequence .= $firstBase;
	                					$match = 0;
	                				} else {
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
	                $count++;
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
	
	# Populate the sequence hash
	# Ensure that $sequence exists, otherwise, just indicate the lack of sequence as "" in the hash
	if ($sequence) {
		$sequence{$fileName}{$category}{$geneName} = $sequence;
	} else {
		$sequence{$fileName}{$category}{$geneName} = "";
	}
	
	# Populate the presence/absence hash
	# Using a cut-off of 70% right now - this may be increased later
	if ($avgID > 70) {
		$targetPresence{$fileName}{$category}{$geneName} = "+";
	} else {
		$targetPresence{$fileName}{$category}{$geneName} = "-";
	}
	
	#Return the hashes
	return(\%rawStats, \%targetPresence, \%sequence);
}

##########################################################
# This subroutine allows the printing of the local time each time it is invoked - it allows for the user to see how much time has passed
sub running_time {
	my $time =  localtime;
	my $hms = "[" . $time->hms . "] @_\n";
	print $hms;
}

#################################################################
sub commas {
	my $sepchar = grep(/,/ => @_) ? ";" : ",";
	(@_ == 0) ? ''			:
	(@_ == 1) ? $_[0]		:
	(@_ == 2) ? join(", ", @_)	:
		    join("$sepchar ", @_[0 .. ($#_-1)], "$_[-1]");
}
