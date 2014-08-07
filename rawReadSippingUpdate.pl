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

# This start time will be used in calculating the total time of the run
my $start_time = time;

# Initialize variables
my ($sequenceName, $forwardReverse, $geneName, @folders, $fastaTitle, $reads, $investigator, $flowCell, $project);
my (%results, %sampleSheet, %sippr, %targets);

my $lane = 1;

# Determine the number of threads present in the system
my @cpus = `awk '/^processor/ { N++} END { print N }' /proc/cpuinfo`;
chomp @cpus;

# This section finds the most recent run folder on the MiSeq and determines which cycle the run is currently on
# NB: This script must be run after the run has been initialised, or the wrong folder will be identified as the current folder, and this will not work!
my $miSeqPath = "/media/miseq/MiSeqOutput";
chdir($miSeqPath);

# Grab all the folders
my @miSeqFolders = glob("*");

# Since the folders are all named starting with the date, they can be sorted (from highest to lowest) and the first folder will be the current run
@miSeqFolders = sort{$b cmp $a}(@miSeqFolders);
my $folder = $miSeqFolders[0];

# Get the flowcell ID from the end of the folder name
$folder =~ /.+_.+_.+_(\S+)/;
$flowCell = $1;

# As this run will not be performed connected to the NAS, the path will have to be modified
#my $folderPath = "/media/nas/akoziol/WGS/RawReadSipping/$folder";
my $folderPath = "/home/blais/Bioinformatics/RawReadSipping/$folder";



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
		chop $investigator;
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
			$reads = $_;
			# Get rid of \n
			chomp $reads;
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
			#my $index = $line[5];
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

close SHEET;
close OUTPUT;

# Find out how many cycles have been completed
chdir($miSeqPath . "/" . $folder . "/Thumbnail_images/L001");

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
my $baseMask = "Y" . $numReads . ",I8,I8,n*";

# Call configureBclToFastq.pl - in order to prevent the compression of the fastq files, I had to manually edit the Config.mk file in /usr/local/share/bcl2fastq-1.8.3/makefiles to not include the compression and compression suffix (commented out lines 174 and 175)
unless(-e("$folderPath/Unaligned")) {
	print "Calling script\n";
	system("configureBclToFastq.pl --input-dir $miSeqPath/$folder/Data/Intensities/BaseCalls/ --output-dir $folderPath/Unaligned --force --sample-sheet $folderPath/SampleSheet_modified.csv --mismatches 1 --no-eamss --fastq-cluster-count 0 --compression none --use-bases-mask $baseMask");
#	system("configureBclToFastqNoZIP.pl --input-dir $miSeqPath/$folder/Data/Intensities/BaseCalls/ --output-dir $folderPath/Unaligned --force --sample-sheet $folderPath/SampleSheet_modified.csv --mismatches 1 --no-eamss --fastq-cluster-count 0 --compression none --use-bases-mask $baseMask");
	chdir("$folderPath/Unaligned");
	system("nohup make -j 16 r1");
}

# Chdir to the working directory
#$project = "000000000-AAERG";
print "$folderPath/Unaligned/Project_$flowcell";
#chdir("$folderPath/Unaligned/Project_$project") or die $!;


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
			print "fastq files\n";
			print join("\n", @fastqFiles);
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
my @targetGenes = glob("*.fa");

foreach my $gene (@targetGenes) {
	# Perform smalt indexing of the target genes (if necessary) - it might eventually be advisable to play around with the k-mer and step sizes in this command
	(my $geneName = $gene) =~ s/.fa//;
	unless (-e ("$geneName.smi")) {
		system ("smalt index -k 5 -s 1 $geneName $gene");
	}
	# Create the index .fai file for use in mapping
	unless (-e ("$geneName.fai")) {
		system ("samtools faidx $gene");
	}
	# This creates a hash that will be used in determining the presence/absence of individual genes at the end of the analysis
	open(TARGETS, "<", "$gene");
	while (<TARGETS>) {
		if (/^\>/) {
			(my $targetTitle = $_) =~ s/>//;
			chomp $targetTitle;
			$targets{$geneName}{$targetTitle} = "+";
		}
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
				foreach my $gene (@targetGenes) {
				# Make a new thread using the blast subroutine - pass appropriate values to subroutine
				my $r = threads->new(\&rawMapping, $path, $filesPath, $file, $gene, \%targets);
				# Output data to @threads. Right now, all this does is increase the size for the scalar(@threads) portion of the while loop
				push(@threads,$r);
			}
		}
	}
	#my $time =  localtime;
	#my $date = $time->ymd;
	#my $hms = $time->hms("_");
	#$dateTime = $date.$hms;

	# This loop ensures that each thread is complete before terminating
	foreach (@threads) {
		my $num = $_->join;
		# As the values returned from the blast subroutine are a hash of hashes, $num must be treated as such - %$num
		# $results{$file}{$geneName}{$fastaTitle} = "+"
		foreach my $geneTarget (sort keys %$num) {
			foreach my $strain (sort keys %{ $$num{$geneTarget} }) {
				foreach my $fastaTitle (sort keys %{ $$num{$geneTarget}{$strain} }){
					chomp $fastaTitle;
					$sippr{$geneTarget}{$strain}{$fastaTitle} = $$num{$geneTarget}{$strain}{$fastaTitle}
				}
			}
		}
	}
	chdir ("$path");
	foreach my $a (sort keys %sippr) {
		foreach my $b (sort keys %{$sippr{$a}}) {
			foreach my $c (sort keys %{$sippr{$a}{$b}}) {
				print "$a\t$b\t$c\t$sippr{$a}{$b}{$c}\n";
			}
		}
	}

#print Dumper(%sippr);

#print Dumper(%targets);

#make_path("$path/reports");

# Return the run time
my $end_time = time;
my $total_time = $end_time - $start_time;
my $legible_time = sprintf("%.1d", $total_time);

print "-------------------------------------------\n";
print "The total run time was $legible_time seconds.\n\n";

exit;

###############################################################################
sub rawMapping {
	my ($path, $filesPath, $rawFile, $gene, $targets) = @_;
	(my $geneName = $gene) =~ s/.fa//;
	(my $file = $rawFile) =~ s/_R1_001.fastq//;
	my %targetGenes = %$targets;
	make_path("$path/results/tmp");
	chdir ("$path/query");
	my $name = "$file" . "_$geneName";
	# Perform smalt mapping
	unless (-e ("$path/results/tmp/$name.sam")) {
		system("smalt map -f samsoft -o $path/results/tmp/$name.sam -n 24 $filesPath/$geneName $rawFile");
	}
	chdir ("$path/results/tmp");
	# Use samtools view to convert sam to bam
	unless (-e("$name.bam")){
		system ("samtools view -b -S $name.sam -o $name.bam");
	}
	# Sort the bam file
	unless (-e("$name" . "_sorted.bam")) {
		system ("samtools sort $name.bam $name" . "_sorted");
	}
	# Produce a fastq file from the sorted bam file using the reference and these piped tools
	unless (-e("$name" . "_consensus.fastq")) {
		system ("samtools mpileup -uf $filesPath/$gene $name" . "_sorted.bam | bcftools view -cg - | vcfutils.pl vcf2fq > $name" . "_consensus.fastq");
	}
	# Use BioPerl to convert the fastq to fasta
	unless (-e("$name" . "_consensus.fa")) {
		my $seqin = Bio::SeqIO -> new (-format => 'fastq',-file => "$name" . "_consensus.fastq");
		my $seqout = Bio::SeqIO -> new (-format => 'fasta',-file => ">$name" . "_consensus.fa");
		while (my $seq_obj = $seqin -> next_seq) {
			$seqout -> write_seq($seq_obj);
		}
	}

	copy("$name" . "_consensus.fa", "$path/results/$name" . "_consensus.fa");
	open(FASTA, "<", "$name" . "_consensus.fa");
	while (<FASTA>) {
		if (/^\>/) {
			($fastaTitle = $_) =~ s/>//;
			chomp $fastaTitle;
			$results{$geneName}{$file}{$fastaTitle} = "+";
		}
		
	}


	foreach my $gName (sort keys %targetGenes) {
		foreach my $gTitle (sort keys %{ $targetGenes{$gName} }) {
			if (exists $results{$gName} and exists $results{$gName}{$file}{$gTitle}) {
				# %results -> strain name; length of reads; target gene series name; target gene name; presence
				$resultsOut{$gName}{$file}{$gTitle} = "+";
				#print "$gName $file $gTitle Found key!\n";
			} elsif (not exists $results{$gName}) {
				next;

			} else {
				$resultsOut{$gName}{$file}{$gTitle} = "-";
				#print "$gName $file $gTitle Key absent!\n";
			}
		}
	}


#print Dumper(%resultsOut);
	#chdir ("$path");
	return \%resultsOut;
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
