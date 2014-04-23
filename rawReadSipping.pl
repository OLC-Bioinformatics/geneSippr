#usr/bin/perl

use warnings;
use strict;
use Cwd;
use Time::Piece;
use threads;
use File::Path qw(make_path remove_tree);
use File::Find;
use Bio::SeqIO;
use File::Copy;
use Data::Dumper qw(Dumper);

my $path = getcwd;

# This start time will be used in calculating the total time of the run
my $start_time = time;

# Initialize variables
my ($sequenceName, $forwardReverse, $geneName, @folders, $fastaTitle, $reads, $investigator, $flowCell);
my (%results, %sampleSheet);
my @lengths = ("50", "45", "40", "35", "30", "25", "20");

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

my $folderPath = "/media/nas/akoziol/WGS/RawReadSipping/$folder";

make_path($folderPath);

chdir($miSeqPath . "/" . $folder);
open(OUTPUT, ">", "$folderPath/SampleSheet_modified.csv") or die $!;
print OUTPUT "FCID\tLane\tSampleID\tSampleRef\tIndex\tDescription\tControl\tRecipe\tOperator\tSampleProject\n";
# Open the parameters file to determine the run length - this is necessary, as we will be pulling data once the first index has been processed (just over halfway through the run)
open(SHEET, "<", "SampleSheet.csv");
while (<SHEET>) {
	# Grab the investigator name
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
			chop $line[5]; chop $line[7]; chop $line[9]; chop $line[9];
			my $index = $line[5] . "-" . $line[7];
			print OUTPUT "$flowCell\t", 	# flowcell
						 "$lane\t",			# lane 1
						 "$line[0]\t",  	# sample ID
						 "no_ref\t",    	# no_ref sampleRef
						 "$index\t",    	# index
						 "$line[9]\t",  	# description
						 "N\t",      		# N control
						 "NA\t",		 	# NA recipe
						 "$investigator\t", # operator
						 "$line[8]\n";		# sample_project

		}
	}
}

=com
# The files must be in the "query" subfolder. This subfolder must only have sequences that you wish to examine, or the program won't be able to find them
chdir ("$path/query");


#Grab the fastq files for manipulations
my @sequenceFiles = glob("*.fastq");

###Note, this piece of code will only have to be used when generating truncated reads
# Use fastq_trimmer to trim the fastq files to the appropriate lengths
foreach my $sequenceFile (@sequenceFiles) {
	# Determine the forward and reverse fastq files based on name (eg. OLC797_R1_001.fastq vs OLC797_R2_001.fastq)
	if ($sequenceFile =~ /R1/) {
		$forwardReverse = "F";
	} elsif ($sequenceFile =~ /R2/) {
		$forwardReverse = "R";
	}
	($sequenceName = $sequenceFile) =~ s/_R\d.+//;
	foreach my $length (@lengths) {
		make_path("$path/query/$length/$sequenceName");
		unless (-e("$length/$sequenceName/$sequenceName" . "$forwardReverse" . "_$length.fastq")) {
			print "Trimming reads from $sequenceName direction $forwardReverse to $length nucleotides in length\n";
			system ("fastx_trimmer -Q33 -l $length -i $sequenceFile -o $length/$sequenceName/$sequenceName" . "$forwardReverse" . "_$length.fastq");
		}
	}
}


chdir ("$path/target");
my @targetGenes = glob("*.fa");

foreach my $gene (@targetGenes) {
	(my $geneName = $gene) =~ s/.fa//;
	unless (-e ("$geneName.smi")) {
		system ("smalt index -k 5 -s 1 $geneName $gene");
	}
	# Create the index .fai file for use in mapping
	unless (-e ("$geneName.fai")) {
		system ("samtools faidx $gene");
	}
}
##### smalt index -k 5 -s 1 eae eae.fa

foreach my $length (@lengths) {
	chdir ("$path/query/$length");
	opendir(DIR, "$path/query/$length") or die "can't find $path/query/$length: $!";
	while (defined(my $folder = readdir(DIR))) {
		# Ignore special files
		next if $folder =~ /^\.\.?$/;
		if (-d "$path/query/$length/$folder"){
			#print "$length/$folder\n";
			chdir ("$path/query/$length/$folder");
			my @fastqF = glob("*F_$length.fastq");
			my @fastqR = glob("*R_$length.fastq");
			foreach my $gene (@targetGenes) {
				chdir ("$path/query/$length/$folder");
				(my $geneName = $gene) =~ s/.fa//;
				make_path("$path/results/$geneName/$folder/$length");
				my $name = "$folder" . "_F_$geneName";
				# Perform smalt mapping
				unless (-e ("$path/results/$geneName/$folder/$length/$name.sam")) {
					system("smalt map -f samsoft -o $path/results/$geneName/$folder/$length/$name.sam -n 24 $path/target/$geneName @fastqF");
				}
				chdir ("$path/results/$geneName/$folder/$length");
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
					system ("samtools mpileup -uf $path/target/$gene $name" . "_sorted.bam | bcftools view -cg - | vcfutils.pl vcf2fq > $name" . "_consensus.fastq");
				}
				# Use BioPerl to convert the fastq to fasta
				unless (-e("$name" . "_consensus.fa")) {
					my $seqin = Bio::SeqIO -> new (-format => 'fastq',-file => "$name" . "_consensus.fastq");
					my $seqout = Bio::SeqIO -> new (-format => 'fasta',-file => ">$name" . "_consensus.fa");
					while (my $seq_obj = $seqin -> next_seq) {
	   					$seqout -> write_seq($seq_obj);
					}
				}
				make_path("$path/results/reports");
				copy("$name" . "_consensus.fa", "$path/results/reports/$name" . "_$length" . "_consensus.fa");
				open(FASTA, "<", "$name" . "_consensus.fa");
				while (<FASTA>) {
					if (/^\>/) {
						($fastaTitle = $_) =~ s/>//;
						chomp $fastaTitle;
						# %results -> strain name; length of reads; target gene series name; target gene name; presence
						$results{$folder}{$length}{$geneName}{$fastaTitle} = "+"
					}
				}
			}
		}
	}
}

#print Dumper(%results);

foreach my $strain (sort keys %results) {
	print "$strain\n";
	foreach my $nt (sort keys %{ $results{$strain} }) {
		print "$nt\t";
		foreach my $geneTarget (sort keys %{ $results{$strain}{$nt} }){
			foreach my $marker (sort keys %{ $results{$strain}{$nt}{$geneTarget} }) {
				print "$marker\t";
			}
		}
		print "\n";
	}
}
# Return the run time
my $end_time = time;
my $total_time = $end_time - $start_time;
my $legible_time = sprintf("%.1d", $total_time);

print "-------------------------------------------\n";
print "The total run time was $legible_time seconds.\n\n";

exit;


exit;

#################################################################
sub commas {
	my $sepchar = grep(/,/ => @_) ? ";" : ",";
	(@_ == 0) ? ''			:
	(@_ == 1) ? $_[0]		:
	(@_ == 2) ? join(", ", @_)	:
		    join("$sepchar ", @_[0 .. ($#_-1)], "$_[-1]");


}
