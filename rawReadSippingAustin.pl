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

my $path = getcwd;

# This start time will be used in calculating the total time of the run
my $start_time = time;

# Initialize variables
my ($sequenceName, $forwardReverse, $geneName, @folders, $fastaTitle);
my (%results);
#my @lengths = ("50", "45", "40", "35", "30", "25", "20");

# Determine the number of threads present in the system
my @cpus = `awk '/^processor/ { N++} END { print N }' /proc/cpuinfo`;
chomp @cpus;

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


# The files must be in the "query" subfolder. This subfolder must only have sequences that you wish to examine, or the program won't be able to find them
chdir ("$path/query");

#Grab the fastq files for manipulations
my @sequenceFiles = glob("*R1_001.fastq");

##### smalt index -k 5 -s 1 eae eae.fa

foreach my $file (@sequenceFiles) {
	foreach my $gene (@targetGenes) {
		chdir ("$path/query");
		(my $geneName = $gene) =~ s/.fa//;
		(my $folder = $file) =~ s/.fastq//;
		make_path("$path/results/$folder/$geneName");
		my $name = "$folder" . "_$geneName";
		# Perform smalt mapping
		unless (-e ("$path/results/$folder/$geneName/$name.sam")) {
			system("smalt map -f samsoft -o $path/results/$folder/$geneName/$name.sam -n 24 $path/target/$geneName $file");
		}
		chdir ("$path/results/$folder/$geneName");
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
		copy("$name" . "_consensus.fa", "$path/results/reports/$name" . "_consensus.fa");
		open(FASTA, "<", "$name" . "_consensus.fa");
		while (<FASTA>) {
			if (/^\>/) {
				($fastaTitle = $_) =~ s/>//;
				chomp $fastaTitle;
				# %results -> strain name; length of reads; target gene series name; target gene name; presence
				$results{$folder}{$geneName}{$fastaTitle} = "+"
			}
		}
	}
}

#print Dumper(%results);

foreach my $strain (sort keys %results) {
	print "$strain\n";
	foreach my $geneTarget (sort keys %{ $results{$strain} }) {
		print "$geneTarget\t";
		foreach my $marker (sort keys %{ $results{$strain}{$geneTarget} }){
			print "$marker\t";
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
