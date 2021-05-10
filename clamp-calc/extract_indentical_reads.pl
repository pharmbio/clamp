#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

# Read command line arguments
my ($in_file, $filetype, $min_qual, $silent) = get_args_and_error_check();

my %reads;
my %quals;

my %sequence_counts;

if($filetype eq "fastq"){
	my ($h1, $h2) = read_fastq_file($in_file);
	%reads = %$h1;
	%quals = %$h2;
}
if($filetype eq "fasta"){ 
	%reads=read_fasta_file($in_file);
}


foreach my $key (keys %reads) {
	my $header=$key;
	my $sequence=$reads{$key};
	my $qual="undefined";

	if($filetype eq "fastq"){
		$qual=$quals{$key};
	}

	if(!defined($unique_sequences{$sequnce})){
		$sequence_counts{$sequnce}=1;
	}
	else{
		$sequence_counts{$sequnce}++;
	}
	
}



sub reverse_complement_IUPAC {

	my $dna = shift;

	# reverse the DNA sequence
	my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

	return $revcomp;
}


# Function for reading a FASTA sequence file
sub read_fasta_file{
	my $fasta_file=shift;

	my %id2seq = ();
	my $id = '';

	open(IN_FILE, "<", $fasta_file) or die "cannot open: $!";

	while(<IN_FILE>){
		chomp;
		if($_ =~ /^>(.+)/){
			$id = $1;
		}
		else{
			if($_ =~ /^\s*([ACGT]+)\s*$/){
				$id2seq{$id} .= $1;
			}
			else{
				die "Input file not FASTA formatted!! Error on line:\n$_\n";
			}
		}
	}

	close(IN_FILE);

	return(%id2seq);
}


# Function for reading a FASTQ sequence file: TODO, fix this function!! Allow for either fastq or fasta as input!!
sub read_fastq_file {
	my $fastq_file=shift;

	my %id2seq = ();
	my %id2qual = ();
	my $id = '';
	my $line;

	
	open(IN_FILE, "<", $fastq_file) or die "cannot open: $!";

	while(<IN_FILE>){
		chomp;
		if($_ =~ /^@(.+)$/){
			$id = $1;

			$line=<IN_FILE>;
			if($line =~ /^\s*([ACGT]+)\s*$/){
				$id2seq{$id}=$line;
			}
			else{
				die "Input file not FASTQ formatted!! Expected DNA sequence on line:\n$1\n";
			}
			
			$line=<IN_FILE>;
			if(!($line =~ /^\+.*$/)){
				die "Input file not FASTQ formatted!! Expected '+' separator on line:\n$_\n";
			}

			$line=<IN_FILE>;
			if($line =~ /^\s*(.*)\s*$/){
				$id2qual{$id}=$line;
			}

			if( !(length($id2seq{$id}) == length($id2qual{$id}) )){
				die "Error in FASTQ file. Sequence and quality strings of different lengths for read id:\n$id\n";
			}
		}
		else{
			die "Input file not FASTQ formatted!! Expected read header on line:\n$_\n";
		}
	}

	close(IN_FILE);

	return(\%id2seq, \%id2qual);
}




# Argument and error handling

sub get_args_and_error_check{

	if (@ARGV == 0) {pod2usage(-exitval => 2, -verbose => 0);}

	my ($ccs_file, $filetype, $fwd_primer, $rev_primer, $silent);

	my $result = GetOptions("--help"           => sub{local *_=\$_[1];
							                      pod2usage(-exitval =>2, -verbose => 1)},
		                    "-f=s"               =>\$ccs_file,
							"-t=s"               =>\$filetype,
							"-min_qual=i"        =>\$min_qual,
							"-silent!"           =>\$silent)|| pod2usage(-exitval => 2, -verbose => 1);

	my $error_to_print;

	unless(defined($ccs_file)) {
		$error_to_print .= "\tNo ccs file specified.\n"
	}

	unless(defined($filetype)) {
		$filetype="fastq";
	}
	else{
		if(!( ($filetype eq "fastq") || ($filetype eq "fasta")) ){
			$error_to_print .= "\tFile type (-t) not correctly formatted, must be either 'fastq' of 'fasta'.\n";
		}
		
	}
	
	if(defined $error_to_print) {
		my $error_msg="ERROR(s):\n$error_to_print\n";
		pod2usage(-message => $error_msg, -exitval => 2, -verbose => 0);
	}

	else{
		return ($ccs_file, $filetype, $min_qual, $silent);
	}
}

__END__

=head1 NAME

count_indentical_reads.pl

=head1 SYNOPSIS

./count_identical_reads.pl [options] B<--help> B<-f> B<-t> B<-min_qual> B<--silent>

=head1 OPTIONS

=over 8

=item [REQUIRED]

=item B<-f>

CCS read input file (FASTA or FASTQ format).

=item [OPTIONAL]

=item B<-t>

Type of input file, either 'fastq' or 'fasta' (default = 'fastq').

=item B<-min_qual>

Minimum quality score used for filtering (optional)

=item B<--silent>

Do not print status to stdout.

=back

=head1 DESCRIPTION

B<This program> will filter all reads that match the specified primer sequences.

=cut


