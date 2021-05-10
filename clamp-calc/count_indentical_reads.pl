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
my %sequence_counts_rc;

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
	my $sequence_rc=reverse_complement_IUPAC($sequence);
	my $qual="undefined";
	my $qual_ok=1;
	
	my $sequence_match=0;

	if($filetype eq "fastq"){
		if($min_qual >= 1){
			$qual=$quals{$key};
			$qual_ok = check_min_qual($qual, $min_qual);
		}
	}

	if($qual_ok == 1){

		if(defined($sequence_counts{$sequence})){
			if(!defined($sequence_counts_rc{$sequence_rc})){
				die "Error1!!!";
			}                        
			$sequence_counts{$sequence}++;
			$sequence_match=1;
		}
                
		if(defined($sequence_counts_rc{$sequence})){
			if(!defined($sequence_counts{$sequence_rc})){
				die "Error2!!!";
			}
			$sequence_counts_rc{$sequence}++;
			$sequence_match=1;
		}
                
		if($sequence_match == 0){ ## New sequence, not seen before
			$sequence_counts{$sequence}=1;
			$sequence_counts_rc{$sequence_rc}=0;
		}
		
	}
}


for my $seq (sort { $sequence_counts{$b} <=> $sequence_counts{$a} } keys %sequence_counts) {
	my $seq_rc = reverse_complement_IUPAC($seq);
	my $count = $sequence_counts{$seq};
	my $count_rc= $sequence_counts_rc{$seq_rc};
	if( ($count >= 1) && ($count_rc >= 1) ){
		my $count_total = $count+$count_rc;
		print ">count:$count_total fwd:$count revcomp:$count_rc\n$seq\n";
	}
}



# Returns 1 if all quality values in a string is above some threshold, otherwise 0.
sub check_min_qual{
	my $qual_string=shift;
	my $min_qual_value=shift;
	my $qual_score_offset=33;

	my @qualVec = unpack "C*", $qual_string;

	my $qualOk=1;

	for my $i (0 .. $#qualVec-1){
		my $q = $qualVec[$i]-$qual_score_offset;
		if($q < $min_qual_value){
			$qualOk=0;
		}
	}

	return($qualOk);
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
				$id2seq{$id}=$1;
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
				$id2qual{$id}=$1;
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

	unless(defined($min_qual)) {
		$min_qual=0;
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

count_identical_reads.pl

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

B<This program> will find all reads that are identical and count their occurrences.

=cut


