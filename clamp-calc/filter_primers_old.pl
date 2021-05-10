#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

# Read command line arguments
my ($in_file, $filetype, $fwd_primer, $rev_primer, $silent) = get_args_and_error_check();

## trim N bases from primer sequences...
my $nbases = 3;

## trimmed primer seqeunce must be find within fist/last M bases in reads
my $mbases = 20;

my $fwd_primer_trimmed = substr $fwd_primer, $nbases, length($fwd_primer);
my $rev_primer_trimmed = substr $rev_primer, 0, length($rev_primer)-$nbases;

my $fwd_primerRC = reverse_complement_IUPAC($fwd_primer);
my $rev_primerRC = reverse_complement_IUPAC($rev_primer);

my $fwd_primerRC_trimmed = reverse_complement_IUPAC($fwd_primer_trimmed);
my $rev_primerRC_trimmed = reverse_complement_IUPAC($rev_primer_trimmed);

my $fwd_primer_len=length($fwd_primer);
my $rev_primer_len=length($rev_primer);

my %reads;
my %quals;

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

	if($filetype eq "fasta"){
		## Make sure all matching reads starts and ends with primer sequences!
		if($sequence =~ /^[ACGT]{0,$mbases}$fwd_primer_trimmed([ACGT]+)$rev_primer_trimmed[ACGT]{0,$mbases}$/){
			my $trimmed_sequence = $1;
			print ">$header\n";
			print "$fwd_primer$trimmed_sequence$rev_primer\n";
		}
		if($sequence =~ /^[ACGT]{0,$mbases}$rev_primerRC_trimmed([ACGT]+)$fwd_primerRC_trimmed([ACGT]{0,$mbases})$/){
			my $trimmed_sequence = $1;
			print ">$header\n";
			print "$rev_primerRC$trimmed_sequence$fwd_primerRC\n";
		}
	}

	if($filetype eq "fastq"){
		$qual=$quals{$key};

		if(length($qual) != length($sequence)){
			die "Error!\n";
		}
		

		if($sequence =~ /^[ACGT]{0,$mbases}$fwd_primer_trimmed([ACGT]+)$rev_primer_trimmed[ACGT]{0,$mbases}$/){
			my $trimmed_sequence = $1;
			my $match_pos = -1;
			my $tmpseq=$sequence;
			my $tmpqual=$qual;
			if($tmpseq =~ /^([ACGT]{0,$mbases}$fwd_primer_trimmed)/g){
				$match_pos = pos($tmpseq);
			}
			my $match_len = length($trimmed_sequence);
			my $sequence_len = length($tmpseq);
			my $match_substr = substr($tmpseq, $match_pos, $match_len);
			my $qual_substr = substr($tmpqual, $match_pos, $match_len);
			my $first_qual = substr($qual_substr, 0, 1);
			my $last_qual = substr($qual_substr,-1);
			my $dummy_qual_start = $first_qual x $fwd_primer_len;
			my $dummy_qual_end = $last_qual x $rev_primer_len;
			print "\@$header\n";
			print "$fwd_primer$trimmed_sequence$rev_primer\n";
			print "+\n";
		    print "$dummy_qual_start$qual_substr$dummy_qual_end\n";
		}
		if($sequence =~ /^[ACGT]{0,$mbases}$rev_primerRC_trimmed([ACGT]+)$fwd_primerRC_trimmed([ACGT]{0,$mbases})$/g){
			my $trimmed_sequence = $1;
			my $match_pos = -1;
			my $tmpseq=$sequence;
			my $tmpqual=$qual;
			if($tmpseq =~ /^([ACGT]{0,$mbases}$rev_primerRC_trimmed)/g){
				$match_pos = pos($tmpseq);
			}
			my $match_len = length($trimmed_sequence);
			my $sequence_len = length($tmpseq);
			my $match_substr = substr($tmpseq, $match_pos, $match_len);
			my $qual_substr = substr($qual, $match_pos, $match_len);
			my $first_qual = substr($qual_substr, 0, 1);
			my $last_qual = substr($qual_substr,-1);
			my $dummy_qual_start = $first_qual x $rev_primer_len;
			my $dummy_qual_end = $last_qual x $fwd_primer_len;
			print "\@$header\n";
			print "$rev_primerRC$trimmed_sequence$fwd_primerRC\n";
			print "+\n";
		    print "$dummy_qual_start$qual_substr$dummy_qual_end\n";
		}
		
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
							"-fwd=s"             =>\$fwd_primer,
							"-rev=s"             =>\$rev_primer,
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
	
    unless(defined($fwd_primer)) {
		$fwd_primer = "TGACCAACTC";
	}

	unless(defined($rev_primer)) {
		$rev_primer = "ACGAAGTGGA";
	}

	if(defined $error_to_print) {
		my $error_msg="ERROR(s):\n$error_to_print\n";
		pod2usage(-message => $error_msg, -exitval => 2, -verbose => 0);
	}

	else{
		return ($ccs_file, $filetype, $fwd_primer, $rev_primer, $silent);
	}
}

__END__

=head1 NAME

filter_primers.pl

=head1 SYNOPSIS

./filter_primers.pl [options] B<--help> B<-f> B<-t> B<-fwd> B<-rev> B<--silent>

=head1 OPTIONS

=over 8

=item [REQUIRED]

=item B<-f>

CCS read input file (FASTA or FASTQ format).

=item [OPTIONAL]

=item B<-t>

Type of input file, either 'fastq' or 'fasta' (default = 'fastq').

=item B<-fwd>

Forward primer sequence (default: TGACCAACTC)

=item B<-rev>

Reverse primer sequence (default: ACGAAGTGGA)

=item B<--silent>

Do not print status to stdout.

=back

=head1 DESCRIPTION

B<This program> will filter all reads that match the specified primer sequences.

=cut


