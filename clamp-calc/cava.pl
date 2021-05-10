#!/usr/bin/perl -w

# Copyright (C) 2015 Adam Ameur

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#########################################################################################
##
## File: cava.pl
##
## CAVA: Count-based analysis of variants in long amplicons
##
## Author: Adam Ameur
##
## Input: Two parameters are required
##
## 1) FASTQ or FASTA file containg full-length amplicon sequence data,
##    produced for example by PacBio sequencing
##
## 2) A pattern string containing the DNA sequence, with mutation indicated by '[]'
##    Example: TGGAACGCACGGACATCACC[A/G]TGAAGCACAAGCTGGGCGGG
##
## The number of allowed mismatches (sequence errors) in the pattern string is passed as
## an optional parameter.
##
## Output: All matching hits of wild-type and mutation sequences and their orientation
## within the reads are reported in the output. Optionally, detailed information about
## for each individual read is printed.
##
##########################################################################################

################
##
## Main program
##
################

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use cava;

# Read command line arguments
my ($in_file, $query, $filetype, $max_mismatch, $minq, $minqn, $output_read_ids, $qual_off, $hash_off, $silent) = get_args_and_error_check();

my %reads;
my %quals;

if($filetype eq "fastq"){
	my ($h1, $h2) = cava::read_fastq_file($in_file);
	%reads = %$h1;
	%quals = %$h2;
}
if($filetype eq "fasta"){
	%reads = cava::read_fasta_file($in_file);
}

my $use_qual=1;

if(defined($qual_off)){
	$use_qual=0;
}

if($filetype eq "fasta"){
	$use_qual=0;
}


# Execute counting of number of occurrences of query
if($query =~/^\s*[ACGT]+\s*$/){
	my %res = cava::count_matches(\%reads, \%quals, $query, $max_mismatch, $minq, $minqn, $output_read_ids, $use_qual, $hash_off, $silent);

	if(!$output_read_ids){
		print "$res{'query'}\t$res{'match_fwd'}\t$res{'match_rev'}\n";
	}
}

# Execute mutation matching for single base substitutions
if($query =~ /^\s*([ACGT]+)\[([ACGT])\/([ACGT])\]([ACGT]+)\s*$/){
	my %res = cava::screen_mutation(\%reads, \%quals, $query, $max_mismatch, $minq, $minqn, $output_read_ids, $use_qual, $hash_off, $silent);

	if(!$output_read_ids){
		print "$res{'mut'}\t$res{'ref_fwd'}\t$res{'alt_fwd'}\t$res{'other_fwd'}\t$res{'freq_fwd'}\t$res{'ref_rev'}\t$res{'alt_rev'}\t$res{'other_rev'}\t$res{'freq_rev'}\t$res{'ref_total'}\t$res{'alt_total'}\t$res{'other_total'}\t$res{'freq_total'}\n";
	}
}
else{
	# Execute mutation matching for indel substitutions
	if($query =~ /^\s*([ACGT]+)\[([ACGT-]*)\/([ACGT-]*)\]([ACGT]+)\s*$/){

		## TODO: Allow for quality score filtering and mismatches for indels
		print STDERR "Mismatch and quality score filtering not yet implemented for indels!\n";
		$max_mismatch=0;

		my $seq_up=$1;
		my $var1=$2;
		my $var2=$3;
		my $seq_down=$4;

		my $query1 = "$seq_up$var1$seq_down";
		my $query2 = "$seq_up$var2$seq_down";

		$query1 =~ s/-//g;
		$query2 =~ s/-//g;

		my $match_name_q1="ref";
		my $match_name_q2="alt";

		## TODO: send in an argument here to say if the read is reference(for query1) or alternative (query2)
		my %res1 = cava::count_matches(\%reads, \%quals, $query1, $max_mismatch, $minq, $minqn, $output_read_ids, $use_qual, $hash_off, $silent, $match_name_q1);
		my %res2 = cava::count_matches(\%reads, \%quals, $query2, $max_mismatch, $minq, $minqn, $output_read_ids, $use_qual, $hash_off, $silent, $match_name_q2);

		my $freq_fwd = 0;
		my $freq_rev = 0;
		my $freq = 0;

		if( ($res1{'match_fwd'}+$res2{'match_fwd'}) > 0 ){
			$freq_fwd=$res2{'match_fwd'}/($res1{'match_fwd'}+$res2{'match_fwd'});
		}
		if( ($res1{'match_rev'}+$res2{'match_rev'}) > 0 ){
			$freq_rev=$res2{'match_rev'}/($res1{'match_rev'}+$res2{'match_rev'});
		}
		if( ($res1{'match_total'}+$res2{'match_total'}) > 0){
			$freq=$res2{'match_total'}/($res1{'match_total'}+$res2{'match_total'});
		}

		if(!$output_read_ids){
			print "$query\t$res1{'match_fwd'}\t$res2{'match_fwd'}\t-1\t$freq_fwd\t$res1{'match_rev'}\t$res2{'match_rev'}\t-1\t$freq_rev\t$res1{'match_total'}\t$res2{'match_total'}\t-1\t$freq\n";
		}

	}
}



# Function for reading arguments and error handling
sub get_args_and_error_check{

	if (@ARGV == 0) {pod2usage(-exitval => 2, -verbose => 0);}

	my ($in_file, $mutation, $filetype, $max_mismatch, $minq, $minqn, $output_read_ids, $qual_off, $hash_off, $silent);

	my $result = GetOptions("--help"           => sub{local *_=\$_[1];
							                      pod2usage(-exitval =>2, -verbose => 1)},
		                    "-f=s"             =>\$in_file,
							"-q=s"             =>\$query,
							"-max_mm=i"        =>\$max_mismatch,
							"-t=s"             =>\$filetype,
							"-minq=i"          =>\$minq,
							"-minqn=i"         =>\$minqn,
							"-d!"              =>\$output_read_ids,
							"-nofilter!"       =>\$qual_off,
							"-h_off!"          =>\$hash_off,
		                    "-silent!"         =>\$silent)|| pod2usage(-exitval => 2, -verbose => 1);

	my $error_to_print;

	unless(defined($in_file)) {
		$error_to_print .= "\tNo input FASTA file specified.\n";
	}

    unless(defined($query)) {
		$error_to_print .= "\tNo mutation specified.\n";
	}
	else{
		my $query_ok=0;

		if($query =~ /^\s*([ACGT]+)\[([ACGT])\/([ACGT])\]([ACGT]+)\s*$/){
			$query_ok=1;
		}
		if($query =~ /^\s*([ACGT]+)\s*$/){
			$query_ok=1;
		}
		if($query =~ /^\s*([ACGT]+)\[([ACGT-]+)\/([ACGT-]+)\]([ACGT]+)\s*$/){
			$query_ok=1;
		}
		if($query_ok == 0){
			$error_to_print .= "\tQuery string (-q) not correctly formatted.\n";
		}
	}

	unless(defined($filetype)) {
		$filetype="fastq";
	}
	else{
		if(!( ($filetype eq "fastq") || ($filetype eq "fasta")) ){
			$error_to_print .= "\tFile type (-t) not correctly formatted, must be either 'fastq' of 'fasta'.\n";
		}

	}

	unless(defined($max_mismatch)) {
		$max_mismatch=5;
	}

	unless(defined($minq)) {
		$minq=0;
	}

	unless(defined($minqn)) {
		$minqn=0;
	}

	if(defined $error_to_print) {
		my $error_msg="ERROR(s):\n$error_to_print\n";
		pod2usage(-message => $error_msg, -exitval => 2, -verbose => 0);
	}

	else{
		return ($in_file, $query, $filetype, $max_mismatch, $minq, $minqn, $output_read_ids, $qual_off, $hash_off, $silent);
	}
}


__END__

=head1 NAME


cava.pl

=head1 SYNOPSIS

./cava.pl [options] B<--help> B<-f> B<-q> B<-t> B<-max_mm> B<-minq> B<-minqn> B<-d> B<-nofilter> B<-h_off> B<--silent>

=head1 OPTIONS

=over 8

=item [REQUIRED]

=item B<-f>

FASTA or FASTQ formatted file containing reads.

=item B<-q>

Query string to use for screening FASTA/FASTQ file.

=item [OPTIONAL]

=item B<-t>

Type of input file, either 'fastq' or 'fasta' (default = 'fastq').

=item B<-max_mm>

Maximum number of mismatches allowed in query string (default = 5).

=item B<-minq>

Minimum quality score at variant position. Only used for FASTQ input (default = 0).

=item B<-minqn>

Minimum quality score at neighboring positions, adjacent to variant. Only used for FASTQ input (default = 0).

=item B<-d>

Print detailed matching info for each read (default = no).

=item B<-nofilter>

If set, filtering on quality values is disabled.

=item B<-h_off>

If set, use of hash tables is reduced in the matching. This slows down the program but reduces memory consumption (default = hash on).

=item B<--silent>

Do not print status to stdout.

=back

=head1 DESCRIPTION

B<This program> performs a matching of mutations in amplicon sequence data.

=cut


