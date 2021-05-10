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


############################################################################
##
## File: cava_snpcaller.pl
##
## CAVA-SNPcaller: A method for de novo mutation screening in an amplicon,
##  based on the CAVA algortithm.
##
## Author: Adam Ameur
##
## Input: Two parameters are required
##
## 1) FASTQ or FASTA file containg full-length amplicon sequence data,
##    produced for example by PacBio sequencing
##
## 2) A reference sequence for the amplicon
##
#############################################################################

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
my ($in_file, $reference_file, $filetype, $max_mismatch, $minq, $minqn, $nreads, $offset, $output_read_ids, $qual_off, $hash_off, $silent) = get_args_and_error_check();

my %reads_orig;
my %quals_orig;

## Store subset of reads in new hash (at most N)
my %reads_subset;
my %quals_subset;

if($nreads <= 0){
	%reads_subset=%reads_orig;
	%quals_subset=%quals_orig;
}
else{
	## Extract top N highest quality reads from fastq file
	if($filetype eq "fastq"){
		my ($h1, $h2) = cava::read_fastq_file_topN($in_file, $nreads);

		%reads_subset = %$h1;
		%quals_subset = %$h2;
	}

	## Extract N random reads from fasta file (since qual is not available)
	if($filetype eq "fasta"){
		%reads_orig = cava::read_fasta_file($in_file);

		my @read_ids = keys %reads_orig;
		my $curr_readnr = 0;

		while (@read_ids && $curr_readnr < $nreads) {
			my $curr_readid = pop(@read_ids);
			$reads_subset{$curr_readid} = $reads_orig{$curr_readid};
			if(defined($quals_orig{$curr_readid})){
				$quals_subset{$curr_readid} = $quals_orig{$curr_readid};
			}
			$curr_readnr++;
		}
	}
}


my $use_qual=1;

if(defined($qual_off)){
	$use_qual=0;
}

if($filetype eq "fasta"){
	$use_qual=0;
}


my %reference = cava::read_fasta_file($reference_file);

my $ref_seq = "undefined";

## Use first sequence in FASTA file as reference
foreach my $key (keys %reference){
	if($ref_seq ne "undefined"){
		die "Only 1 entry allowed in FASTA reference file!";
	}
	$ref_seq = $reference{$key};
}

print "pos\tcontext\tcoverage\tref\talt1\talt2\talt3\tcov_fwd\tref_fwd\talt1_fwd\talt2_fwd\talt3_fwd\tcov_rev\tref_rev\talt1_rev\talt2_rev\talt3_rev\n";

my $reflen=length $ref_seq;

for(my $i=$offset; $i<($reflen-$offset); $i++){

	if(!$silent){
		my $lastpos = ($reflen-$offset-1);
		print STDERR "\r - Analyzing de novo variants (position $i of $lastpos)";
	}


	my $seq_up = substr($ref_seq, $i-$offset, $offset);
	my $base = substr($ref_seq, $i, 1);
	my $seq_dn = substr($ref_seq, $i+1, $offset);

	my @subst_vec = substitute_vector($base);

	my $wt_base = $base;
	my $alt1 = $subst_vec[0];
	my $alt2 = $subst_vec[1];
	my $alt3 = $subst_vec[2];

	my $reference_string = "$seq_up\[N\]$seq_dn";
	my $mutation_str1 = "$seq_up\[$wt_base\/$alt1\]$seq_dn";
	my $mutation_str2 = "$seq_up\[$alt2\/$alt3\]$seq_dn";

	# Execute mutation matching
	my %res1 = cava::screen_mutation(\%reads_subset, \%quals_subset, $mutation_str1, $max_mismatch, $minq, $minqn, $output_read_ids, $use_qual, $hash_off, $silent);

	my $wt_fwd = $res1{'ref_fwd'};
	my $alt1_fwd = $res1{'alt_fwd'};
	my $alt2_fwd = 0;
	my $alt3_fwd = 0;

	my $wt_rev = $res1{'ref_rev'};
	my $alt1_rev = $res1{'alt_rev'};
	my $alt2_rev = 0;
	my $alt3_rev = 0;

	if( ($res1{'other_fwd'}>0) || ($res1{'other_rev'}>0) ){

		my %res2 = cava::screen_mutation(\%reads_subset, \%quals_subset, $mutation_str2, $max_mismatch, $minq, $minqn, $output_read_ids, $use_qual, $hash_off, $silent);

		$alt2_fwd = $res2{'ref_fwd'};
		$alt3_fwd = $res2{'alt_fwd'};

		$alt2_rev = $res2{'ref_rev'};
		$alt3_rev = $res2{'alt_rev'};
	}

	my $wt_tot = $wt_fwd+$wt_rev;
	my $alt1_tot = $alt1_fwd+$alt1_rev;
	my $alt2_tot = $alt2_fwd+$alt2_rev;
	my $alt3_tot = $alt3_fwd+$alt3_rev;

	my $cov_fwd = $wt_fwd+$alt1_fwd+$alt2_fwd+$alt3_fwd;
	my $cov_rev = $wt_rev+$alt1_rev+$alt2_rev+$alt3_rev;

	my $cov_tot = $cov_fwd+$cov_rev;

	print "$i\t$reference_string\t$cov_tot\t$wt_base\t$alt1\t$alt2\t$alt3\t$cov_fwd\t$wt_fwd\t$alt1_fwd\t$alt2_fwd\t$alt3_fwd\t$cov_rev\t$wt_rev\t$alt1_rev\t$alt2_rev\t$alt3_rev\n";

}

if(!$silent){
	print STDERR " done\n";
}


sub substitute_vector{

	my $base = shift;

	my @subst_vec;

	if($base eq 'A'){
		@subst_vec = ('C','G','T');
	}
	if($base eq 'C'){
		@subst_vec = ('A','G','T');
	}
	if($base eq 'G'){
		@subst_vec = ('A','C','T');
	}
	if($base eq 'T'){
		@subst_vec = ('A','C','G');
	}

	return(@subst_vec);
}


# Function for reading arguments and error handling
sub get_args_and_error_check{

	if (@ARGV == 0) {pod2usage(-exitval => 2, -verbose => 0);}

	my ($in_file, $reference_file, $filetype, $max_mismatch, $minq, $minqn, $nreads, $offset, $output_read_ids, $qual_off, $hash_off, $silent);

	my $result = GetOptions("--help"           => sub{local *_=\$_[1];
							                      pod2usage(-exitval =>2, -verbose => 1)},
		                    "-f=s"             =>\$in_file,
							"-r=s"             =>\$reference_file,
							"-max_mm=i"        =>\$max_mismatch,
							"-t=s"             =>\$filetype,
							"-minq=i"          =>\$minq,
							"-minqn=i"         =>\$minqn,
							"-d!"              =>\$output_read_ids,
							"-n=i"             =>\$nreads,
							"-w=i"             =>\$offset,
							"-nofilter!"       =>\$qual_off,
							"-h_off!"          =>\$hash_off,
		                    "-silent!"         =>\$silent)|| pod2usage(-exitval => 2, -verbose => 1);

	my $error_to_print;

	unless(defined($in_file)) {
		$error_to_print .= "\tNo input FASTA/FASTQ file specified.\n";
	}

    unless(defined($reference_file)) {
		$error_to_print .= "\tNo input FASTA reference file specified.\n";
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

	unless(defined($nreads)) {
		$nreads=-1;
	}

	unless(defined($offset)) {
		$offset=20;
	}

	if(defined $error_to_print) {
		my $error_msg="ERROR(s):\n$error_to_print\n";
		pod2usage(-message => $error_msg, -exitval => 2, -verbose => 0);
	}

	else{
		return ($in_file, $reference_file, $filetype, $max_mismatch, $minq, $minqn, $nreads, $offset, $output_read_ids, $qual_off, $hash_off, $silent);
	}
}


__END__

=head1 NAME


cava.pl

=head1 SYNOPSIS

./cava_snpcaller.pl [options] B<--help> B<-f> B<-m> B<-t> B<-max_mm> B<-minq> B<-minqn> B<-n> B<-w> B<-nofilter> B<-h_off> B<--silent>

=head1 OPTIONS

=over 8

=item [REQUIRED]

=item B<-f>

FASTA or FASTQ formatted file containing reads.

=item B<-r>

FASTA file containing reference sequence.

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

=item B<-n>

Maximum number of reads to used for snp calling (default = all reads).

=item B<-w>

Window size of number of bases up- and downstream of position to be used in snp calling (default = 20).

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


