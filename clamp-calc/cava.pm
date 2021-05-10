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


######################################################
##
## File: cava.pm
##
## Contents: Routines for the cava mutation caller
##
######################################################

package cava;

# Global hash variables. These are used for speeding up the matching code.
my %seen_substrings;  # Hash with number of mismatches in substrings of same length as pattern
my %unmatched_reads;  # Hash with all reads that do not contain any match to pattern

# Function that counts the number of matches to a specific query, allowing for a number of mismatches
# TODO: Add qual values in this function!
sub count_matches{
	my $reads_ref = shift;
	my $quals_ref = shift;
	my $query = shift;
	my $max_mismatch = shift;
	my $minq = shift;
	my $minqn = shift;
	my $output_read_ids = shift;
	my $use_qual = shift;
	my $hash_off = shift;
	my $silent = shift;
	my $match_name = shift;

	my $match_prefix="match";

	if(defined($match_name)){
		$match_prefix=$match_name;
	}

	my $query_rc = reverse_complement_IUPAC($query);

	my %reads=%$reads_ref;
	my %quals=%$quals_ref;

	my %updated_reads=%reads; # Hash to keep track of reads that remain to be screened

	%seen_substrings = (); # Make sure hashes are empty...
	%unmatched_reads = ();

	my %counts = ('match_fwd',0,'match_rev',0);

	## First count all perfect matches...
	foreach my $key (keys %reads) {
		my $mutation_result="undefined";
		my $match_sequence="undefined";
		my $qual="undefined";
		my $qual_sequence="undefined";
		my $nr_mismatches=-1;
		my $sequence = $reads{$key};

		if($use_qual == 1){
			$qual = $quals{$key};
		}

		if($sequence =~/($query)/g){
			$match_sequence="$1";
			my $match_pos = pos($sequence);
			my $match_len = length($match_sequence);
			if($use_qual == 1){
				$qual_sequence = substr($qual, $match_pos-$match_len, $match_len);
			}
			$mutation_result=$match_prefix."_fwd";
			$nr_mismatches=0;
		}

		if($sequence =~/($query_rc)/g){
			$match_sequence="$1";
			my $match_pos = pos($sequence);
			my $match_len = length($match_sequence);
			if($use_qual == 1){
				$qual_sequence = substr($qual, $match_pos-$match_len, $match_len);
			}
			$mutation_result=$match_prefix."_rev";
			$nr_mismatches=0;
		}


		## Report results and remove from hash if there is a match
		if(!($mutation_result eq "undefined")){
			## TODO: Here use average qual in region instead!!
			##if($use_qual == 1){
			## Extract quality scores at variable position and neighboring bases
			##my ($qUp, $qPos, $qDn) = extract_qual_values($qual_sequence, $mutation_startpos, $mutation_endpos);

			## Only report reads with quality values above thresholds
			##if(($qPos >= $minq) && ($qUp >= $minqn) && ($qDn >= $minqn)){

			#	$counts{$mutation_result}++;

			#	if($output_read_ids){
			#		print "$key\t$match_sequence\t$mutation_result\t$nr_mismatches\t$qual_sequence\t$qUp\t$qPos\t$qDn\n";
			#	}
			#}
			#}
			## Report results with no quality filtering
			#else{
			$counts{$mutation_result}++;

			if($output_read_ids){
				print "$key\t$match_sequence\t$mutation_result\t$nr_mismatches\n";
			}
			#}

			delete $updated_reads{$key};
		}
	}

	## Go through all remaing reads and count occurences of ref/alt by fuzzy matching...
	if($max_mismatch > 0){
		while (my ($key) = each(%updated_reads)) {

			my $mutation_result="undefined";
			my $match_sequence="undefined";
			my $qual_sequence="undefined";
			my $orientation="undefined";
			my $qual="undefined";
			my $nr_mismatches=-1;
			my $sequence = $updated_reads{$key};
			my $alignment_pos=-1;
			my $mutation_startpos = -1;
			my $mutation_endpos = -1;

			if($use_qual == 1){
				$qual = $quals{$key};
			}

			## Check if read is in a global hash of 'unmatching' reads. Execute the commands
			## only if this is not true.
			if( !(defined $unmatched_reads{$sequence}) ){
                ## make a fuzzy matching on fwd strand
				my %best_alignment_fwd = get_best_alignment($sequence, $query, $max_mismatch, $hash_off);
				my @mismatch_vec_fwd = (keys %best_alignment_fwd);
				my $min_mismatch = $mismatch_vec_fwd[0];

				if($min_mismatch <= $max_mismatch){
					$match_sequence = $best_alignment_fwd{$min_mismatch};
					$orientation = "fwd";
					$nr_mismatches = $min_mismatch;
				}

				## make a fuzzy matching on rev strand
				if($min_mismatch > $max_mismatch){
					my $sequence_rev = reverse_complement_IUPAC($sequence);

					my %best_alignment_rev = get_best_alignment($sequence_rev, $query, $max_mismatch, $hash_off, $qual);

					my @mismatch_vec_rev = (keys %best_alignment_rev);
					my $min_mismatch_rev = $mismatch_vec_rev[0];

					if($min_mismatch_rev <= $max_mismatch){
						$match_sequence = $best_alignment_rev{$min_mismatch_rev};
						$orientation = "rev";
						$min_mismatch = $min_mismatch_rev;
						$nr_mismatches = $min_mismatch;
					}
				}

				## If a match is found either on fwd or rev, go through hash and report all other reads
				## having the exact same matching sequence
				if(($nr_mismatches >= 0) && ($nr_mismatches <= $max_mismatch)){

					if($orientation eq "fwd"){
						$mutation_result = $match_prefix."_fwd";;
					}

					if($orientation eq "rev"){
						$mutation_result = $match_prefix."_rev";;
					}

					if($orientation eq "rev"){
						$match_sequence = reverse_complement_IUPAC($match_sequence);
					}


					## Loop through, report all matching reads and remove from hash
					foreach my $k (keys %updated_reads) {
						my $s = $reads{$k};

						#print "$match_sequence\n";

						if($s =~/$match_sequence/g){
							## TODO: add qual filtering here!!!

							$counts{$mutation_result}++;

							if($output_read_ids){
								print "$k\t$match_sequence\t$mutation_result\t$nr_mismatches\n";
							}

							delete $updated_reads{$k};
						}
					}
				}
				else{
					if(!(defined $hash_off)){
						$unmatched_reads{$sequence}=$nr_mismatches;
					}
					delete $updated_reads{$key};
				}
			}

			if(defined $updated_reads{$key}){
				delete $updated_reads{$key};
			}
		}
	}

	my $match_total=$counts{'match_fwd'}+$counts{'match_rev'};

	return('query',$query,'match_fwd',$counts{'match_fwd'},'match_rev',$counts{'match_rev'},'match_total',$match_total);

}



# Function that performs the mutation screening
sub screen_mutation{

	my $reads_ref = shift;
	my $quals_ref = shift;
	my $mutation = shift;
	my $max_mismatch = shift;
	my $minq = shift;
	my $minqn = shift;
	my $output_read_ids = shift;
	my $use_qual = shift;
	my $hash_off = shift;
	my $silent = shift;

	my %reads=%$reads_ref;
	my %quals=%$quals_ref;

	my %updated_reads=%reads; # Hash to keep track of reads that remain to be screened

	%seen_substrings = (); # Make sure hashes are empty...
	%unmatched_reads = ();

	my $seq_up = "undefined";
	my $seq_down = "undefined";
	my $mut = "undefined";
	my $wt = "undefined";
	my $seq_up_rc = "undefined";
	my $seq_down_rc = "undefined";
	my $mut_rc = "undefined";
	my $wt_rc = "undefined";
	my $len_up=-1;
	my $len_down=-1;
	my $len_wt=-1;
	my $len_mut=-1;
	my %counts = ('ref_fwd',0,'ref_rev',0,'alt_fwd',0,'alt_rev',0,'other_fwd',0,'other_rev',0);

	if($mutation =~ /^\s*([ACGT]+)\[([ACGT]*)\/([ACGT]*)\]([ACGT]+)\s*$/){
		$seq_up = $1;
		$wt = $2;
		$mut = $3;
		$seq_down = $4;
		$seq_up_rc = reverse_complement_IUPAC($seq_up);
		$seq_down_rc = reverse_complement_IUPAC($seq_down);
		$mut_rc = reverse_complement_IUPAC($mut);
		$wt_rc = reverse_complement_IUPAC($wt);
		$len_up=length $seq_up;
		$len_down=length $seq_down;
		$len_wt=length $wt;
		$len_mut=length $mut;
	}

	my $pat = "$seq_up$wt$seq_down";

	## First count all identical matches to wt/mut...
	foreach my $key (keys %reads) {
		my $mutation_result="undefined";
		my $match_sequence="undefined";
		my $qual="undefined";
		my $qual_sequence="undefined";
		my $nr_mismatches=-1;
		my $sequence = $reads{$key};
		my $mutation_startpos = -1;
		my $mutation_endpos = -1;

		if($use_qual == 1){
			$qual = $quals{$key};
		}

		if($sequence =~/($seq_up)([^$mut])($seq_down)/g){
			$match_sequence="$1$2$3";
			my $match_pos = pos($sequence);
			my $match_len = length($match_sequence);
			$mutation_startpos = length($1);
			$mutation_endpos = $match_len-length($3);
			if($use_qual == 1){
				$qual_sequence = substr($qual, $match_pos-$match_len, $match_len);
			}
			if($sequence =~/$seq_up$wt$seq_down/){
				$nr_mismatches = 0;
				$mutation_result = "ref_fwd";
			}
			else{
				$nr_mismatches = 1;
				$mutation_result = "other_fwd";
			}
		}
		if($sequence =~/($seq_down_rc)([^$mut_rc])($seq_up_rc)/g){
			$match_sequence="$1$2$3";
			my $match_pos = pos($sequence);
			my $match_len = length($match_sequence);
			$mutation_startpos = length($1);
			$mutation_endpos = $match_len-length($3);
			if($use_qual == 1){
				$qual_sequence = substr($qual, $match_pos-$match_len, $match_len);
			}
			if($sequence =~/$seq_down_rc$wt_rc$seq_up_rc/){
				$nr_mismatches = 0;
				$mutation_result = "ref_rev";
			}
			else{
				$nr_mismatches = 1;
				$mutation_result = "other_rev";
			}
		}
		if($sequence =~/($seq_up)($mut)($seq_down)/g){
			$match_sequence="$1$2$3";
			my $match_pos = pos($sequence);
			my $match_len = length($match_sequence);
			$mutation_startpos = length($1);
			$mutation_endpos = $match_len-length($3);
			if($use_qual == 1){
				$qual_sequence = substr($qual, $match_pos-$match_len, $match_len);
			}
			$mutation_result = "alt_fwd";
			$nr_mismatches = 0;
		}
		if($sequence =~/($seq_down_rc)($mut_rc)($seq_up_rc)/g){
			$match_sequence="$1$2$3";
			my $match_pos = pos($sequence);
			my $match_len = length($match_sequence);
			$mutation_startpos = length($1);
			$mutation_endpos = $match_len-length($3);
			if($use_qual == 1){
				$qual_sequence = substr($qual, $match_pos-$match_len, $match_len);
			}
			$mutation_result = "alt_rev";
			$nr_mismatches = 0;
		}

		## Report results and remove from hash if there is a match
		if(!($mutation_result eq "undefined")){

			## Perform quality filtering and report results
			if($use_qual == 1){
				## Extract quality scores at variable position and neighboring bases
				my ($qUp, $qPos, $qDn) = extract_qual_values($qual_sequence, $mutation_startpos, $mutation_endpos);

				## Only report reads with quality values above thresholds
				if(($qPos >= $minq) && ($qUp >= $minqn) && ($qDn >= $minqn)){

					$counts{$mutation_result}++;

					if($output_read_ids){
						print "$key\t$match_sequence\t$mutation_result\t$nr_mismatches\t$qual_sequence\t$qUp\t$qPos\t$qDn\n";
					}
				}
			}
			## Report results with no quality filtering
			else{
				$counts{$mutation_result}++;

				if($output_read_ids){
					print "$key\t$match_sequence\t$mutation_result\t$nr_mismatches\n";
				}
			}

			delete $updated_reads{$key};
		}
	}


	## Go through all remaing reads and count occurences of ref/alt by fuzzy matching...
	if($max_mismatch > 0){
		while (my ($key) = each(%updated_reads)) {

			my $mutation_result="undefined";
			my $match_sequence="undefined";
			my $qual_sequence="undefined";
			my $orientation="undefined";
			my $qual="undefined";
			my $nr_mismatches=-1;
			my $sequence = $updated_reads{$key};
			my $alignment_pos=-1;
			my $mutation_startpos = -1;
			my $mutation_endpos = -1;

			if($use_qual == 1){
				$qual = $quals{$key};
			}

			## Check if read is in a global hash of 'unmatching' reads. Execute the commands
			## only if this is not true.
			if( !(defined $unmatched_reads{$sequence}) ){
                ## make a fuzzy matching on fwd strand
				my %best_alignment_fwd = get_best_alignment($sequence, $pat, $max_mismatch, $hash_off);

				my @mismatch_vec_fwd = (keys %best_alignment_fwd);
				my $min_mismatch = $mismatch_vec_fwd[0];

				if($min_mismatch <= $max_mismatch){
					$match_sequence = $best_alignment_fwd{$min_mismatch};
					$orientation = "fwd";
					$nr_mismatches = $min_mismatch;
				}

				## make a fuzzy matching on rev strand
				if($min_mismatch > $max_mismatch){
					my $sequence_rev = reverse_complement_IUPAC($sequence);

					my %best_alignment_rev = get_best_alignment($sequence_rev, $pat, $max_mismatch, $hash_off, $qual);

					my @mismatch_vec_rev = (keys %best_alignment_rev);
					my $min_mismatch_rev = $mismatch_vec_rev[0];

					if($min_mismatch_rev <= $max_mismatch){
						$match_sequence = $best_alignment_rev{$min_mismatch_rev};
						$orientation = "rev";
						$min_mismatch = $min_mismatch_rev;
						$nr_mismatches = $min_mismatch;
					}
				}

				## If a match is found either on fwd or rev, go through hash and report all other reads
				## having the exact same matching sequence
				if(($nr_mismatches >= 0) && ($nr_mismatches <= $max_mismatch)){

					my $tmp_trimmed = substr $match_sequence, $len_up; # trim bases upstream
					my $trimmed_len = length($tmp_trimmed);
					my $variable_seq = substr $tmp_trimmed, 0, $trimmed_len-$len_down; # trim bases downstream

					if($orientation eq "fwd"){
						$mutation_startpos=$len_up;
						$mutation_endpos=$mutation_startpos+$len_wt;

						if($variable_seq eq $wt){
							$mutation_result = "ref_fwd";
						}
						if($variable_seq eq $mut){
							$mutation_result = "alt_fwd";
							$mutation_endpos=$mutation_startpos+$len_mut;
						}
						if($mutation_result eq "undefined"){
							$mutation_result = "other_fwd";
						}
					}
					if($orientation eq "rev"){
						$mutation_startpos=$len_down;
						$mutation_endpos=$mutation_startpos+$len_wt;

						if($variable_seq eq $wt){
							$mutation_result = "ref_rev";
						}
						if($variable_seq eq $mut){
							$mutation_result = "alt_rev";
							$mutation_endpos=$mutation_startpos+$len_mut;
						}
						if($mutation_result eq "undefined"){
							$mutation_result = "other_rev";
						}
					}

					if($orientation eq "rev"){
						$match_sequence = reverse_complement_IUPAC($match_sequence);
					}


					## Loop through, report all matching reads and remove from hash
					foreach my $k (keys %updated_reads) {
						my $s = $reads{$k};

						if($s =~/$match_sequence/g){

							if($use_qual == 1){ # Report matches, quality filtering on

								my $q = $quals{$k};
								my $m_pos = pos($s);
								my $m_len = length($match_sequence);
								my $q_sequence = substr($q, $m_pos-$m_len, $m_len);

								my ($qUp, $qPos, $qDn) = extract_qual_values($q_sequence, $mutation_startpos, $mutation_endpos);

								## Only report reads with quality values above thresholds
								if(($qPos >= $minq) && ($qUp >= $minqn) && ($qDn >= $minqn)){

									$counts{$mutation_result}++;

									if($output_read_ids){
										print "$k\t$match_sequence\t$mutation_result\t$nr_mismatches\t$q_sequence\t$qUp\t$qPos\t$qDn\n";
									}
								}
							}
							else{ # Report matches, quality filtering off
								$counts{$mutation_result}++;

								if($output_read_ids){
									print "$k\t$match_sequence\t$mutation_result\t$nr_mismatches\n";
								}
							}
							delete $updated_reads{$k};
						}
					}
				}
				else{
					if(!(defined $hash_off)){
						$unmatched_reads{$sequence}=$nr_mismatches;
					}
					delete $updated_reads{$key};
				}
			}

			if(defined $updated_reads{$key}){
				delete $updated_reads{$key};
			}
		}
	}

	my $freq_fwd = 0;
	if(($counts{'ref_fwd'})+$counts{'alt_fwd'}>0){
		$freq_fwd = $counts{'alt_fwd'}/($counts{'ref_fwd'}+$counts{'alt_fwd'}+$counts{'other_fwd'});
	}
	my $freq_rev = 0;
	if(($counts{'ref_rev'})+$counts{'alt_rev'}>0){
		$freq_rev = $counts{'alt_rev'}/($counts{'ref_rev'}+$counts{'alt_rev'}+$counts{'other_rev'});
	}
	my $ref_count_total=$counts{'ref_fwd'}+$counts{'ref_rev'};
	my $alt_count_total=$counts{'alt_fwd'}+$counts{'alt_rev'};
	my $other_count_total=$counts{'other_fwd'}+$counts{'other_rev'};
	my $freq = 0;

	if(($ref_count_total+$alt_count_total+$other_count_total)>0){
		$freq = $alt_count_total/($ref_count_total+$alt_count_total+$other_count_total);
	}


	return('mut',$mutation,'ref_fwd',$counts{'ref_fwd'},'alt_fwd',$counts{'alt_fwd'},'other_fwd',$counts{'other_fwd'},'freq_fwd',$freq_fwd,'ref_rev',$counts{'ref_rev'},'alt_rev',$counts{'alt_rev'},'other_rev',$counts{'other_rev'},'freq_rev',$freq_rev,'ref_total',$ref_count_total,'alt_total',$alt_count_total,'other_total',$other_count_total,'freq_total',$freq);
}



# Function for extracting quatilty values from FASTQ coded string at a site of a mutation
sub extract_qual_values{
	my $qual_string=shift;
	my $mut_startpos=shift;
	my $mut_endpos=shift;
	my $qual_score_offset=33;

	my @qualVec = unpack "C*", $qual_string;

	my $qUp = $qualVec[$mut_startpos-1];
	my $qPos = $qualVec[$mut_startpos];
	my $qDn = $qualVec[$mut_endpos];

	## Subtract -33 from QUAL values (standard QUAL encoding for PacBio CCS reads)
	$qUp = $qUp-$qual_score_offset;
	$qPos = $qPos-$qual_score_offset;
	$qDn = $qDn-$qual_score_offset;

	return($qUp, $qPos, $qDn);
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

# Function for reading a FASTQ sequence file
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



# Function for reading a FASTQ sequence file, and return top N reads with highest quality
sub read_fastq_file_topN {
	my $fastq_file=shift;
	my $topN=shift;

	my %id2seq = ();
	my %id2qualstring = ();
	my %id2qualvalue = ();
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
				$id2qualstring{$id}=$1;
				my @qual_vec = unpack "C*", $1;
				my $total_qual=0;
				foreach (@qual_vec) {
					$total_qual += $_;
				}

				my $avg_qual = $total_qual/@qual_vec;
				$id2qualvalue{$id} = $avg_qual;
			}

			if( !(length($id2seq{$id}) == length($id2qualstring{$id}) )){
				die "Error in FASTQ file. Sequence and quality strings of different lengths for read id:\n$id\n";
			}
		}
		else{
			die "Input file not FASTQ formatted!! Expected read header on line:\n$_\n";
		}
	}

	close(IN_FILE);

	#my @sorted_ids = sort ({$id2qualvalue{$a} <=> $id2qualvalue{$b}} keys %id2qualvalue);

	my @sorted_ids = keys (%id2qualvalue); #select "random" reads (top n in) instead of the ones with highest quality
		
	my %id2seq_topN=();
	my %id2qual_topN=();

	my $i=0;

	while (@sorted_ids && $i < $topN) {
		my $curr_readid = pop(@sorted_ids);
		$id2seq_topN{$curr_readid} = $id2seq{$curr_readid};
		$id2qual_topN{$curr_readid} = $id2qualstring{$curr_readid};
		$i++;
		#print STDERR "$id2seq_topN{$curr_readid}\n";
		#print STDERR "$id2qual_topN{$curr_readid}\n";
		#print STDERR "$id2qualvalue{$curr_readid}\n";
		$id2qual_topN{$curr_readid} = $id2qualstring{$curr_readid};
	}


	return(\%id2seq_topN, \%id2qual_topN);
}



# Function that reports the first position in a sequence that contains a match to a given
# pattern, allowing for a given number of mismatches.
sub get_best_alignment{
	my $seq = shift;
	my $pat = shift;
	my $max_mismatches = shift;
	my $hash_off = shift;
	my $qual = shift;

	my $wt_str="undefined";
	my $mut_str="undefided";

	my $len_in  = length $seq;
	my $len_pat = length $pat;

	my %best_match_sequence;
	my $best_alignment = "undefined";
	my $best_mm = $len_pat;

	for my $i ( 0 .. $len_in - $len_pat ) {

		#Extract k-mer of length = $len_pat from position $i
		my $substr = substr($seq, $i, $len_pat);

		my $current_mm = $len_pat;

		## Only make a new matching if substring has not been seen before!
		## Store all seen substrings in a global hash
		if(!(defined $seen_substrings{$substr})){
			$current_mm = ($pat ^ $substr) =~ tr/\0//c;
			if(!(defined $hash_off)){
				$seen_substrings{$substr} = $current_mm;
			}
		}
		else{
			$current_mm = $seen_substrings{$substr};
		}

		if($current_mm < $best_mm){
			$best_alignment = $substr;
			$best_mm = $current_mm;
			if($best_mm <= $max_mismatches){
				$best_match_sequence{$best_mm}=$best_alignment;
				return(%best_match_sequence);
			}
		}

	}

	$best_match_sequence{$best_mm}=$best_alignment;

	return(%best_match_sequence);
}


# Function to perform a reverse complement of an IUPAC sequence string
sub reverse_complement_IUPAC{

	my $dna = shift;

	# reverse the DNA sequence
	my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

	return $revcomp;
}


return 1;
