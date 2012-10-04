#!/usr/bin/env perl
use strict;
$/=undef;

## Given a length 15 string, assumed to be codon, this will output all of the appropriate probabilities
## for each pair of codons (dependencies) for both neutral and dependencies, the probability of the
## triplet in first order markov dependence, and the rate away from the string.

if (@ARGV != 1) { die "Usage: perl check_string_tau-ij.pl <string_of_codons>\n"; }

my (@CCDS, @GRCh37);

my @array;
for my $i (0 .. 63) {
	$CCDS[$i] = [ @array ];
	$GRCh37[$i] = [ @array ];
	for my $j (0 .. 63) {
		$CCDS[$i][$j] = 0;
		$GRCh37[$i][$j] = 0;
	}
}

open IN, "CCDS_nucleotide-1OMm.current.dep.norm";
my $slurp = <IN>;
close IN;

my @tmp = split /\n/, $slurp;
for my $i (0 .. @tmp-1) {
	my @tmp2 = split /\s+/, $tmp[$i];
	# First element is triplets. Not interested...
	for my $j (1 .. @tmp2-1) {
		$CCDS[$i][$j-1] = $tmp2[$j];
	}
}

open IN, "GRCh37_non-coding-1OMm.dep.norm";
my $slurp = <IN>;
close IN;

my @tmp = split /\n/, $slurp;
for my $i (0 .. @tmp-1) {
	my @tmp2 = split /\s+/, $tmp[$i];
	# First element is triplets. Not interested...
	for my $j (1 .. @tmp2-1) {
		$GRCh37[$i][$j-1] = $tmp2[$j];
	}
}

#for my $i (0 .. 63) {
#	print "$i: " ;
#	for my $j (0 .. 63) {
#		print " $j $GRCh37[$i][$j]\n";
#	}
#	print "\n";
#}


my $pasted_seq = shift;
print "$pasted_seq\n";

my @env;
my @sequence = split //, $pasted_seq;
my (@codons, @codon_index);

## Create this to keep track of rates away per position.
my (@rate_away);
for my $i (0 .. @sequence-1) { $rate_away[$i] = [ @array ];	for my $j (0 .. 3) { $rate_away[$i][$j] = 0; } }
my @nucleotides = ( "A", "C", "G", "T" );
my ($j_codon, $i_codon);


for (my $i = 0; $i < @sequence; $i+=3) {
	$codons[$i/3] = $sequence[$i] . $sequence[$i+1] . $sequence[$i+2];
	print "codon " . ($i/3) . ": $codons[$i/3]\n";
	$codon_index[$i/3] = &get_index($codons[$i/3]);
	$env[$i] = $env[$i+1] = $env[$i+2] = 1;
}
$env[0] = $env[1] = $env[2] = 0; $env[@sequence-3] = $env[@sequence-2] = $env[@sequence-1] = 2;

## Sequence and env set up just like iSG.
for my $i (0 .. @sequence-1) {
	print "$sequence[$i] $env[$i]\n";
}


## See if we can make the rates match with iSG...
my ($codon_minus1, $codon_plus1);
for (my $i = 3; $i < @sequence; $i+= 3) {
	$i_codon = $sequence[$i] . $sequence[$i+1] . $sequence[$i+2];
	$codon_minus1 = $sequence[$i-3] . $sequence[$i-2] . $sequence[$i-1];
	if ($i+4 < @sequence) { $codon_plus1 = $sequence[$i+3] . $sequence[$i+4] . $sequence[$i+5]; }
	else { $codon_plus1 = ""; }
	for my $j (0 .. 3) {
		if ($nucleotides[$j] ne $sequence[$i]) {
			$j_codon = $nucleotides[$j] . $sequence[$i+1] . $sequence[$i+2];
			if (!&stop_codon($j_codon)) {
				$rate_away[$i][$j] = 0.25 * &tau_ij_eq($env[$i], $i_codon, $j_codon, $codon_minus1, $codon_plus1, $env[$i]);
			}
		}
	}

	for my $j (0 .. 3) {
		if ($nucleotides[$j] ne $sequence[$i+1]) {
			$j_codon = $sequence[$i] . $nucleotides[$j] .  $sequence[$i+2];
			if (!&stop_codon($j_codon)) {
				$rate_away[$i+1][$j] = 0.25 * &tau_ij_eq($env[$i], $i_codon, $j_codon, $codon_minus1, $codon_plus1, $env[$i]);
			}
		}
	}

	for my $j (0 .. 3) {
		if ($nucleotides[$j] ne $sequence[$i+2]) {
			$j_codon = $sequence[$i] . $sequence[$i+1] . $nucleotides[$j] ;
			if (!&stop_codon($j_codon)) {
				$rate_away[$i+2][$j] = 0.25 * &tau_ij_eq($env[$i], $i_codon, $j_codon, $codon_minus1, $codon_plus1, $env[$i]);
			}
		}
	}

	print "|$codon_minus1|*$i_codon*|$codon_plus1|\n";
	print "$sequence[$i]: ";
	for my $j (1 .. 3) {
		$rate_away[$i][$j] += $rate_away[$i][$j-1];
		$rate_away[$i+1][$j] += $rate_away[$i+1][$j-1];
		$rate_away[$i+2][$j] += $rate_away[$i+2][$j-1];
	}
	for my $j (0 .. 3) {
		print " $rate_away[$i][$j]";
	}
	print "\n";
	print "$sequence[$i+1]: ";
	for my $j (0 .. 3) {
		print " $rate_away[$i+1][$j]";
	}
	print "\n";
	print "$sequence[$i+2]: ";
	for my $j (0 .. 3) {
		print " $rate_away[$i+2][$j]";
	}
	print "\n";

}

my $sequence_rate_away = 0;
for my $i (0 .. @sequence-1) {
	$sequence_rate_away += $rate_away[$i][3];
}
print "SEQ_RATE_AWAY: $sequence_rate_away\n";

exit(0);

sub stop_codon
{
	my ($seq) = @_;
	
	if ($seq eq "TAG") { return 1; }
	if ($seq eq "TGA") { return 1; }
	if ($seq eq "TAA") { return 1; }
	return 0;
}

sub tau_ij_eq
{
	my ($env, $i_codon, $j_codon, $codon_left, $codon_right) = @_;
	my $tau_ij;
	my $diffPji;
	my $diffP0ji;

	my ($idx_codon_i, $idx_codon_j, $idx_codon_left, $idx_codon_right);
	$idx_codon_i = &get_index($i_codon);
	$idx_codon_j = &get_index($j_codon);
	$idx_codon_left = &get_index($codon_left);
	$idx_codon_right = &get_index($codon_right);

#	print "Pij\n";
#	print " i-1 i   |$codon_left|$i_codon| = $CCDS[$idx_codon_left][$idx_codon_i]\n";
#	if($env == 1) { print " i   i+1 |$i_codon|$codon_right| = $CCDS[$idx_codon_i][$idx_codon_right]\n"; }
#	print " j-1 j   |$codon_left|$j_codon| = $CCDS[$idx_codon_left][$idx_codon_j]\n";
#	if($env == 1) { print " j   j+1 |$j_codon|$codon_right| = $CCDS[$idx_codon_j][$idx_codon_right]\n"; }

#	print "P0ij\n";
#	print " i-1 i   |$codon_left|$i_codon| = $GRCh37[$idx_codon_left][$idx_codon_i]\n";
#	if($env == 1) { print " i   i+1 |$i_codon|$codon_right| = $GRCh37[$idx_codon_i][$idx_codon_right]\n"; }
#	print " j-1 j   |$codon_left|$j_codon| = $GRCh37[$idx_codon_left][$idx_codon_j]\n";
#	if($env == 1) { print " j   j+1 |$j_codon|$codon_right| = $GRCh37[$idx_codon_j][$idx_codon_right]\n"; }
	
	$diffPji = markov_ratio($env, $idx_codon_left, $idx_codon_right, $idx_codon_i, $idx_codon_j, \@CCDS);
	$diffP0ji = markov_ratio($env, $idx_codon_left, $idx_codon_right, $idx_codon_i, $idx_codon_j, \@GRCh37);

#	print "diffPji:  $diffPji\n";
#	print "diffP0ji: $diffP0ji\n";

	$tau_ij = $diffPji / $diffP0ji;

#	print "tau_ij = diffPji/diffP0ji = $tau_ij\n";
	
	return (log($tau_ij) / (1-1/$tau_ij));
}



for my $i (1 .. @codons-1) {
	print "doublet: $codons[$i-1]|$codons[$i] ----> $codon_index[$i-1]*64+$codon_index[$i]\n";
	print "P($codons[$i-1]|$codons[$i]) = $CCDS[$codon_index[$i-1]][$codon_index[$i]]";
	print "  $GRCh37[$codon_index[$i-1]][$codon_index[$i]]\n";
}

sub markov_ratio
{
	my ($env, $left, $right, $i, $j, $matrix_ref) = @_;
	my $ratio;
	
	if ($env == 1) {
		$ratio = ( $matrix_ref->[$left][$j] * $matrix_ref->[$j][$right] ) / ( $matrix_ref->[$left][$i] * $matrix_ref->[$i][$right] );
	} else {
		## This is environment at end of the sequence. We assumed that beginning of sequence is all 0 rates away, so
		## this is safe.
		$ratio = $matrix_ref->[$left][$j] / $matrix_ref->[$left][$i];
	}
	return $ratio;
}

sub get_index
{
	my $codon = shift;
	my @c = split //, $codon;
	my $index = 16 * &nucl_val($c[0]) + 4 * &nucl_val($c[1]) + &nucl_val($c[2]);
}

sub nucl_val
{
	my $nucl = shift;
	
	if ($nucl eq "A") { return 0; }
	if ($nucl eq "C") { return 1; }
	if ($nucl eq "G") { return 2; }
	if ($nucl eq "T") { return 3; }
}