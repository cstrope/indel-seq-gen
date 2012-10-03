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
my @nucleotides = { "A", "C", "G", "T" };
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
for (my $i = 3; $i < @sequence; $i+= 3) {
	$i_codon = $sequence[$i] . $sequence[$i+1] . $sequence[$i+2];
	for my $j (0 .. 3) {
		if ($nucleotides[$j] ne $sequence[$i]) {
			$j_codon = $nucleotides[$j] . $sequence[$i+1] . $sequence[$i+2];
			$rate_away[$i][$j] = 0.25 * &tau_ij_eq($env[$i], $i_codon, $j_codon);
		}

exit(0);

		if ($nucleotides[$j] ne $sequence[$i+1]) {
			$j_codon = $sequence[$i] . $nucleotides[$j] .  $sequence[$i+2];
		}
#		$rate_away[$i+1][$j] =

		if ($nucleotides[$j] ne $sequence[$i+2]) {
			$j_codon = $sequence[$i+2] . $sequence[$i+1] . $nucleotides[$j] ;
		}
#		$rate_away[$i+2][$j] =
	}
}

exit(0);


sub tau_ij_eq
{
	my ($env, $i_codon, $j_codon) = @_;
	my $tau_ij;
	my $diffPji;
	my $diffP0ji;

### NEED THE SEQUENCE INDICES!!!!!!!!!!!!!!!!!!!!!!!!!!! Pass them into here. ###

	$diffPji = markov_ratio($env, $i, $j);
	
	return (log(tau_ij) / (1-1/$tau_ij));
}

for my $i (1 .. @codons-1) {
	print "doublet: $codons[$i-1]|$codons[$i] ----> $codon_index[$i-1]*64+$codon_index[$i]\n";
	print "P($codons[$i-1]|$codons[$i]) = $CCDS[$codon_index[$i-1]][$codon_index[$i]]";
	print "  $GRCh37[$codon_index[$i-1]][$codon_index[$i]]\n";

	
}

### Calculate the rate away for each sequence chunk.
my $tau_ij = 0;
for my $i (1 .. @codons-1) {
	$tau_ij += 
		  	   ( $GRCh37[$codon_index[$i-1]][$codon_index[$i]] 
			     * $GRCh37[$codon_index[$i-1]][$codon_index[$i]] 
			   ) /
	 		   (
	 		     $CCDS[$codon_index[$i-1]][$codon_index[$i]] 
	 		     * $CCDS[$codon_index[$i-1]][$codon_index[$i]]
	 		   );
}

my $Rij = 0.25 * (log($tau_ij)/(1-1/$tau_ij));
print "$Rij\n";

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