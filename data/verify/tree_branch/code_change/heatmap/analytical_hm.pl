#!/usr/bin/env perl
use strict;

##########
### Get values for heatmaps analytically, print them to file.
##########

my $L = 1000;		### Length of sequences ###
my ($L_I);			### Number of identical sites ###
my ($Pik, $Pjk);	### Transition probabilities ###
my ($PDjk, $PIjk);	### Transition probabilities for different/identical sites ###
my $nt = 1000;		### The number of time slices.
my $begin = 0	;	### Begin time.
my $Qidot_k;		### Value of interest, for full sequence.

my @rates_away;
my @proportion_different;
my @array;

for my $i ($begin .. $nt) {
	$rates_away[$i] = [ @array ];
	$proportion_different[$i] = [ @array ];
	for my $j (0 .. $L) {
		$rates_away[$i][$j] = 0;
		$proportion_different[$i][$j] = 0;
	}
}

### For all time slices.
for my $dt ($begin .. $nt-1) {
	my $time_to_go = (($nt-$dt)/$nt)*0.75;
	my $exponential = exp( -(4.0/3.0) * $time_to_go );
	### Calculate PDjk, PIjk
	$PDjk = 0.25 - 0.25*$exponential;
	$PIjk = 0.25 + 0.75*$exponential;

#	print "$PDjk $PIjk\n";
	### For all possible numbers of differences
	for my $L_D (0 .. $L) {
		$L_I = $L-$L_D;
		$Qidot_k 
		=   0.25 * $L_D * (2 + $PIjk/$PDjk)		## Differing sites changing to either one of 2 differing states or ID state
		  + 0.25 * $L_I * (3*$PDjk/$PIjk);		## Identical sites changing to any of the 3 differing states

		$rates_away[$dt][$L_D] = $Qidot_k;

		$proportion_different[$dt][$L_D]
		= (0.25 * $L_D * (2 + $PIjk/$PDjk)) / $Qidot_k;
		
#		if (
#			$L_D == 750 && $dt == 1 or
#			$L_D == 750 && $dt == 999 or
#			$L_D == 250 && $dt == 1 or
#			$L_D == 250 && $dt == 999 or
#			$L_D == 500 && $dt == 500 or
#			$L_D == 750 && $dt == 1 or
#			$L_D == 750 && $dt == 999 or
#			$L_D ==  50 && $dt ==  50
#		   ) {
#			print "" 
#				  . ($dt/1000) 
#				  . " $Qidot_k\n";	
#				  . " ld: $L_D $proportion_different[$dt][$L_D]  length_diff: " . ($proportion_different[$dt][$L_D]/$L_D) ."\n"; 
#		}

	}
}

open OUT, ">Pdiff";
print OUT "time\t";
print OUT "seq0";
for my $i (1 .. 1000) { print OUT "\tseq$i"; }
print OUT "\n";
for my $i (0 .. 1000) {
	print OUT "t$i\t";
	print OUT $proportion_different[$i][0];
	for my $j (1 .. 1000) {
		print OUT "\t" . $proportion_different[$i][$j];
	}
	print OUT "\n";
}
close OUT;

open OUT, ">seqR";
print OUT "time\t";
print OUT "seq0";
for my $i (1 .. 1000) { print OUT "\tseq$i"; }
print OUT "\n";
for my $i (0 .. 1000) {
	print OUT "t$i\t";
	print OUT $rates_away[$i][0];
	for my $j (0 .. 1000) {
		print OUT "\t" . $rates_away[$i][$j];
	}
	print OUT "\n";
}
close OUT;

open OUT, ">propdiff1.dat";
for my $i (0 .. 999) {
	print OUT "" . ($i/1000) . " $proportion_different[$i][1]\n";
}

open OUT, ">propdiff250.dat";
for my $i (0 .. 999) {
	print OUT "" . ($i/1000) . " $proportion_different[$i][250]\n";
}
open OUT, ">propdiff500.dat";
for my $i (0 .. 999) {
	print OUT "" . ($i/1000) . " $proportion_different[$i][500]\n";
}
open OUT, ">propdiff750.dat";
for my $i (0 .. 999) {
	print OUT "" . ($i/1000) . " $proportion_different[$i][750]\n";
}
open OUT, ">propdiff999.dat";
for my $i (0 .. 999) {
	print OUT "" . ($i/1000) . " $proportion_different[$i][999]\n";
}

