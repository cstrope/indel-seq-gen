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
my $Qidot_k;		### Value of interest, for full sequence.

my @rates_away;
my @proportion_different;
my @array;

for my $i (0 .. $nt) {
	$rates_away[$i] = [ @array ];
	$proportion_different[$i] = [ @array ];
	for my $j (0 .. $L) {
		$rates_away[$i][$j] = 0;
		$proportion_different[$i][$j] = 0;
	}
}

### For all time slices.
for my $dt (0 .. $nt-1) {
	my $time_to_go = ($nt-$dt)/$nt;
	my $exponential = exp( -(4.0/3.0) * $time_to_go );
	### Calculate PDjk, PIjk
	$PDjk = 0.25 - 0.25*$exponential;
	$PIjk = 0.25 + 0.75*$exponential;

	#print "$PDjk $PIjk\n";
	### For all possible numbers of differences
	for my $L_D (0 .. $L) {
		$L_I = $L-$L_D;
		$Qidot_k 
		=   0.25 * $L_D * (2 + $PIjk/$PDjk)		## Differing sites changing to either one of 2 differing states or ID state
		  + 0.25 * $L_I * (3*$PDjk/$PIjk);		## Identical sites changing to any of the 3 differing states

		$rates_away[$dt][$L_D] = $Qidot_k;

		$proportion_different[$dt][$L_D]
		= (0.25 * $L_D * (2 + $PIjk/$PDjk)) / $Qidot_k;
		
#		if ($L_D == 0) {
			print "Time: $time_to_go   #diff: $L_D     val = $Qidot_k    %diff = $proportion_different[$dt][$L_D]\n"; 
#		}


	}
}
