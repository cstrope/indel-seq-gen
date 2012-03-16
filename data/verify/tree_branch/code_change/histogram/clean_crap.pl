#!/usr/bin/env perl
$/=undef;

my $burn_in = 50000;
my $num_bins = 1000;

open IN, "one_round_fwdepc-1-1.results";
$slurp = <IN>;
close IN;

### This removes extra output from instances where 2,4 occurs.
@tmp = split /pt(.|\n)+?\n(A|R)/, $slurp;
print "temp is size " . (scalar @tmp) . " after split.\n";
open OUT, ">one_round_fwdepc-1-1_clean.results";
for my $i (0..@tmp) {
	print OUT $tmp[$i];
}

