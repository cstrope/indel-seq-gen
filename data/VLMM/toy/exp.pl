#!/usr/bin/env perl

my $total = 1;
for my $i (0 .. @ARGV-1) {
	$total *= $ARGV[$i];
}

print "$total " . (log($total)) . "\n";
