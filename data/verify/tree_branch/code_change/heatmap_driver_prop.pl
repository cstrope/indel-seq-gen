#!/usr/bin/env perl
use strict;

my $seqlen = 1000;

for my $i (0 .. $seqlen-1) {
	print "-----------------------$i differences-------------------------\n"; sleep(1);
	system("perl heat_map_prop.pl 9967 $seqlen $i");
}
