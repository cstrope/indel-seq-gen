#!/usr/bin/env perl

if (@ARGV != 1) {
	print STDERR "Usage: perl plot_bins.pl <bin_file_prefix_to_plot>.\n";
	exit(0);
}

$bin_file = shift;

open OUT, ">out.gnu";

print OUT "set term postscript eps color\n";

print OUT "set output \"$bin_file.eps\"\n";

print OUT "plot \"$bin_file.bins\" usi 1:2 w points ti \"EPC\",";
print OUT "\"$bin_file.bins\" usi 1:3 w points ti \"FWD\"\n";

system("gnuplot out.gnu");
system("open $bin_file.eps");
