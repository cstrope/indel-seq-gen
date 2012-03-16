#!/usr/bin/env perl

if (@ARGV != 1) {
	print STDERR "Usage: perl plot_bins.pl <bin_file_prefix_to_plot>.\n";
	exit(0);
}

$bin_file = shift;

open OUT, ">out.gnu";

print OUT "set term postscript eps color\n";
print OUT "set style line 1 lt 1 lw 3\n";
print OUT "set style line 2 lt 3 lw 3\n";
print OUT "set style line 3 lt 4 lw 3\n";

print OUT "set output \"$bin_file.eps\"\n";

print OUT "plot \"$bin_file.Qcmp\" usi 1:2 w points ls 1 ti \"EPC Qi.\",";
print OUT "\"$bin_file.Qcmp\" usi 1:3 w points ls 2 ti \"EPC Qi.|k(t)\",";
print OUT "\"$bin_file.Qcmp\" usi 1:4 w points ls 3 ti \"FWD Qi.\"\n";

system("gnuplot out.gnu");
system("open $bin_file.eps");

