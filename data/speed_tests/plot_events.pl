#!/usr/bin/env perl

if (@ARGV != 1) {
	print STDERR "Usage: perl plot_bins.pl <event_file_prefix_to_plot>.\n";
	exit(0);
}

$event_file = shift;

open OUT, ">out.gnu";

print OUT "set term postscript eps color\n";

print OUT "set output \"$event_file.eps\"\n";
print OUT "set style line 2 lt 1 lw 3\n";
print OUT "set style line 1 lt 4 lw 3\n";

print OUT "plot \"$event_file.avg_hits\" usi (\$1-0.015):2:3 w yerrorbars ti \"EPC\" ls 2,";
print OUT "\"$event_file.avg_hits\" usi (\$1+.015):4:5 w yerrorbars ti \"FWD\" ls 1\n";

system("gnuplot out.gnu");
system("open $event_file.eps");
