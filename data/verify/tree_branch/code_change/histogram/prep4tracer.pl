#!/usr/bin/env perl
$/=undef;
use strict;

if (@ARGV != 2) {
	print STDERR "Usage: perl prep4tracer.pl <mcmc.outfile> <sample_each>.\n";
	exit(0);
}

my $infile = shift;
my $sample_each = shift;

my @tmp = split /\./, $infile;

my $outfile_prefix = $tmp[0];
for my $i (1..@tmp-2) {
	$outfile_prefix .= "." . $tmp[$i];
}

$outfile_prefix .= "-$sample_each" . "c";

open IN, "$infile" or die "WTF is $infile????\n";
my $slurp = <IN>;
close IN;

my @dat = split /\n/, $slurp;

open OUT, ">$outfile_prefix.p";

print OUT "[Data for EPC run $infile]\n";
print OUT "Gen\tLnL\tNS\n";

for (my $i = 0; $i < @dat-1; $i += $sample_each) {
	my @line_split = split /\s+/, $dat[$i];
	print OUT "$line_split[0]\t$line_split[9]\t$line_split[8]\n";
}

close OUT;

print "Printed to file $outfile_prefix.p.\n";