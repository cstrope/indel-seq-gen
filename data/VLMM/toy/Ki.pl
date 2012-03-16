#!/usr/bin/env perl
use strict;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $file = "run_" . $hour . "." . $min . "." . $sec;
my $num_bins = 1000;
my $low_val = 1000000000;
my $high_val = -10000000;
my $num_runs = 1000;

my (@events,@values);
my $return_val;
my @ret_split;

my @ctr;

for my $i (0 .. $num_runs-1) {
	print "################### $i/$num_runs #####################\n";
	my $forward_val = `./indel-seq-gen -m JC69 -D -e toy -o f < toy.tree`;

	@ret_split = split /\n/, $forward_val;
	my $forward_log_prob = $ret_split[0];

	for my $j (0 .. $num_runs-1) {
		print "       ################### $i.$j/$num_runs #####################\n";
		my $return_val = `./indel-seq-gen -m JC69 -D -E toy.sim.ma < toy.tree`;
	
		@ret_split = split /\n/, $return_val;
		my $epc_log_prob = $ret_split[0];

		if ($epc_log_prob > $forward_log_prob) {
			$ctr[$i] += 1.0/$num_runs;
		}
	}
	print "Forward probability: $forward_val, Number of more probable EPC paths: " . ($ctr[$i]*$num_runs) . "\n";
}

my @bins;
for my $i (0 .. $num_runs-1) {
	$bins[int($ctr[$i]*$num_runs)]++;
}

open OUT, ">$file.k_i.dat";
for my $i (0 .. $num_runs) {
	print OUT "" . ($i/$num_runs) . " ";
	if ($bins[$i] != 0) { print OUT "$bins[$i]\n"; }
	else { print OUT "0\n"; }
}
close OUT;

my ($mu, $stdev) = &mu_rho(\@ctr);
print "average: $mu     standard deviation: $stdev\n";

open OUT, ">$file.k_i.gnu";
print OUT "set term postscript eps color\n";

print OUT "set xrange [0:1]\n";
print OUT "set yrange [0:$num_runs]\n";
print OUT "set style data histogram\n";
print OUT "set style histogram cluster gap 1\n";
print OUT "set style fill solid border -1\n";
print OUT "set boxwidth " . (1/$num_runs) . "\n";
print OUT "set output \"$file.k_i.eps\"\n";
print OUT "plot \"$file.k_i.dat\" w boxes\n";
close OUT;

exit(0);

print "min: $low_val     max: $high_val\n";

my (@bin);
for my $i (0 .. $num_runs) {
	chop($values[$i]);
	$bin[&binNo($values[$i], $low_val, $high_val, $num_bins)]++;
}

open OUT, ">$file.1";
for my $i (0 .. $num_bins) {
	print OUT "" . &hash_mark_val($i, $low_val, $high_val, $num_bins) . " " . int($bin[$i]) . "\n";
}
close OUT;

my $boxwidth = ($high_val - $low_val) / $num_bins;

print "$file.1.eps\n";
open OUT, ">gnuplot1.gnu";
print OUT "set term postscript eps color\n";
print OUT "set auto x\n";
#print OUT "set yrange [0:300000]";
print OUT "set style data histogram\n";
print OUT "set style histogram cluster gap 1\n";
print OUT "set style fill solid border -1\n";
print OUT "set boxwidth $boxwidth\n";
#print OUT "set xtic rotate by -45 scale 0\n";
#set bmargin 10
print OUT "set output \"$file.1.eps\"\n";
print OUT "plot \"$file.1\" w boxes\n";
close OUT;
system("gnuplot gnuplot1.gnu");

open OUT, ">$file.2" or die "fuck.\n";
for my $i (0 .. $num_runs) {
	print OUT "$values[$i] $events[$i]\n";
	print "$values[$i] $events[$i]\n";
}
close OUT;

print "$file.2.eps\n";
open OUT, ">gnuplot2.gnu";
print OUT "set term postscript eps color\n";
print OUT "set auto x\n";
print OUT "set output \"$file.2.eps\"\n";
print OUT "plot \"$file.2\" w points\n";
close OUT;
system("gnuplot gnuplot2.gnu");

sub hash_mark_val {
	my ($value, $low, $high, $num_bins) = @_;
	
	my $val = (($high-$low)*$value)/$num_bins + $low;
#	print "$value corresponds to bin $val\n";

	return $val;
}

sub binNo {
	my ($value, $low, $high, $num_bins) = @_;

	my $bin_no = int((($value - $low) / ($high-$low)) * $num_bins);
#	print "$value goes into bin: $bin_no\n";

	return $bin_no;
}

sub mu_rho
{
	my ($ref) = @_;
	my ($average, $stdev);
	for my $i (0 .. @{ $ref }-1) {
		$average += $ref->[$i];
	}
	$average /= scalar(@{ $ref });#-1;
	for my $i (0 .. @{ $ref }-1) {
		$stdev += ($average - $ref->[$i]) * ($average - $ref->[$i]);
	}
	$stdev /= scalar(@{ $ref });#-1;
	return ($average, sqrt($stdev));
}
