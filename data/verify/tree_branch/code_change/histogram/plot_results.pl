#!/usr/bin/env perl
$/=undef;
use strict;

my $infile = shift;

my $burn_in = 50000;
my $num_bins = 1000;

open IN, "$infile";
my $slurp = <IN>;
close IN;

my (@tmp, @tmp2, @tmp3);
### This removes extra output from instances where 2,4 occurs.
#@tmp = split /pt(.|\n)+?\n(A|R)/, $slurp;
#print "temp is size " . (scalar @tmp) . " after split.\n";
#open OUT, ">2-1_clean.results";
#for my $i (0..@tmp) {
#	print OUT $tmp[$i];
#}

### This splits the accept and reject from the number of differences between sequences.
#open OUT, ">2-1.results";
#@tmp = split /\n/, $slurp;
#for my $i (0..@tmp) {
#	@tmp2 = split /\s+/, $tmp[$i];
#	for my $j (0..5) { print OUT "$tmp2[$j] "; }
#	@tmp3 = split //, $tmp2[6];
#	$tmp3_end = scalar @tmp3;
#	for my $j (0..($tmp3_end-2)) { print OUT "$tmp3[$j]"; }
#	print OUT " $tmp3[$tmp3_end-1]";
#	print OUT " $tmp2[7]";
#	print OUT " $tmp2[8]";
#	print OUT "\n";
#}
#close OUT;

##########
### Binning routine for accepted paths versus number of paths simulated, by length
##########
my ($diff, $subst, $fwd, @change_in_probability, @accepted_path_length);
my (@bin_array, @accepted_array, @tmpprev);
@change_in_probability = ();
@tmp = split /\n/, $slurp;
for my $i ($burn_in..@tmp-1) {
	@tmp2 = split /\s+/, $tmp[$i];
	@tmpprev = split /\s+/, $tmp[$i-1];
	$diff = (int(($tmp2[2] - $tmp2[1])*$num_bins))/$num_bins;
	$subst = $tmp2[8];
	$fwd = $tmp2[9];
	$bin_array[$diff*$num_bins]++;
	if (&accept_path($tmp2[7])) { 
		$accepted_array[$diff*$num_bins]++; 
		push (@change_in_probability, $fwd - $tmpprev[9]);
		push (@accepted_path_length, $diff);
	}
}

##########
### Create arrays holding the bin array reference number and acceptance data.
##########
my (@x_data, @y_data);
for my $i (0..$num_bins-1) {
	if ($bin_array[$i]) {
		$x_data[$i] = $accepted_array[$i]/$bin_array[$i];
	} else {
		$x_data[$i] = 0;
	}
	$y_data[$i] = $i/$num_bins;
}


&gnuplot(
		 $infile . "_accepted_histo",
		 "Subpath length (t_E - t_B)",
		 "Percent M-H acceptance",
		 "Percent accepted path lengths, $infile",
		 \@x_data,
		 \@y_data
	    );

@x_data = ();
@y_data = ();
for my $i (0..@change_in_probability-2) {
	$x_data[$i] = $change_in_probability[$i];
	$y_data[$i] = $accepted_path_length[$i];
}

&gnuplot(
		 $infile . "_pvl_plot",
		 "Subpath length (t_z - t_{z-1})",
		 "Forward probability change",
		 "Change in Probability versus path length, $infile",
		 \@x_data,
		 \@y_data
	    );

sub gnuplot {
	my ($outfile, $xlabel, $ylabel, $title, $x_data_ref, $y_data_ref) = @_;

	open OUT, ">$outfile.dat";
	for my $i (0..@{ $y_data_ref }-1) {
		print OUT "" . ( $y_data_ref->[$i]) . " " .  ( $x_data_ref->[$i] ) . "\n";
	}

	open OUT, ">$outfile.gnu";
	print OUT "set term postscript eps color\n";
	print OUT "set title \"$title\"\n";
	print OUT "set xlabel \"$xlabel\"\n";
	print OUT "set ylabel \"$ylabel\"\n";
	print OUT "set output \"$outfile.eps\"\n";
	print OUT "plot \"$outfile.dat\" w points ti \"$infile\"\n";
	close OUT;
	
	system("gnuplot $outfile.gnu");
}

sub accept_path {
	my ($char) = @_;
	
	if ($char eq "R") { return 0; }
	else { return 1; }
}