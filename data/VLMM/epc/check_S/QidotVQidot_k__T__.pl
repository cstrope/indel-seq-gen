#!/usr/bin/env perl
use strict;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $file = "run_" . $hour . "." . $min . "." . $sec;
my $num_runs = 100;
my $num_bins = $num_runs;

my $INDEPENDENT_SITES = 0;

my $return_val_fwd;
my $return_val_epc;
my @ret_split;

my $j = 0;
my @array;

my @S = (
#		1.0,
#		0.1,
#		0.3,
		0.42,
#		0.7,
#		2.0,
#		100.0,
		);

for my $s (0 .. @S-1) {
	open OUT, ">$file-$S[$s].dat";
	open OUT2, ">$file-$S[$s].num_diff";
	for my $i (0 .. $num_runs) {
		print "################### $i/$num_runs, S = $S[$s] #####################\n";
		if ($INDEPENDENT_SITES) {
			$return_val_fwd = `./indel-seq-gen -m JC69 -e toy -o f < toy.tree`;
			$return_val_epc = `./indel-seq-gen -m JC69 -E toy.sim.ma -S $S[$s] < toy.tree`;
			print OUT "$return_val_epc";
		} else {
			$return_val_fwd = `./indel-seq-gen -m JC69 -D -e toy -o f < toy.tree`;
			$return_val_epc = `./indel-seq-gen -m JC69 -D -E toy.sim.ma -S $S[$s] < toy.tree`;
			@ret_split = split(">>", $return_val_epc);
			print OUT "$ret_split[0]";
			for my $j (1 .. @ret_split) {
				my @ret_split2 = split /\n/, $ret_split[$j];
				print OUT2 "$ret_split2[0]\n";
				print OUT "$ret_split2[1]\n";
			}
		}
	}
	close OUT;

	open OUT, ">$file-$S[$s].0.gnu";
	print OUT "set term postscript eps color\n";
	print OUT "set log y\n";
	print OUT "set xrange [0:1.5]\n";
	print OUT "set output \"$file-$S[$s].0.eps\"\n";
	print OUT "plot \"$file-$S[$s].dat\" usi 1:2 w points ti \"$S[$s] Qi.\",";
	print OUT " \"$file-$S[$s].dat\" usi 1:3 w points ti \"$S[$s] Qi.|k(t)\"\n";
	close OUT;
	
	system("gnuplot $file-$S[$s].0.gnu");
	system("open $file-$S[$s].0.eps");

	open OUT, ">$file-$S[$s].1.5.gnu";
	print OUT "set term postscript eps color\n";
	print OUT "set log y\n";
	print OUT "set xrange [1.5:2.75]\n";
	print OUT "set output \"$file-$S[$s].1.5.eps\"\n";
	print OUT "plot \"$file-$S[$s].dat\" usi 1:2 w points ti \"$S[$s] Qi.\",";
	print OUT " \"$file-$S[$s].dat\" usi 1:3 w points ti \"$S[$s] Qi.|k(t)\"\n";
	close OUT;
	
	system("gnuplot $file-$S[$s].1.5.gnu");
	system("open $file-$S[$s].1.5.eps");

	open OUT, ">$file-$S[$s].2.75.gnu";
	print OUT "set term postscript eps color\n";
	print OUT "set log y\n";
	print OUT "set xrange [2.75:3]\n";
	print OUT "set output \"$file-$S[$s].2.75.eps\"\n";
	print OUT "plot \"$file-$S[$s].dat\" usi 1:2 w points ti \"$S[$s] Qi.\",";
	print OUT " \"$file-$S[$s].dat\" usi 1:3 w points ti \"$S[$s] Qi.|k(t)\"\n";
	close OUT;
	
	system("gnuplot $file-$S[$s].2.75.gnu");
	system("open $file-$S[$s].2.75.eps");

	open OUT, ">$file-$S[$s].gnu2";
	print OUT "set term postscript eps color\n";
	print OUT "set log y\n";
	print OUT "set xrange [0:3]\n";
	print OUT "set ylabel \"Rate away\"\n";
	print OUT "set y2label \"Number of differences\"\n";
	print OUT "set y2range [0:1500]\n";
	print OUT "set y2tic 0,100\n";
	print OUT "set xlabel \"time on branch\"\n";

	print OUT "set output \"$file-$S[$s].num_diff.eps\"\n";
	print OUT "plot \"$file-$S[$s].num_diff\" axes x1y2 noti, ";
#    print OUT "\"$file-$S[$s].dat\" usi 1:3 axes x1y1 noti,";
    print OUT "\"$file-$S[$s].dat\" usi 1:4 axes x1y1 noti";
#    print OUT "\"$file-$S[$s].dat\" usi 1:2 axes x1y1 noti";
	close OUT;

	system("gnuplot $file-$S[$s].gnu2");
	system("open $file-$S[$s].num_diff.eps");
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