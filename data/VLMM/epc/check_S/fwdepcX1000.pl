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
my (@forward_value, @forward_event, @epc_values, @epc_events, @event_time_bins);

my @S = (
#		1.0,
#		0.1,
#		0.3,
		0.42,
#		0.7,
#		2.0,
#		100.0,
		);

open OUT, ">$file.dat";
open OUT3, ">$file.events_table";
for my $s (0 .. @S-1) {
	@event_time_bins = ();
	@forward_value = ();
	@forward_event = ();
	@epc_values = ();
	@epc_events = ();
	for my $i (0 .. $num_runs) {
		print "################### $i/$num_runs, S = $S[$s] #####################\n";
		if ($INDEPENDENT_SITES) {
			$return_val_fwd = `./indel-seq-gen -m JC69 -e toy -o f < toy.tree`;
			@ret_split = split /\n/, $return_val_fwd;
			#$forward_value[$i] = $ret_split[0];
			$forward_event[$i] = $ret_split[0];
			$return_val_epc = `./indel-seq-gen -m JC69 -E toy.sim.ma -S $S[$s] < toy.tree`;
		} else {
			$return_val_fwd = `./indel-seq-gen -m JC69 -D -e toy -o f < toy.tree`;
			print "return_val_fwd: $return_val_fwd\n";
			@ret_split = split /\n/, $return_val_fwd;
			#$forward_value[$i] = $ret_split[0];
			$forward_event[$i] = $ret_split[0];
			$return_val_epc = `./indel-seq-gen -m JC69 -D -E toy.sim.ma -S $S[$s] < toy.tree`;
			print "return_val_epc: $return_val_epc\n";
		}
		my @tmp_split = split /\n/, $return_val_epc;
		@ret_split = split /\n/, $tmp_split[1];
		#$epc_values[$i] = $ret_split[0];
		$epc_events[$i] = $ret_split[0];

#		my @event_split = split /\n/, $tmp_split[0];
#		for my $a (0 .. @event_split-1) {
#			my @space_split = split /\s+/, $event_split[$a];
#			my $bin_no = int($space_split[0]*$num_bins);
#			$event_time_bins[$bin_no]++;
#		}

	}

	print "EPC:\n";
	my ($mu, $stdev);
#	($mu, $stdev) = &mu_rho(\@epc_values);
#	print "average logPr: $mu     standard deviation: $stdev\n";
#	print OUT "$S[$s] $mu $stdev ";
	($mu, $stdev) = &mu_rho(\@epc_events);
	print "average subst: $mu     standard deviation: $stdev\n";
	print OUT "$mu $stdev\n";
	print OUT3 "$S[$s] $mu $stdev    ";
	print "FWD:\n";
#	($mu, $stdev) = &mu_rho(\@forward_value);
#	print "average logPr: $mu     standard deviation: $stdev\n";
	($mu, $stdev) = &mu_rho(\@forward_event);
	print "average subst: $mu     standard deviation: $stdev\n";
	print OUT3 " $mu $stdev\n";

#	open OUT2, ">$file.$S[$s].dat";
#	for my $a (0 .. $num_bins) {
#		if ($event_time_bins[$a] != 0) {
#			print OUT2 "" . ($a/$num_bins) . " " . ($event_time_bins[$a]/$num_bins) . "\n";
#		} else {
#			print OUT2 "" . ($a/$num_bins) . " 0\n";
#		}
#	}
#	close OUT2;
}
close OUT;
close OUT3;

print "output is in: $file.dat\n";

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