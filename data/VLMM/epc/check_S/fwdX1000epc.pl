#!/usr/bin/env perl
use strict;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $file = "run_" . $hour . "." . $min . "." . $sec;
my $num_runs = 1000;
my $num_bins = $num_runs;

my $INDEPENDENT_SITES = 1;

my $return_val;
my @ret_split;

my $j = 0;
my @array;
my ($forward_value, $forward_event, @epc_values, @epc_events, @event_time_bins);

if ($INDEPENDENT_SITES) {
	$return_val = `./indel-seq-gen -m JC69 -e toy -o f < toy.tree`;
} else {
	$return_val = `./indel-seq-gen -m JC69 -D -e toy -o f < toy.tree`;
}
@ret_split = split /\n/, $return_val;
$forward_value = $ret_split[0];
$forward_event = $ret_split[1];
open OUT, ">$file.dat";
open OUT3, ">$file.events_table";
print OUT3 "$forward_event events in forward simulation.\n";
for (my $S =0.1; $S <= 2; $S += 0.1, $j++) {
	@event_time_bins = ();
	for my $i (0 .. $num_runs) {
		print "################### $i/$num_runs, S = $S #####################\n";
		if ($INDEPENDENT_SITES) {
			$return_val = `./indel-seq-gen -m JC69 -E toy.sim.ma -S $S < toy.tree`;
		} else {
			$return_val = `./indel-seq-gen -m JC69 -D -E toy.sim.ma -S $S < toy.tree`;
		}
		my @tmp_split = split /\n\n/, $return_val;
		@ret_split = split /\n/, $tmp_split[1];
		$epc_values[$i] = $ret_split[0];
		$epc_events[$i] = $ret_split[1];

		my @event_split = split /\n/, $tmp_split[0];
		for my $a (0 .. @event_split-1) {
			my @space_split = split /\s+/, $event_split[$a];
			my $bin_no = int($space_split[0]*$num_bins);
			$event_time_bins[$bin_no]++;
		}

	}

	my ($mu, $stdev) = &mu_rho(\@epc_values);
	print "average logPr: $mu     standard deviation: $stdev\n";
	print OUT "$S $mu $stdev ";
	($mu, $stdev) = &mu_rho(\@epc_events);
	print "average subst: $mu     standard deviation: $stdev\n";
	print OUT "$mu $stdev\n";

	print OUT3 "$S $mu $stdev\n";

	open OUT2, ">$file.$S.dat";
	for my $a (0 .. $num_bins) {
		if ($event_time_bins[$a] != 0) {
			print OUT2 "" . ($a/1000) . " " . ($event_time_bins[$a]/$num_runs) . "\n";
		} else {
			print OUT2 "" . ($a/1000) . " 0\n";
		}
	}
	close OUT2;
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