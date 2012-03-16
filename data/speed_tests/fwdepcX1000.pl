#!/usr/bin/env perl
use strict;
use POSIX;
srand(100111);

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $file = "run_" . $hour . "." . $min . "." . $sec;
my $num_runs = 100;
my $num_bins = 1000;
my $sequence_length = 1000;

my $indel_seq_gen = "isg-test";
my $freqs = "-f 0.1,0.2,0.3,0.4";
$freqs = " ";

my $seed;

my $INDEPENDENT_SITES = 0;
my $EMULATE = 1;

my $return_val_fwd;
my $return_val_epc;
my @ret_split;

my $filename = "$file-fwdepc";
if ($INDEPENDENT_SITES) { $filename .= "I"; $file .= "I"; }
$filename = $file . "/" . $filename;

mkdir $file;

my @branch_length = (
#		0.01,
#		0.25,
#		0.5,
#		0.75,
#		1.0,
#		1.25,
#		1.5,
#		1.75,
#		2.0,
#		5.0, 
		10.0,
#		20.0,
#		50.0,
		);

my @dependence_superscript =
		(
#		 0.0,
#		 0.01,
#		 0.1,
#		 0.5,	#sqrt
		 1,
#		 2,		#squared
#		 5,
		);


my @bin_offset =
	(
	0.0,
	0.1,
	0.2,
	0.3,
	0.4,
	0.5,
	0.6,
	0.7,
	0.8,
	0.9,
	);

my $offset = 0.5;

my @array;

my (@bin_first_events_AoA);
@bin_first_events_AoA = ();
for my $x (0 .. $num_bins-1) {
	$bin_first_events_AoA[$x] = [ @array ];
	for my $y (0 .. $num_runs-1) {
		$bin_first_events_AoA[$x][$y] = 0;
	}
}

my (@test_binning, @test_binningsq) ;

my (
	@forward_event_times,
	@forward_event_times_bin_array,
	@epc_event_times, 
	@epc_event_times_bin_array,
	@forward_num_events, 
	@epc_num_events,
	@forward_Qidot,
	@epc_Qidot,
	@forward_Qidot_k_,
	@epc_Qidot_k_,
	@forward_Qidot_bin_array,
	@epc_Qidot_bin_array,
	@forward_Qidot_k_bin_array,
	@epc_Qidot_k_bin_array,
   );

##########
### Use these arrays to replace the above arrays. Use record_data_2D function.
##########
my (@forward_AoA_event_array, @forward_AoA_bin_array);
my (@epc_AoA_event_array, @epc_AoA_bin_array);
my $NUM_STATS = 5;
for my $i (0 .. $NUM_STATS) {
	$forward_AoA_bin_array[$i] = [ @array ];
	$epc_AoA_bin_array[$i] = [ @array ];
}

my (
	$DT,
	$LAMBDA_T,
	$SITE_RATE_AWAY,
	$EVENTS,
	$Qidot,
	$Qidot_k_,
   ) = (0, 1, 2, 3, 4, 5);

##########
### Data structure explanations.
##########
## These arrays capture one value per run, the total number of events that occur.
##  * Coincide with the output "EVENTS:"
@forward_num_events = ();
@epc_num_events = ();
## These arrays capture many events per run, capture each time point that an event occurs,
## to be used for binning the event times.
##  * Coincide with the output "DT:"
@forward_event_times = ();
@epc_event_times = ();
## These arrays are the binning arrays for the previous. Add all DT times to the appropriate bin.
@forward_event_times_bin_array = ();
@epc_event_times_bin_array = ();
## These arrays deal with Qi. and Qi.|k over time.

open OUT, ">$file/avg_hits";
for my $BL (0 .. @branch_length-1) {
	@forward_AoA_event_array = ();
	@epc_AoA_event_array = ();
	for my $i (0 .. $NUM_STATS) {
		$forward_AoA_event_array[$i] = [ @array ];
		$epc_AoA_event_array[$i] = [ @array ];
	}
	@forward_num_events = ();
	@epc_num_events = ();
	@forward_event_times = ();
	@epc_event_times = ();
	@forward_event_times_bin_array = ();
	@epc_event_times_bin_array = ();
	@forward_Qidot_bin_array = ();
	@forward_Qidot_k_bin_array = ();
	@epc_Qidot_bin_array = ();
	@epc_Qidot_k_bin_array = ();
	##########
	### Output the guide tree for this branch length.
	##########
	open OUTt, ">$filename-$branch_length[$BL].tree";
	print OUTt "[$sequence_length](T_1:$branch_length[$BL],T_2:0);\n";
	close OUTt;
	for my $i (0 .. $num_runs-1) {
		print "################### $i/$num_runs #####################\n";
		@forward_AoA_event_array = ();
		@epc_AoA_event_array = ();
		for my $i (0 .. $NUM_STATS) {
			$forward_AoA_event_array[$i] = [ @array ];
			$epc_AoA_event_array[$i] = [ @array ];
		}
		@forward_event_times = ();
		@epc_event_times = ();
		@forward_Qidot = ();
		@forward_Qidot_k_ = ();
		@epc_Qidot = ();
		@epc_Qidot_k_ = ();
		my $emulation_option = "";
		if ($EMULATE) { $emulation_option = "-M $filename-$branch_length[$BL].sim.trace "; }
		$seed = "-z " . int(rand()*1000) . "," . int(rand()*1000) . "," . int(rand()*1000) . "," . int(rand()*1000) . " ";
		if ($INDEPENDENT_SITES) {
			$return_val_fwd = `./$indel_seq_gen -m JC69 $seed -e $filename-$branch_length[$BL] -o f < $filename-$branch_length[$BL].tree`;
			&record_data($return_val_fwd, \@forward_num_events, \@forward_event_times, \@forward_Qidot, \@forward_Qidot_k_);
			&record_data_2D($return_val_fwd, \@forward_AoA_event_array);
			$return_val_epc = `./$indel_seq_gen -m JC69 -seed $emulation_option -E $filename-$branch_length[$BL].sim.ma $emulation_option < $filename-$branch_length[$BL].tree`;
			&record_data($return_val_epc, \@epc_num_events, \@epc_event_times, \@epc_Qidot, \@epc_Qidot_k_);
			&record_data_2D($return_val_epc, \@epc_AoA_event_array);
		} else {
			my $command = "./$indel_seq_gen -m JC69 -D $dependence_superscript[0] $seed -e $filename-$branch_length[$BL] -o f < $filename-$branch_length[$BL].tree";
			$return_val_fwd = `$command`;
			&record_data($return_val_fwd, \@forward_num_events, \@forward_event_times, \@forward_Qidot, \@forward_Qidot_k_);
			&record_data_2D($return_val_fwd, \@forward_AoA_event_array);
			$command = "./$indel_seq_gen -m JC69 -D $filename-$branch_length[$BL].sim.dep $seed $emulation_option -E $filename-$branch_length[$BL].sim.ma < $filename-$branch_length[$BL].tree";
			print "$command\n";
			$return_val_epc = `$command`;
			&record_data($return_val_epc, \@epc_num_events, \@epc_event_times, \@epc_Qidot, \@epc_Qidot_k_);
			&record_data_2D($return_val_epc, \@epc_AoA_event_array);
		}

		##########
		### Increment the bins for 
		##########
		my $bin;
		for my $a (0 .. @forward_event_times-1) {
			$bin = &bin_no($forward_event_times[$a], $branch_length[$BL], 0);
			$forward_event_times_bin_array[$bin]++;
			$forward_Qidot_bin_array[$bin] += $forward_Qidot[$a];
			$forward_Qidot_k_bin_array[$bin] += $forward_Qidot_k_[$a];

			$bin = &bin_no($forward_AoA_event_array[$DT][$a], $branch_length[$BL]);
			$forward_AoA_bin_array[$DT][$bin]++;
			$forward_AoA_bin_array[$Qidot][$bin] += $forward_AoA_event_array[$Qidot][$a];
			$forward_AoA_bin_array[$Qidot_k_][$bin] += $forward_AoA_event_array[$Qidot_k_][$a];
			$forward_AoA_bin_array[$SITE_RATE_AWAY][$bin] += $forward_AoA_event_array[$SITE_RATE_AWAY][$a];

#			print "$forward_event_times_bin_array[$bin] $forward_AoA_bin_array[$DT][$bin]\n";
#			print "$forward_Qidot_bin_array[$bin] $forward_AoA_bin_array[$Qidot][$bin]\n";
#			print "$forward_Qidot_k_bin_array[$bin] $forward_AoA_bin_array[$Qidot_k_][$bin]\n";
		}
		for my $a (0 .. @epc_event_times-1) {
			$bin = &bin_no($epc_event_times[$a], $branch_length[$BL], 0);
			$epc_event_times_bin_array[$bin]++;
			$epc_Qidot_bin_array[$bin] += $epc_Qidot[$a];
			$epc_Qidot_k_bin_array[$bin] += $epc_Qidot_k_[$a];

			$bin = &bin_no($epc_AoA_event_array[$DT][$a], $branch_length[$BL]);
			$epc_AoA_bin_array[$DT][$bin]++;
			$epc_AoA_bin_array[$Qidot][$bin] += $epc_AoA_event_array[$Qidot][$a];
			$epc_AoA_bin_array[$Qidot_k_][$bin] += $epc_AoA_event_array[$Qidot_k_][$a];
			$epc_AoA_bin_array[$SITE_RATE_AWAY][$bin] += $epc_AoA_event_array[$SITE_RATE_AWAY][$a];

#			print "$epc_event_times_bin_array[$bin] $epc_AoA_bin_array[$DT][$bin]\n";
#			print "$epc_Qidot_bin_array[$bin] $epc_AoA_bin_array[$Qidot][$bin]\n";
#			print "$epc_Qidot_k_bin_array[$bin] $epc_AoA_bin_array[$Qidot_k_][$bin]\n";
		}

		##########
		### First, fill in the bins with ONLY their first hit, based on the offset, of course.
		### * This may skip bins. That will be dealt with later.
		##########
		for my $a (0 .. @epc_event_times-1) {
			$bin = &bin_no($epc_event_times[$a], $branch_length[$BL], ($branch_length[$BL]*$offset) / $num_bins);
			if ($bin_first_events_AoA[$bin][$i] == 0) {
				$bin_first_events_AoA[$bin][$i] = $epc_Qidot_k_[$a];
			}
		}

		##########
		### All the bins that are first hit are now filled. However, some events might skip bins.
		### Those bins that have been skipped will still hold the value of 0.0. Fill those bins in
		### with the appropriate value.
		##########
		for my $a (1 .. $num_bins-1) {
			if ($bin_first_events_AoA[$a][$i] == 0) {
				#### First bin is empty. Need to fill it with the value of the immediately next bin.
				if ($bin_first_events_AoA[$a-1][$i] == 0) {
					my $x = $a;
					while ($bin_first_events_AoA[$x][$i] == 0) { 
						$x++; 
						if ($x > $num_bins) { print STDERR "Something wonky happening...\n"; exit(0);  }
					}
					$bin_first_events_AoA[$a-1][$i] = $bin_first_events_AoA[$x][$i];
					$bin_first_events_AoA[$a][$i] = $bin_first_events_AoA[$x][$i];
				} else {
					$bin_first_events_AoA[$a][$i] = $bin_first_events_AoA[$a-1][$i]; 
				}
			}
		}

		#for my $a (0 .. $num_bins-1) {
		#	print "bin $a: $bin_first_events_AoA[$a][$i]\n";
		#}

		#&bin_first_events($i, $branch_length[$BL], $num_bins, \@epc_AoA_event_array[$DT], \@epc_AoA_event_array[$Qidot_k_] )
	}

	my ($mu, $stdev);

	open OUT, ">$file/$branch_length[$BL]-alt.Qcmp";
	for my $a (0 .. $num_bins-1) {
		#for my $b (0 .. $num_runs-1) { print "$bin_first_events_AoA[$a][$b] "; }
		($mu, $stdev) = &mu_rho(\@{ $bin_first_events_AoA[$a] });
		my $bin = ($a + $offset) / ($num_bins*$branch_length[$BL]);
		print OUT "$bin $mu $stdev\n";
	}

	#print "EPC:\n";
	($mu, $stdev) = &mu_rho(\@epc_num_events);
	#print "average subst: $mu     standard deviation: $stdev\n";
	print OUT "$branch_length[$BL] ";
	print OUT "$mu $stdev ";
	#print "FWD:\n";
	($mu, $stdev) = &mu_rho(\@forward_num_events);
	#print "average subst: $mu     standard deviation: $stdev\n";
	print OUT "$mu $stdev\n";

	open OUT2, ">$file/$branch_length[$BL].bins";
	for my $a (0 .. $num_bins-1) {
		print OUT2 "" . (($a/$num_bins)*$branch_length[$BL]). " " . (int($epc_event_times_bin_array[$a])/$num_bins);
		print OUT2 " " . (int($forward_event_times_bin_array[$a])/$num_bins) . "\n";
	}
	close OUT2;
	print "Binned times in file \"$file/$branch_length[$BL].bins\".\n";

	open OUT2, ">$file/$branch_length[$BL].Qcmp";
	for my $a (0 .. $num_bins-1) {
		print OUT2 (($a/$num_bins)*$branch_length[$BL]);
		if ( $epc_event_times_bin_array[$a] != 0 ) {
			print OUT2 " " . (int($epc_Qidot_bin_array[$a])/$epc_event_times_bin_array[$a]);
			print OUT2 " " . (int($epc_Qidot_k_bin_array[$a])/$epc_event_times_bin_array[$a]);
		} else { print OUT2 " " . 0 . " " . 0; }
		if ( $forward_event_times_bin_array[$a] != 0 ) {
			print OUT2 " " . (int($forward_Qidot_bin_array[$a])/$forward_event_times_bin_array[$a]);
			print OUT2 " " . (int($forward_Qidot_k_bin_array[$a])/$forward_event_times_bin_array[$a]);
		} else { print OUT2 " " . 0 . " " . 0; }
			print OUT2 "\n";
	}
	close OUT2;

	open OUT2, ">$file/$branch_length[$BL]-2.Qcmp";
	for my $a (0 .. $num_bins-1) {
		print OUT2 (($a/$num_bins)*$branch_length[$BL]);
		if ( $epc_AoA_bin_array[$DT][$a] != 0 ) {
			print OUT2 " " . (int($epc_AoA_bin_array[$Qidot][$a])/$epc_AoA_bin_array[$DT][$a]);
			print OUT2 " " . (int($epc_AoA_bin_array[$Qidot_k_][$a])/$epc_AoA_bin_array[$DT][$a]);
		} else { print OUT2 " " . 0 . " " . 0; }
		if ( $forward_AoA_bin_array[$DT][$a] != 0 ) {
			print OUT2 " " . (int($forward_AoA_bin_array[$Qidot][$a])/$forward_AoA_bin_array[$DT][$a]);
			print OUT2 " " . (int($forward_AoA_bin_array[$Qidot_k_][$a])/$forward_AoA_bin_array[$DT][$a]);
		} else { print OUT2 " " . 0 . " " . 0; }
		print OUT2 "\n";
	}
	close OUT2;

#	for my $a (0 .. $num_bins-1) {
#		print " " . (int($epc_Qidot_bin_array[$a])/$epc_event_times_bin_array[$a]) . " " . (int($epc_AoA_bin_array[$Qidot][$a])/$epc_AoA_bin_array[$DT][$a]) . "\n";
#		print " " . (int($epc_Qidot_k_bin_array[$a])/$epc_event_times_bin_array[$a]) . " " . (int($epc_AoA_bin_array[$Qidot_k_][$a])/$epc_AoA_bin_array[$DT][$a]) . "\n";
#		print " " . (int($forward_Qidot_bin_array[$a])/$forward_event_times_bin_array[$a]) . " " . (int($forward_AoA_bin_array[$Qidot][$a])/$forward_AoA_bin_array[$DT][$a]) . "\n";
#		print " " . (int($forward_Qidot_k_bin_array[$a])/$forward_event_times_bin_array[$a]) . " " . (int($forward_AoA_bin_array[$Qidot_k_][$a])/$forward_AoA_bin_array[$DT][$a]) . "\n";
#		if ((int($epc_Qidot_bin_array[$a])/$epc_event_times_bin_array[$a]) != (int($epc_AoA_bin_array[$Qidot][$a])/$epc_AoA_bin_array[$DT][$a]) ) {
#			print "#1 not equal.\n";
#		}
#		if ((int($epc_Qidot_k_bin_array[$a])/$epc_event_times_bin_array[$a]) != (int($epc_AoA_bin_array[$Qidot_k_][$a])/$epc_AoA_bin_array[$DT][$a]) ) {
#			print "#2 not equal.\n";
#		}
#		if ((int($forward_Qidot_bin_array[$a])/$forward_event_times_bin_array[$a]) != (int($forward_AoA_bin_array[$Qidot][$a])/$forward_AoA_bin_array[$DT][$a]) ) {
#			print "#3 not equal.\n";
#		}
#		if ((int($forward_Qidot_k_bin_array[$a])/$forward_event_times_bin_array[$a]) != (int($forward_AoA_bin_array[$Qidot_k_][$a])/$forward_AoA_bin_array[$DT][$a]) ) {
#			print "#4 not equal.\n";
#		}
#	}


}
close OUT;

#sub bin_first_events
#{
#	my ($replicate, $BL, $num_bins, $event_times_ref, $Qij_k_ref ) = @_;
#	my ($bin_size);
#
	##########
	### Two arrays are global in scope: $bin_first_events_AoA[val][rep] and $bin_offset[], 
	### that specifies where in the bin we are collecting the first event.
	##########

#	$bin_size = $BL / $num_bins;

#	for my $offset (0 .. @bin_offset-1) {
#		my $in_bin_cutoff = $bin_size * $bin_offset[$b];
#		my $cutoff = $in_bin_cutoff;

		### Cycle through all of the events that occurred... this index works for all ptrs to arrays. ###
#		for my $a (0 .. @{ @event_times_ref }-1) {
#			if ( $cutoff > @{ $event_times_ref }->[$a] ) {
#				while ( $cutoff > @{ $event_times_ref }->[$a] && $a < $num_bins ) {
#					$bin_first_events_AoA[$offset][$a][$replicate] = @{ $Qij_k_ref }->[$a];
#					$a++;
#				}
#				$cutoff+=$bin_size;
#			}
#		}

#		exit(0);
#	}
#}

sub bin_no
{
	my ($val, $branch_length, $offset) = @_;

#	print "val: $val   offset: $offset   BL: $branch_length   num_bins: $num_bins\n";
#	print " (($val-$offset)*$num_bins)/$branch_length\n";

	### if $val is 0.9005, offset is 0.0004, BL = 1, then this will return 900. If offset=0.0006, then 899, which is appropriate.
	my $return_val = (($val-$offset)* $num_bins) / $branch_length;

#	print " " . int($return_val) . "\n";

	return int($return_val);
}

sub record_data_2D
{
	my ($data, $AoA_array_ref) = @_;
	my (@ret_split);

	@ret_split = split/\n\n+/, $data;

	for my $i (0 .. @ret_split-1) {
		my @line_split = split /\n/, $ret_split[$i];
		
		for my $j (0 .. @line_split-1) {
			my @colon_split = split(":", $line_split[$j]);
			if ($colon_split[0] eq "DT") {
				push(@{ $AoA_array_ref->[$DT] }, $colon_split[1]);
			} elsif ($colon_split[0] eq "LAMBDA_T") {
				push(@{ $AoA_array_ref->[$LAMBDA_T] }, $colon_split[1]);
			} elsif ($colon_split[0] eq "SITE_RATE_AWAY") {
				push(@{ $AoA_array_ref->[$SITE_RATE_AWAY] }, $colon_split[1]);
			} elsif ($colon_split[0] eq "EVENTS") {
				push(@{ $AoA_array_ref->[$EVENTS] }, $colon_split[1]);
			} elsif ($colon_split[0] eq "Qi.") {
				push(@{ $AoA_array_ref->[$Qidot] }, $colon_split[1]);
			} elsif ($colon_split[0] eq "JC_Rij") {
			} elsif ($colon_split[0] eq "Qi.|k") {
				push(@{ $AoA_array_ref->[$Qidot_k_] }, $colon_split[1]);
			} else {
				print $data;
				print STDERR "I do not know what the event for \"$colon_split[0]\" is.\n";
				exit(0);
			}
		}
	}
}

sub record_data
{
	my ($data, $num_event_array_ref, $event_time_array_ref, $Qidot_array_ref, $Qidot_k_array_ref) = @_;
	my (@ret_split);

	@ret_split = split/\n\n+/, $data;

	for my $i (0 .. @ret_split-1) {
		my @line_split = split /\n/, $ret_split[$i];
		
		for my $j (0 .. @line_split-1) {
			my @colon_split = split(":", $line_split[$j]);
			if ($colon_split[0] eq "DT") {
				push(@{ $event_time_array_ref }, $colon_split[1]);
			} elsif ($colon_split[0] eq "LAMBDA_T") {
			} elsif ($colon_split[0] eq "SITE_RATE_AWAY") {
				push(@{ $num_event_array_ref }, $colon_split[1]);
			} elsif ($colon_split[0] eq "EVENTS") {
				push(@{ $num_event_array_ref }, $colon_split[1]);
			} elsif ($colon_split[0] eq "Qi.") {
				push(@{ $Qidot_array_ref }, $colon_split[1]);
			} elsif ($colon_split[0] eq "JC_Rij") {
			} elsif ($colon_split[0] eq "Qi.|k") {
				push(@{ $Qidot_k_array_ref }, $colon_split[1]);
			} else {
				print $data;
				print STDERR "I do not know what the event for \"$colon_split[0]\" is.\n";
				exit(0);
			}
		}

	}
}

sub mu_rho
{
	my ($ref) = @_;
	my ($average, $stdev);
	if (scalar (@{ $ref }) != 0) { 
		for my $i (0 .. @{ $ref }-1) {
			$average += $ref->[$i];
		}
		$average /= scalar(@{ $ref });#-1;
		for my $i (0 .. @{ $ref }-1) {
			$stdev += ($average - $ref->[$i]) * ($average - $ref->[$i]);
		}
		$stdev /= scalar(@{ $ref });#-1;
		return ($average, sqrt($stdev));
	} else { return (0, 0); }
}