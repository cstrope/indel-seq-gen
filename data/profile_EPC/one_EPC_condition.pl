#!/usr/bin/env perl
use strict;
use POSIX;
srand(100111);

if (@ARGV != 6) {
	print STDERR "Usage: perl one_EPC_condition.pl <BL> <dep_superscript> <sequence_length> <#reps> <#bins> <emulation?>\n";
	exit(0);
}

my ($BL, $dep_superscript, $sequence_length, $num_runs, $num_bins, $EMULATE) = @ARGV;

print "::$BL $dep_superscript $sequence_length $num_runs $num_bins $EMULATE ::\n"; 

##########
### Create a unique filename.
##########
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $file = "run_" . $hour . "." . $min . "." . $sec . "." . $BL . "-" . $dep_superscript;
if ($EMULATE) { $file .= "emu"; }

my $indel_seq_gen = "isg-test";
my $freqs = " ";
my $seed;
my $INDEPENDENT_SITES = 0;
if ($dep_superscript < 0) { $INDEPENDENT_SITES = 1; }

my $return_val_fwd;
my $return_val_epc;
my @ret_split;
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


##########
### Make work directory.
##########
my $filename = "$file-fwdepc";
if ($INDEPENDENT_SITES) { $filename .= "I"; $file .= "I"; }
$filename = $file . "/" . $filename;
mkdir $file;

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
open OUTt, ">$filename-$BL.tree";
print OUTt "[$sequence_length](T_1:$BL,T_2:0);\n";
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
	if ($EMULATE) { $emulation_option = "-M $filename-$BL.sim.trace "; }
	$seed = "-z " . int(rand()*1000) . "," . int(rand()*1000) . "," . int(rand()*1000) . "," . int(rand()*1000) . " ";
	if ($INDEPENDENT_SITES) {
		$return_val_fwd = `./$indel_seq_gen -m JC69 $seed -e $filename-$BL -o f < $filename-$BL.tree`;
		&record_data($return_val_fwd, \@forward_num_events, \@forward_event_times, \@forward_Qidot, \@forward_Qidot_k_);
		&record_data_2D($return_val_fwd, \@forward_AoA_event_array);
		$return_val_epc = `./$indel_seq_gen -m JC69 -seed $emulation_option -E $filename-$BL.sim.ma $emulation_option < $filename-$BL.tree`;
		&record_data($return_val_epc, \@epc_num_events, \@epc_event_times, \@epc_Qidot, \@epc_Qidot_k_);
		&record_data_2D($return_val_epc, \@epc_AoA_event_array);
	} else {
		my $command = "./$indel_seq_gen -m JC69 -D $dep_superscript $seed -e $filename-$BL -o f < $filename-$BL.tree";
		$return_val_fwd = `$command`;
		&record_data($return_val_fwd, \@forward_num_events, \@forward_event_times, \@forward_Qidot, \@forward_Qidot_k_);
		&record_data_2D($return_val_fwd, \@forward_AoA_event_array);
		$command = "./$indel_seq_gen -m JC69 -D $filename-$BL.sim.dep $seed $emulation_option -E $filename-$BL.sim.ma < $filename-$BL.tree";
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
		$bin = &bin_no($forward_event_times[$a], $BL, 0);
		$forward_event_times_bin_array[$bin]++;
		$forward_Qidot_bin_array[$bin] += $forward_Qidot[$a];
		$forward_Qidot_k_bin_array[$bin] += $forward_Qidot_k_[$a];
		$bin = &bin_no($forward_AoA_event_array[$DT][$a], $BL);
		$forward_AoA_bin_array[$DT][$bin]++;
		$forward_AoA_bin_array[$Qidot][$bin] += $forward_AoA_event_array[$Qidot][$a];
		$forward_AoA_bin_array[$Qidot_k_][$bin] += $forward_AoA_event_array[$Qidot_k_][$a];
		$forward_AoA_bin_array[$SITE_RATE_AWAY][$bin] += $forward_AoA_event_array[$SITE_RATE_AWAY][$a];
	}
	for my $a (0 .. @epc_event_times-1) {
		$bin = &bin_no($epc_event_times[$a], $BL, 0);
		$epc_event_times_bin_array[$bin]++;
		$epc_Qidot_bin_array[$bin] += $epc_Qidot[$a];
		$epc_Qidot_k_bin_array[$bin] += $epc_Qidot_k_[$a];

		$bin = &bin_no($epc_AoA_event_array[$DT][$a], $BL);
		$epc_AoA_bin_array[$DT][$bin]++;
		$epc_AoA_bin_array[$Qidot][$bin] += $epc_AoA_event_array[$Qidot][$a];
		$epc_AoA_bin_array[$Qidot_k_][$bin] += $epc_AoA_event_array[$Qidot_k_][$a];
		$epc_AoA_bin_array[$SITE_RATE_AWAY][$bin] += $epc_AoA_event_array[$SITE_RATE_AWAY][$a];
	}

	##########
	### First, fill in the bins with ONLY their first hit, based on the offset, of course.
	### * This may skip bins. That will be dealt with later.
	##########
	for my $a (0 .. @epc_event_times-1) {
		$bin = &bin_no($epc_event_times[$a], $BL, ($BL*$offset) / $num_bins);
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

}

my ($mu, $stdev);
open OUT, ">$file/$BL-alt.Qcmp";
for my $a (0 .. $num_bins-5) {
	($mu, $stdev) = &mu_rho(\@{ $bin_first_events_AoA[$a] });
	my $bin = (($a+$offset)/$num_bins)*$BL; #($a + $offset) / ($num_bins*$BL);
	print OUT "$bin $mu $stdev\n";
}
close OUT;

#print "EPC:\n";
#($mu, $stdev) = &mu_rho(\@epc_num_events);
#print "average subst: $mu     standard deviation: $stdev\n";
#print OUT "$BL ";
#print OUT "$mu $stdev ";
#print "FWD:\n";
#($mu, $stdev) = &mu_rho(\@forward_num_events);
#print "average subst: $mu     standard deviation: $stdev\n";
#print OUT "$mu $stdev\n";

open OUT2, ">$file/$BL.bins";
for my $a (0 .. $num_bins-5) {
	print OUT2 "" . (($a/$num_bins)*$BL). " " . (int($epc_event_times_bin_array[$a])/$num_bins);
	print OUT2 " " . (int($forward_event_times_bin_array[$a])/$num_bins) . "\n";
}
close OUT2;
print "Binned times in file \"$file/$BL.bins\".\n";
open OUT2, ">$file/$BL.Qcmp";
for my $a (0 .. $num_bins-5) {
	print OUT2 (($a/$num_bins)*$BL);
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

open OUT2, ">$file/$BL-2.Qcmp";
for my $a (0 .. $num_bins-5) {
	print OUT2 (($a/$num_bins)*$BL);
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

&plot("$file/$BL");

sub plot
{
	my ($file) = @_;
	open OUT, ">$file.gnu";
	print OUT "set term postscript eps color\n";
	print OUT "set style line 1 lt 1 lw 3\n";
	print OUT "set style line 2 lt 3 lw 3\n";
	print OUT "set style line 3 lt 4 lw 3\n";
	print OUT "set style line 4 lt 5 lw 3\n";
	print OUT "set style line 5 lt 6 lw 3\n";
	print OUT "set title \"$file: seq_len:$sequence_length reps:$num_runs bins:$num_bins\"\n";
	print OUT "set xlabel \"Time on Branch\"\n";
	print OUT "set ylabel \"Rate Away\"\n";
	print OUT "set output \"$file.eps\"\n";
	print OUT "plot \"$file.Qcmp\" usi 1:2 w points ls 1 ti \"EPC Qi.\",";
	print OUT "\"$file.Qcmp\" usi 1:3 w points ls 2 ti \"EPC Qi.|k(t)\",";
	print OUT "\"$file.Qcmp\" usi 1:4 w points ls 3 ti \"FWD Qi.\",";
	print OUT "\"$file-alt.Qcmp\" usi 1:2 w points ls 4 ti \"EPC Qi.|k(t)-alt\"\n";
	close OUT;
	system("gnuplot $file.gnu");
	system("open $file.eps");
}

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