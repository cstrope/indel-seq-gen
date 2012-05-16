#/usr/bin/env perl
use strict;

##########
## Assess paths from iSG runs.
##########

{	######### MAIN #########
	if (@ARGV != 4) { die "Usage: perl assess_paths <delta_t> <path_filename> <first_cycle_to_sample> <sample_each_x_cycles>\n"; }
	my ($delta_t, $filename, $first_cycle_to_sample, $sample_each_x_cycles) = @ARGV;
	my ($numStates);
	$numStates = 4; 	## For now, assume nucleotides. Can change to codons or higher order contexts later. ##

	open my $infile, '<', $filename or die "Cannot open file $filename: $?\n";
	open IN, "$filename";
	my $tree = <$infile>;
	my $sequence_data = "";
	my $initial_cycle;
	print "$tree\n";

	my ($t_0, $T) = (0, 1);		## Known, for the time being. ##
	my $time_bins = 1000;
	my $time_increment = ($T - $t_0) / $time_bins;

	($sequence_data, $initial_cycle) = &get_next_path($infile);

	my ($i_0, $k_0, @work, @each_seq);
	@each_seq = split />.+\n/, $sequence_data;

	$i_0 = $each_seq[1];
	$k_0 = $each_seq[2];

	my ($sequence_size);
	@work = split //, $i_0;
	$sequence_size = scalar @work;

	##########
	## Build AVERAGE PATH array. 3D array, $average_path[D1][D2][D3], where
	## D1 = time_bin
	## D2 = Sequence_site
	## D3 = state at site (A,C,G,T) for 0-th order MM, (A|A, A|C, ..., T|G, T|T) for 1-st order MM, etc.
	##########
	my (@average_path, @array);
	for my $i ( 0 .. $time_bins ) {
		$average_path[$i] = [ @array ];
		for my $j (0 .. $sequence_size-1) {
			$average_path[$i][$j] = [ @array ];
			for my $k (0 .. $numStates) {
				$average_path[$i][$j][$k] = 0;
			}
		}
	}

	##########
	### Have basic data structures, now need to put together the MCMC runs.
	##########
	my (
		$nextpath_cycle,## When < current_cycle, stored path is invalid. New path needed.
		@tmp, 			## Exactly what it says it is.
		$current_cycle, ## Current cycle, value will be (($cc-$first_cycle_to_sample)% sample_x_cycles == 0)
		$path, 			## Path previous to $next_cycle path.
	   );
	$current_cycle = $nextpath_cycle = $initial_cycle;

	## Read paths until reaching the first potential cycle to sample.
	while($nextpath_cycle < $first_cycle_to_sample and !eof($infile)) {
		($path, $nextpath_cycle) = &get_next_path($infile);
	}
	## In the end, set current_cycle to the first cycle to sample. This is because there may be
	## non-accepted paths before the first true path to sample, which could potentially skip a
	## number of possible samples of the "current" path state.
	$current_cycle = $first_cycle_to_sample;
	($path, $nextpath_cycle) = &get_next_sample($infile, $current_cycle, $sample_each_x_cycles);

	## Next sample is set to -1 when there are no more paths to be read from MCMC
	## After reading to the next path, we will have the next pathID (cycle for which it was accepted),
	## and the previous path (since next pathID will be larger than the current cycle).
	while ($nextpath_cycle != -1) {
		### Set the work sequence to the initial state ###
		@work = split //, $i_0;
		### Get the next path ###
		($path, $nextpath_cycle, $current_cycle) 
		= &get_next_sample(
						   $infile, 
						   $current_cycle, 
						   $sample_each_x_cycles, 
						   $path, 
						   $nextpath_cycle
						  );
		## Should have path to profile && cycle of interest.
		print "$current_cycle $nextpath_cycle\n"; 
		&add_to_profile($i_0, $t_0, $T, $path, \@average_path, $time_increment, $time_bins);
	}
	
	&normalize_profile(\@average_path, $time_bins, $sequence_size, $numStates);
	&print_profile(\@average_path);
}	### End main ###

sub print_profile
{
	my ($prof_ref, $time_bins, $sequence_size, $numStates) = @_;
	
	open OUT, ">avg_path.stat";

	for my $i ( 0 .. $time_bins ) {
		for my $j (0 .. $sequence_size-1) {
			my $num_hits = 0;
			for my $k (0 .. $numStates) {
				print OUT "$i $j $k $prof_ref->[$i][$j][$k]\n";
			}
		}
	}
	close OUT;
}

sub normalize_profile
{
	my ($prof_ref, $time_bins, $sequence_size, $numStates) = @_;

	for my $i ( 0 .. $time_bins ) {
		for my $j (0 .. $sequence_size-1) {
			my $num_hits = 0;
			for my $k (0 .. $numStates) {
				$num_hits += $prof_ref->[$i][$j][$k];
			}

			if ($num_hits != 0) {
				for my $k (0 .. $numStates) {
					$prof_ref->[$i][$j][$k] /= $num_hits;
				}
			}
		}
	}
}

sub add_to_profile
{
	my (
		$i_0,				## Begin sequence.
		$t_0, 				## Start time of i_0 (begin time before time_bin increments).
		$T,					## End time.
		$path_to_add, 		## Path from $i_0--->$k_0
		$profile_array_ref, ## 3D array, from above.
		$dt,				## Amount of time to increment at each step.
		$total_time_bins	## To increment the array at the proper point.
	   ) = @_;
	
	my @path_events = split /\n/, $path_to_add;
	pop @path_events;		## Branch ending event. ##
	my @work_sequence = split //, $i_0;
	my $last_event_ptr = 0;
	my ($prof_site_idx) = (0);

	my ($t, $change, $Qidot, $site, $Qidot_k);
	($t, $site, $change, $Qidot, $Qidot_k) = split(",", $path_events[$last_event_ptr]);
	for (my $time_bin = $t_0; $last_event_ptr < (scalar @path_events); $time_bin += $dt) {
		## Catch work sequence up to the times.
		while ($t < $time_bin and $last_event_ptr < (scalar @path_events) ) {
			## Make change to sequence
			my @nucl_pair = split //, $change;
			if ($work_sequence[$site] ne $nucl_pair[0]) { die "Odd... $work_sequence[$site] != $nucl_pair[0].\n"; } 
			else { $work_sequence[$site] = $nucl_pair[1]; }

			&increment_profile($profile_array_ref, \@work_sequence, $time_bin * $total_time_bins);

			## Finished event, get data of next event.
			$last_event_ptr++;
			($t, $site, $change, $Qidot, $Qidot_k) = split(",", $path_events[$last_event_ptr]);
			#print "$last_event_ptr: $t $site $change $Qidot $Qidot_k\n";
		}
	}
}

sub increment_profile
{
	my ($prof_ref, $seqdata_ref, $time_idx) = @_;
	
	for my $i (0 .. @{ $seqdata_ref }-2) {
		$prof_ref->[$time_idx] [$i] [&pattern_idx($seqdata_ref->[$i])]++;
		#print "prof_ref->[$time_idx][$i][(&pattern_idx($seqdata_ref->[$i]))]++\n";
	}
}

sub pattern_idx
{
	my ($pattern) = @_;
	my @elements = split //, $pattern;
	
	## Naturally, this will have to be better done for higher order Markov models.
	if ($elements[0] eq "A") { return 0; }
	if ($elements[0] eq "C") { return 1; }
	if ($elements[0] eq "G") { return 2; }
	if ($elements[0] eq "T") { return 3; }
	else { die "Unknown element for pattern index: \"$elements[0]\"\n"; }
}

#########
### This will set the file to the appropriate next sample.
#########
sub get_next_sample
{
	my ($fh, $current_cycle, $sample_each_x_cycles, $path, $next_path_cycle) = @_;

	$current_cycle += $sample_each_x_cycles;

	## Sample falls beyond current cycle.
	if ($current_cycle > $next_path_cycle) {
		while ($next_path_cycle < $current_cycle && !eof($fh)) {
			($path, $next_path_cycle) = &get_next_path($fh);
		}
	}

	if (eof($fh)) { $next_path_cycle = -1; }

	return ($path, $next_path_cycle, $current_cycle);
}

sub get_next_path 
{
	my $fh = shift;
	my ($path, $pathID);
	my $line;
	
	my $cont = 1;
	while ($cont and $line = <$fh>) {
		if ($line !~ m/^\#/) {
			$path .= $line;
		} else { $cont = 0; }
	}

	my @tmp = split /\#/, $line;
	$pathID = $tmp[1];
	chomp($pathID);

	return ($path, $pathID);	
}

sub read_to_next_sample