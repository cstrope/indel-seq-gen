#!/usr/bin/env perl
use strict;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $file = "run_" . $hour . "." . $min . "." . $sec;
my $num_runs = 10;
my $num_bins = $num_runs;

my $INDEPENDENT_SITES = 0;

my $return_val_fwd;
my $return_val_epc;
my @ret_split;


my $branch_length = 
#		 0.01,
#		 1.0,
#		 2.0,
		10.0
		;
my $dependence_superscript =
#		 0.01,
#		 0.1,
		 0.5,	#sqrt
#		 1,
#		 2,		#squared
#		 5,
		;
		
my $filename = "one_round_fwdepc-$branch_length-$dependence_superscript";
my $seq_size = 1000;
my $indel_seq_gen = "isg-test-one";

open OUTt, ">$filename.tree";
print OUTt "[$seq_size](T_1:$branch_length,T_2:0);\n";
close OUTt;

if ($INDEPENDENT_SITES) {
	$return_val_fwd = `./$indel_seq_gen -m JC69 -e $filename -o f < $filename.tree`;
	$return_val_epc = `./$indel_seq_gen -m JC69 -E $filename.sim.ma < $filename.tree`;
} else {
	# seq_size = 1000, BL = 1.0
	#  * events.size() = 422
	#  * Forward Probability: -3556.67
	$return_val_fwd = `./$indel_seq_gen -m JC69 -z 5009,2001,3001,4009 -D $dependence_superscript -e $filename -o f < $filename.tree`;
	# seq_size = 1000, BL = 1.0
	#  * events.size() = 444
	#  * EPC->Forward Probability: -3832.44
	$return_val_epc = `./$indel_seq_gen -m JC69 -z 5009,2001,3001,4009 -D $filename.sim.dep -E $filename.sim.ma < $filename.tree`;
}
