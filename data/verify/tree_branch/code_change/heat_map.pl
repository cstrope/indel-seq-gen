#!/usr/bin/env perl
use strict;

##########
### Plan: make heat map of subst. pr. for sequence at specific points during the run when
### fed 2 sequences that are 1 apart. Shows how our method changes the probability of each
### substitution as time progresses. Also, want to print out site-change probability for
### each site and Pij for the time points { begin, mid, end-dt }.
##########

if (@ARGV != 1) {
	die "Usage: perl heat_map.pl <random_number_seed>\n";
}
my $RANDOM_SEED = shift;

#### Branch:
##    i                         j
##  0 |-------------------------| T
####
my $increment = 0.001;
my $begin = 0.000;		### Start at very beginning of the branch
my $end = 0.999;		### Makes no sense to try to get 1 substitution when starting at T
my $branch_length = 1.000;	### T

my $anc = "CGTACGTACG";
my $des = "CGTACGGACG";
my @tmp = split //, $anc;
my $seq_size = (scalar @tmp);

my @heat_map = ();		## Variable to keep the heat map info (the Pr. of sequence Qi.)
my @Qidot = ();			## Just to make sure, get baseline. Should be same for all.

my $indel_seq_gen = "isg-test-one";

my $filename = "results_dir/heat_map_$branch_length-$increment";

########## Need .sim.dep and .sim.ma file for runs. #########
### Output JC69 dependency file. ###
open OUT, ">$filename.sim.dep";
for my $i (0 .. 63) {
	print OUT "" . (1/64) . " 0.25 0.25 0.25 0.25\n";
}
close OUT;
open OUT, ">$filename.sim.ma"; 
print OUT ">T_1\n$anc\n";
print OUT ">T_2\n$des\n";
close OUT;
########## Data for run is done. #########

for my $t ( ($begin * (1/$increment)) .. ($end * (1/$increment)) ) {
	print "----------------------" . ($t*$increment) . "---------------------------\n";
	## Run iSG
	open OUT, ">$filename.tree";
	print OUT "[$seq_size](T_1:" . ($branch_length - $t * $increment) . ",T_2:0);\n";
	close OUT;
	my $epc_command = "";
	$epc_command .= "-m JC69 -z $RANDOM_SEED,2001,3001,4001 -e $filename -O 3 -D $filename.sim.dep -E $filename.sim.ma ";
	$epc_command .= " -I 1";		# EPC sample type. -O order
	$epc_command .= " < $filename.tree";
	print "$epc_command\n ";
	my $return_val_epc = `./$indel_seq_gen $epc_command`;
	print "Done with run.....\n";
	
	## Gather output	
	my @each_line = split /\n/, $return_val_epc;
	$heat_map[$t] = $each_line[0];


}

open OUTPUT, ">$filename.results";
for my $i (0 .. @heat_map - 1) {
	print OUTPUT "" . ($branch_length - $i*$increment) . " $heat_map[$i]\n";
}
close OUTPUT;
