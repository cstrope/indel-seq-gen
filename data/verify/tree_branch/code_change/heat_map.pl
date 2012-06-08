#!/usr/bin/env perl
use strict;

### For this one, need to use iSG-test-one, for reasons pointed in propose_path.cpp.

##########
### Plan: make heat map of subst. pr. for sequence at specific points during the run when
### fed 2 sequences that are 1 apart. Shows how our method changes the probability of each
### substitution as time progresses. Also, want to print out site-change probability for
### each site and Pij for the time points { begin, mid, end-dt }.
##########

if (@ARGV != 3) {
	die "Usage: perl heat_map.pl <random_number_seed> <seqlen> <numdiff>\n";
}
my $RANDOM_SEED = shift;
my $seqlen = shift;
my $numdiff = shift;

my $epc_sample_type = 1;

#### Branch:
##    i                         j
##  0 |-------------------------| T
####
my $increment = 0.001;
my $begin = 0.998;		### Start at very beginning of the branch
my $end = 0.999;		### Makes no sense to try to get 1 substitution when starting at T
my $branch_length = $end+$increment;	### T

my $anc = "";
my $des = "";
for my $L (1 .. $seqlen) { $anc .= "A"; }
for my $L (1 .. $numdiff) { $des .= "C"; }
for my $L ($numdiff+1 .. $seqlen) { $des .= "A"; }
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

open OUTPUT, ">$filename-$seqlen-$numdiff.results";
print OUTPUT "dat-$seqlen-$numdiff";
print OUTPUT "\n";
for my $t ( ($begin * ($branch_length/$increment)) .. ($end * (($branch_length/$increment))) ) {
	print "----------------------$t---------------------------\n";
	## Run iSG
	open OUT, ">$filename.tree";
	print OUT "[$seq_size](T_1:" . ($branch_length - $t*$increment) . ",T_2:0);\n";
	close OUT;
	my $epc_command = "";
	$epc_command .= "-m JC69 -z $RANDOM_SEED,2001,3001,4001 -d 000000 -e $filename -O 3 -D $filename.sim.dep -E $filename.sim.ma ";
	$epc_command .= " -I $epc_sample_type";		# EPC sample type. -O order
	$epc_command .= " < $filename.tree";
	print "$epc_command\n ";
	my $return_val_epc = `./$indel_seq_gen $epc_command`;
	print "Done with run.....\n";
	
	## Gather output	
	my @each_line = split /\n/, $return_val_epc;
	#print OUTPUT "" . ($branch_length - $t*$increment) . ",$each_line[0]";
	print OUTPUT "$each_line[0]\n";
	print "X$each_line[0] X$each_line[1] X$each_line[2] ...\n";
	$heat_map[$t] = $each_line[0];
}

close OUTPUT;

