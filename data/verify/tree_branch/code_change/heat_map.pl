#!/usr/bin/env perl
use strict;

##########
### Plan: make heat map of subst. pr. for sequence at specific points during the run when
### fed 2 sequences that are 1 apart. Shows how our method changes the probability of each
### substitution as time progresses. Also, want to print out site-change probability for
### each site and Pij for the time points { begin, mid, end-dt }.
##########

if (@ARGV != 0) {
	die "Usage: perl heat_map.pl \n";
}

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

my @heat_map = ();		## Variable to keep the heat map info (the Pr. of sequence Qi.)
my @Qidot = ();			## Just to make sure, get baseline. Should be same for all.

my $indel_seq_gen = "isg-heat-map";

########## Need .sim.dep and .sim.ma file for runs. #########
### Output JC69 dependency file. ###
open OUT, ">$filename.sim.dep";
for my $i (0 .. 63) {
	print OUT "" . (1/64) . " 0.25 0.25 0.25 0.25\n";
}
close OUT;
open OUT, ">$filename.sim.ma"; 
print OUT ">T1\n$anc\n";
print OUT ">T2\n$des\n";
close OUT;
########## Data for run is done. #########

	my $epc_command = "";
	$epc_command .= "-m JC69 -z $RANDOM_SEED,2001,3001,4001 -e $filename -O $order -D $filename.sim.dep -E $filename.sim.ma ";
	$epc_command .= " -I $EPC_SAMPLE_TYPE";
	if ($MCMC) {
		$epc_command .= " --mcmc $MCMC ";
		if ($EVEN_SAMPLING) {
			$epc_command .= " --mcmc_sample_evenly ";
		}
	}
	$epc_command .= " < $filename.tree";
	print "$epc_command\n ";
	$return_val_epc = `./$indel_seq_gen $epc_command`;
	open OUT, ">$filename.results";
	print OUT "$return_val_fwd\n$return_val_epc\n";
	close OUT;

for my $t ($begin .. $end) {
	## Run iSG
	
	## Gather output
	
	## Parse output

}