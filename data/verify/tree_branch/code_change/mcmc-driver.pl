#!/usr/bin/env perl
use strict;

if (@ARGV != 7) {
	print "Usage: perl mcmc-driver.pl <dependence_superscript> <num_mcmc_cycles> <sequence_length> <branch_length> <XXX> <sample_evenly> <random_seed>.\n";
	die "XXX:\nrasmus=0\nQdPc=1\nQdP=2\nQPc=3\nQP=4\n";
}

my $dependence_superscript = shift;
my $MCMC = shift;
my $seq_size = shift;
my $branch_length = shift;
my $EPC_SAMPLE_TYPE = shift;
my $EVEN_SAMPLING = shift;
my $RANDOM_SEED = shift;

if (!$RANDOM_SEED) { $RANDOM_SEED = 9966; }

my $num_replicates = 1;
my $order = 3;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $file = "run_" . $hour . "." . $min . "." . $sec;
my $num_runs = 1;

my $INDEPENDENT_SITES = 0;
my $EMULATE = 0;

my $return_val_fwd;
my $return_val_epc;

my @ret_split;

#my $branch_length = 
#		 0.01,
#		 0.5,
#		 0.95,
#		 1.0,
#		 2.0,
#		 5.0,
#		10.0
		;
		
my $filename = "results_dir/mcmc$MCMC-bl$branch_length-ds$dependence_superscript-$EPC_SAMPLE_TYPE$EVEN_SAMPLING-$RANDOM_SEED";
my $indel_seq_gen = "isg-test-one";

my $twoBL = $branch_length/2.0;
my $guide_tree = 
#	"[$seq_size](T_1:$branch_length, T_2:0);\n";		### Standard previous: 1 branch only.
	"[$seq_size](T_1:$twoBL, T_2:$twoBL);\n";	### 2 branches, each with length. Same as prev, but minor improvement.

open OUTt, ">$filename.tree";
print OUTt "$guide_tree";
#print OUTt "[$seq_size](T_1:$branch_length,T_2:0);\n";
close OUTt;

if ($INDEPENDENT_SITES) {
	$return_val_fwd = `./$indel_seq_gen -m JC69 -e $filename -o f < $filename.tree`;
	$return_val_epc = `./$indel_seq_gen -m JC69 -E $filename.sim.ma < $filename.tree`;
} else {
	my $fwd_command = "";
	$fwd_command .= "-m JC69 -z $RANDOM_SEED,2001,3001,4001 -O $order -D $dependence_superscript -e $filename -o f ";
	$fwd_command .= " -n $num_runs ";		## This is only to test if fwd can be run multiple times. Does not yet work with EPC following. ##
	$fwd_command .= "< $filename.tree";
	print "$fwd_command\n ";
	$return_val_fwd = `./$indel_seq_gen $fwd_command`;
	open OUT, ">$filename.fwd";
	print OUT "$return_val_fwd\n";
	close OUT;

	#exit(0);		## While testing new dependency model. ##

	print STDERR "----- END-POINT CONDITIONED RUN -------------\n";
	my $epc_command = "";
	$epc_command .= "-m JC69 -z $RANDOM_SEED,2001,3001,4001 -e $filename -O $order -D $filename.sim.dep -E $filename.sim.ma ";
	if ($EMULATE) {
		$epc_command .= " -M $filename.sim.trace ";
	} else {
		$epc_command .= " -n $num_replicates ";
	}
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
}

print STDERR "True path forward probability: $return_val_fwd\n";
print STDERR "Reps: $num_replicates\n";
print STDERR "Power: $dependence_superscript\n";
print STDERR "MCMC: Last $MCMC lines in file $filename.results.\n";

print " branch length:   $branch_length\n";
print " sequence length: $seq_size\n";
print " depsup:          $dependence_superscript\n";
print " mcmc style:      "; 
if($EPC_SAMPLE_TYPE == 0) { print "rasmus"; } 
else { 
	if ($EPC_SAMPLE_TYPE == 1) { print "QdPc"; }
	elsif ($EPC_SAMPLE_TYPE == 2) { print "QdP"; }
	elsif ($EPC_SAMPLE_TYPE == 3) { print "QPc"; }
	elsif ($EPC_SAMPLE_TYPE == 5) { print "QN"; }
	else { print "QP"; }

	print "-subpath, P(time = T/2) = "; 
	if ($EVEN_SAMPLING) { print "0.333"; } 
	else { print "0.5"; } 
} 
print "\n";
print " mcmc_cycles:     $MCMC\n";
print " seed:            $RANDOM_SEED\n";
