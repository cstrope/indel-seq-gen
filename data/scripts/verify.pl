#!/usr/bin/env perl
use Data::Dumper;

$/=undef;
use strict;

my @params = (
			  "sim_conditions",
			  "bl_choice",
			  "CpG_multiplier",
			  "num_reps",
			  "root_len",
			);

if (@ARGV != @params) {
	print "Usage: ./verify.pl";
	for my $i (0 .. @params-1) {
		print " <$params[$i]>";
	}
	print ".\n";
	exit(0);
}

my $sim_conditions = $ARGV[0];
my $bl_multiplier = $ARGV[1];
my $CpG_multiplier = $ARGV[2];
my $num_reps = $ARGV[3];
my $root_length = $ARGV[4];

##################
### This script will make 2 seq datasets, run path proposal, and gnuplot results.
##################

### Old, no need to use.
my $PLOT = 0;

### Phylogenetic trees. See function outputTreefile at the end of script.
<<<<<<< .mine
<<<<<<< .mine
my $CALC_LIKELIHOOD = 0;
=======
my $CALC_LIKELIHOOD = 4;
=======
my $CALC_LIKELIHOOD = 0;
>>>>>>> .r1895
>>>>>>> .r1886
	my $ZERO_BL_TREE = 0;
		my $ZERO_BL_TREE2 = 0;
my $PRUNING = 0;
my $SMALL = 0;

my $UNSET_NODES = $CALC_LIKELIHOOD + $PRUNING + $SMALL;
### This was originally used in conjunction with PLOT (above). This should always be set to 1.
<<<<<<< .mine
my $TEST = 1000;
=======
my $TEST = $num_reps;
>>>>>>> .r1895
### Root sequence length.
my $test_rootseq_length = $root_length;

### What we want to do. ENDPOINT_COND specifies if we want do to endpoint conditioned runs. 
### FORWARD_SIM does nothing.
my $FORDWARD_SIM = 1;
my $ENDPOINT_COND = 1;
my $OUTPUT_NUCL_STATS = 0;

my $MAX_BIN = 0;

my ($slurp, @tmp, @tmp1, @tmp2);
# print Dumper(\@y_axis_stats);
# print Dumper(\@y_axis_stats[0]);

my $num_bins = 100;

############################################################################
### Below is a listing of the arrays. Only one element in each array should be un-commented
### since I am unsure of the behavior of the programs if more than one is un-commented.
############################################################################


############################################################################
### These arrays correspond with one another. Un-comment the same element of each array.
my @matrix = (
			  "JC69",
			  "JC69",
			  "HKY",
			  "F81",		### Successful: rev. 1816, 4-23-11. 10000 reps, 
						    ###   * NoCpG: 2999.0430
						    ###   * CpG:   
			  "F81",		###
						    ###   * NoCpG: 2999.4942
						    ###   * CpG:   2943.xxxx
			  "F81",
			  "GTR",
			  "SYM",
			 );
my @path_freqs =
		(
		 "0.25 0.25 0.25 0.25",
		 "0.25 0.25 0.25 0.25",
		 "0.3 0.2 0.2 0.3",
		 "0.05 0.05 0.05 0.85",
		 "0.1 0.4 0.4 0.1",
		 "0.05 0.05 0.05 0.85",
		 "0.1 0.2 0.3 0.4",
		 "0.25 0.25 0.25 0.25",
		);
my @iSG_freqs =
		(
		 "0.25,0.25,0.25,0.25",
		 "0.25,0.25,0.25,0.25",
		 "0.3,0.2,0.2,0.3",
		 "0.05,0.05,0.05,0.85",
		 "0.1,0.4,0.4,0.1",
		 "0.05,0.05,0.05,0.85",
		 "0.1,0.2,0.3,0.4",
		 "0.25,0.25,0.25,0.25",
		);
my @iSG_rel_rates = 
		(
		 "",
		 "",
		 "-r 2",
		 "",
		 "",
		 "",
		 "-r 2,3,2,2,3,1",
		 "-r 2,3,2,2,3,1",
		);
my @path_rel_rates = 
		(
		 "",
		 "",
		 "2",
		 "",
		 "",
		 "",
		 "2 3 2 2 3 1",
		 "2 3 2 2 3 1",
		);
my @iSG_gamma_cats =
		(
		 "",
		 "-g 4 -a 0.5",
		 "-g 4 -a 0.5",
		 "",
		 "-g 4 -a 0.5",
		 "-g 4 -a 0.5",
		 "",
		 "-g 4 -a 0.5",
		 "",
		);
my @path_gamma_cats =
		(
		 "",
		 "-g 4 -a 0.5",
		 "-g 4 -a 0.5",
		 "",
		 "4 0.5",
		 "4 0.5",
		 "",
		 "4 0.5",
		 "",
		);

my @junkfile_sim = 
		(
		 "test_0s",
		 "test_01s",
		 "test_1s",
		 "yyys",
		 "aaas",
		 "zzzs",
		 "zzzs",
		 "xxxs",
		 "test_4s",
		);
my @junkfile_pro = 
		(
		 "test_0p",
		 "test_01p",
		 "test_1p",
		 "yyyp",
		 "aaap",
		 "zzzp",
		 "zzzp",
		 "xxxp",
		 "test_4p",
		);
my $treefile;

############################################################################
### These arrays are set separately than the above ones. Self-explanatory names?
my @branch_length = 
		(

		 "0.1",
		 "0.5",
		 "1.0",
		 "2.0",
		 "3.0",
		 "3.9",
		);

my $rootseq_length;
my ($root_seq, $target_seq);
my $pathfile;
my (@simulated_subst, @proposed_subst);
my $num_replicates;
if ($TEST) { $num_replicates = $TEST; $rootseq_length = $test_rootseq_length; } 
else { $num_replicates = 100; $rootseq_length = 1000; }
my @simulated_A;
my @simulated_C;
my @simulated_G;
my @simulated_T;
my @proposed_A;my @proposed_C;my @proposed_G;my @proposed_T;
my @proposed_AA;my @proposed_AC;my @proposed_AG;my @proposed_AT;my @proposed_CA;my @proposed_CC;my @proposed_CG;my @proposed_CT;my @proposed_GA;my @proposed_GC;my @proposed_GG;my @proposed_GT;my @proposed_TA;my @proposed_TC;my @proposed_TG;my @proposed_TT;my @simulated_AA;my @simulated_AC;my @simulated_AG;my @simulated_AT;my @simulated_CA;my @simulated_CC;my @simulated_CG;my @simulated_CT;my @simulated_GA;my @simulated_GC;my @simulated_GG;my @simulated_GT;my @simulated_TA;my @simulated_TC;my @simulated_TG;my @simulated_TT;
my $running_total_simsubs = 0, 
my $running_total_prosubs=0;

my (@simulated_subst_per_model_condition);
my (@sim_CpG, @prop_CpG);

my (@output_data, @array, @bin_hits);
my (@each_subst);
my (@epc_root_rate_away,@epc_root_rate_away2);
my (@sim_cat_mult, @epc_cat_mult);
my (@sim_CpG_per_run, @epc_CpG_per_run);
my (@nucl_freq);	# Nucleotide frequencies. #

for my $i (0 .. @junkfile_sim-1 ) {
	$simulated_subst_per_model_condition[$i] = 0;
}

for my $k (0 .. $num_replicates-1) {
	$simulated_subst[$k] = 0;
	$proposed_subst[$k] = 0;
	$simulated_A[$k] = 0;$simulated_C[$k] = 0;$simulated_G[$k] = 0;$simulated_T[$k] = 0;$proposed_A[$k] = 0;$proposed_C[$k] = 0;$proposed_G[$k] = 0;$proposed_T[$k] = 0;$simulated_AA[$k] = 0;$simulated_AC[$k] = 0;$simulated_AG[$k] = 0;$simulated_AT[$k] = 0;$proposed_AA[$k] = 0;$proposed_AC[$k] = 0;$proposed_AG[$k] = 0;$proposed_AT[$k] = 0;$simulated_CA[$k] = 0;$simulated_CC[$k] = 0;$simulated_CG[$k] = 0;$simulated_CT[$k] = 0;$proposed_CA[$k] = 0;$proposed_CC[$k] = 0;$proposed_CG[$k] = 0;$proposed_CT[$k] = 0;$simulated_GA[$k] = 0;$simulated_GC[$k] = 0;$simulated_GG[$k] = 0;$simulated_GT[$k] = 0;$proposed_GA[$k] = 0;$proposed_GC[$k] = 0;$proposed_GG[$k] = 0;$proposed_GT[$k] = 0;$simulated_TA[$k] = 0;$simulated_TC[$k] = 0;$simulated_TG[$k] = 0;$simulated_TT[$k] = 0;$proposed_TA[$k] = 0;$proposed_TC[$k] = 0;$proposed_TG[$k] = 0;$proposed_TT[$k] = 0;
}

my @filetmp = split /\s+/, $path_gamma_cats[$sim_conditions];
my $file_root = "$matrix[$sim_conditions]" 
				. "_$bl_multiplier" 
				. "_$CpG_multiplier" 
				. "_$filetmp[0]-$filetmp[1]";
$treefile  = $file_root . ".tree";
$pathfile  = $file_root . ".path";
$running_total_simsubs = $running_total_prosubs = 0;
for my $k (0 .. $num_replicates-1) {
	sleep(0.1);
	print STDERR "********$matrix[$sim_conditions], $bl_multiplier******rep $k/$num_replicates*************\n";
	&outputTreefile($treefile, $bl_multiplier);
	my $command_line_options = "";
	$command_line_options .= "-m $matrix[$sim_conditions] ";
	$command_line_options .= "-s $rootseq_length ";
	$command_line_options .= "$iSG_gamma_cats[$sim_conditions] ";
	$command_line_options .= "$iSG_rel_rates[$sim_conditions] ";
	$command_line_options .= "-f $iSG_freqs[$sim_conditions] ";
	$command_line_options .= "-C $CpG_multiplier "; 
	$command_line_options .= "-w a ";
	$command_line_options .= "-d 001010 ";
	$command_line_options .= "-e $file_root ";
	$command_line_options .= "-o f ";
	$command_line_options .= "< $treefile > $junkfile_sim[$sim_conditions]";
	print STDERR "./indel-seq-gen $command_line_options\n";
	system("./indel-seq-gen $command_line_options");
	##########
	### Gather MSA sequences to output
	##########
	my @seqs;

	if ($UNSET_NODES) {		# Basically, there is some sort of pruning going on.
		&getPruningSequences("$file_root.sim.ma", \@seqs);
	} else {
		open IN, "$file_root.sim.ma";
		$slurp = <IN>;
		close IN;
		@tmp = split(">", $slurp);	
		for my $g (1 .. @tmp-1) {
			@tmp1 = split /\n/, $tmp[$g];
			$seqs[$g-1] = ">$tmp1[0]\n$tmp1[1]";
		}
		##########
		### Count the nucleotide frequencies of the descendant sequence.
		##########
		my @tmp_seq = split /\n/, $tmp[2];
		my @sequence = split //, $tmp_seq[1];
		for my $q (0 .. @sequence-1) {
			if ($sequence[$q] eq "A") {
				$nucl_freq[0]++;
			} elsif ($sequence[$q] eq "C") {
				$nucl_freq[1]++;
			} elsif ($sequence[$q] eq "G") {
				$nucl_freq[2]++;
			} elsif ($sequence[$q] eq "T") {
				$nucl_freq[3]++;
			} else {
				print "Skipped $sequence[$q]\n";
				exit(0);
			}
		}
	}
	$simulated_subst[$k]  
	= $simulated_subst_per_model_condition[$sim_conditions] 
	= &getSubstitutions("$file_root.trace", "FORWARD");
	$running_total_simsubs += $simulated_subst[$k];
	print STDERR "Gathered simsubs: $simulated_subst[$k]\n";
	print STDERR "  Running total simsubs: " . ($running_total_simsubs/($k+1)) . "\n";
	open IN, "$junkfile_sim[$sim_conditions]";
	$slurp = <IN>;
	close IN;
	
	my @tmp20 = split /\n/, $slurp;
	$sim_cat_mult[$k] = $tmp20[0];
	$sim_CpG_per_run[$k] = $tmp20[1];

	if ($ENDPOINT_COND) {
		&outputPathFile($pathfile, \@seqs);
		my $command_line_options = "";
		$command_line_options .= "-m $matrix[$sim_conditions] ";
		$command_line_options .= "$iSG_gamma_cats[$sim_conditions] ";
		$command_line_options .= "$iSG_rel_rates[$sim_conditions] ";
		$command_line_options .= "-f $iSG_freqs[$sim_conditions] ";
		$command_line_options .= "-d 001010 ";
		$command_line_options .= "-C $CpG_multiplier "; 
		$command_line_options .= "-e path_proposal ";
		$command_line_options .= "-E $file_root.path ";
		$command_line_options .= "< $treefile > $junkfile_pro[$sim_conditions]";
		print STDERR "./indel-seq-gen $command_line_options\n";
		system("./indel-seq-gen $command_line_options");

		$proposed_subst[$k] = &getSubstitutions("path_proposal.trace", "EPC");
		$running_total_prosubs += $proposed_subst[$k];
		print STDERR "Gathered prosubs: $proposed_subst[$k]\n";
		print STDERR "  Running total prosubs: " . ($running_total_prosubs/($k+1)) . "\n";

		open IN, "$junkfile_pro[$sim_conditions]";
		$slurp = <IN>;
		close IN;

		@tmp1 = split (">>", $slurp);		## This splits for the rate away from each node (from output at end of EvolveStep in propose_path.cpp, printing before entering loop.).
		my @tmp20 = split /\s+/, $tmp1[1];	## Gather the root node rate away.
		$epc_root_rate_away[$k] = $tmp20[1];
		@tmp20 = split /\s+/, $tmp1[4];	## Gather the root node rate away.
		$epc_root_rate_away2[$k] = $tmp20[1];
		$epc_cat_mult[$k] = $tmp1[0];
		@tmp20 = split /\s+/, $tmp1[5];	## Gather avg_cpgs..
		$epc_CpG_per_run[$k] = $tmp20[0];

	}
}

my ($avg_sim, $stdev_sim) = &mu_rho(\@simulated_subst);
my ($avg_sim_cat_mult, $stdev_sim_cat_mult) = &mu_rho(\@sim_cat_mult);
my ($avg_sim_CpG, $stdev_sim_CpG) = &mu_rho(\@sim_CpG_per_run);
print "sim:  ($avg_sim, $stdev_sim)\n";
print "  category_multiplier: ($avg_sim_cat_mult, $stdev_sim_cat_mult)\n";
print "  Average CpG per step: ($avg_sim_CpG, $stdev_sim_CpG)\n";

if ($ENDPOINT_COND) {
	my ($avg_epc, $stdev_epc) = &mu_rho(\@proposed_subst);
	my ($avg_root_lambda_T, $stdev_root_lambda_T) = &mu_rho(\@epc_root_rate_away);
	my ($avg_epc_CpG, $stdev_epc_CpG) = &mu_rho(\@epc_CpG_per_run);
	my ($avg_epc_cat_mult, $stdev_epc_cat_mult) = &mu_rho(\@epc_cat_mult);
	my ($avg_root_lambda_T2, $stdev_root_lambda_T2) = &mu_rho(\@epc_root_rate_away2);
	print "epc:  ($avg_epc, $stdev_epc)\n";
	print "  category_multiplier: ($avg_epc_cat_mult, $stdev_epc_cat_mult)\n";
	print "  Average CpG per step: ($avg_epc_CpG, $stdev_epc_CpG)\n";
	print "lamda_T: ($avg_root_lambda_T, $stdev_root_lambda_T)\n";
	print "lamda_T2: ($avg_root_lambda_T2, $stdev_root_lambda_T2)\n";
}

print "Nucleotide Frequencies: \n";
print "  A: " . ($nucl_freq[0] / $num_reps) . "\n";
print "  C: " . ($nucl_freq[1] / $num_reps) . "\n";
print "  G: " . ($nucl_freq[2] / $num_reps) . "\n";
print "  T: " . ($nucl_freq[3] / $num_reps) . "\n";

sub mu_rho
{
	my ($ref) = @_;
	
	my ($average, $stdev);
	
	for my $i (0 .. @{ $ref }-1) {
		$average += $ref->[$i];
#		print "$i " . ($ref->[$i]) . "\n";
	}

#	print "$average " . (scalar(@{ $ref })) . "\n";

	$average /= scalar(@{ $ref });#-1;

#	print "$average " . (scalar(@{ $ref })) . "\n";

	for my $i (0 .. @{ $ref }-1) {
		$stdev += ($average - $ref->[$i]) * ($average - $ref->[$i]);
	}
	
#	print "STDEV:\n";
#	print "Before dividing by replicates: ";
#	print "$stdev " . (scalar(@{ $ref })) . "\n";

	$stdev /= scalar(@{ $ref });#-1;

#	print "After dividing by replicates:   ";
#	print "$stdev " . (scalar(@{ $ref })) . "\n";
	return ($average, sqrt($stdev));
}

sub outputPathFile
{
	my ($pathfile, $seqs_ref) = @_;
	
	open OUT, ">$pathfile" or die "WTF2??\n";
	for my $i (0 .. @{ $seqs_ref }) {
		print OUT ($seqs_ref->[$i]) . "\n";
	}
	close OUT;

	print $pathfile . "\n";
}

sub outputTreefile
{
	my ($treefile, $branch_length) = @_;

	open OUT, ">$treefile" or die "WTF??\n";
	print OUT "[$rootseq_length]";
	if ($CALC_LIKELIHOOD) {
		##########
		### Gamma category test routine. This works primarily with the code in treefile.cpp,
		### calculateJCLikelihood.
		##########
		my $b2 = $branch_length;
		print STDERR "Total tree length = " . (3*$b2) . "\n";
		if ($ZERO_BL_TREE) {
			if ($ZERO_BL_TREE2) {
				print OUT "((T_1:$b2,T_2:$b2):0, T_3:$b2);\n";
			} else {
				print OUT "((T_1:$b2,T_2:$b2):$b2, T_3:0);\n";
			}
		} else {
			if ($CALC_LIKELIHOOD == 2) {
				print OUT "((T_1:$b2,T_2:$b2):" . ($b2 * 4.0/5) . ", T_3:" . ($b2 * 1.0/5) . ");\n";
			} elsif ($CALC_LIKELIHOOD == 3) {
				print OUT "((T_1:$b2,T_2:$b2):" . ($b2 * 1.0/5) . ", T_3:" . ($b2 * 4.0/5) . ");\n";
			} elsif ($CALC_LIKELIHOOD == 4) {
				print OUT "((T_1:$b2,T_2:$b2):" . ($b2 * 1.0/100) . ", T_3:" . ($b2 * 99.0/100) . ");\n";
			} elsif ($CALC_LIKELIHOOD == 5) {
				print OUT "((T_1:$b2,T_2:$b2):" . ($b2 * 99.0/100) . ", T_3:" . ($b2 * 1.0/100) . ");\n";
			} else {
				print OUT "((T_1:$b2,T_2:$b2):" . ($b2/2) . ", T_3:" . ($b2/2) . ");\n";
			}
		}
	} elsif ($PRUNING == $SMALL and $PRUNING) {
		##########
		### This tree is based on the tree from Ziheng Yangs omputational Molecular Evolution book,
		### in the pruning methodology chapter.	Sum of all branches is 1.3 units * BL.
		##########
		my $b1 = 0.33333*$branch_length;
#		my $b2 = 0.66666*$branch_length;
		my $b2 = $branch_length;
#		print STDERR "Total tree length = " . (1.3*$branch_length) . "\n";
		print OUT "((T_1:$b2,T_2:$b2):" . ($b2/2) . ",T_3:" . ($b2/2) . ");\n";
	} elsif ($PRUNING) {
		##########
		### This tree is based on the tree from Ziheng Yangs omputational Molecular Evolution book,
		### in the pruning methodology chapter.	Sum of all branches is 1.3 units * BL.
		##########
		my $b1 = 0.2*$branch_length;
		my $b2 = 0.1*$branch_length;
		print STDERR "Total tree length = " . (1.3*$branch_length) . "\n";
		print OUT "(((T_1:$b2,T_2:$b2):$b1,T_3:$b2):$b1,(T_4:$b2,T_5:$b2):$b1);\n";
	} else {
		print STDERR "Total tree length = $branch_length\n";
		print OUT "(T_1:$branch_length,T_2:0);\n";
	}
	close OUT;
}

sub getSubstitutions
{
	my ($filename, $method) = @_;

	my (@tmp, @tmp1, $slurp, $num_subst);

	open IN, "$filename";
	$slurp = <IN>;
	close IN;
	if ($method eq "FORWARD") {
		@tmp = split /\n\n/, $slurp;
		@tmp1 = split /\n/, $tmp[1];
	} else {
		@tmp1 = split /\n/, $slurp;
	}
	for my $a (0 .. @tmp1-3) { $num_subst++; }

	return $num_subst;
}

sub getPruningSequences 
{
	my ($infile, $return_seqs) = @_;
	
	my ($rs_index);
	
	open IN, "$infile" or die "No idea what file \"$infile\" is. $!\n";
	my $slurp = <IN>;
	close IN;
	
	open OUT, ">$infile";
	my @tmp = split(">", $slurp);
	$rs_index = 0;
	$return_seqs->[$rs_index++] = ">$tmp[1]";
	for my $i (2 .. @tmp) {
		my @tmp2 = split /\n/, $tmp[$i];
		my @tmp3 = split //, $tmp2[0];
		if ($tmp3[0] eq "T") { $return_seqs->[$rs_index++] = ">$tmp[$i]"; }
	}
}

<<<<<<< .mine
=======






















sub plot_other 
{
#	my ($input_file, $label, $file_root, $mat, $bl, $cat, $CpG_x, $num_reps) = @_;

#	open OUT, ">gnuplot.gnu";
#	print OUT "set term postscript eps color\n";
#	print OUT "set title \"$mat bl:$bl, Gamma:$cat CpG:$CpG_x num_reps:$num_reps\" font \"Helvetica,20\"\n";
#	print OUT "set size 1.25,1.0\n";
#	print OUT "set xlabel \"$statistics[$x_axis_stat-1]\"\n";
#	print OUT "set output \"$file_root.eps\"\n";
#	print OUT "plot ";
#	print OUT " \"$input_file\" usi $x_axis_stat:2 ti \"$label\"";
#	close OUT;

#	system("gnuplot gnuplot.gnu");

#	system("open $file_root.eps");
}

sub plot 
{
#	my ($y_axis_ref, $gnu_opt, $file_root, $mat, $bl, $cat, $CpG_x, $num_reps) = @_;

#	print Dumper(@{ $y_axis_ref });

#	open OUT, ">gnuplot.gnu";
#	print OUT "set term postscript eps color\n";
#	print OUT "set title \"$mat bl:$bl, Gamma:$cat CpG:$CpG_x num_reps:$num_reps\" font \"Helvetica,20\"\n";
#	print OUT "set size 1.25,1.0\n";
#	print OUT "$gnu_opt\n";
#	print OUT "set xlabel \"$statistics[$x_axis_stat-1]\"\n";
#	print OUT "set output \"$file_root.eps\"\n";
#	print OUT "plot ";
	
#	print " " . (scalar @{ $y_axis_ref }) . "\n";
#	for my $i (0 .. @{ $y_axis_ref }-1) {
#		print OUT " \"mu.dat\" usi $x_axis_stat:$y_axis_ref->[$i] ti \"$statistics[$y_axis_ref->[$i]-1]\"";
#		if ($i+1 != (scalar @{ $y_axis_ref })) { print OUT ", "; }
#		print " \"mu.dat\" usi $x_axis_stat:$y_axis_ref->[$i] ti \"$statistics[$y_axis_ref->[$i]-1]\"";
#	}
#	close OUT;

#	system("gnuplot gnuplot.gnu");

#	system("open $file_root.eps");

}

sub outGnuplotFile
{
#	my ($simsub_ptr, $prosub_ptr, $file_root, $model, $branch_length, $gamma) = @_;

#	my $min_x_val = 1000000;
#	my $max_x_val = 0;

#	for my $i (0 .. @{ $simsub_ptr }-1) {
#		if ( $simsub_ptr->[$i] < $min_x_val	) { $min_x_val=$simsub_ptr->[$i]; }
#		if ( $simsub_ptr->[$i] > $max_x_val	) { $max_x_val=$simsub_ptr->[$i]; }
#	}
#	if ($ENDPOINT_COND) {
#		for my $i (0 .. @{ $prosub_ptr }-1) {
#			if ( $prosub_ptr->[$i] < $min_x_val	) { $min_x_val=$prosub_ptr->[$i]; }
#			if ( $prosub_ptr->[$i] > $max_x_val	) { $max_x_val=$prosub_ptr->[$i]; }
#		}
#	}

#	print "($min_x_val,$max_x_val)\n";

#	my (@outsim, @outpro);
#	for my $i (0 .. @{ $simsub_ptr }-1) {
#		$outsim[$simsub_ptr->[$i]]++;
#	}
#	if ($ENDPOINT_COND) {
#		for my $i (0 .. @{ $prosub_ptr }-1) {
#			$outpro[$prosub_ptr->[$i]]++;
#		}
#	}
	
#	open OUT, ">$file_root.simgnu" or die "WTFGNU?\n";
#	for my $i ($min_x_val .. $max_x_val) {
#		if ($outsim[$i] ne "") {
#			print OUT "$i $outsim[$i]\n";
#		} else {
#			print OUT "$i 0\n";
#		}
#	}
#	close OUT;
#	if ($ENDPOINT_COND) {
#		open OUT, ">$file_root.prognu" or die "WTFGNU?\n";
#		for my $i ($min_x_val .. $max_x_val) {
#			if ($outpro[$i] ne "") {
#				print OUT "$i $outpro[$i]\n";
#			} else {
#				print OUT "$i 0\n";
#			}
#		}
#		close OUT;
#	}

#	my ($sim_mu, $sim_sigma, $epc_mu, $epc_sigma);
#	($sim_mu, $sim_sigma) = &mu_sigma($simsub_ptr);
#	if($ENDPOINT_COND) { ($epc_mu, $epc_sigma) = &mu_sigma($prosub_ptr); }

#	open OUT, ">gnuplot.gnu";
#	print OUT "set term postscript eps color\n";
#	print OUT "set title \"Distribution Endpoint and Simulated Substitutions\\n $model matrix with branch length $branch_length, Gamma $gamma\" font \"Helvetica,24\"\n";
#	print OUT "set style line 1 lt 1 lw 0 lc rgb \"red\"\n";
#	print OUT "set style line 1 lt 1 lw 0 lc rgb \"black\"\n";
#	print OUT "set boxwidth 0.5\n";
#	print OUT "set size 1.25,1.0\n";
#	print OUT "set xlabel \"Number of Substitutions\"\n";
#	print OUT "set ylabel \"Number of Replicate Datasets\"\n";
#	print OUT "set output \"$file_root.eps\"\n";
#	print OUT "plot \"$file_root.simgnu\" usi (\$1-0.25):2 w boxes fs solid 0.75 ls 1 lc rgb \"red\" ti \"Simulated " . sprintf("%.2f",$sim_mu) . " +/- " . sprintf("%.2f",$sim_sigma) . "\"";
#	if ($ENDPOINT_COND) { print OUT ", \\\n     \"$file_root.prognu\" usi (\$1+0.25):2 w boxes fs solid 0.75 ls 1 lc rgb \"black\" ti \"Proposed " . sprintf("%.2f",$epc_mu) . " +/- " . sprintf("%.2f",$epc_sigma) . "\"\n"; }

#	system("gnuplot gnuplot.gnu");

#	system("open $file_root.eps")
}

>>>>>>> .r1895
