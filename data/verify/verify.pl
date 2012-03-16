#!/usr/bin/env perl
use Data::Dumper;

$/=undef;
use strict;

##################
### This script will make 2 seq datasets, run path proposal, and gnuplot results.
##################

my $PLOT = 1;
my $TEST = 100;
my $FORDWARD_SIM = 1;
my $ENDPOINT_COND = 1;
my $OUTPUT_NUCL_STATS = 0;

my $MAX_BIN = 0;

my ($slurp, @tmp, @tmp1, @tmp2);

my $x_axis_stat = 1;
my @y_axis_stats = (
					[ 2, 3 ],
					[ 4, 5 ],
					[ 6, 7 ],
					[ 8, 9 ],
					[ 11, 12, 13, 14 ],
					[ 15 ],
					[ 16 ],
					[ 17 ],
					[ 10 ],
				   );

my @gnuplot_options = (
					"",
					"",
					"",
					"",
					"set log y",
					"",
					"",
					"",
					"",
				   );

# print Dumper(\@y_axis_stats);
# print Dumper(\@y_axis_stats[0]);

my $num_bins = 100;

my @statistics = (
				  "t",										# 1
				  "dQij sum Pjk/Pij ratios, same sites",	# 2
				  "dQij sum Pjk/Pij ratios, diff sites",	# 3
				  "dQij sum, same sites",					# 4
				  "dQij sum, diff sites",					# 5
				  "Qij sum Pjk/Pij ratios, same sites",		# 6
				  "Qij sum Pjk/Pij ratios, diff sites",		# 7
				  "Qij sum, same sites",					# 8
				  "Qij sum, diff sites",					# 9
				  "Number of different sites",				# 10
				  "dQij, dt = 0.01",						# 11
				  "dQij, dt = 0.001",						# 12
				  "dQij, dt = 0.0001",						# 13
				  "dQij, dt = 0.00001",						# 14
				  "Number of substitutions",				# 15
				  "dt",										# 16
				  "rate_away",								# 17
				  "",
				  "",
				  "",
				  "",
				  "",
				 );

my @matrix = (
			  "JC69",
#			  "JC69",
#			  "HKY",
#			  "F81",		### Successful: rev. 1816, 4-23-11. 10000 reps, 
						    ###   * NoCpG: 2999.0430
						    ###   * CpG:   
#			  "F81",		###
						    ###   * NoCpG: 2999.4942
						    ###   * CpG:   2943.xxxx
#			  "F81",
#			  "GTR",
#			  "SYM",
			 );
my @path_freqs =
		(
		 "0.25 0.25 0.25 0.25",
#		 "0.25 0.25 0.25 0.25",
#		 "0.3 0.2 0.2 0.3",
#		 "0.05 0.05 0.05 0.85",
#		 "0.1 0.4 0.4 0.1",
#		 "0.05 0.05 0.05 0.85",
#		 "0.1 0.2 0.3 0.4"
#		 "0.25 0.25 0.25 0.25",
		);
my @iSG_freqs =
		(
		 "0.25,0.25,0.25,0.25",
#		 "0.25,0.25,0.25,0.25",
#		 "0.3,0.2,0.2,0.3",
#		 "0.05,0.05,0.05,0.85",
#		 "0.1,0.4,0.4,0.1",
#		 "0.05,0.05,0.05,0.85",
#		 "0.1,0.2,0.3,0.4"
#		 "0.25,0.25,0.25,0.25",
		);
my @iSG_rel_rates = 
		(
		 "",
#		 "",
#		 "-r 2",
#		 "",
#		 "",
#		 "",
#		 "-r 2,3,2,2,3,1"
#		 "-r 2,3,2,2,3,1"
		);
my @path_rel_rates = 
		(
		 "",
#		 "",
#		 "2",
#		 "",
#		 "",
#		 "",
#		 "2 3 2 2 3 1"
#		 "2 3 2 2 3 1"
		);
my @iSG_gamma_cats =
		(
		 "",
#		 "-g 4 -a 0.5",
#		 "",
#		 "-g 4 -a 0.5",
#		 "-g 4 -a 0.5",
#		 "",
#		 "-g 4 -a 0.5",
#		 ""
		);
my @path_gamma_cats =
		(
		 "",
#		 "-g 4 -a 0.5",
#		 "",
#		 "4 0.5",
#		 "4 0.5",
#		 "",
#		 "4 0.5",
#		 ""
		);

my @junkfile_sim = 
		(
		 "test_0s",
#		 "test_01s",
#		 "test_1s",
#		 "yyys",
#		 "aaas",
#		 "zzzs",
#		 "zzzs",
#		 "xxxs",
#		 "test_4s",
		);
my @junkfile_pro = 
		(
		 "test_0p",
#		 "test_01p",
#		 "test_1p",
#		 "yyyp",
#		 "aaap",
#		 "zzzp",
#		 "zzzp",
#		 "xxxp",
#		 "test_4p",
		);
my $treefile;
my @branch_length = 
		(
#		 "0.1",
#		 "0.5",
#		 "1.0",
#		 "2.0",
		 "3.0"
		);

my @CpG_multiplier =
		(
#		 "1",
#		 "2",
#		 "5",
		 "10",
#		 "50",
		);

my $rootseq_length;
my ($root_seq, $target_seq);
my $pathfile;
my (@simulated_subst, @proposed_subst);
my $num_replicates;
if ($TEST) { $num_replicates = $TEST; $rootseq_length = 1000; } 
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

##########
### Set up an array to grab all output from iSG, used to sum up over runs and plot.
##########
for my $i (0 .. 20) {
	$output_data[$i] = [ @array ];
	for my $j (0 .. 10000) {
		$output_data[$i][$j] = 0;
		$bin_hits[$j] = 0;
		$each_subst[$j] = 0;
	}
}

for my $i (0 .. @junkfile_sim-1 ) {
	$simulated_subst_per_model_condition[$i] = 0;
}

for my $i (0 .. @matrix-1) {
	for my $j (0 .. @branch_length-1) {
		for my $c (0 .. @CpG_multiplier-1) {
			for my $k (0 .. $num_replicates-1) {
				$simulated_subst[$k] = 0;
				$proposed_subst[$k] = 0;
				$simulated_A[$k] = 0;$simulated_C[$k] = 0;$simulated_G[$k] = 0;$simulated_T[$k] = 0;$proposed_A[$k] = 0;$proposed_C[$k] = 0;$proposed_G[$k] = 0;$proposed_T[$k] = 0;$simulated_AA[$k] = 0;$simulated_AC[$k] = 0;$simulated_AG[$k] = 0;$simulated_AT[$k] = 0;$proposed_AA[$k] = 0;$proposed_AC[$k] = 0;$proposed_AG[$k] = 0;$proposed_AT[$k] = 0;$simulated_CA[$k] = 0;$simulated_CC[$k] = 0;$simulated_CG[$k] = 0;$simulated_CT[$k] = 0;$proposed_CA[$k] = 0;$proposed_CC[$k] = 0;$proposed_CG[$k] = 0;$proposed_CT[$k] = 0;$simulated_GA[$k] = 0;$simulated_GC[$k] = 0;$simulated_GG[$k] = 0;$simulated_GT[$k] = 0;$proposed_GA[$k] = 0;$proposed_GC[$k] = 0;$proposed_GG[$k] = 0;$proposed_GT[$k] = 0;$simulated_TA[$k] = 0;$simulated_TC[$k] = 0;$simulated_TG[$k] = 0;$simulated_TT[$k] = 0;$proposed_TA[$k] = 0;$proposed_TC[$k] = 0;$proposed_TG[$k] = 0;$proposed_TT[$k] = 0;
			}
			for my $i (0 .. 20) {
				for my $j (0 .. $branch_length[$j]*$num_bins) {
					$output_data[$i][$j] = 0;
					$bin_hits[$j] = 0;
					$each_subst[$j] = 0;
				}
			}
			my @filetmp = split /\s+/, $path_gamma_cats[$i];
			my $file_root = "$matrix[$i]" . "_$branch_length[$j]" . "_$CpG_multiplier[$c]" . "_$filetmp[0]-$filetmp[1]";
			$treefile  = $file_root . ".tree";
			$pathfile  = $file_root . ".path";
			$running_total_simsubs = $running_total_prosubs = 0;
			for my $k (0 .. $num_replicates-1) {
				sleep(0.1);
				print STDERR "********$matrix[$i], $branch_length[$j]******rep $k/$num_replicates*************\n";
	
				&outputTreefile($treefile, $branch_length[$j]);
				my $command_line_options = "";
				$command_line_options .= "-m $matrix[$i] ";
				$command_line_options .= "-s $rootseq_length ";
				$command_line_options .= "$iSG_gamma_cats[$i] ";
				$command_line_options .= "$iSG_rel_rates[$i] ";
				$command_line_options .= "-f $iSG_freqs[$i] ";
				$command_line_options .= "-C $CpG_multiplier[$c] "; 
				$command_line_options .= "-w a ";
				$command_line_options .= "-d 001010 ";
				$command_line_options .= "-e $file_root ";
				$command_line_options .= "-o f ";
				$command_line_options .= "< $treefile > $junkfile_sim[$i]";
				print "./indel-seq-gen $command_line_options\n";
				system("./indel-seq-gen $command_line_options");
				##########
				### Gather MSA sequences to output
				##########
				open IN, "$file_root.ma";
				$slurp = <IN>;
				close IN;
				@tmp = split(">", $slurp);	

				my (@seqs);
				for my $g (1 .. @tmp) {
					@tmp1 = split /\n/, $tmp[$g];
					$seqs[$g-1] = $tmp1[1];
				}

				
				$simulated_subst[$k] = 
				$simulated_subst_per_model_condition[$i] 
				= &getSubstitutions("$file_root.trace", "FORWARD");
				$running_total_simsubs += $simulated_subst[$k];
				print STDERR "Gathered simsubs: $simulated_subst[$k]\n";
				print STDERR "  Running total simsubs: " . ($running_total_simsubs/($k+1)) . "\n";

				open IN, "$junkfile_sim[$i]";
				$slurp = <IN>;
				close IN;
			
				@tmp = split("<<", $slurp);
				@tmp1 = split /\n/, $tmp[0];
				for my $i ( 0 .. @tmp1-1 ) {
					@tmp = split /\s+/, $tmp1[$i];
					$sim_CpG[$i] += $tmp[1];
				}

				if ($ENDPOINT_COND) {
					&outputPathFile($pathfile, \@seqs);
					my $command_line_options = "";
					$command_line_options .= "-m $matrix[$i] ";
					$command_line_options .= "$iSG_gamma_cats[$i] ";
					$command_line_options .= "$iSG_rel_rates[$i] ";
					$command_line_options .= "-f $iSG_freqs[$i] ";
					$command_line_options .= "-d 000010 ";
					$command_line_options .= "-C $CpG_multiplier[$c] "; 
					$command_line_options .= "-e path_proposal ";
					$command_line_options .= "-E $file_root.ma ";
					$command_line_options .= "< $treefile > $junkfile_pro[$i]";
					print "./indel-seq-gen $command_line_options\n";
					system("./indel-seq-gen $command_line_options");

					$proposed_subst[$k] = &getSubstitutions("path_proposal.trace", "EPC");
					$running_total_prosubs += $proposed_subst[$k];
					print STDERR "Gathered prosubs: $proposed_subst[$k]\n";
					print STDERR "  Running total prosubs: " . ($running_total_prosubs/($k+1)) . "\n";

					open IN, "$junkfile_pro[$i]";
					$slurp = <IN>;
					close IN;

					@tmp1 = split /\n/, $slurp;
					for my $j ( 0 .. @tmp1-7 ) {
						@tmp = split /\s+/, $tmp1[$j];
						my $bin_no = int($tmp[0]*$num_bins);
						$bin_hits[$bin_no]++;
						if ($bin_no > $MAX_BIN) { $MAX_BIN = $bin_no; }
						$each_subst[$bin_no]++;
						for my $i (1 .. @tmp) {
							$output_data[$i][$bin_no] += $tmp[$i];
						}
					}
				}
			}
			if (!$TEST) { &outGnuplotFile(\@simulated_subst, \@proposed_subst, $file_root, $matrix[$i], $branch_length[$j], $path_gamma_cats[$i]); }
			if ($PLOT) {
				open OUT, ">mu.dat";
				my $printed = 0;
				for my $j (0 .. $num_bins*$branch_length[$j]) {
					my $out = ($j/$num_bins) . " ";
					$printed = 0;
					for my $i (1 .. 19) {
						if ($bin_hits[$j]) {
							$printed = 1;
							$out .= $output_data[$i][$j] . " ";# / $bin_hits[$j] . " ";
						} else { $out .= "0 "; }
					}
					if ($printed) { print OUT "$out\n"; }
				}
				close OUT;


				for my $y (0 .. @y_axis_stats-1) {
					print Dumper(\@y_axis_stats[$y]);
					&plot(
						  @y_axis_stats[$y],
						  $gnuplot_options[$y],
						  $file_root . "_$y",
						  $matrix[$i], 
						  $branch_length[$j], 
						  $path_gamma_cats[$i],
						  $CpG_multiplier[$c],
						  $num_replicates
						 );
				}

			}

			open OUT, ">binhits.dat";
			for my $j (0 .. $num_bins*$branch_length[$j]) {
				print OUT ($j/$num_bins) . " ";
				if ($bin_hits[$j]) {
					print OUT $bin_hits[$j] . " ";
				} else { print OUT "0 "; }
				print OUT "\n";
			}
			close OUT;

			&plot_other(
				  "binhits.dat",
				  "data points in bin",
				  $file_root . "_binhits",
				  $matrix[$i], 
				  $branch_length[$j], 
				  $path_gamma_cats[$i],
				  $CpG_multiplier[$c],
				  $num_replicates
				 );

			open OUT, ">subst.dat";

			for my $j (0 .. $num_bins*$branch_length[$j]) {
				print OUT ($j/$num_bins) . " ";
				if ($bin_hits[$j]) {
					print OUT $bin_hits[$j] . " ";
				} else { print OUT "0 "; }
				print OUT "\n";
			}
			close OUT;
			&plot_other(
				  "subst.dat",
				  "number of substitutions in bin",
				  $file_root . "_subst",
				  $matrix[$i], 
				  $branch_length[$j], 
				  $path_gamma_cats[$i],
				  $CpG_multiplier[$c],
				  $num_replicates
				 );
		}
	}
}

sub plot_other {
	my ($input_file, $label, $file_root, $mat, $bl, $cat, $CpG_x, $num_reps) = @_;

	open OUT, ">gnuplot.gnu";
	print OUT "set term postscript eps color\n";
	print OUT "set title \"$mat bl:$bl, Gamma:$cat CpG:$CpG_x num_reps:$num_reps\" font \"Helvetica,20\"\n";
	print OUT "set size 1.25,1.0\n";
	print OUT "set xlabel \"$statistics[$x_axis_stat-1]\"\n";
	print OUT "set output \"$file_root.eps\"\n";
	print OUT "plot ";
	print OUT " \"$input_file\" usi $x_axis_stat:2 ti \"$label\"";
	close OUT;

	system("gnuplot gnuplot.gnu");

	system("open $file_root.eps");
}

sub plot {
	my ($y_axis_ref, $gnu_opt, $file_root, $mat, $bl, $cat, $CpG_x, $num_reps) = @_;

	print Dumper(@{ $y_axis_ref });

	open OUT, ">gnuplot.gnu";
	print OUT "set term postscript eps color\n";
	print OUT "set title \"$mat bl:$bl, Gamma:$cat CpG:$CpG_x num_reps:$num_reps\" font \"Helvetica,20\"\n";
	print OUT "set size 1.25,1.0\n";
	print OUT "$gnu_opt\n";
	print OUT "set xlabel \"$statistics[$x_axis_stat-1]\"\n";
	print OUT "set output \"$file_root.eps\"\n";
	print OUT "plot ";
	
	print " " . (scalar @{ $y_axis_ref }) . "\n";
	for my $i (0 .. @{ $y_axis_ref }-1) {
		print OUT " \"mu.dat\" usi $x_axis_stat:$y_axis_ref->[$i] ti \"$statistics[$y_axis_ref->[$i]-1]\"";
		if ($i+1 != (scalar @{ $y_axis_ref })) { print OUT ", "; }
		print " \"mu.dat\" usi $x_axis_stat:$y_axis_ref->[$i] ti \"$statistics[$y_axis_ref->[$i]-1]\"";
	}
	close OUT;

	system("gnuplot gnuplot.gnu");

	system("open $file_root.eps");

}

for my $i (0 .. @junkfile_sim-1 ) {
	print STDERR "********$matrix[$i]****$path_freqs[$i]****$path_gamma_cats[$i]********\n";
	print "$simulated_subst_per_model_condition[$i] / $num_replicates)\n";
	print "" . ($simulated_subst_per_model_condition[$i] / $num_replicates) . "\n";
}

open OUT, ">forward_CpG.gnu";
for my $i ( 0 .. 59999 ) {
	print OUT "$i " . ($sim_CpG[$i] / $num_replicates) . "\n";
}
close OUT;

open OUT, ">epc_CpG.gnu";
for my $i ( 0 .. 59999 ) {
	print OUT "$i " . ($prop_CpG[$i] / $num_replicates) . "\n";
}
close OUT;

sub outGnuplotFile
{
	my ($simsub_ptr, $prosub_ptr, $file_root, $model, $branch_length, $gamma) = @_;

	my $min_x_val = 1000000;
	my $max_x_val = 0;

	for my $i (0 .. @{ $simsub_ptr }-1) {
		if ( $simsub_ptr->[$i] < $min_x_val	) { $min_x_val=$simsub_ptr->[$i]; }
		if ( $simsub_ptr->[$i] > $max_x_val	) { $max_x_val=$simsub_ptr->[$i]; }
	}
	if ($ENDPOINT_COND) {
		for my $i (0 .. @{ $prosub_ptr }-1) {
			if ( $prosub_ptr->[$i] < $min_x_val	) { $min_x_val=$prosub_ptr->[$i]; }
			if ( $prosub_ptr->[$i] > $max_x_val	) { $max_x_val=$prosub_ptr->[$i]; }
		}
	}

	print "($min_x_val,$max_x_val)\n";

	my (@outsim, @outpro);
	for my $i (0 .. @{ $simsub_ptr }-1) {
		$outsim[$simsub_ptr->[$i]]++;
	}
	if ($ENDPOINT_COND) {
		for my $i (0 .. @{ $prosub_ptr }-1) {
			$outpro[$prosub_ptr->[$i]]++;
		}
	}
	
	open OUT, ">$file_root.simgnu" or die "WTFGNU?\n";
	for my $i ($min_x_val .. $max_x_val) {
		if ($outsim[$i] ne "") {
			print OUT "$i $outsim[$i]\n";
		} else {
			print OUT "$i 0\n";
		}
	}
	close OUT;
	if ($ENDPOINT_COND) {
		open OUT, ">$file_root.prognu" or die "WTFGNU?\n";
		for my $i ($min_x_val .. $max_x_val) {
			if ($outpro[$i] ne "") {
				print OUT "$i $outpro[$i]\n";
			} else {
				print OUT "$i 0\n";
			}
		}
		close OUT;
	}

	my ($sim_mu, $sim_sigma, $epc_mu, $epc_sigma);
	($sim_mu, $sim_sigma) = &mu_sigma($simsub_ptr);
	if($ENDPOINT_COND) { ($epc_mu, $epc_sigma) = &mu_sigma($prosub_ptr); }

	open OUT, ">gnuplot.gnu";
	print OUT "set term postscript eps color\n";
	print OUT "set title \"Distribution Endpoint and Simulated Substitutions\\n $model matrix with branch length $branch_length, Gamma $gamma\" font \"Helvetica,24\"\n";
	print OUT "set style line 1 lt 1 lw 0 lc rgb \"red\"\n";
	print OUT "set style line 1 lt 1 lw 0 lc rgb \"black\"\n";
	print OUT "set boxwidth 0.5\n";
	print OUT "set size 1.25,1.0\n";
	print OUT "set xlabel \"Number of Substitutions\"\n";
	print OUT "set ylabel \"Number of Replicate Datasets\"\n";
	print OUT "set output \"$file_root.eps\"\n";
	print OUT "plot \"$file_root.simgnu\" usi (\$1-0.25):2 w boxes fs solid 0.75 ls 1 lc rgb \"red\" ti \"Simulated " . sprintf("%.2f",$sim_mu) . " +/- " . sprintf("%.2f",$sim_sigma) . "\"";
	if ($ENDPOINT_COND) { print OUT ", \\\n     \"$file_root.prognu\" usi (\$1+0.25):2 w boxes fs solid 0.75 ls 1 lc rgb \"black\" ti \"Proposed " . sprintf("%.2f",$epc_mu) . " +/- " . sprintf("%.2f",$epc_sigma) . "\"\n"; }

	system("gnuplot gnuplot.gnu");

	system("open $file_root.eps")
}

sub mu_sigma
{
	my ($ref) = @_;
	
	my ($average, $stdev);
	
	for my $i (0 .. @{ $ref }-1) {
		$average += $ref->[$i];
	}

	$average /= scalar(@{ $ref });

	for my $i (0 .. @{ $ref }-1) {
		$stdev += ($average - $ref->[$i]) * ($average - $ref->[$i]);
	}
	
	$stdev /= scalar(@{ $ref });

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
}

sub outputTreefile
{
	my ($treefile, $branch_length) = @_;

	open OUT, ">$treefile" or die "WTF??\n";
	print OUT "[$rootseq_length](t1:$branch_length,t2:0);\n";
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

