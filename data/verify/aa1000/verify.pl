#!/usr/bin/env perl
$/=undef;
use strict;

##################
### This script will make 2 seq datasets, run path proposal, and gnuplot results.
##################

my $GIL = 1;

my ($slurp, @tmp, @tmp1, @tmp2);
my @matrix = (
#			  "",
#			  ",
			  "WAG",	# Running, terminal 4.
#			  "GTR",	# Running, terminal 2.
#			  "",
			 );
my @path_freqs =
		(
#		 "0.25 0.25 0.25 0.25",
#		 "0.3 0.2 0.2 0.3",
		 "0.0866 0.0440 0.0391 0.0570 0.0193 0.0367 0.0581 0.0833 0.0244 0.0485 0.0862 0.0620 0.0195 0.0384 0.0458 0.0695 0.0610 0.0144 0.0353 0.0709",
#		 "0.1 0.2 0.3 0.4"
#		 "0.25 0.25 0.25 0.25",
		);
my @iSG_freqs =
		(
#		 "0.25,0.25,0.25,0.25",
#		 "0.3,0.2,0.2,0.3",
		 "0.0866,0.0440,0.0391,0.0570,0.0193,0.0367,0.0581,0.0833,0.0244,0.0485,0.0862,0.0620,0.0195,0.0384,0.0458,0.0695,0.0610,0.0144,0.0353,0.0709",
#		 "0.1,0.2,0.3,0.4"
#		 "0.25,0.25,0.25,0.25",
		);
my @iSG_rel_rates = 
		(
#		 "",
#		 "-r 2",
		 "",
#		 "-r 2,3,2,2,3,1"
#		 "-r 2,3,2,2,3,1"
		);
my @path_rel_rates = 
		(
#		 "",
#		 "2",
		 "",
#		 "2 3 2 2 3 1"
#		 "2 3 2 2 3 1"
		);
my @iSG_gamma_cats =
		(
#		 "",
#		 "",
		 "",
#		 "-g 4 -a 0.5",
#		 ""
		);
my @path_gamma_cats =
		(
#		 "",
#		 "",
		 "",
#		 "4 0.5",
#		 ""
		);

my @junkfile = 
		(
#		 "test_0",
#		 "test_1",
		 "yyy",
#		 "xxx",
#		 "test_4",
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
my $rootseq_length = 1000;
my ($root_seq, $target_seq);
my $pathfile;
my (@simulated_subst, @proposed_subst);
my $num_replicates = 1000;
my @simulated_A;
my @simulated_C;
my @simulated_G;
my @simulated_T;
my @proposed_A;my @proposed_C;my @proposed_G;my @proposed_T;
my @proposed_AA;my @proposed_AC;my @proposed_AG;my @proposed_AT;my @proposed_CA;my @proposed_CC;my @proposed_CG;my @proposed_CT;my @proposed_GA;my @proposed_GC;my @proposed_GG;my @proposed_GT;my @proposed_TA;my @proposed_TC;my @proposed_TG;my @proposed_TT;my @simulated_AA;my @simulated_AC;my @simulated_AG;my @simulated_AT;my @simulated_CA;my @simulated_CC;my @simulated_CG;my @simulated_CT;my @simulated_GA;my @simulated_GC;my @simulated_GG;my @simulated_GT;my @simulated_TA;my @simulated_TC;my @simulated_TG;my @simulated_TT;


for my $i (0 .. @matrix-1) {
	for my $j (0 .. @branch_length-1) {
		for my $k (0 .. $num_replicates) {
			$simulated_subst[$k] = 0;
			$proposed_subst[$k] = 0;
			$simulated_A[$k] = 0;$simulated_C[$k] = 0;$simulated_G[$k] = 0;$simulated_T[$k] = 0;$proposed_A[$k] = 0;$proposed_C[$k] = 0;$proposed_G[$k] = 0;$proposed_T[$k] = 0;$simulated_AA[$k] = 0;$simulated_AC[$k] = 0;$simulated_AG[$k] = 0;$simulated_AT[$k] = 0;$proposed_AA[$k] = 0;$proposed_AC[$k] = 0;$proposed_AG[$k] = 0;$proposed_AT[$k] = 0;$simulated_CA[$k] = 0;$simulated_CC[$k] = 0;$simulated_CG[$k] = 0;$simulated_CT[$k] = 0;$proposed_CA[$k] = 0;$proposed_CC[$k] = 0;$proposed_CG[$k] = 0;$proposed_CT[$k] = 0;$simulated_GA[$k] = 0;$simulated_GC[$k] = 0;$simulated_GG[$k] = 0;$simulated_GT[$k] = 0;$proposed_GA[$k] = 0;$proposed_GC[$k] = 0;$proposed_GG[$k] = 0;$proposed_GT[$k] = 0;$simulated_TA[$k] = 0;$simulated_TC[$k] = 0;$simulated_TG[$k] = 0;$simulated_TT[$k] = 0;$proposed_TA[$k] = 0;$proposed_TC[$k] = 0;$proposed_TG[$k] = 0;$proposed_TT[$k] = 0;
		}
		my @filetmp = split /\s+/, $path_gamma_cats[$i];
		my $file_root = "$matrix[$i]" . "_$branch_length[$j]" . "_$filetmp[0]-$filetmp[1]";
		$treefile  = $file_root . ".tree";
		$pathfile  = $file_root . ".path";
		for my $k (0 .. $num_replicates) {
			print STDERR "********$matrix[$i], $branch_length[$j]******rep $k/$num_replicates*************\n";
			&outputTreefile($treefile, $branch_length[$j]);
			my $command_line_options = "";
			$command_line_options .= "-m $matrix[$i] ";
			$command_line_options .= "-s $rootseq_length ";
			$command_line_options .= "$iSG_gamma_cats[$i] ";
			$command_line_options .= "$iSG_rel_rates[$i] ";
			$command_line_options .= "-f $iSG_freqs[$i] ";
			if ($GIL) { $command_line_options .= "-j gil "; }
			else { $command_line_options .= "-j trs "; }
			$command_line_options .= "-d 001010 ";
			$command_line_options .= "-e $file_root ";
			$command_line_options .= "-o f ";
			$command_line_options .= "< $treefile > $junkfile[$i]";
			print "./indel-seq-gen $command_line_options\n";
			system("./indel-seq-gen $command_line_options");
			##########
			### Gather MSA sequences to output
			##########
			open IN, "$file_root.ma";
			$slurp = <IN>;
			close IN;
			@tmp = split(">", $slurp);
			@tmp1 = split /\n/, $tmp[1];
			$target_seq = $tmp1[1];
			@tmp1 = split /\n/, $tmp[2];
			$root_seq = $tmp1[1];

			open IN, "$file_root.trace";
			$slurp = <IN>;
			close IN;
			@tmp = split /\n\n/, $slurp;
			@tmp1 = split /\n/, $tmp[1];
			for my $a (0 .. @tmp1-3) {
				$simulated_subst[$k]++;
			}

			open IN, "$junkfile[$i]";
			$slurp = <IN>;
			close IN;
			
			@tmp = split("<<", $slurp);
			@tmp1 = split /\n/, $tmp[1];
			$simulated_A[$k] = $tmp1[0];
			$simulated_C[$k] = $tmp1[1];
			$simulated_G[$k] = $tmp1[2];
			$simulated_T[$k] = $tmp1[3];

			@tmp1 = split /\n/, $tmp[2];
			$simulated_AA[$k] = $tmp1[0];
			$simulated_AC[$k] = $tmp1[1];
			$simulated_AG[$k] = $tmp1[2];
			$simulated_AT[$k] = $tmp1[3];
			$simulated_CA[$k] = $tmp1[4];
			$simulated_CC[$k] = $tmp1[5];
			$simulated_CG[$k] = $tmp1[6];
			$simulated_CT[$k] = $tmp1[7];
			$simulated_GA[$k] = $tmp1[8];
			$simulated_GC[$k] = $tmp1[9];
			$simulated_GG[$k] = $tmp1[10];
			$simulated_GT[$k] = $tmp1[11];
			$simulated_TA[$k] = $tmp1[12];
			$simulated_TC[$k] = $tmp1[13];
			$simulated_TG[$k] = $tmp1[14];
			$simulated_TT[$k] = $tmp1[15];
	
			&outputPathFile($pathfile, $root_seq, $target_seq, $matrix[$i], $branch_length[$j], $path_rel_rates[$i], $path_freqs[$i], $path_gamma_cats[$i]);
			system("./indel-seq-gen -E $pathfile > $junkfile[$i]");
			open IN, "path_proposal.trace" or die "Did not generate path proposal file.\n";
			$slurp = <IN>;
			close IN;
			@tmp = split /\n/, $slurp;
			for my $a (1 .. @tmp-1) { 
				$proposed_subst[$k]++;
			}

			open IN, "$junkfile[$i]";
			$slurp = <IN>;
			close IN;
			
			@tmp = split("<<", $slurp);
			@tmp1 = split /\n/, $tmp[1];
			$proposed_A[$k] = $tmp1[0];
			$proposed_C[$k] = $tmp1[1];
			$proposed_G[$k] = $tmp1[2];
			$proposed_T[$k] = $tmp1[3];

			@tmp1 = split /\n/, $tmp[2];
			$proposed_AA[$k] = $tmp1[0];
			$proposed_AC[$k] = $tmp1[1];
			$proposed_AG[$k] = $tmp1[2];
			$proposed_AT[$k] = $tmp1[3];
			$proposed_CA[$k] = $tmp1[4];
			$proposed_CC[$k] = $tmp1[5];
			$proposed_CG[$k] = $tmp1[6];
			$proposed_CT[$k] = $tmp1[7];
			$proposed_GA[$k] = $tmp1[8];
			$proposed_GC[$k] = $tmp1[9];
			$proposed_GG[$k] = $tmp1[10];
			$proposed_GT[$k] = $tmp1[11];
			$proposed_TA[$k] = $tmp1[12];
			$proposed_TC[$k] = $tmp1[13];
			$proposed_TG[$k] = $tmp1[14];
			$proposed_TT[$k] = $tmp1[15];

#			print "$simulated_A[$k] $proposed_A[$k]\n"; 
#			print "$simulated_C[$k] $proposed_C[$k]\n"; 
#			print "$simulated_G[$k] $proposed_G[$k]\n"; 
#			print "$simulated_T[$k] $proposed_T[$k]\n"; 
#
#			print "$simulated_AA[$k] $proposed_AA[$k]\n"; 
#			print "$simulated_AC[$k] $proposed_AC[$k]\n"; 
#			print "$simulated_AG[$k] $proposed_AG[$k]\n"; 
#			print "$simulated_AT[$k] $proposed_AT[$k]\n"; 
#			print "$simulated_CA[$k] $proposed_CA[$k]\n"; 
#			print "$simulated_CC[$k] $proposed_CC[$k]\n"; 
#			print "$simulated_CG[$k] $proposed_CG[$k]\n"; 
#			print "$simulated_CT[$k] $proposed_CT[$k]\n"; 
#			print "$simulated_GA[$k] $proposed_GA[$k]\n"; 
#			print "$simulated_GC[$k] $proposed_GC[$k]\n"; 
#			print "$simulated_GG[$k] $proposed_GG[$k]\n"; 
#			print "$simulated_GT[$k] $proposed_GT[$k]\n"; 
#			print "$simulated_TA[$k] $proposed_TA[$k]\n"; 
#			print "$simulated_TC[$k] $proposed_TC[$k]\n"; 
#			print "$simulated_TG[$k] $proposed_TG[$k]\n"; 
#			print "$simulated_TT[$k] $proposed_TT[$k]\n"; 
#			print STDERR "$simulated_subst[$k] $proposed_subst[$k] \n";
		}

#		for my $k (0 .. $num_replicates) {
#			print "REPLICATE $k:\n";
#			print "$simulated_A[$k] $proposed_A[$k]\n"; 
#			print "$simulated_C[$k] $proposed_C[$k]\n"; 
#			print "$simulated_G[$k] $proposed_G[$k]\n"; 
#			print "$simulated_T[$k] $proposed_T[$k]\n"; 
#			print "\n";

#			print "$simulated_AA[$k] $proposed_AA[$k]\n"; 
#			print "$simulated_AC[$k] $proposed_AC[$k]\n"; 
#			print "$simulated_AG[$k] $proposed_AG[$k]\n"; 
#			print "$simulated_AT[$k] $proposed_AT[$k]\n"; 
#			print "$simulated_CA[$k] $proposed_CA[$k]\n"; 
#			print "$simulated_CC[$k] $proposed_CC[$k]\n"; 
#			print "$simulated_CG[$k] $proposed_CG[$k]\n"; 
#			print "$simulated_CT[$k] $proposed_CT[$k]\n"; 
#			print "$simulated_GA[$k] $proposed_GA[$k]\n"; 
#			print "$simulated_GC[$k] $proposed_GC[$k]\n"; 
#			print "$simulated_GG[$k] $proposed_GG[$k]\n"; 
#			print "$simulated_GT[$k] $proposed_GT[$k]\n"; 
#			print "$simulated_TA[$k] $proposed_TA[$k]\n"; 
#			print "$simulated_TC[$k] $proposed_TC[$k]\n"; 
#			print "$simulated_TG[$k] $proposed_TG[$k]\n"; 
#			print "$simulated_TT[$k] $proposed_TT[$k]\n"; 
#
#			print "\n";
#		}

		my ($simulated_mu, $proposed_mu, $simulated_sigma, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_A);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_A );
		printf("A: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_C);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_C );
		printf("C: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);

		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_G);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_G );
		printf("G: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);

		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_T);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_T );
		printf("T: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		print "\n";

		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_AA);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_AA );
		printf("AA: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_AC);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_AC );
		printf("AC: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_AG);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_AG );
		printf("AG: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_AT);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_AT );
		printf("AT: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		print "\n";

		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_CA);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_CA );
		printf("CA: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_CC);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_CC );
		printf("CC: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_CG);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_CG );
		printf("CG: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_CT);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_CT );
		printf("CT: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		print "\n";

		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_GA);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_GA );
		printf("GA: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_GC);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_GC );
		printf("GC: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_GG);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_GG );
		printf("GG: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_GT);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_GT );
		printf("GT: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		print "\n";

		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_TA);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_TA );
		printf("TA: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_TC);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_TC );
		printf("TC: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_TG);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_TG );
		printf("TG: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		($simulated_mu, $simulated_sigma) = &mu_sigma(\@simulated_TT);
		($proposed_mu,  $proposed_sigma) =  &mu_sigma(\@proposed_TT );
		printf("TT: %10.4f %8.4f    %10.4f %8.4f\n", $simulated_mu, $simulated_sigma, $proposed_mu, $proposed_sigma);
		print "\n";

		&outGnuplotFile(\@simulated_subst, \@proposed_subst, $file_root, $matrix[$i], $branch_length[$j], $path_gamma_cats[$i]);
	}
}


sub outGnuplotFile
{
	my ($simsub_ptr, $prosub_ptr, $file_root, $model, $branch_length, $gamma) = @_;

	my $min_x_val = 1000000;
	my $max_x_val = 0;

	for my $i (0 .. @{ $simsub_ptr }-1) {
		if ( $simsub_ptr->[$i] < $min_x_val	) { $min_x_val=$simsub_ptr->[$i]; }
		if ( $simsub_ptr->[$i] > $max_x_val	) { $max_x_val=$simsub_ptr->[$i]; }
	}
	for my $i (0 .. @{ $prosub_ptr }-1) {
		if ( $prosub_ptr->[$i] < $min_x_val	) { $min_x_val=$prosub_ptr->[$i]; }
		if ( $prosub_ptr->[$i] > $max_x_val	) { $max_x_val=$prosub_ptr->[$i]; }
	}

	print "($min_x_val,$max_x_val)\n";

	my (@outsim, @outpro);
	for my $i (0 .. @{ $simsub_ptr }-1) {
		$outsim[$simsub_ptr->[$i]]++;
	}
	for my $i (0 .. @{ $prosub_ptr }-1) {
		$outpro[$prosub_ptr->[$i]]++;
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
	open OUT, ">$file_root.prognu" or die "WTFGNU?\n";
	for my $i ($min_x_val .. $max_x_val) {
		if ($outpro[$i] ne "") {
			print OUT "$i $outpro[$i]\n";
		} else {
			print OUT "$i 0\n";
		}
	}
	close OUT;

	my ($sim_mu, $sim_sigma, $epc_mu, $epc_sigma);
	($sim_mu, $sim_sigma) = &mu_sigma($simsub_ptr);
	($epc_mu, $epc_sigma) = &mu_sigma($prosub_ptr);

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
	print OUT "plot \"$file_root.simgnu\" usi (\$1-0.25):2 w boxes fs solid 0.75 ls 1 lc rgb \"red\" ti \"Simulated " . sprintf("%.2f",$sim_mu) . " +/- " . sprintf("%.2f",$sim_sigma) . "\", \\\n";
	print OUT "     \"$file_root.prognu\" usi (\$1+0.25):2 w boxes fs solid 0.75 ls 1 lc rgb \"black\" ti \"Proposed " . sprintf("%.2f",$epc_mu) . " +/- " . sprintf("%.2f",$epc_sigma) . "\"\n";

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
	my ($pathfile, $root, $target, $Q, $BL, $rel_rates, $path_freqs, $gammastuff) = @_;
	
	open OUT, ">$pathfile" or die "WTF2??\n";
	print OUT "$Q\n";
	print OUT "$rel_rates\n";
	print OUT "$gammastuff\n";
	print OUT "$path_freqs\n";
	print OUT "$BL\n";
	print OUT "$root\n";
	print OUT "$target\n";
	close OUT;
}

sub outputTreefile
{
	my ($treefile, $branch_length) = @_;

	open OUT, ">$treefile" or die "WTF??\n";
	print OUT "[$rootseq_length](t1:$branch_length,t2:0);\n";
	close OUT;
}
