#!/usr/bin/env perl
use strict;

if (@ARGV != 1) { die "Usage: ./assign_markov_order_rates.pl <order>.\n"; }

my $SIMPLE_MODEL = 1;
my $SIMPLEST_MODEL = 1;
my @simple_model_transition_probs;

my $order = shift;

print "Order of Markov chain: $order.\n";

my @nucl = ( "A", "C", "G", "T");

srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);

my $order_rows = 4**$order-1;
my $random_number;
my (@xtuplet, @transition);
my $xtuplet_total = 0;
for my $i (0 .. $order_rows) {
	$xtuplet[$i] = rand();
	$xtuplet_total += $xtuplet[$i];
}

my $xtuplet_sum = 0;
for my $i (0 .. $order_rows - 1) {
	$xtuplet[$i] /= $xtuplet_total;
	$xtuplet_sum += $xtuplet[$i];
}

my $total = 0;

$xtuplet[$order_rows] = 1.0 - $xtuplet_sum;
$xtuplet_sum += $xtuplet[$order_rows];

for my $i (0 .. $order_rows) {
	my $rand_sum = 0;

	print "\t//" . &triplet($i) . "\n";
	for my $j (0 .. 3) {
		$transition[$j] = rand();
		$rand_sum += $transition[$j];
	}

	my $trans_sum = 0;
	for my $j (0 .. 2) {
		$transition[$j] /= $rand_sum;
		$trans_sum += $transition[$j];
	}

	$transition[3] = 1.0 - $trans_sum;
	$trans_sum += $transition[3];

	if ($i == 0) {
		$simple_model_transition_probs[0] = $transition[0];
		$simple_model_transition_probs[1] = $transition[1];
		$simple_model_transition_probs[2] = $transition[2];
		$simple_model_transition_probs[3] = $transition[3];
		if ($SIMPLEST_MODEL) {
			$simple_model_transition_probs[0] = 0.251;
			$simple_model_transition_probs[1] = 0.252;
			$simple_model_transition_probs[2] = 0.249;
			$simple_model_transition_probs[3] = 0.248;
		}
	}

	if ($SIMPLE_MODEL or $SIMPLEST_MODEL) {
		my ($first, $second, $third) = &trip($i);
		$xtuplet[$i] = $simple_model_transition_probs[$first] * $simple_model_transition_probs[$second] * $simple_model_transition_probs[$third];
		$total += $xtuplet[$i];
	}
	print "\ttriplet_pi.push_back($xtuplet[$i]);\n";
	print "\ttriplet_pi_inv.push_back(" . (1.0/$xtuplet[$i]) . ");\n";
	for my $j (0 .. 3) {
		if ($SIMPLE_MODEL) { print "\tpi_tmp.push_back($simple_model_transition_probs[$j]);\n"; }
		else { print "\tpi_tmp.push_back($transition[$j]);\n"; }
	}

	print "\tpi.push_back(pi_tmp);\n";
	print "\tfor (vector<double>::iterator it = pi_tmp.begin(); it != pi_tmp.end(); ++it) (*it) = 1.0/(*it);\n";
	print "\tpi_inv.push_back(pi_tmp); pi_tmp.clear();\n";

}

#print "Total triplet freq: $total\n";

sub trip
{
	my ($ival) = @_;

	my $first = int($ival / 16);
	$ival -= $first * 16;
	my $second = int($ival / 4);
	$ival -= $second * 4;
	my $third = $ival;
	
	return ($first, $second, $third);
}

sub triplet 
{
	my ($ival) = @_;

	my $first = int($ival / 16);
	$ival -= $first * 16;
	my $second = int($ival / 4);
	$ival -= $second * 4;
	my $third = $ival;
	
	return "$nucl[$first]$nucl[$second]$nucl[$third]";
}