#!/usr/bin/env perl
$/ = undef;


open IN, "calc_ur";
$slurp = <IN>;
close IN;

@each_val = split /\n/, $slurp;
$num_reps = scalar @each_val;
$sum = 0;

for my $i (0 .. @each_val-1) {
	@tmp = split /\s+/, $each_val[$i];
	$sum += $tmp[1];
}

$u = $sum/$num_reps;
$var = 0;

for my $i (0 .. @each_val-1) {
	@tmp = split /\s+/, $each_val[$i];
	$var += ($u - $tmp[1])**2;
}

$var /= $num_reps;


print "$u $var\n";
