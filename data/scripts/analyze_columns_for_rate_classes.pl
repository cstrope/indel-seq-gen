#!/usr/bin/env perl
$/=undef;
use strict;

my ($file_root, $treefile, $command_line_options);
my (@branch_length, $j, $slurp, @tmp, @seq1, @seq2, @seq3);
my $b2 = 0.5;
my $rootseq_length = 100000;

my $file_root = "JC69_0.5_1_4-0.5";
$treefile  = $file_root . ".tree";

open OUT, ">$treefile" or die "WTF??\n";
print OUT "[$rootseq_length]";
#print OUT "((T_1:$b2,T_2:$b2):$b2, T_3:0);\n";
#print OUT "((T_1:$b2,T_2:$b2):0, T_3:$b2);\n";
print OUT "((T_1:$b2,T_2:$b2):" . ($b2/2) . ",T_3:" . ($b2/2) . ");\n";
close OUT;

my $command_line_options = "";
$command_line_options .= "-m JC69 ";
$command_line_options .= "-s $rootseq_length ";
$command_line_options .= "-g 4 -a 0.5 ";
$command_line_options .= "-C 1 "; 
$command_line_options .= "-d 001010 ";
$command_line_options .= "-e $file_root ";
$command_line_options .= "-o f ";
$command_line_options .= "< $treefile > acfrc.junk";
print STDERR "./indel-seq-gen $command_line_options\n";
system("./indel-seq-gen $command_line_options");

open IN, "$file_root.ma";
$slurp = <IN>;
close IN;

@tmp = split /\n/, $slurp;

@seq1 = split //, $tmp[1];
@seq2 = split //, $tmp[3];
@seq3 = split //, $tmp[5];

my @count;

for (my $i = 0; $i < 64; $i++) {
	$count[$i] = 0;
}

my $total = 0;

for (my $i = 0; $i < $rootseq_length; $i++) {
	$count[&pos($seq1[$i], $seq2[$i], $seq3[$i])]++;
}

my @n = ("A", "C", "G", "T");

for (my $i = 0; $i < 64; $i++) {
	my $n1 = int($i / 16);
	my $n2 = int(($i - ($n1*16)) / 4);
	my $n3 = ($i - $n1*16 - $n2*4);
	print " $n[$n1]$n[$n2]$n[$n3] $count[$i]\n";
	$total += $count[$i];
}

print "Total: $total\n";

sub pos
{
	my ($a, $b, $c) = @_;
	my $position = 0;

	my ($ta, $tb, $tc);

	if ($a eq "A") { $ta = 0; }
	if ($a eq "C") { $ta = 1; }
	if ($a eq "G") { $ta = 2; }
	if ($a eq "T") { $ta = 3; }

	if ($b eq "A") { $tb = 0; }
	if ($b eq "C") { $tb = 1; }
	if ($b eq "G") { $tb = 2; }
	if ($b eq "T") { $tb = 3; }

	if ($c eq "A") { $tc = 0; }
	if ($c eq "C") { $tc = 1; }
	if ($c eq "G") { $tc = 2; }
	if ($c eq "T") { $tc = 3; }

	return $ta*16 + $tb*4 + $tc;

}
