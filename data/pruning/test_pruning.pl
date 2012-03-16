#!/usr/bin/env perl;
$/=undef;

#print "TEST 1: No ancestral sequence.\n";
#print "./indel-seq-gen -m K80 -r 2 -w a -o f -e out -d 001010 < ZY_Comp-Mol-Evol.tree > ZY_CME_sim.junk\n";
#system("./indel-seq-gen -m K80  -r 2 -w a -o f -e out -d 001010 < ZY_Comp-Mol-Evol.tree > ZY_CME_sim.junk");
#open IN, "out.ma";
#$slurp = <IN>;
#close IN;
#open OUT, ">out.ma";
#@tmp = split(">", $slurp);
#print OUT ">$tmp[1]";
#for my $i (2 .. @tmp) {
#	@tmp2 = split /\n/, $tmp[$i];
#	@tmp3 = split //, $tmp2[0];
#	if ($tmp3[0] eq "T") { print OUT ">$tmp[$i]"; }
#}
#print "./indel-seq-gen -m K80  -r 2 -d 000010 -E out.ma -e path_proposal < ZY_comp-mol-Evol.tree > ZY_CME_epc.junk\n";
#system("./indel-seq-gen -m K80  -r 2 -d 000010 -E out.ma -e path_proposal < ZY_comp-mol-Evol.tree > ZY_CME_epc.junk");
print "./indel-seq-gen -m K80 -g 4 -a 0.5 -r 2 -d 000010 -E out.ma -e path_proposal < ZY_comp-mol-Evol.tree > ZY_CME_epc.junk\n";
system("./indel-seq-gen -m K80  -g 4 -a 0.5  -r 2 -d 000010 -E out.ma -e path_proposal < ZY_comp-mol-Evol.tree > ZY_CME_epc.junk");

