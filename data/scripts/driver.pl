#!/usr/bin/env perl

for ($i = 0; $i < 3; $i++) {
	if ($i == 0) {
		$model = "HKY";
		$treefile = "ORFclade_labelled.tree";
		$freq_opt = "3,4,5,2";
		$lineage_file = "pseudogene_nuclineage.spec";
	}
	if ($i == 2) {
		$model = "JTT";
		$treefile = "clade_labelled.tree";
		$freq_opt = "3,4,5,1,4,2,5,6,7,2,3,4,7,4,5,3,2,4,6,7";
#		$lineage_file = "lineage.spec";
		$lineage_file = "pseudogene_lineage.spec";
	}
	if ($i == 1) {
		$model = "HKY";
		$treefile = "nucclade_labelled.tree";
		$freq_opt = "3,4,5,2";
#		$lineage_file = "nuclineage.spec";
		$lineage_file = "pseudogene_nuclineage.spec";
	}



	for ($j = 0; $j < 2; $j++) {
		if ($j == 0) {
			$evolutionary_step_option = "des";
		} else {
			$evolutionary_step_option = "trs";
		}
		&run("./indel-seq-gen -m $model -j $evolutionary_step_option < $treefile");
		&run("./indel-seq-gen -m $model -k $lineage_file -j $evolutionary_step_option < $treefile");
		&run("./indel-seq-gen -m $model -k $lineage_file -j $evolutionary_step_option -i 0.2 -w a < $treefile");
		&run("./indel-seq-gen -m $model -k $lineage_file -i 0.2 -j $evolutionary_step_option -w a -a 2.3 < $treefile");
		&run("./indel-seq-gen -m $model -k $lineage_file -i 0.2 -w a -a 2 -g 10 -j $evolutionary_step_option < $treefile");
		&run("./indel-seq-gen -m $model -k $lineage_file -j $evolutionary_step_option -i 0.2 -w a -f $freq_opt < $treefile");
	}
}

sub run {
	system($_[0]);
	print $_[0] . "\n";
	$wait = <STDIN>;
}
