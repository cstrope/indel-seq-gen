// Copyright (C) 2008 Cory Strope <corystrope@gmail.com>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

#include "seqGenOptions.h"
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include "inTree.h"
#include "inClade.h"

using namespace std;

const int seqGenOptions::JTT             = 0;
const int seqGenOptions::WAG             = 1;
const int seqGenOptions::PAM             = 2;
const int seqGenOptions::BLOSUM          = 3;
const int seqGenOptions::MTREV           = 4;
const int seqGenOptions::CPREV           = 5;
const int seqGenOptions::GENERAL         = 6;
const char *seqGenOptions::modelTitles[] = {
                                            "JTT: Jones, Taylor & Thornton (1992) CABIOS  8:275-282\n              DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
                                            "WAG: Whelan & Goldman (2001) Mol Biol Evol 18:691-699",
                                            "PAM: Dayhoff, Schwartz & Orcutt (1978)\n              DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
                                            "BLOSUM62: Henikoff & Henikoff (1992) PNAS USA 89:10915-10919",
                                            "MTREV24: Adachi & Hasegawa (1996) J Mol Evol 42:459-468",
                                            "CPREV45: Adachi et al. (2000) J Mol Evol 50:348-358",
                                            "GENERAL: General time reversible (amino acids)"
                                           };

seqGenOptions::seqGenOptions():max_subseqs(50),debug(false)
{
	init();
}

seqGenOptions::seqGenOptions(int max, int prec, bool deb):max_subseqs(max),debug(deb)
{
	init();
}

void seqGenOptions::SpoolWarnings(string warn, bool program_end) 
{
	if (program_end) {
		if ( !warnings.empty() ) cerr << "WARNINGS: " << endl;
		for (vector<string>::iterator it = warnings.begin(); it != warnings.end(); it++)
			cerr << "  " << (*it) << endl;
	} else {
		bool identical = false;
		for (vector<string>::iterator it = warnings.begin(); it != warnings.end(); it++)
			if ( (*it).compare(warn) == 0) identical = true;
		if (!identical) warnings.push_back(warn);
	}
}

void seqGenOptions::postProcess()
{
	if (output_files == DEFAULT_STRING) {
		for (vector<bool>::iterator it = output_file_flags.begin(); it != output_file_flags.end(); it++)
			(*it) = false;
		output_file_flags_unset = true;
	} else output_file_flags_unset = false;

	if(global_model == -1 && !treelength_check) {
		cerr << "You must specify the substitution model (-m)." << endl << endl;
		exit(EXIT_FAILURE);
	}
	if(default_rateHetero == DiscreteGammaRates && alpha == -1) {
		cerr << "You must set the alpha value (-a <float>) when using discrete Gamma categories." << endl;
		exit(EXIT_FAILURE);
	}

	if(!userSeed) {
		srand(time(NULL));
		randomSeed = CreateSeed();
	} else srand(1);	// Make run reproducible.

	if (isNucModel && random_sequence_proportion_motif) {
		SpoolWarnings("PROSITE motifs can only be used for amino acid sequence runs. Deactivating option.", false);
		random_sequence_proportion_motif = 0;
	}
	if (help) printHelp();

}

void seqGenOptions::init()
{
	output_width                            = 60;
	alpha                                   = -1;
	output_files                            = "!";
	output_indel_gaps						= false;
	length                                  = 1000;
	indel_fill_model                        = "random";
	codon_only_indels						= false;
	number_of_datasets_per_tree             = 1;
	default_branchScale                    	= 1;
	num_discrete_gamma_categories           = -1;
	default_proportion_invariable           = 0;
	default_AAorNUC_frequencies             = "!";
	random_sequence_proportion_motif		= 0.0;
	events_to_track[0]						= true;
	events_to_track[1]						= true;
	events_to_track[2]						= false;
	treelength_check						= false;
	transition_vs_transversion_ratio        = -1;
	ancestral_sequence_number               = -1;
	output_file_format                      = 'p';
	fileFormat								= PHYLIPFormat;
	writeAncestors 							= false;
	write_site_rates_or_ancestral_sequences = '!';
	quiet                                   = false;
	help                                    = false;
	inputRoot                               = false;
	lineageSpecificFile                     = "";
	default_rateHetero						= NoRates;
	invariableSites							= false;
	simulation_step_type					= DISCRETE_EVOLUTIONARY_STEPS;
	global_model							= -1;
	userSeed								= false;
	filePath								= "./";
	for (int i = 0; i < MAX_RATE_CATS; i++) default_catRate[i] = 0;
	output_file_flags 					    = vector<bool>(6, true);
	output_file_flags[verb]					= false;
	output_file_flags_unset					= true;
	perturbTree								= 1.0;
	warnings.clear();
}

void seqGenOptions::InitGlobals() 
{
	isNucModel = -1;
	equalFreqs = 1;
	equalTstv = 1;
	tstv=0.50002;
	for (int i = 0; i < NUM_AA_REL_RATES; i++) {
		aaRelativeRate[i] = 1.0;
	}
	for (int i = 0; i < NUM_AA; i++) {
		aaFreq[i] = 1.0;
	}
	aaFreqSet = 0;
}

void seqGenOptions::readOptions(int argc, char  *argv[])
{
	int c;
	size_t i, j;
	list<string> arg_split;
	string rate_matrix_values;
	string codon_position_rates;
	string events2track;
	string str;
	
	InitGlobals();

	while (1) {
	
		static struct option long_options[] =
		{
            {"alpha",     		required_argument,  0, 'a'},
            {"branch_scale",  	required_argument, 	0, 'b'},
            {"codon_rates",  	required_argument, 	0, 'c'},
            {"select_outputs", 	required_argument, 	0, 'd'},
            {"outfile",  		required_argument, 	0, 'e'},
            {"frequencies",  	required_argument, 	0, 'f'},
            {"num_gamma_cats",  required_argument, 	0, 'g'},
            {"help",  			no_argument, 		0, 'h'},
            {"invar",  			required_argument, 	0, 'i'},
            {"step_type",		required_argument,	0, 'j'},
            {"lineage",			required_argument,  0, 'k'},
            {"length",  		required_argument, 	0, 'l'},
            {"treelength",  	no_argument, 		0, 'L'},
            {"matrix",  		required_argument, 	0, 'm'},
            {"num_runs",  		required_argument, 	0, 'n'},
            {"outfile_format",  required_argument, 	0, 'o'},
            {"path",			required_argument,	0, 'p'},
            {"quiet",  			no_argument,		0, 'q'},
            {"rel_rates",  		required_argument, 	0, 'r'},
            {"output_width",  	required_argument,	0, 's'},
            {"tstv",  			required_argument, 	0, 't'},
            {"perturb_tree",	required_argument,	0, 'T'},
            {"indel_fill",  	required_argument, 	0, 'u'},
            {"indel_gaps",		no_argument,		0, 'v'},
            {"write_anc", 		no_argument, 		0, 'w'},
            {"full_codon_indels", no_argument,		0, 'x'},
            {"track_events",	required_argument,	0, 'y'},
            {"rng_seed",  		required_argument, 	0, 'z'},
            {"proportion_motif",required_argument,  0, '1'},
            {0, 0, 0, 0}
		};

        // getopt_long stores the option index here.
        int option_index = 0;
     
        c = getopt_long (argc, argv, "a:b:c:d:e:f:g:hi:j:k:l:Lm:n:o:p:qr:s:t:u:vwxy:z:1:T:",long_options, &option_index);
     
        // Detect the end of the options.
        if (c == -1) break;
	
		switch(c) {
			case 'a':
				if(default_rateHetero==CodonRates) {
					cerr << "You can only have codon rates or gamma rates, not both" << endl << endl;
					exit(EXIT_FAILURE);
				}
				if(default_rateHetero==NoRates) default_rateHetero=GammaRates;
				alpha = atof(optarg);
				if(alpha <= 0.0) {
					cerr << "Bad Gamma Shape: " << alpha << endl << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'b':
				default_branchScale = atof(optarg);
				if (default_branchScale < 0.0) {
					cerr << "Scaling average root-to-tip path length to " << (-default_branchScale) << "." << endl;
					treelength_scale = true;
				} else if(default_branchScale == 0.0) {
					cerr << "Bad branch length scale: " << default_branchScale << endl << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'c':
				if(!isNucModel) {
					cerr << "You can only have codon rates when using nucleotide models" << endl << endl;
					exit(EXIT_FAILURE);
				}
				if(default_rateHetero==GammaRates) {
					cerr << "You can only have codon rates or gamma rates, not both" << endl << endl;
					exit(EXIT_FAILURE);
				}
				num_discrete_gamma_categories = 3;
				default_rateHetero=CodonRates;
				codon_position_rates = optarg;
				arg_split = split(codon_position_rates, ",");
				if(arg_split.size() != 3) {
					cerr << "You must specify rates for the 3 categories of codon positions when using the -c option" << endl;
					exit(EXIT_FAILURE);
				}
				j = 0;
				for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); it++) {
					default_catRate[j] = atof((*it).c_str());
					if(default_catRate[j] <= 0) {
						cerr << "Bad Category Rates: " << codon_position_rates << endl;
						exit(EXIT_FAILURE);
					}
					j++;
				}
				break;
			case 'd':
				str = optarg;
				if (str.size() != 6) {
					cerr << "Need to specify 6 binary values of output file selection (-d)." << endl;
					exit(EXIT_FAILURE);
				}
				output_file_flags_unset = false;
				for (int k = 0; k < 6; k++)
					if (str.at(k) != '0') output_file_flags[k] = true;
					else output_file_flags[k] = false;
				break;
			case 'e':
				output_files = optarg;
				break;
			case 'f':
				SetFrequencies(optarg);
				break;
			case 'g':
				if(default_rateHetero==CodonRates) {
					cerr << "You can only have codon rates or gamma rates, not both" << endl << endl;
					exit(EXIT_FAILURE);
				}
				num_discrete_gamma_categories = atoi(optarg);
				default_rateHetero = DiscreteGammaRates;
				if(num_discrete_gamma_categories < 2 || num_discrete_gamma_categories > MAX_RATE_CATS) {
					cerr << "The number of Gamma categories specified \"" << num_discrete_gamma_categories << "\" must be between 2..32 " << endl << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'h':
				printHelp();
				break;
			case 'i':
				default_proportion_invariable = atof(optarg);
				if(default_proportion_invariable < 0.0 || default_proportion_invariable >= 1.0) {
					cerr << "Bad proportion invariable sites: " << default_proportion_invariable << endl;
					exit(EXIT_FAILURE);
				} 
				invariableSites = true;
				break;
			case 'j':
				if (strcmp(optarg,"trs") == 0) {
					simulation_step_type = TIME_RELATIVE_STEPS;
				} else if (strcmp(optarg, "des") == 0) {
					simulation_step_type = DISCRETE_EVOLUTIONARY_STEPS;
				} else if (strcmp(optarg, "gil") == 0) {
					simulation_step_type = GILLESPIE;
				} else {
					cerr << "Invalid option for simulation type. Value must be \"trs\" or \"des\"" << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'k':
				lineageSpecificFile = optarg;
				break;
			case 'l':
				SpoolWarnings("Length is no longer a necessary options.");
				break;
			case 'L':
				treelength_check = true;
				quiet = true;
				break;
			case 'm':
				global_model = FindModel(optarg);
				break;
			case 'n':
				number_of_datasets_per_tree = atoi(optarg);
				break;
			case 'o':
				output_file_format = toupper(optarg[0]);
				if(output_file_format == 'P') {
					fileFormat=PHYLIPFormat;
				} else if (output_file_format == 'R') {
					fileFormat=PHYLIPFormat;
				} else if (output_file_format == 'F') {
					fileFormat=FASTAFormat;
				} else if (output_file_format == 'N') {
					fileFormat=NEXUSFormat;
				} else {
					cerr << "Invalid output file format \"-o" << output_file_format << "\"\n";
					exit(EXIT_FAILURE);
				}
				break;
			case 'p':
				filePath = optarg;
				break;
			case 'q':
				quiet = true;
				break;
			case 'r':
				rate_matrix_values = optarg;
				if(global_model == GTR) {
					arg_split = split(rate_matrix_values, ",");
					j = 0;
					for(list<string>::iterator it = arg_split.begin(); it != arg_split.end(); it++) {
						nucRelativeRates[j] = atof((*it).c_str());
						j++;
					}
					if(j != NUM_NUC_REL_RATES) {
						cerr << "Bad general nucleotide rate matrix: " << rate_matrix_values << endl << endl;
						exit(EXIT_FAILURE);
					}
					if(nucRelativeRates[NUM_NUC_REL_RATES-1] != 1.0) {
						for(j=0; j<NUM_NUC_REL_RATES-1; j++) 
							nucRelativeRates[j] /= nucRelativeRates[NUM_NUC_REL_RATES-1];
						nucRelativeRates[NUM_NUC_REL_RATES-1] = 1.0;							
					}
				} else if(global_model == GENERAL) {
					arg_split = split(rate_matrix_values, ",");
					j=0;
					for(list<string>::iterator it = arg_split.begin(); it != arg_split.end(); it++) {
						aaRelativeRate[j] = atof((*it).c_str());
						j++;
					}
					if(j != NUM_AA_REL_RATES) {
						cerr << "Bad General amino acid rate matrix: " << rate_matrix_values << endl << endl;
						exit(EXIT_FAILURE);
					}
				} else {
					cerr << "You can only have a general rate matrix when using GTR or GENERAL models" << endl << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 's':
				output_width = atoi(optarg);
				break;
			case 't':
				if(global_model != HKY && global_model != F84) {
					cerr << "You can only have a transition/transversion ratio when using HKY or F84 models" << endl << endl;
					exit(EXIT_FAILURE);
				}
				equalTstv = 0;
				transition_vs_transversion_ratio = atof(optarg);
				if (transition_vs_transversion_ratio < 0.0 || transition_vs_transversion_ratio > 1.0) {
					cerr << "Transition versus transversion ratio \"" << transition_vs_transversion_ratio << "\" not within acceptable range of [0:1]." << endl;
					exit(EXIT_FAILURE);
				}
				tstv = transition_vs_transversion_ratio;
				break;
			case 'T':
				perturbTree = atof(optarg);
				if (perturbTree <= 0.0) {
					cerr << "Cannot scale branch lengths by a multiplier less than or equal to zero." << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'u':
				indel_fill_model = optarg;
				if (indel_fill_model == "random" || indel_fill_model == "ran") {
					indel_fill_aas = RANDOM_FILL;
				} else if (indel_fill_model == "xia") {
					indel_fill_aas = XIA_FILL;
				} else if (indel_fill_model == "sp") {
					indel_fill_aas = SP_FILL;
				} else {
					cerr << "Invalid indel sequence fill model \"-u " << indel_fill_model << "\" "
						 << "(\"random\", \"ran\", \"xia\", and \"sp\" accepted)." << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'v':
				output_indel_gaps = true;
				break;
			case 'w':
				write_site_rates_or_ancestral_sequences = 'A';
				switch(write_site_rates_or_ancestral_sequences) {
				case 'A': writeAncestors = 1; break;
				default:
					cerr << "Unknown write mode: " << write_site_rates_or_ancestral_sequences << endl << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'x':
				codon_only_indels = true;
				break;
			case 'y':
				events2track = optarg;
				arg_split = split(events2track, ",");
				if (arg_split.size() != 3) {
					cerr << "Error parsing flag -y (--track_events). Argument requires a comma delimited list (without spaces) of size 3. Ex: a,b,c" << endl;
					exit(EXIT_FAILURE);
				}
				
				i = 0;
				for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); it++, i++) {
					if ((*it).c_str() == "1") events_to_track[i] = true;
					else if ((*it).c_str() == "0") events_to_track[i]= false;
					else {
						cerr << "Error parsing flag -y (--track_events). Arguments in list must be either a '1' or a '0'. Order of arguments is insertion, deletion, substitution." << endl;
						exit(EXIT_FAILURE);
					}
				}
				if (events_to_track[0]) cerr << "Insertions" << endl;
				if (events_to_track[1]) cerr << "Deletions" << endl;
				if (events_to_track[2]) cerr << "Substitutions" << endl;
				break;
			case 'z':
				userSeed = true;
				randomSeed = strtoul(optarg, NULL, 0);
				if (randomSeed == 0) {
					userSeed = false;
					SpoolWarnings("Could not convert input random seed. Performed simulation with computed seed.");
				} else if (randomSeed == ULONG_MAX) {
					cerr << "Input random seed number " << optarg << " is outside of unsigned long integer maximum value. Exiting..." << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case '1':
				random_sequence_proportion_motif = atof(optarg);
				if (random_sequence_proportion_motif > 1.0 || random_sequence_proportion_motif < 0.0) {
					cerr << "--proportion_motif (-1) must be a floating point value between 0.0 and 1.0.\n";
					exit(EXIT_FAILURE);
				}
				break;
			case '?':
				printHelp();
				exit(EXIT_FAILURE);
				break;
			default:
				abort();
		}
	}
}

void seqGenOptions::printHelp()
{

	cout << PROGRAM_NAME << VERSION_NUMBER << endl;
	cout << "\n\tUsage: ./indel-seq-gen [-bdefghilmnoqsuwz] < [tree_file]\n\n";
	
	cout << "Options (short long)  Input     Description" << endl;
	cout << "--------------------  -----     -----------\n";
	cout << "-a --alpha            <float>   Shape (alpha) for gamma rate heterogeneity.\n";
	cout << "-b --option_width     <int>     The number of residues per line on the multiple\n"; 
	cout << "                                alignment output [default = 60].\n";
	cout << "-c --codon_rates      <list>    #1,#2,#3 = rates for codon position heterogeneity.\n";
	cout <<	"-d --tree_scale       <float>   Total tree scale [default = use branch lengths].\n";
	cout <<	"-e --outfile          <string>  Output file name for the .root, .seq, and .ma files.\n";
	cout <<	"-f --frequencies      <list>    4 or 20 floats separated by commas. Represent\n"; 
	cout << "                                nucleotide or amino acid frequencies.\n";
	cout << "                                order: ACGT or ARNDCQEGHILKMFPSTWYV, resp., or\n";
	cout <<	"                                use \"e\" for equal freqs. [default = matrix freqs].\n"; 
	cout <<	"-g --num_gamma_cats   <float>   Shape for gamma rate heterogeneity [default = none].\n";
	cout <<	"-h --help             none      Give this help message.\n";
	cout <<	"-i --invar            <float>   Proportion invariable [default: none].\n";
	cout <<	"-l --length           <int>     Sequence length [default = 1000].\n"; 
	cout <<	"-m --matrix           <string>  MODEL = HKY, F84, GTR, PAM, JTT, MTREV, GENERAL\n"; 
	cout << "                                HKY, F84 and GTR are for nucleotides and the rest\n";
	cout << "                                are for amino acids\n";
	cout <<	"-n --num_runs         <int>     Number of simulated datasets per tree [default = 1].\n"; 
	cout <<	"-o --outfile_format   <char>    Output file format (p = phylip, n = nexus, f = fasta)\n";
	cout <<	"                                [default = p].\n";
	cout <<	"-q --quiet            none      Quiet mode.\n";
	cout << "-r --rel_rates        <list>    #1,#2,#3,#4,#5,#6 = general rate matrix.\n";
	cout <<	"-s --branch_scale     <float>   Branch length scaling factor [default = 1.0].\n"; 
	cout << "-t --tstv             <float>   transition-transversion ratio.\n";
	cout << "-T --perturb_tree     <float>   Randomly rescale each branch by value between [1/value, value].\n";
	cout <<	"-u --indel_fill       <string>  INDEL FILL MODEL, based on neighbor effects:\n";
	cout <<	"                                  'xia' = E. coli K-12 proteins,\n";
	cout <<	"                                  'sp'  = all Swiss-Prot proteins,\n"; 
	cout <<	"                                  'ran' = no neighbor effect\n";
	cout <<	"-w --write_anc        none      Write ancestral sequences.\n";
	cout <<	"-z --rng_seed         <int>     Manually set the random number seed.\n";
	cout << "-1 --proportion_motif <float>   Set maximum proportion of sequence to set as PROSITE motif.\n";
	cout << "                                Valid only for random root sequences.\n\n";
	cout <<	"tree_file             <string>  Name of the specification file.\n";
	cout << "                                (See manual for format).\n\n";
	
	exit(EXIT_SUCCESS);
}
