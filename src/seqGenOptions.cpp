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
#include <sstream>
#include <unistd.h>
#include <getopt.h>
#include "inTree.h"
#include "inClade.h"

using namespace std;

extern bool MCMC_sample_evenly;
extern bool Pc, Qd, nij;

const char *seqGenOptions::modelTitles[] = {
                                            "JTT: Jones, Taylor & Thornton (1992) CABIOS  8:275-282\n              DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
                                            "WAG: Whelan & Goldman (2001) Mol Biol Evol 18:691-699",
                                            "PAM: Dayhoff, Schwartz & Orcutt (1978)\n              DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
                                            "BLOSUM62: Henikoff & Henikoff (1992) PNAS USA 89:10915-10919",
                                            "MTREV24: Adachi & Hasegawa (1996) J Mol Evol 42:459-468",
                                            "CPREV45: Adachi et al. (2000) J Mol Evol 50:348-358",
                                            "GENERAL: General time reversible (amino acids)"
                                           };

seqGenOptions::seqGenOptions():max_subseqs(50)
{
	init();
}

seqGenOptions::seqGenOptions(int max, int prec, bool deb):max_subseqs(max)
{
	init();
}

void seqGenOptions::SpoolWarnings(string warn, bool program_end) 
{
	if (program_end) {
		if ( !warnings.empty() ) cerr << "WARNINGS: " << endl;
		for (vector<string>::iterator it = warnings.begin(); it != warnings.end(); ++it)
			cerr << "  " << (*it) << endl;
	} else {
		bool identical = false;
		for (vector<string>::iterator it = warnings.begin(); it != warnings.end(); ++it)
			if ( (*it).compare(warn) == 0) identical = true;
		if (!identical) warnings.push_back(warn);
	}
}

void seqGenOptions::postProcess()
{
	stringstream ss (stringstream::in | stringstream::out);
	RateMatrix *temp_rates = new RateMatrix();
	//////////
	/// When we propose a path, we read in the necessary values from the input file rather than cmd-line.
	//////////
	if(!userSeed) {
		srand(time(NULL));
		randomSeed = CreateSeed();
	} else srand(1);	// Make run reproducible.

	if (path_proposal) {
		simulation_step_type = TIME_RELATIVE_STEPS;
		for (vector<bool>::iterator it = output_file_flags.begin(); it != output_file_flags.end(); ++it)
			(*it) = false;
		output_file_flags.at(ma) = true;
		output_file_flags.at(trace) = true;
		//////////
		/// With a path proposition, we need an output file name regardless of whether it is specified
		/// or not. Name to default "path_proposal".
		//////////
		if (output_files == DEFAULT_STRING) output_files = "path_proposal";
		return;
	}

	//////////
	/// Checks to see if model exists, sets numStates.
	//////////
	global_model = temp_rates->FindModel(inputModel);
	if (global_pi.empty()) {
		for (int i = 0; i < numStates; i++) {
			ss << 1.0/numStates;
			if (i != numStates-1) ss << ",";
			global_pi = ss.str();
		}
	} else ;//convertDelims(global_pi, ",", " ");
	if ( !rate_matrix_values.empty() ) ;//convertDelims(rate_matrix_values, ",", " ");
	if(default_rateHetero == DiscreteGammaRates) {
 		if (alpha == -1) {
 			cerr << "You must set the alpha value (-a <float>) when using discrete Gamma categories." << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (output_files == DEFAULT_STRING) 
		for (vector<bool>::iterator it = output_file_flags.begin(); it != output_file_flags.end(); ++it)
			(*it) = false;
	if (isNucModel && random_sequence_proportion_motif) {
		SpoolWarnings("PROSITE motifs can only be used for amino acid sequence runs. Deactivating option.", false);
		random_sequence_proportion_motif = 0;
	}
	
	delete temp_rates;
}

string convertDelims(
					 string& str, 
					 string fromDelim, 
					 string toDelim
					)
{
	stringstream ss (stringstream::in | stringstream::out);
	list<string> arg_split;

	arg_split = split(str,fromDelim);
	for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); ++it)
		ss << (*it) << toDelim;
	str = ss.str();
	str.erase(str.end()-1);	// Get rid of the extra delimiter at the end. //

	return str;
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
	random_sequence_proportion_motif		= 0.0;
	events_to_track[track_fossil]			= true;
	events_to_track[track_substitution]		= true;
	events_to_track[track_deletion]			= true;
	events_to_track[track_insertion]		= true;
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
	paleo_root_age = fossil_deposition_rate = 0.0;
	deposit_fossils							= false;
	lineageSpecificFile                     = "";
	default_rateHetero						= NoRates;
	invariableSites							= false;
	simulation_step_type					= TIME_RELATIVE_STEPS;
	global_model							= -1;
	userSeed								= false;
	filePath								= "./";
	for (int i = 0; i < MAX_RATE_CATS; i++) default_catRate[i] = 0;
	output_file_flags 					    = vector<bool>(6, true);
	output_file_flags[verb]					= false;
	perturbTree								= 1.0;
	end_point_condition_file				= "";
	path_proposal							= false;
	warnings.clear();
	kappa.clear();
	global_pi.clear();
	inputModel								= "";
	dependency_file.clear();
	dependence_superscript					= "1.0";		// In creating dependence model (3rd-order Markov), raise each random number to this power.
	forward_sim_event_file_to_emulate 		= "";
	epc_emulate_forward_simulation 			= false;
	num_mcmc_steps							= 0;
	context_order 							= 0;
	rasmus_independent_proposals			= false;
}

void seqGenOptions::InitGlobals() 
{
	isNucModel = -1;
}

void seqGenOptions::readOptions(
								int argc, 
								char  *argv[]
							   )
{
	int c;
	size_t i, j;
	list<string> arg_split;
	string codon_position_rates;
	string str;
	list<string> seeds;
	vector<unsigned int>::iterator jt;
	InitGlobals();
	list<string>::iterator it;
	enum proposal_style_t { independent_sites, QdPc, QdP, QPc, QP, Nij };
	size_t acc;

	while (1) {
	
		static struct option long_options[] =
		{
            {"alpha",     			required_argument,  0, 'a'},
            {"branch_scale",  		required_argument, 	0, 'b'},
            {"codon_rates",  		required_argument, 	0, 'c'},
			{"test_alternate",		no_argument,  		0, 'C'},
            {"select_outputs", 		required_argument, 	0, 'd'},
            {"dependency_input", 	required_argument, 	0, 'D'},
            {"outfile",  			required_argument, 	0, 'e'},
            {"path_proposal",  		required_argument, 	0, 'E'},
            {"frequencies",  		required_argument, 	0, 'f'},
			{"epc_Pij_frequencies", required_argument,  0, 'F'},
            {"num_gamma_cats",  	required_argument, 	0, 'g'},
            {"help",  				no_argument, 		0, 'h'},
			{"mcmc_proposals",		required_argument,	0, 'I'},
            {"invar",  				required_argument, 	0, 'i'},
            {"step_type",			required_argument,	0, 'j'},
            {"lineage",				required_argument,  0, 'k'},
            {"kappa",				required_argument,  0, 'K'},
            {"length",  			required_argument, 	0, 'l'},
            {"treelength",  		no_argument, 		0, 'L'},
            {"matrix",  			required_argument, 	0, 'm'},
			{"emulate_fwd_events",	required_argument,  0, 'M'},
            {"num_runs",  			required_argument, 	0, 'n'},
            {"outfile_format",  	required_argument, 	0, 'o'},
			{"context_order",		required_argument,  0, 'O'},
            {"path",				required_argument,	0, 'p'},
            {"paleo",           	required_argument,  0, 'P'},
            {"quiet",  				no_argument,		0, 'q'},
            {"rel_rates",  			required_argument, 	0, 'r'},
            {"output_width",  		required_argument,	0, 's'},
            {"tstv",  				required_argument, 	0, 't'},
            {"perturb_tree",		required_argument,	0, 'T'},
            {"indel_fill",  		required_argument, 	0, 'u'},
            {"mcmc_sample_evenly",	no_argument,		0, 'U'},
            {"all_gaps_out",		no_argument,		0, 'v'},
            {"write_anc", 			no_argument, 		0, 'w'},
            {"full_codon_indels", 	no_argument,		0, 'x'},
            {"mcmc",				required_argument,	0, 'y'},
            {"rng_seed",  			required_argument, 	0, 'z'},
            {"proportion_motif",	required_argument,  0, '1'},
            {0, 0, 0, 0}
		};

        // getopt_long stores the option index here.
        int option_index = 0;
     
        c = getopt_long (argc, argv, "a:b:c:Cd:D:e:E:f:F:g:hi:I:j:k:K:l:Lm:M:n:o:O:p:P:qr:s:t:u:Uvwxy:z:1:T:",long_options, &option_index);
     
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
				for (list<string>::iterator it = arg_split.begin(); it != arg_split.end(); ++it) {
					default_catRate[j] = atof((*it).c_str());
					if(default_catRate[j] <= 0) {
						cerr << "Bad Category Rates: " << codon_position_rates << endl;
						exit(EXIT_FAILURE);
					}
					j++;
				}
				break;
			case 'C':
				test_alternate_representations = true;
				break;
			case 'd':
				str = optarg;
				if (str.size() != 6) {
					cerr << "Need to specify 6 binary values of output file selection (-d)." << endl;
					exit(EXIT_FAILURE);
				}
				for (int k = 0; k < 6; k++)
					if (str.at(k) != '0') output_file_flags[k] = true;
					else output_file_flags[k] = false;
				//setEventTracking(optarg);
				break;
			case 'D':
				order_3_markov = true;
				dependence_superscript = optarg;
				break;
			case 'e':
				output_files = optarg;
				break;
			case 'E':
				end_point_condition_file = optarg;
				path_proposal = true;
				break;
			case 'f':
				global_pi = optarg;
				//SetFrequencies(optarg);
				break;
			case 'F':
				epc_Pij_pi = optarg;
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
				printNucHelp();
				exit(EXIT_SUCCESS);
//			case 'H':
//				break;
			case 'i':
				default_proportion_invariable = atof(optarg);
				if(default_proportion_invariable < 0.0 || default_proportion_invariable >= 1.0) {
					cerr << "Bad proportion invariable sites: " << default_proportion_invariable << endl;
					exit(EXIT_FAILURE);
				} 
				invariableSites = true;
				break;
			case 'I':
				acc = atoi(optarg);
				rasmus_independent_proposals = false; Qd = false; Pc=false;
				if (acc == independent_sites) rasmus_independent_proposals = true;
				else if (acc == QdPc) { Qd = Pc = true; nij = false; }
				else if (acc == QdP) { Qd = true; Pc = false; nij = false; }
				else if (acc == QPc) { Qd = false; Pc = true; nij = false; }
				else if (acc == QP) { Qd = false; Pc = false; nij = false; }
				else if (acc == Nij) { Qd = true; Pc = false; nij = true; }
				else {
					cerr << "Input option -I: Unknown proposal style \"" << acc << "\"." << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'j':
				if (strcmp(optarg,"trs") == 0) {
					simulation_step_type = TIME_RELATIVE_STEPS;
				} else if (strcmp(optarg, "des") == 0) {
					simulation_step_type = DISCRETE_EVOLUTIONARY_STEPS;
				} else if (strcmp(optarg, "gil") == 0) {
					simulation_step_type = UNIFORMIZATION;
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
				inputModel = optarg;
				break;
			case 'M':
				forward_sim_event_file_to_emulate = optarg;
				epc_emulate_forward_simulation = true;
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
			case 'O':
				context_order = atof(optarg);
				break;
			case 'p':
				filePath = optarg;
				break;
			case 'P':
				arg_split = split(optarg, ",");
				if (arg_split.size() != 2) {
					cerr << "Option --paleo (-P) expects two values." << endl;
					exit(EXIT_FAILURE);
				}
				deposit_fossils = true;
				paleo_root_age = atof(arg_split.front().c_str());			// In Mya
				fossil_deposition_rate = atof(arg_split.back().c_str());	// In Mya
				break;
			case 'q':
				quiet = true;
				break;
			case 'r':
				rate_matrix_values = optarg;
				break;
			case 's':
				output_width = atoi(optarg);
				break;
//			case 'S':
//				break;
			case 't':
				//equalTstv = 0;
				transition_vs_transversion_ratio = atof(optarg);
				if (transition_vs_transversion_ratio < 0.0) {
					cerr << "Transition versus transversion ratio \"" << transition_vs_transversion_ratio << "\" not within acceptable range of [0:1]." << endl;
					exit(EXIT_FAILURE);
				}
				//tstv = transition_vs_transversion_ratio;
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
			case 'U':
				MCMC_sample_evenly = true;
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
				num_mcmc_steps = atoi(optarg);
				break;
			case 'z':
				seeds = split(optarg,",");
				gil_Seed.assign(4, 0);
				userSeed = true;
				if (seeds.size() != 4) {
					cerr << "Option -z requires a comma-delimited list of 4 integers for the seed." << endl;
					exit(EXIT_FAILURE);
				}
				randomSeed = strtoul(seeds.front().c_str(), NULL, 0);
				jt = gil_Seed.begin();
				for (list<string>::iterator it = seeds.begin(); it != seeds.end(); ++it, ++jt) {
					(*jt) = atoi((*it).c_str());
				}
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
	
	postProcess();
}

void seqGenOptions::setEventTracking( 
									 string eventsToTrack 
									)
{
	if ( eventsToTrack.find_first_of("Ff") != string::npos ) events_to_track[track_fossil] = true;
	else events_to_track[track_fossil] = false;
	if ( eventsToTrack.find_first_of("Ss") != string::npos ) events_to_track[track_substitution] = true;
	else events_to_track[track_substitution] = false;
	if ( eventsToTrack.find_first_of("Ii") != string::npos ) events_to_track[track_insertion] = true;
	else events_to_track[track_insertion] = false;
	if ( eventsToTrack.find_first_of("Dd") != string::npos ) events_to_track[track_deletion] = true;
	else events_to_track[track_deletion] = false;

	if ( eventsToTrack.find_first_not_of("FfSsIiDd") != string::npos ) {
		cerr << "Unknown events included in option -y " << eventsToTrack << "." << endl;
		cerr << " Acceptable values are to track one or more of: f (fossil), s (substitution), i (insertion), and d (deletions)." << endl;
		exit(EXIT_FAILURE);
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
	cout <<	"-d --spec_outfiles    123456    Specify which files to output, where each is either a 0 or a 1:\n";
	cout << "                                  1: .root \n";
	cout << "                                  2: .seq  \n";
	cout << "                                  3: .ma   \n";
	cout << "                                  4: .tree \n";
	cout << "                                  5: .trace\n";
	cout << "                                  6: .verb \n";
	cout <<	"-e --outfile          <string>  Output file name for the .root, .seq, and .ma files.\n";
	cout <<	"-f --frequencies      <list>    4 or 20 floats separated by commas. Represent\n"; 
	cout << "                                nucleotide or amino acid frequencies.\n";
	cout << "                                order: ACGT or ARNDCQEGHILKMFPSTWYV, resp., or\n";
	cout <<	"                                use \"e\" for equal freqs. [default = matrix freqs].\n"; 
	cout <<	"-g --num_gamma_cats   <float>   Shape for gamma rate heterogeneity [default = none].\n";
	cout <<	"-h --help             none      Give this help message.\n";
	cout << "-H --Help             none      Print help message for nucleotide models (see --rel_rates).\n";
	cout <<	"-i --invar            <float>   Proportion invariable [default: none].\n";
	cout << "-j --step_type        <string>  Simulation step type:\n";
	cout << "                                 des = discrete evolutionary steps,\n";
	cout << "                                 trs = time relative steps,\n";
	cout << "                                 gil = time relative, using gillespie algorithm\n";
	cout << "-k --lineage          <string>  Filename that lists the lineage and motif specifications.\n";
	cout <<	"-l --length           <int>     Sequence length [default = 1000].\n"; 
	cout << "-L --treelength       <float>   Scale trees, where average root-to-tip path length equals input value.\n"; 
	cout <<	"-m --matrix           <string>  MODEL = JC69, K80, K81, K81ne, F81, HKY, F84, T92, TN93, TN93eq, TIM, TIMeq, TVM, TVMeq, SYM, GTR, PAM,\n"; 
	cout << "								 JTT, MTREV, GENERAL\n"; 
	cout << "                                Models listed before model GTR are for nucleotides and\n";
	cout << "                                the rest are for amino acids\n";
	cout <<	"-n --num_runs         <int>     Number of simulated datasets per tree [default = 1].\n"; 
	cout <<	"-o --outfile_format   <char>    Output file format (p = phylip, n = nexus, f = fasta)\n";
	cout <<	"                                [default = p].\n";
	cout << "-p --path                       Used only for suiteMSA run if indel-seq-gen.\n";
	cout << "-P --paleo            2x<float> Input tree branch lengths interpreted as time, rather than substitutions per\n";
	cout << "                                site. Two values are specified: The amount of time at the root of the input\n";
	cout << "                                tree (in millions of years), and a fossil deposition model?\n";
	cout <<	"-q --quiet            none      Quiet mode.\n";
	cout << "-r --rel_rates        <list>    Passes parameters to nucleotide rate matrices:\n";
	cout << "                                * a,b,c,d,e,f = general rate matrix (GTR and SYM).\n";
	cout <<	"                                * a,b,c,d     = TVM rate matrix only.\n";
	cout << "                                * a,b,c       = TIM rate matrix only.\n";
	cout << "                                * a,b         = TN93 (a=K1,b=K2) and K81 (a=K,b=alpha) rate matrices.\n";
	cout << "                                * a           = K80, HKY, and F84 (a=kappa for all).\n";
	cout <<	"-s --branch_scale     <float>   Branch length scaling factor [default = 1.0].\n"; 
	cout << "-t --tstv             <float>   transition-transversion ratio.\n";
	cout << "-T --perturb_tree     <float>   Randomly rescale each branch by value between [1/value, value].\n";
	cout <<	"-u --indel_fill       <string>  INDEL FILL MODEL, based on neighbor effects:\n";
	cout <<	"                                  'xia' = E. coli K-12 proteins,\n";
	cout <<	"                                  'sp'  = all Swiss-Prot proteins,\n"; 
	cout <<	"                                  'ran' = no neighbor effect\n";
	cout << "-v --all_gaps_out     none      Output all-gap columns in the true alignment.\n";
	cout <<	"-w --write_anc        none      Write ancestral sequences.\n";
	cout <<	"-z --rng_seed         <int>     Manually set the random number seed.\n";
	cout << "-1 --proportion_motif <float>   Set maximum proportion of sequence to set as PROSITE motif.\n";
	cout << "                                Valid only for random root sequences.\n\n";
	cout <<	"tree_file             <string>  Name of the specification file.\n";
	cout << "                                (See manual for format).\n\n";
}

void seqGenOptions::printNucHelp()
{

	cout << PROGRAM_NAME << VERSION_NUMBER << endl << endl;
	cout << "This help message shows where each value for the input relative rates (--rel_rates, -r)\n";
	cout << "are applied for nucleotide models.\n\n";
	cout << "-r <list> or --rel_rates <list>\n\n";
	cout << " <list> = a -----------------------------------\n";
	cout << "     model K80/HKY/F84\n";
	cout << "           A C G T  \n";
	cout << "         | *   a   |\n";
	cout << "     Q = |   *   a |\n";
	cout << "         | a   *   |\n";
	cout << "         |   a   * |\n";
	cout << "\n";
	cout << " <list> = a,b ---------------------------------\n";
	cout << "     model K81/K81ne    model TN93/TN93eq\n";
	cout << "           A C G T            A C G T\n";
	cout << "         | *   a b |        | *   a   |\n";
	cout << "     Q = |   * b a |    Q = |   *   b |\n";
    cout << "         | a b *   |        | a   *   |\n";
    cout << "         | b a   * |        |   b   * |\n";
	cout << "\n";
	cout << " <list> = a,b,c -------------------------------\n";
	cout << "     model TIM/TIMeq\n";
	cout << "           A C G T   \n";
	cout << "         | *   a b | \n";
	cout << "     Q = |   * b c | \n";
    cout << "         | a b *   | \n";
    cout << "         | b c   * | \n";
	cout << "\n";
	cout << " <list> = a,b,c,d -----------------------------\n";
	cout << "     model TVM/TVMeq\n";
	cout << "           A C G T   \n";
	cout << "         | *   a b | \n";
	cout << "     Q = |   * c b | \n";
    cout << "         | a c * d | \n";
    cout << "         | b b d * | \n";
	cout << "\n";
	cout << " <list> = a,b,c,d,e,f--------------------------\n";
	cout << "     model GTR/SYM\n";
	cout << "           A C G T   \n";
	cout << "         | * a b c | \n";
	cout << "     Q = | a * d e | \n";
    cout << "         | b d * f | \n";
    cout << "         | c e f * | \n";
	cout << "\n";
	cout << " ----------------------------------------------\n";
	cout << "\n\n";

	cout << "If you want general help, use option \"-h\" or \"--help\".\n";
}

bool seqGenOptions::tracking_defined()
{
	if (
		events_to_track[track_fossil] ||
		events_to_track[track_substitution] ||
		events_to_track[track_insertion] ||
		events_to_track[track_deletion]
	   ) return true;
	return false;
}