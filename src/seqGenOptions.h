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

#ifndef SEQGENOPTIONS_H_
#define SEQGENOPTIONS_H_

#if HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef ULONG_MAX
#define ULONG_MAX ((unsigned long)(0xffffffff))
#endif

using namespace std;

#include <string>
#include<list>
#include "model.h"
#include "insert_freqs.h"
#include "propose_path.h"

//////////
/// Definitions
//////////
#define DISCRETE_EVOLUTIONARY_STEPS 0
#define TIME_RELATIVE_STEPS 1
#define UNIFORMIZATION -1
#define ENDPOINT_CONDITIONED 2
#define INDEPENDENT_ENDPOINT_CONDITIONED 3
#define PROGRAM_NAME "indel-Seq-Gen "
#define VERSION_NUMBER "v2.1.03"
#define DEFAULT_STRING "!"

class motifInfo;

extern double CpG_multiplier;
extern bool	  order_3_markov;
extern bool   Human_Data_simulation;

enum {
	PHYLIPFormat,
	RelaxedFormat,
	NEXUSFormat,
	FASTAFormat
};

string convertDelims(string& str, string fromDelim, string toDelim);


class seqGenOptions : private Counter<seqGenOptions>
{
private:
	void init();
public:
	using Counter<seqGenOptions>::howMany;
	static const char *modelTitles[];
	
	bool	deposit_fossils;
	bool	treelength_check;
	bool 	treelength_scale;
	bool 	userSeed;
	bool 	invariableSites;
	bool	writeAncestors;
	bool 	quiet;
	bool 	help;
	bool 	inputRoot;
	bool	events_to_track[4];		// Ins, Del, Sub, Fossil
	bool	codon_only_indels;
	bool	output_indel_gaps;
	bool	path_proposal;
	bool	rasmus_independent_proposals;
	vector<bool> output_file_flags;
	char 	output_file_format;
	char 	write_site_rates_or_ancestral_sequences;
	double 	alpha;
	double 	default_branchScale;
	double 	default_proportion_invariable;
	double	perturbTree;
	double	random_sequence_proportion_motif;
	double 	transition_vs_transversion_ratio;
	double  paleo_root_age, fossil_deposition_rate;
	double  default_catRate[MAX_RATE_CATS];
	double epc_transition_scale;
	double context_order;
	string dependence_superscript;
	vector<double>	kappa;	// For matrices such as HKY85, F84, and such.
	int		default_rateHetero;
	int 	max_subseqs;
	int 	output_width;
	int 	fileFormat;
	int 	length;
	int 	number_of_datasets_per_tree;
	int 	ancestral_sequence_number;
	int 	num_discrete_gamma_categories;
	int		simulation_step_type;
	int		global_model;
	int 	num_mcmc_steps;
	string	end_point_condition_file;
	string	inputModel;
	string 	rate_matrix_values;
	string  indel_fill_model;
	string 	output_files;
	string 	lineageSpecificFile;
	string	filePath;
	string	global_pi;
	string	epc_Pij_global_pi;
	string	epc_Pij_pi;
	string 	dependency_file;
	string	forward_sim_event_file_to_emulate;	// These two variables are used to emulate forward simulation
	string	dependence_model_counts;
	string 	neutral_model_counts;
	bool	epc_emulate_forward_simulation;  	// events, should be helpful to judge Qi.|k(t) vs Qi.

	string end_point_branch_length;
	string end_point_start_sequence;
	string end_point_target_sequence;

	unsigned long randomSeed;
	vector<unsigned int> gil_Seed;
	vector<string>	warnings;
	
	enum {
		track_insertion,
		track_deletion,
		track_substitution,
		track_fossil
	};

	enum {
		root,
		seq,
		ma,
		tree,
		trace,
		verb
	};

	bool tracking_defined();
	list<motifInfo*> motif_list;
	seqGenOptions();
	seqGenOptions(int, int, bool);
	void setEventTracking(string eventsToTrack);
	void SpoolWarnings(string warn, bool program_end = false);
	void setModelConditions(string rate_matrix, string rates);
	void readOptions(int argc, char *argv[]);
	void printHelp();
	void printNucHelp();
	void InitGlobals();
	void postProcess();
};

#endif /*SEQGENOPTIONS_H_*/
