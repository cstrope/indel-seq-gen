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

//////////
/// Definitions
//////////
#define DISCRETE_EVOLUTIONARY_STEPS 0
#define TIME_RELATIVE_STEPS 1
#define GILLESPIE 2
#define PROGRAM_NAME "indel-Seq-Gen "
#define VERSION_NUMBER "v2.1.0"
#define DEFAULT_STRING "!"

#include <string>
#include "model.h"
#include "aamodels.h"
#include "nucmodels.h"
#include "insert_freqs.h"
#include<list>

using namespace std;

class motifInfo;

enum {
	PHYLIPFormat,
	RelaxedFormat,
	NEXUSFormat,
	FASTAFormat
};

enum {
	root,
	seq,
	ma,
	tree,
	trace,
	verb
};

class seqGenOptions : private Counter<seqGenOptions>
{
private:
	void init();
public:
	using Counter<seqGenOptions>::howMany;
	static const int JTT;
	static const int WAG;
	static const int PAM;
	static const int BLOSUM;
	static const int MTREV;
	static const int CPREV;
	static const int GENERAL;
	static const char *modelTitles[];
	
	bool	treelength_check;
	bool 	treelength_scale;
	bool 	userSeed;
	bool 	invariableSites;
	bool	writeAncestors;
	bool 	quiet;
	bool 	help;
	bool 	inputRoot;
	bool 	debug;
	bool	events_to_track[3];		// Ins, Del, ?Sub
	bool	codon_only_indels;
	bool	output_indel_gaps;
	vector<bool> output_file_flags;
	bool	output_file_flags_unset; // If trace is suppressed, need to know if any of the other flags are set.
	char 	output_file_format;
	char 	write_site_rates_or_ancestral_sequences;
	double 	alpha;
	double 	default_branchScale;
	double 	default_proportion_invariable;
	double	perturbTree;
	double	random_sequence_proportion_motif;
	double 	transition_vs_transversion_ratio;
	double  default_catRate[MAX_RATE_CATS];
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
	string  indel_fill_model;
	string 	output_files;
	string 	default_AAorNUC_frequencies;
	string 	lineageSpecificFile;
	string	filePath;
	unsigned long randomSeed;
	vector<string>	warnings;
	
	list<motifInfo*> motif_list;
	seqGenOptions();
	seqGenOptions(int, int, bool);
	void SpoolWarnings(string warn, bool program_end = false);
	void readOptions(int argc, char *argv[]);
	void printHelp();
	void InitGlobals();
	void postProcess();
};

#endif /*SEQGENOPTIONS_H_*/
