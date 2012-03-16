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
#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <list>
#include <cstring>
#include "tree.h"
#include "seqGenOptions.h"
#include "inTree.h"
#include "treefile.h"
#include "evolve.h"
#include "model.h"
#include "nucmodels.h"
#include "aamodels.h"
#include "progress.h"
#include "twister.h"
#include "random.h"
#include "trace.h"

bool	output_trs_example = false;
int		run_no = 0;
bool	empty_column_warning_spooled = false;

int numDatasets, numTrees;
int ancestorSeq;
// Motif accuracy checking:
int total_motif_positions = 0;
int correct_motif_positions = 0;
int num_template_violations = 0;
ofstream verbose_output;
extern ofstream verbose_output;
vector<int> num_siteProperty_copies (4,0);
int num_copy_ctor = 0, num_dtor = 0, num_ctor = 0;
int num_insert = 0, num_delete = 0, len_insert = 0, len_delete = 0;
size_t TNode_InheritMotifSites = 0;
vector<int> motifSite_setActiveProps (3,0);
int num_ms_siteProps = 0;
int num_inserted_positions = 0;
int num_deleted_positions = 0;

using namespace std;

// MAIN prototypes

void Print_Root(list<inTree*> inputTrees, ostream& root_out, seqGenOptions *options);
void Print_Trace(list<inTree*> inputTrees, list<eventTrack*> events, ostream& trace_out, int step_type, int treeNo, int setNo);
void Print_Seq(list<inTree*> inputTrees, ostream& seq_out, seqGenOptions *options, int dataset_num);
void Print_MA(list<inTree*> inputTrees, ostream& ma_out, seqGenOptions *options, int dataset_num, list<eventTrack*> events);
void Print_Motifs(list<inTree*> inputTrees, ostream& motif_out);
void WriteAncestralMA(vector<vector<char> >&, int *print_array_row, TNode *branch, int *n, TTree *tree, char **anc_names, seqGenOptions *options);
void WriteAncestralSequencesNode(vector<vector<char> >& print, int *print_array_row, TTree *tree, int *nodeNo, TNode *des, char **anc_names);
void PrintTitle();
void Print_NEXUS_Header(list<inTree*> inputTrees, seqGenOptions *options, ostream& out);
void Reset_Run(list<inTree*>& inTrees, seqGenOptions *options, bool reset);
void readInputTrees(list<inTree*>& inputTrees, 	seqGenOptions *options, inClade *global_environment, 
					vector<string>& OTU_names, int *numTrees, int numTaxa);
inline unsigned int rand_seed();
void setTRS(list<inTree*>& inputTrees, seqGenOptions *options, double max_path_length, double min_branch);
void CheckTrees(list<inTree*>& inputTrees, seqGenOptions *options);
size_t trimMSA(vector<vector<char> >& print, vector<bool>& empty_columns, seqGenOptions *options);
long print_memory(string message, size_t howMany, int object_size);
void Leakage();
void fiddle_with_trees(list<inTree*>& inTrees, seqGenOptions *options);

int main(int argc, char *argv[])
{	
	int numTaxa = 0;		// To check if all trees have the same number of taxa.
	clock_t totalStart;
	double totalSecs;
	vector<string> OTU_names;
	inClade *global_environment;
	vector<unsigned int> gil_Seed;
	gil_Seed.clear();

	totalStart = clock();

	// Options parsing, error-checking //
	seqGenOptions *options = new seqGenOptions();
	options->debug = true;
	options->readOptions(argc, argv);
	options->postProcess();
	if (gil_Seed.empty() ) {
		for (unsigned int i = 0; i < 4; i++) {
			gil_Seed.push_back(rand_seed());
		}
	}
	mt_srand(&gil_Seed[0], gil_Seed.size());
	SetSeed(options->randomSeed);
	model = options->global_model;
	SetModel(model, NULL);
	global_environment = new inClade("Global Environment", options);	

	//////////
	/// Gathering inTrees, checking sanity 
	//////////
	list<inTree*> inputTrees;
	int numTrees = 0;
	readInputTrees(inputTrees,options,global_environment,OTU_names,&numTrees,numTaxa);
	CheckTrees(inputTrees, options);

	if(options->writeAncestors) {
		string node_names = "";
		inputTrees.front()->CheckPhylogeneticAncestralNodes(node_names, SET);
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
			if(!(*it)->CheckPhylogeneticAncestralNodes(node_names, CHECK)) {
				options->writeAncestors = false;
				options->SpoolWarnings("Cannot write ancestral sequences unless all partition tree topologies are equivalent.");
				break;
			}
		}
	}

	string verbose_outfile;
	ofstream root_out, seq_out, ma_out, trace_out, tree_out;
	bool anc_tree_not_printed = true;
	string tree_outfile, root_outfile, seq_outfile, ma_outfile, trace_outfile;
	if(options->output_file_flags[tree]) {
		tree_outfile = options->output_files + ".anc_tree";
		tree_out.open(tree_outfile.c_str(), ios::trunc | ios::out);
	}
	if(options->output_file_flags[root]) {
		root_outfile = options->output_files + ".root";	
		root_out.open(root_outfile.c_str(), ios::trunc | ios::out);
	}
	if(options->output_file_flags[seq]) {
		seq_outfile = options->output_files + ".seq";
		seq_out.open(seq_outfile.c_str(), ios::trunc | ios::out);
	}
	if(options->output_file_flags[ma]) {
		ma_outfile = options->output_files + ".ma";
		ma_out.open(ma_outfile.c_str(), ios::trunc | ios::out);
	}
	if(options->output_file_flags[trace]) {
		trace_outfile = options->output_files + ".trace";
		trace_out.open(trace_outfile.c_str(), ios::trunc | ios::out);
	}
	if (options->output_file_flags[verb]) {
		verbose_outfile = options->output_files + ".verb";
		verbose_output.open(verbose_outfile.c_str(), ios::trunc | ios::out);
	}

	if (!options->quiet) {
		if (options->fileFormat == NEXUSFormat) 
			Print_NEXUS_Header(inputTrees,options,((options->output_file_flags[ma])?ma_out:cout));

		stringstream header;
		header << "Taxa bit representation:" << endl;
		inputTrees.front()->Print_Trace_Header(header, OTU_names, options);

		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
			header << "Clade List (<taxon members> = <clade name>): " << endl;
			for (list<TNode*>::iterator itn = (*it)->my_tree->nodeList.begin(); itn != (*it)->my_tree->nodeList.end(); itn++) {
				if ( !((*itn)->clade_label.empty()) ) {
					for (int j = 0; j < (*it)->my_tree->numTips; j++) header << (*itn)->bipartition.at(j);
					header << " = " << (*itn)->clade_label << endl;
				}
			}
		}

		if (!options->output_file_flags[trace]) cout << header.str() << endl;
		else trace_out << header.str() << endl;
	}

	list<eventTrack*> events;
	for (int i = 0; i < options->number_of_datasets_per_tree; i++, run_no++) {
		double sum = 0;
		int numSites = 0;
		indelNo = 0;
		events.clear();
		int size_prev_msa_positions = 0;

		//////////
		/// Set guide tree branch lengths.
		//////////
		//
		// (1) Apply all tree scaling options (perturbTree, scaleBranches || scaleTree)
		//
		fiddle_with_trees(inputTrees, options);
		//
		// (2) Calculate the min/max path lengths
		//
		double max_path_length = 0.0, thisTree_max_path, thisTree_min_branch1, thisTree_min_branch2;
		double min_branch = numeric_limits<double>::max();
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
			if ( (*it)->my_tree->rooted )  {
				thisTree_max_path = (*it)->CalculatePathLengths();
				if(thisTree_max_path > max_path_length) 
					max_path_length = thisTree_max_path;
				thisTree_min_branch1= (*it)->SearchMinBranch((*it)->my_tree->root->branch1,min_branch);
				thisTree_min_branch2= (*it)->SearchMinBranch((*it)->my_tree->root->branch2,min_branch);
				if(thisTree_min_branch2 < min_branch && thisTree_min_branch2 != 0.0) 
					min_branch = thisTree_min_branch2;
				if(thisTree_min_branch1 < min_branch && thisTree_min_branch1 != 0.0) 
					min_branch = thisTree_min_branch1;
			} 
		}
		//
		// (3) Finally, set TRS steps
		//
		setTRS(inputTrees, options, max_path_length, min_branch);
		//
		// Print out tree with perturbation (Perturbed tree will print out different trees for each
		// dataset) if the user wants the anc_trees.
		//
		if (options->output_file_flags[tree]) {
			if (options->simulation_step_type == TIME_RELATIVE_STEPS || options->simulation_step_type == GILLESPIE) {
				if (options->perturbTree != 1 || anc_tree_not_printed) {
					for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++)
						(*it)->Print_Time_Rel_Tree(tree_out, options->writeAncestors);
					anc_tree_not_printed = false;
				}
			} else {
				if (options->perturbTree != 1 || anc_tree_not_printed) {
					for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++)
						(*it)->Print_Tree(tree_out, false);
					anc_tree_not_printed = false;
				}
			}
		}
		//
		//////////
		//////////
		/// Setup of sequences for this run.
		//////////
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
			(*it)->Setup_Tree(options);
			sum += (*it)->partitionRate * (*it)->partitionLength;
			numSites += (*it)->partitionLength;
			if ((*it)->partitionLength < 2) {
				cerr << "Tree " << (*it)->treeNum << ": Partition length must be 2 sites or larger." << endl;
				exit(EXIT_FAILURE);
			}
			if ( isNucModel ) {
				if ( (*it)->my_tree->root->nodeEnv->rateHetero == CodonRates) {
					list<inTree*>::iterator rit = inputTrees.begin();
					while ( (*rit) != (*it) ) {
						if ( (*rit)->my_tree->root->nodeEnv->rateHetero == CodonRates ) {
							(*it)->codon_offset += (*rit)->partitionLength;
							(*it)->codon_offset = ((*it)->codon_offset) % 3;
						}
						rit++;
					}
				}
			}
		}
		CreateMatrix();

		//////////
		/// One final check to see if all coding regions add up to 0 (mod 3)
		//////////
		if (isNucModel) {
			for(list<inTree*>::reverse_iterator rit = inputTrees.rbegin(); rit != inputTrees.rend(); rit++) {
				if ( (*rit)->my_tree->root->nodeEnv->rateHetero == CodonRates ) {
					if ( ((*rit)->partitionLength + (*rit)->partitionLength) % 3 != 0 )
						options->SpoolWarnings("The lengths of all codon sequences do not add up to a factor of three.");
				}
			}
		}

		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) 
			(*it)->partitionRate *= numSites / sum;
		
		size_t x = 0;
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++, x++) {
			EvolveSequences(*it, &events, options);
			// Sort out the positions of the sequences.
			int msa_size = (*it)->my_tree->root_numSites;
			if (traceEvents) {
				for (list<eventTrack*>::iterator it2 = events.begin(); it2 != events.end(); it2++) {
					msa_size = (*it2)->Compute_MSA_Positions((*it)->my_tree, size_prev_msa_positions);
				}
				size_prev_msa_positions += msa_size;
			}
			if (options->simulation_step_type == DISCRETE_EVOLUTIONARY_STEPS && traceEvents) {
				Print_Trace(
							inputTrees, 
							events, 
							((options->output_file_flags[trace]) ? trace_out : cout), 
							options->simulation_step_type, 
							(*it)->treeNum, 
							i
						   );
				events.clear();		// DES prints each partition separately, since they cannot be related.
			}
		}

		if (options->output_file_flags[verb]) verbose_output << "DATASET " << i << ":" << endl;
		if (traceEvents) {
			if (options->simulation_step_type == TIME_RELATIVE_STEPS || options->simulation_step_type == GILLESPIE) {
				Print_Trace(inputTrees, 
							events, 
							((options->output_file_flags[trace]) ? trace_out : cout), 
							options->simulation_step_type, 
							0, 
							i
						   );
			} else if (options->output_file_flags[trace]) trace_out << endl;
			else cout << endl;
		}

		Print_Root(inputTrees, ((options->output_file_flags[root]) ? root_out : cout), options );
		Print_Seq(inputTrees, ((options->output_file_flags[seq]) ? seq_out : cout), options, i);
		Print_MA(inputTrees, ((options->output_file_flags[ma]) ? ma_out : cout), options, i, events);
		//Print_Motifs(inputTrees, ((options->output_files == "!") ? cout : cout)/*, options, i*/);

		Reset_Run(inputTrees, options, (i+1 != options->number_of_datasets_per_tree));

//		cerr << "Run " << i << " leakage:" << endl;
//		Leakage();
	}

	if (options->output_file_flags[tree]) tree_out.close();
	if (options->output_file_flags[root]) root_out.close();
	if (options->output_file_flags[seq]) seq_out.close();
	if (options->output_file_flags[ma]) ma_out.close();
	if (options->output_file_flags[verb]) verbose_output.close();
	size_t x = 1;
	if (options->output_file_flags[trace]) {
		for (list<inTree*>::iterator xt = inputTrees.begin(); xt != inputTrees.end(); xt++, x++)
			trace_out << "GUIDE_TREE_PARTITION" << x << "=" << (*xt)->tree << ";" << endl;
		trace_out.close();
	}

	options->SpoolWarnings("",true);

	totalSecs = (double)(clock() - totalStart) / CLOCKS_PER_SEC;
	if(!options->quiet) fprintf(stderr,"Time taken: %G seconds\n", totalSecs);

	//////////
	/// Cleanup all objects created in run.
	//////////
	delete global_environment;
	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
		for (list<TNode*>::iterator it2 = (*it)->my_tree->nodeList.begin(); it2 != (*it)->my_tree->nodeList.end(); it2++) {
			delete (*it2)->branch;
			delete (*it2);
		}
		for (list<inClade*>::iterator it2 = (*it)->my_tree->treeEnv.begin(); it2 != (*it)->my_tree->treeEnv.end(); it2++)
			delete (*it2);
		for (list<inMotif*>::iterator it2 = (*it)->motif_specs.begin(); it2 != (*it)->motif_specs.end(); it2++)
			delete (*it2);
		delete (*it)->my_tree;
		delete (*it);
	}
	delete options;

//	Leakage();
	
	return EXIT_SUCCESS;
}

//////////\\\\\\\\\\
//// FUNCTIONS \\\\\
//////////\\\\\\\\\\

//////////
/// This applies tree perturbations and treeScales to tree.
//////////
void fiddle_with_trees(list<inTree*>& inTrees, seqGenOptions *options) 
{
	double all_paths_sum, mean, std_dev;

	//////////
	/// Scale and perturb tree
	//////////
	//
	// If option is set, output tree stats:
	// * Average root-to-tip path lengths (+ std. dev).
	size_t x = 0;
	for (list<inTree*>::iterator it = inTrees.begin(); it != inTrees.end(); it++, x++) {
		if (options->treelength_check || options->treelength_scale || (*it)->scaleTree) {
			(*it)->calcDistancesFromRoot();
			all_paths_sum = mean = std_dev = 0.0;
			for (vector<TNode*>::iterator jt = (*it)->my_tree->tips.begin(); jt != (*it)->my_tree->tips.end(); jt++)
				all_paths_sum += (*jt)->DistanceFromRoot;
			mean = all_paths_sum / (*it)->my_tree->tips.size();
			for (vector<TNode*>::iterator jt = (*it)->my_tree->tips.begin(); jt != (*it)->my_tree->tips.end(); jt++)
				std_dev += ((*jt)->DistanceFromRoot-mean)*((*jt)->DistanceFromRoot-mean);
			std_dev /= (*it)->my_tree->tips.size();
			if ((*it)->scaleTree) {
				//////////
				/// Apply the branch scaling
				//////////
				//
				(*it)->apply_branchScales(-((*it)->branchScale/mean));
				//
				// Recalculate the distances from the root
				(*it)->calcDistancesFromRoot();
				//
				// Test if it is correct.
				all_paths_sum = 0.0;
				mean = std_dev = 0.0;
				for (vector<TNode*>::iterator jt = (*it)->my_tree->tips.begin(); jt != (*it)->my_tree->tips.end(); jt++)
					all_paths_sum += (*jt)->DistanceFromRoot;
				mean = all_paths_sum / (*it)->my_tree->tips.size();
				for (vector<TNode*>::iterator jt = (*it)->my_tree->tips.begin(); jt != (*it)->my_tree->tips.end(); jt++)
					std_dev += ((*jt)->DistanceFromRoot-mean)*((*jt)->DistanceFromRoot-mean);
				std_dev /= (*it)->my_tree->tips.size();
				cout << mean << " " << sqrt(std_dev) << endl;
			} else {
				cout << mean << " ... " << sqrt(std_dev) << endl;
				exit(EXIT_SUCCESS);
			}
		}
	}
	//
	// Set the perturbation values for the trees.
	//
	for (list<inTree*>::iterator it = inTrees.begin(); it != inTrees.end(); it++)
		(*it)->perturbTree(options->perturbTree);


	ofstream treescale_out;
	string treescale_outfile;
	treescale_outfile = options->output_files + ".scale_tree";
	treescale_out.open(treescale_outfile.c_str(), ios::trunc | ios::out);
	for (list<inTree*>::iterator it = inTrees.begin(); it != inTrees.end(); it++) 
		(*it)->Print_Tree(treescale_out, false);

	return;
}

//////////
/// Remove single-replicate specific objects, which will be recreated in next run
//////////
void Reset_Run(list<inTree*>& inTrees, seqGenOptions *options, bool reset)
{
	//////////
	/// This small chunk of code applies ONLY to random sequences using the PROSITE library
	/// (option -1). 
	//////////
	for (list<inTree*>::iterator it = inTrees.begin(); it != inTrees.end(); it++) {
		if ((*it)->rootSeqType == RANDOM && options->random_sequence_proportion_motif) {
			(*it)->motif_specs.clear();
			// Destroy inMotifs, will reset them for next run.
			for (list<inMotif*>::iterator jt = (*it)->my_tree->treeEnv.front()->my_motifs.begin(); 
										  jt != (*it)->my_tree->treeEnv.front()->my_motifs.end(); 
										  jt++) 
			{
				delete *jt;
			}
			(*it)->my_tree->treeEnv.front()->my_motifs.clear();
		}

		for (list<TNode*>::iterator it2 = (*it)->my_tree->nodeList.begin(); it2 != (*it)->my_tree->nodeList.end(); it2++) {
			(*it2)->evolvingSequence->evolutionaryAttributes.clear();
			delete (*it2)->evolvingSequence;
			// Variable region lists.
			if (reset) (*it2)->addGeneral_varSites();
		}
	}
	DestroyMatrix();
}

//////////
/// Read the trees in the partition file, parse them, get statistics of trees, and perform
/// necessary bookkeeping.
//////////
void readInputTrees(list<inTree*>& inputTrees, 	seqGenOptions *options, inClade *global_environment, 
					vector<string>& OTU_names, int *numTrees, int numTaxa)
{
	string treeInput;
	string tmpInput;
	bool foundOneTree;

	while (cin.good()) {
		tmpInput.clear();
		treeInput.clear();
		foundOneTree = false;
		do {
			tmpInput.clear();
			getline(cin, tmpInput);
			// Checking for errors in the tree file. All trees need [], so if more than one tree is found, this will catch it.
			if (tmpInput.find("[") != tmpInput.npos) {
				if (foundOneTree) {
					cerr << "Treefile error: Treefiles must contain the root sequence option \"[]\" and trees must end with a ';'." << endl  
					     << "Input tree: " << treeInput << endl << endl;
					exit(EXIT_FAILURE);
				}
				foundOneTree = true;
			}
			treeInput += trim(tmpInput);
		} while (treeInput.find(";") == treeInput.npos && cin.good());
		
		if(!treeInput.empty() && !cin.good()) {
			cerr << "Treefile error: Trailing characters after last tree. (Last tree may not end with ';')" << endl << treeInput << endl;
			exit(EXIT_FAILURE);
		}

		if(treeInput == "")
			continue;
		(*numTrees)++;
		inTree *tempTree = new inTree(options, *numTrees, global_environment);
		tempTree->perform_checks(treeInput, options);
		tempTree->parseTree(treeInput, options);
		if(options->lineageSpecificFile != "") tempTree->parseLineages(options);
		if(!(numTaxa)) {
			numTaxa = tempTree->my_tree->numTips;
			for (int i = 0; i < numTaxa; i++) 
				OTU_names.push_back(tempTree->my_tree->names.at(i));
		}

		tempTree->Define_Ancestors();
		tempTree->propagateTreePointers();
		tempTree->Define_Bipartitions(OTU_names);
		tempTree->SortTaxonTips(OTU_names);
		//////////
		/// If scaling branches, do now. If, however, scaling tree, this will be done
		/// later in the main loop.
		//////////
		if (!tempTree->scaleTree) tempTree->apply_branchScales(tempTree->branchScale);
		inputTrees.push_back(tempTree);
	}
}

//////////
/// This function checks to see if all lineages have marker matches, as well as other 
/// items of interest.
//////////
void CheckTrees(list<inTree*>& inputTrees, seqGenOptions *options)
{
	stringstream warning;
	bool set = false;
	// Now check each inMotif, find out if each marker defined is used or not. //
	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
		if ((*it)->rootSeqType != RANDOM) {
			for (list<inMotif*>::iterator jt = (*it)->motif_specs.begin(); jt != (*it)->motif_specs.end(); jt++) {
				if ( (*jt)->sitemap.empty() ) {
					warning << "Motif " << (*jt)->marker << " is not used in any of the input root sequences." << endl;
					set = true;
				}
			}
		}
	}
	if (set) options->SpoolWarnings(warning.str());
}

void Print_Motifs(list<inTree*> inputTrees, ostream& motif_out) 
{
	size_t curr_tree = 0;

	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
		motif_out << curr_tree << ":" << endl;
		for (list<inMotif*>::iterator jt = (*it)->motif_specs.begin(); jt != (*it)->motif_specs.end(); jt++) {
			motif_out << "Name: " << (*jt)->name << endl;
			motif_out << "Regular Expression: " << (*jt)->regex << endl;
			motif_out << "Placement on Root Sequence: " << (*jt)->sitemap << endl;
		}
	}
}

void Print_Trace(list<inTree*> inputTrees, list<eventTrack*> events, ostream& trace_out, int step_type, int treeNo, int setNo)
{

	trace_out << ">Dataset_" << setNo;

	if(step_type == TIME_RELATIVE_STEPS || step_type == GILLESPIE) {
		trace_out << endl;
		vector<eventTrack*> thisTip_event;
		thisTip_event.clear();
		for (list<eventTrack*>::iterator it = events.begin(); it != events.end(); it++) {
			thisTip_event.push_back((*it));
		}

		if (thisTip_event.size() > 0) {
			double min_time;
			int insert_sort_index;
			eventTrack *tmp;
			string out_event;
			for (size_t j = 0; j < thisTip_event.size()-1; j++) {
				insert_sort_index = j;
				min_time = thisTip_event.at(j)->eventTime;
				for (size_t k = j+1; k < thisTip_event.size(); k++) {
					if (thisTip_event.at(k)->eventTime < min_time) {
						min_time = thisTip_event.at(k)->eventTime;
						insert_sort_index = k;
					}
				}
				tmp = thisTip_event.at(j);
				thisTip_event.at(j) = thisTip_event.at(insert_sort_index);
				thisTip_event.at(insert_sort_index) = tmp;
			}

			for (vector<eventTrack*>::iterator it = thisTip_event.begin(); it != thisTip_event.end(); it++)
				trace_out << (*it)->Print_Event();

			thisTip_event.clear();
		}
	} else {	// DISCRETE_EVOLUTIONARY_STEPS
		// Locate the current tree:
		string out_event;
		inTree *curr_tree = NULL;
		if(events.size() > 0) {
			for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
				if ((*it)->treeNum == treeNo) curr_tree = (*it);
			}
	
			if(curr_tree == NULL) {
				cerr << "Print_Trace: curr_tree has null value." << endl;
				exit(EXIT_FAILURE);
			}

			trace_out << "__partition_" << treeNo << endl;
			for (list<eventTrack*>::iterator it = events.begin(); it != events.end(); it++) {
				out_event = (*it)->Print_Event();
				trace_out << out_event;
			}
		}  else {
			trace_out << "__partition_" << treeNo << endl << "-" << endl;
		}
	}
	
	// Free all of the event tracking.
	for (list<eventTrack*>::iterator it = events.begin(); it != events.end(); it++) {
		delete (*it);
	}
}

void Print_Root(list<inTree*> inputTrees, ostream& root_out, seqGenOptions *options) 
{
	if (options->output_file_flags[verb]) verbose_output << "Root: ";

	for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
		for(int i = 0; i < (*it)->my_tree->root_numSites; i++) {
			root_out << stateCharacters[(*it)->my_tree->root->seq_evo.at(i).returnState()];
			if (options->output_file_flags[verb]) verbose_output << stateCharacters[(*it)->my_tree->root->seq_evo.at(i).returnState()];
		}
	}
	root_out << endl;
	if (options->output_file_flags[verb]) verbose_output << endl << endl;
}

void Print_NEXUS_Header(list<inTree*> inputTrees, seqGenOptions *options, ostream& out) 
{
	int root_length = 0;
	int partitions = 0;
	bool gene_flag = false;

	for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
		root_length += (*it)->partitionLength;
		if((*it)->my_tree->treeEnv.front()->rateHetero == CodonRates) gene_flag = true;
		partitions++;
	}

	out << "#NEXUS" << endl << "[" << endl;
	out << "Generated by " << PROGRAM_NAME << " " << VERSION_NUMBER << endl << endl;
	out << "Simulation of " << inputTrees.front()->my_tree->numTips << " sequences of " << root_length;
	out << " " << ((isNucModel) ? "nucleotides" : "amino acids") << endl;
	out << "   for " << options->number_of_datasets_per_tree << " datasets and " << partitions;
	out << " partitions per dataset" << endl << endl;

	// Rates etc?
		
	out << "Partition-specific parameters: " << endl;
	partitions = 1;
	for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
		out << "(Partition " << partitions << ") " << (*it)->label << endl;
		out << " -Tree = " << (*it)->tree << endl;
		if((*it)->scaleBranches) {
			out << " -Branch lengths of tree multiplied by factor of " << (*it)->branchScale << endl;
		} else {
			out << " -Branch lengths assumed to be the number of substitutions per site." << endl;
		}
		out << " -Invariable sites model = ";
		if((*it)->my_tree->treeEnv.front()->invariableSites) {
			if((*it)->randomInvariableAssignment) {
				out << (*it)->proportion_invariable * 100 << "% invariable." << endl;
			} else {
				out << "User-supplied root sequence" << endl;
				out << "\tConserved sites: ";
				bool no_indel_site = false;
				string no_indel_run;
				no_indel_run.clear();

				for (int i = 0; i < (*it)->partitionLength; i++) {
					if (atoi((*it)->invariable_by_partition.substr(i,1).c_str()) == 1) {
						if(no_indel_site) {
							out << no_indel_run << " ";
							no_indel_run.clear();
						}
						out << (*it)->rootseq_by_partition.at(i) << " ";
						no_indel_site = false;
					} else if (atoi((*it)->invariable_by_partition.substr(i,1).c_str()) == 2) {
						no_indel_run.push_back('X');
						no_indel_site = true;
					} else if (atoi((*it)->invariable_by_partition.substr(i,1).c_str()) == 3) {
						no_indel_run.push_back((*it)->rootseq_by_partition.at(i));
						no_indel_site = true;
					} else {
						if(no_indel_site) {
							out << no_indel_run << " ";
							no_indel_run.clear();
						}
						no_indel_site = false;
					}
				}
				if(no_indel_site) out << no_indel_run;
				out << endl;
			}
		} else {
			out << "No invariable sites" << endl;
		}

		if(gene_flag) {
			out << " -Gene sequence: ";
			if((*it)->my_tree->treeEnv.front()->rateHetero == CodonRates) out << "Exon" << endl;
			else out << "Intron" << endl;
		}

		out << " -Model = " << modelNames[options->global_model] << endl;
		out << " -Insertion/Deletion model: " << endl;
		if((*it)->my_tree->treeEnv.front()->indelFlag) {
			out << "\tMaximum indel size = " << (*it)->my_tree->treeEnv.front()->maxIndel << endl;
			if((*it)->my_tree->treeEnv.front()->P_ins_ == 0 && (*it)->my_tree->treeEnv.front()->P_del_ == 0 && (*it)->my_tree->treeEnv.front()->indelFlag == true) {
				out << "\tChang, M.S.S, and Benner, S. (2004). J. Mol. Biol. 341:617-631." << endl;
			} else {
				if((*it)->my_tree->treeEnv.front()->P_ins_) {
					out << "\tP(ins) = " << (*it)->my_tree->treeEnv.front()->P_ins_ << "\t";
					for (int i = 1; i <= (*it)->my_tree->treeEnv.front()->maxIndel; i++) 
						out << i << ":" << (*it)->my_tree->treeEnv.front()->insert_lengthDistribution.at(i) << " ";
					out << endl;
				}
				if((*it)->my_tree->treeEnv.front()->P_del_) {
					out << "\tP(del) = " << (*it)->my_tree->treeEnv.front()->P_del_ << "\t";
					for (int i = 1; i <= (*it)->my_tree->treeEnv.front()->maxIndel; i++) 
						out << i << ":" << (*it)->my_tree->treeEnv.front()->delete_lengthDistribution.at(i) << " ";
					out << endl;
				}
			}
		} else out << "None." << endl;

		out << " -Character Frequencies:";
		for(int i = 0; i < numStates; i++) {
			if(i % 10 == 0) { out << "\n\t"; }
			out << stateCharacters[i] << ":";
			if(isNucModel)
				out << setprecision(3) << setw(6) << nucFreq[i] << " ";
			else 
				out << setprecision(3) << setw(6) << aaFreq[i] << " ";
		}
		out << endl;
		if((*it)->my_tree->treeEnv.front()->rateHetero == CodonRates) {
			out << " -Codon Relative Rates: " << endl;
			out << "\t*Pos. 1: " << (*it)->my_tree->treeEnv.front()->catRate[0];
			out << "\t*Pos. 2: " << (*it)->my_tree->treeEnv.front()->catRate[1];
			out << "\t*Pos. 3: " << (*it)->my_tree->treeEnv.front()->catRate[2];
		}
		out << endl << endl;
		
		partitions++;
	}
	out << "]" << endl << endl;
}

void Print_Seq(
			   list<inTree*> inputTrees, 
			   ostream& seq_out, 
			   seqGenOptions *options, 
			   int dataset_num
			  ) 
{   
	int	numTips = inputTrees.front()->my_tree->numTips;
	if(options->writeAncestors) {
		if(inputTrees.front()->my_tree->rooted) 
			numTips = (2 * (inputTrees.front()->my_tree->numTips + 1)) - 3;
		else
			numTips = (2 * (inputTrees.front()->my_tree->numTips + 1)) - 4;
	} else {
		numTips = inputTrees.front()->my_tree->numTips;
	}
	int n;
	int print_array_row;
	vector<char> v;
	vector<vector<char> > print(numTips,v);
	char **anc_names = NULL;;
	if(options->writeAncestors) {
		anc_names = new char* [numTips+1];
		for(int x = 0; x < numTips+1; x++) {
			anc_names[x] = new char [MAX_NAME_LEN];
			for(int y = 0; y < MAX_NAME_LEN; y++) {
				anc_names[x][y] = '\0';
			}
		}
	}

	if (options->output_file_flags[verb]) verbose_output << "Sequences:" << endl;

	if(options->writeAncestors) {
		for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
			print_array_row = 0;
			n = inputTrees.front()->my_tree->numTips + 1;
			sprintf(anc_names[0], "%d", n);
			for (int j = 0; j<(*it)->my_tree->root_numSites; j++)
				print.at(0).push_back(stateCharacters[(*it)->my_tree->root->seq_evo.at(j).returnState()]);
            print_array_row++;
    		if (!(*it)->my_tree->rooted) 
    			WriteAncestralSequencesNode(print, &print_array_row, (*it)->my_tree, &n, (*it)->my_tree->root->branch0, anc_names);
    		WriteAncestralSequencesNode(print, &print_array_row, (*it)->my_tree, &n, (*it)->my_tree->root->branch1, anc_names);
    		WriteAncestralSequencesNode(print, &print_array_row, (*it)->my_tree, &n, (*it)->my_tree->root->branch2, anc_names);
		}

		for(size_t i = 0; i < print.size(); i++) {
			seq_out << ">" << anc_names[i] << endl;

			if (options->output_file_flags[verb]) {
				int j;
				for(j = 0; j < 9 && anc_names[i][j]; j++)
					verbose_output << anc_names[i][j];
				while(j++ < 10) verbose_output << " ";
			}
			
			for(size_t j = 0; j < print.at(i).size(); j++) {
				seq_out << print.at(i).at(j);
				// If the second condition is true, there is a double-endline in the output, which is bad.
				if((j+1) % options->output_width == 0 && (j+1) != print.at(i).size()) 
					seq_out << endl;
				if (options->output_file_flags[verb]) verbose_output << print.at(i).at(j);
			}
			seq_out << endl;
			if (options->output_file_flags[verb]) verbose_output << endl;
		}
		seq_out << endl;
		if (options->output_file_flags[verb]) verbose_output << endl;
	} else {
		// tips were correctly sorted in the beginning, so using the front tree only is OK.
		for (int i = 0; i < numTips; i++) {
			seq_out << ">" << inputTrees.front()->my_tree->names.at(i) << endl;
			if (options->output_file_flags[verb]) {
				int j;
				for(j = 0; j < 9 && j < inputTrees.front()->my_tree->names.at(i).size(); j++) {
					verbose_output << inputTrees.front()->my_tree->names.at(i).at(j);
				}
				while(j++ < 10) verbose_output << " ";
			}

			int aggregate = 0;
			for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
				for(int j = 0; j < (*it)->my_tree->tips.at(i)->seq_evo.size(); aggregate++, j++) {
					seq_out << stateCharacters[(*it)->my_tree->tips.at(i)->seq_evo.at(j).returnState()];
					if(options->fileFormat == FASTAFormat) {
						if((aggregate+1) % options->output_width == 0 && (aggregate+1) != (*it)->my_tree->tips.at(i)->seq_evo.size()) 
							seq_out << endl;
					}
					if (options->output_file_flags[verb]) verbose_output << stateCharacters[(*it)->my_tree->tips.at(i)->seq_evo.at(j).returnState()];
				}
			}
			seq_out << endl;
			if (options->output_file_flags[verb]) verbose_output << endl;
		}
		seq_out << endl;
		if (options->output_file_flags[verb]) verbose_output << endl;
	}
}      

void WriteAncestralSequencesNode(vector<vector<char> >& print, int *print_array_row, TTree *tree, 
								 int *nodeNo, TNode *des, char **anc_names)
{ 
    int j;       
     
    if (des->tipNo==-1) { 
        (*nodeNo)++; 
		sprintf(anc_names[*print_array_row], "%d", *nodeNo);
	} else {
		strcpy(anc_names[*print_array_row],(tree->names.at(des->tipNo)).c_str());
	}
                         
    for (j=0; j<des->seq_evo.size(); j++) { 
		print.at(*print_array_row).push_back(stateCharacters[des->seq_evo.at(j).returnState()]);
    } 

	(*print_array_row)++;
	if(des->tipNo==-1) {         
        WriteAncestralSequencesNode(print, print_array_row, tree, nodeNo, des->branch1, anc_names);
        WriteAncestralSequencesNode(print, print_array_row, tree, nodeNo, des->branch2, anc_names);
    }   
} 

void Print_MA(list<inTree*> inputTrees, ostream& ma_out, seqGenOptions *options, int dataset_num,
			  list<eventTrack*> events) 
{
	int numTips;
	
	if(options->writeAncestors) {
		if(inputTrees.front()->my_tree->rooted) 
			numTips = (2 * (inputTrees.front()->my_tree->numTips + 1)) - 3;
		else
			numTips = (2 * (inputTrees.front()->my_tree->numTips + 1)) - 4;
	} else {
		numTips = inputTrees.front()->my_tree->numTips;
	}
	int i;	// Goes with the global arrays
	bool gene_flag = false;
	vector<char> v;
	vector<vector<char> > print(numTips,v);
	vector<int> invar;
	string label_print;
	label_print.clear();
	int n = inputTrees.front()->my_tree->numTips + 1;
	int print_array_row, taxa_ptr;
	char **anc_names = NULL;
	if(options->writeAncestors) {
		anc_names = new char* [numTips+1];
		for(int x = 0; x < numTips+1; x++) {
			anc_names[x] = new char [MAX_NAME_LEN];
			for(int y = 0; y < MAX_NAME_LEN; y++) {
				anc_names[x][y] = '\0';
			}
		}
	}

	int last_iter_size = 0;
	size_t this_iter_size;
	char overflow_label = 'A';
	list<string> overflows;
	for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
		if((*it)->my_tree->treeEnv.front()->rateHetero == CodonRates) gene_flag = true;
		if(options->writeAncestors) {
			print_array_row = i = taxa_ptr = 0;
			n = inputTrees.front()->my_tree->numTips+1;
			sprintf(anc_names[0], "%d", n);
			if((*it)->my_tree->global_arrays.at(0).action == 'd') { i=0; } else { i=1; }
			while((*it)->my_tree->global_arrays.at(i).action != 'i' && i < (*it)->my_tree->global_arrays.size()) 
				i++;			

			while(i < (*it)->my_tree->global_arrays.size()) {
				if((*it)->my_tree->global_arrays.at(i).fromAnc == -1) {
					label_print += " ";
					if ((!(*it)->randomInvariableAssignment && (*it)->rootSeqType == RANDOM) ||
						(!(*it)->my_tree->treeEnv.front()->invariableSites)) {
						invar.push_back(0);
					} else {
						invar.push_back((*it)->my_tree->root->seq_evo.at(taxa_ptr).returnInvariableState());
					}
					print.at(print_array_row).push_back(stateCharacters[(*it)->my_tree->root->seq_evo.at(taxa_ptr).returnState()]);
					taxa_ptr++;
				} else {
					if(options->output_indel_gaps) print.at(print_array_row).push_back('+');
					else print.at(print_array_row).push_back('-');
					label_print += " ";
					invar.push_back(0);
				}
				i++;
				if (i < (*it)->my_tree->global_arrays.size())
					while((*it)->my_tree->global_arrays.at(i).action != 'i' && i < (*it)->my_tree->global_arrays.size())
						i++;
			}

			label_print.at(last_iter_size) = '|';
			this_iter_size = label_print.size() - last_iter_size;
			if((*it)->label.size() > (this_iter_size-2)) {
				label_print.at(last_iter_size+this_iter_size/2-1) = ':';
				label_print.at(last_iter_size+this_iter_size/2) = overflow_label;
				overflow_label++;
				overflows.push_back((*it)->label);
			} else {
				int placement = last_iter_size+(this_iter_size/2-(*it)->label.size()/2)+1;
				for (size_t x = 0; x < (*it)->label.size(); x++, placement++)
					label_print.at(placement) = (*it)->label.at(x);
			}
			last_iter_size = label_print.size();
			label_print.at(label_print.size()-1) = '|';

			print_array_row++;
			if(!(*it)->my_tree->rooted)
				WriteAncestralMA(print, &print_array_row, (*it)->my_tree->root->branch0, &n, (*it)->my_tree, anc_names, options);
			WriteAncestralMA(print, &print_array_row, (*it)->my_tree->root->branch1, &n, (*it)->my_tree, anc_names, options);
			WriteAncestralMA(print, &print_array_row, (*it)->my_tree->root->branch2, &n, (*it)->my_tree, anc_names, options);
		} else {
			vector<int> taxa_ptrs(numTips);
			for(int x = 0; x < numTips; x++) { taxa_ptrs.at(x) = 0; }
			i=0;
			if((*it)->my_tree->global_arrays.at(0).action == 'd') i=0; else i=1;
			bool push_back, label_push_back;
			while(i < (*it)->my_tree->global_arrays.size()) {
				push_back = label_push_back = false;
				while((*it)->my_tree->global_arrays.at(i).action != 'i' && i < (*it)->my_tree->global_arrays.size()) 
					i++;
				for(int j = 0; j < numTips; j++) {
					if(!isAnc((*it)->my_tree->tips.at(j),(*it)->my_tree->global_arrays.at(i).fromAnc)) {
						if(options->output_indel_gaps) print.at(j).push_back('+');
						else print.at(j).push_back('-');
						if(!label_push_back) { label_print += " "; label_push_back = true; }
					} else {	
						// If it is an ancestor, and prev is insert, then this is an insert
						if((*it)->my_tree->global_arrays.at(i-1).action == 'i') {
							print.at(j).push_back(stateCharacters[(*it)->my_tree->tips.at(j)->seq_evo.at(taxa_ptrs.at(j)).returnState()]);
							if ((*it)->my_tree->treeEnv.front()->invariableSites && !push_back &&
								(*it)->my_tree->tips.at(j)->nodeEnv->invariableSites) {
								invar.push_back((*it)->my_tree->tips.at(j)->seq_evo.at(taxa_ptrs.at(j)).returnInvariableState());
								push_back = true;
							}
							if(!label_push_back) { label_print += " "; label_push_back = true; }
							taxa_ptrs.at(j)++;
						} else {
							// Take care of deletions...
							int keep = i;
							i--;
							while ((*it)->my_tree->global_arrays.at(i).action != 'i' && i >=0)
								i--;
							i++;
							if (keep == i) {
								if ((*it)->my_tree->global_arrays.at(i).fromAnc == (*it)->my_tree->tips.at(j)->tipNo) {
									if (!label_push_back) { label_print += " "; label_push_back = true; }
									print.at(j).push_back(stateCharacters[(*it)->my_tree->tips.at(j)->seq_evo.at(taxa_ptrs.at(j)).returnState()]);
									if ((*it)->my_tree->treeEnv.front()->invariableSites && !push_back &&
										(*it)->my_tree->tips.at(j)->nodeEnv->invariableSites) {
										invar.push_back((*it)->my_tree->tips.at(j)->seq_evo.at(taxa_ptrs.at(j)).returnInvariableState());
										push_back = true;
									}
									taxa_ptrs.at(j)++;
								} else {
									if (!label_push_back) { label_print += " "; label_push_back = true; }
									print.at(j).push_back('-');
								}
							} else {
								int deletion = 0;
								while ((*it)->my_tree->global_arrays.at(i).action != 'i' && i < (*it)->my_tree->global_arrays.size()) {
									if (isAnc((*it)->my_tree->tips.at(j),(*it)->my_tree->global_arrays.at(i).fromAnc)) {
										deletion = 1;
									}
									i++;
								}
								if (deletion) {
									if (!label_push_back) { label_print += " "; label_push_back = true; }
									print.at(j).push_back('-');
								} else {
									print.at(j).push_back(stateCharacters[(*it)->my_tree->tips.at(j)->seq_evo.at(taxa_ptrs.at(j)).returnState()]);
									if (!label_push_back) { label_print += " "; label_push_back = true; }
									if ((*it)->my_tree->treeEnv.front()->invariableSites && !push_back &&
										(*it)->my_tree->tips.at(j)->nodeEnv->invariableSites) {
										invar.push_back((*it)->my_tree->tips.at(j)->seq_evo.at(taxa_ptrs.at(j)).returnInvariableState());
										push_back = true;
									}
									taxa_ptrs.at(j)++;
								}
							}
						}
					}
				}
				if (!push_back) invar.push_back(0);
				i++;
			}
			label_print.at(last_iter_size) = '|';
			this_iter_size = label_print.size() - last_iter_size;
			if((*it)->label.size() > (this_iter_size-2)) {
				label_print.at(last_iter_size+this_iter_size/2-1) = ':';
				label_print.at(last_iter_size+this_iter_size/2) = overflow_label;
				overflow_label++;
				overflows.push_back((*it)->label);
			} else {
				int placement = last_iter_size+(this_iter_size/2-(*it)->label.size()/2)+1;
				for (size_t x = 0; x < (*it)->label.size(); x++, placement++)
					label_print.at(placement) = (*it)->label.at(x);
			}
			last_iter_size = label_print.size();
			label_print.at(label_print.size()-1) = '|';
		}
	}

	//////////
	/// Detect empty columns
	//////////
	vector<bool> empty_columns(print.front().size(), false);
	int num_empty_col = 0;
	num_empty_col = trimMSA(print, empty_columns, options);
	size_t maxNameLen = 0;
	int MA_length = print.at(0).size() - num_empty_col;
	
	//////////
	/// MSA header
	//////////
	if(options->fileFormat == NEXUSFormat) {
		size_t len;
		if(options->writeAncestors) {
			for(int i = 0; i < numTips; i++) {
				len = strlen(anc_names[i]);
				if(len > maxNameLen) maxNameLen = strlen(anc_names[i]);
			}
		} else {
			maxNameLen = (inputTrees.front()->my_tree->names.at(0)).length();
			for (int i = 1; i < numTips; i++) {
				len = (inputTrees.front()->my_tree->names.at(i)).length();
				if (len > maxNameLen) maxNameLen = len;
			}
		}
		ma_out << "Begin DATA;\t[Dataset " << dataset_num << "]" << endl;
		ma_out << "\tDimensions NTAX=" << numTips << " NCHAR=" << MA_length << ";" << endl;
		ma_out << "\tFormat MISSING=? GAP=- DATATYPE=" << ((isNucModel)?((gene_flag)?"CODON":"DNA"):"PROTEIN") << ";" << endl;
		ma_out << "\tMatrix" << endl;
	} else if (options->fileFormat == PHYLIPFormat) {
		maxNameLen = 9;
		ma_out << " " << numTips << " " << MA_length << endl;
	}

	//////////
	/// MSA output
	//////////
	if(options->fileFormat == PHYLIPFormat || options->fileFormat == NEXUSFormat) {
		int k = 0;
		vector<char>::size_type m = 0;		
		while(m < print.at(0).size()) {	
			for(size_t i = 0; i < print.size(); i++) {
				int j=0;
				//////////
				/// Phylip only prints labels on first block. Rest of blocks omit labels.
				//////////
				if (m == 0 || options->fileFormat == NEXUSFormat) {
					if(options->writeAncestors) 
						for(j = 0; j < maxNameLen && anc_names[i][j]; j++)
							ma_out << anc_names[i][j];
					else
						for(j = 0; j < maxNameLen && j < inputTrees.front()->my_tree->names.at(i).size(); j++) 
							ma_out << inputTrees.front()->my_tree->names.at(i).at(j);
				}
				while(j++ <= maxNameLen) ma_out << " ";			

				//////////
				/// Sequence output.
				//////////
				for(k = 0; (m + k) < print.at(i).size() && k < options->output_width; k++)
					if ( !empty_columns.at(m+k) )
						ma_out << print.at(i).at(m+k);
				ma_out << endl;
			}
			m += options->output_width;
			ma_out << endl;
		}
	} else if (options->fileFormat == FASTAFormat) {
		for(size_t i = 0; i < print.size(); i++) {
			if(options->writeAncestors) ma_out << ">" << anc_names[i] << endl;
			else ma_out << ">" << inputTrees.front()->my_tree->names.at(i) << endl;

			//////////
			/// Sequence output
			//////////
			size_t num_empty = 0;
			for(size_t j = 0; j < print.at(i).size(); j++) {
				if (empty_columns.at(j)) num_empty++;
				else ma_out << print.at(i).at(j);
				if((j+1-num_empty) % options->output_width == 0 && (j+1) != print.at(i).size()) 
					ma_out << endl;
			}
			ma_out << endl;
		}
	} else { cerr << "What format is this??? " << options->fileFormat << endl; exit(EXIT_FAILURE); }

	//////////
	/// VERBOSE-SPECIFIC OUTPUT
	//////////
	if (options->output_file_flags[verb]) {
		int k = 0;
		size_t m = 0;
		verbose_output << "Multiple Alignment:" << endl;
		while(m < print.at(0).size()) {
			//////////
			/// Labels
			//////////
			for(int j = 0; j < 10; j++) verbose_output << " ";
			for(k = 0; (m + k) < print.at(0).size() && k < options->output_width; k++)
				verbose_output << label_print.at(m+k);
			verbose_output << endl;

			//////////
			/// Invariable
			//////////
			for(int j = 0; j < 10; j++) verbose_output << " ";
			for(k = 0; (m + k) < print.at(0).size() && k < options->output_width; k++) 
				verbose_output << invar.at(m+k);
			verbose_output << endl;

			for(size_t i = 0; i < print.size(); i++) {
				int j;
				if(options->writeAncestors)
					for(j = 0; j < 9 && anc_names[i][j]; j++)
						verbose_output << anc_names[i][j];
				else
					for(j = 0; j < 9 && j < inputTrees.front()->my_tree->names.at(i).size(); j++) 
						verbose_output << inputTrees.front()->my_tree->names.at(i).at(j);
				while(j++ < 10) verbose_output << " ";
				for(k = 0; (m + k) < print.at(i).size() && k < options->output_width; k++) 
					verbose_output << print.at(i).at(m+k);
				verbose_output << endl;
			}
			m += options->output_width;
			verbose_output << endl;
		}
	
		if(overflows.size() > 0) {
			verbose_output << "Labels for Short Sequences:" << endl;
			char ch = 'A';
			for(list<string>::iterator str = overflows.begin(); str != overflows.end(); str++) {
				verbose_output << ":" << ch << ".\t" << *str << endl;
			}
		}
		verbose_output << endl;
	}

	//////////
	/// MSA format footer data.
	//////////
	if(options->fileFormat == NEXUSFormat) ma_out << "\t;" << endl << "END;" << endl << endl;
	if(options->fileFormat == FASTAFormat) ma_out << endl;
}

void WriteAncestralMA(vector<vector<char> >& print, int *print_array_row, TNode *des, int *nodeNo, 
					  TTree *tree, char** anc_names, seqGenOptions *options)
{
	int i, taxa_ptr, keep, found;
	i = taxa_ptr = 0;
	
	if(des->tipNo == -1) {
		(*nodeNo)++;
		sprintf(anc_names[*print_array_row], "%d", *nodeNo);
	} else {
		strcpy(anc_names[*print_array_row], (tree->names.at(des->tipNo).c_str()));
	}

	if(tree->global_arrays.at(0).action == 'd') { i=0; } else { i=1; }
    while(tree->global_arrays.at(i).action != 'i' && i < tree->global_arrays.size()) i++;
	
    while(i < tree->global_arrays.size()) {
		if(!isAnc(des,tree->global_arrays.at(i).fromAnc)) {
			if(options->output_indel_gaps) print.at(*print_array_row).push_back('+');
	    	else print.at(*print_array_row).push_back('-');
		} else {
	    	keep = i;
	    	while(tree->global_arrays.at(--i).action == 'd' && i > 0) ;
	    	i++;
	    	if(keep == i) {
	    		print.at(*print_array_row).push_back(stateCharacters[des->seq_evo.at(taxa_ptr).returnState()]);
				taxa_ptr++;
	    	} else {
				found = 0;
				while(i < keep && !found) {
		    		if(isAnc(des,tree->global_arrays.at(i).fromAnc)) {
	    				print.at(*print_array_row).push_back('-');
						found = 1;
		    		} else {
		    			i++;
		    		}
				}
				i=keep;
				if(!found) {
		    		print.at(*print_array_row).push_back(stateCharacters[des->seq_evo.at(taxa_ptr).returnState()]);
		    		taxa_ptr++;
				}
	    	}
		}

		i++;
		if (i < tree->global_arrays.size())
			while(tree->global_arrays.at(i).action != 'i' && i <= tree->global_arrays.size()) 
				i++;
    }

    (*print_array_row)++;	// Go to the next row.

    if(des->tipNo==-1) {
		WriteAncestralMA(print, print_array_row, des->branch1, nodeNo, tree, anc_names, options);
		WriteAncestralMA(print, print_array_row, des->branch2, nodeNo, tree, anc_names, options);
    }
}

//////////
/// This function will detect the columns that are empty, and set those columns to true for use
/// in Print_MSA. If output_file_flags[trace] is unset, this function will not set empty_columns
/// to true in any place, since each site is important for the trace file. 
//////////
size_t trimMSA(vector<vector<char> >& print, vector<bool>& empty_columns, seqGenOptions *options)
{
	size_t num_columns = print.front().size();	// Number of sequences
	size_t num_rows = print.size(); // Number of columns.
	size_t empty_column_detected = 0;
	
	for (size_t i = 0; i < num_columns; i++) {
		bool empty = true;
		for (size_t j = 0; j < num_rows; j++) {
			if (print.at(j).at(i) != '-') empty = false;
		}
		if (empty) {
			empty_columns.at(i) = true;
			empty_column_detected++;
		}
	}

	if (empty_column_detected && !empty_column_warning_spooled) {
		stringstream warning;
		warning << "Output MSAs contain columns consisting of all gaps";
		if (options->output_file_flags[trace]) {
			warning << ".";
			empty_columns.assign(empty_columns.size(), false);
			empty_column_detected = 0;
		} else warning << ", which have been removed from the .ma output file. The .verb file, "
					   << "if created, contains the MSA with all-gap columns present.";
		options->SpoolWarnings(warning.str());
		empty_column_warning_spooled = true;
	}
	
	return empty_column_detected;
}

void PrintTitle()
{
	cerr << "Sequence Generator - " << PROGRAM_NAME << VERSION_NUMBER << endl;
	cerr << "Substitution Engine modified from Seq-Gen v1.3.2:" << endl;
	cerr << "(c) 2005--present Cory Strope" << endl;
}

//////////
/// Fast and dirty random seed algorithm; based off of NR's randqd1
//////////
inline unsigned int rand_seed()
{
	static unsigned int u = (unsigned int)(time(NULL)+3*getpid());
	return (u = u*1664525u + 1013904223u);
}

void setTRS(list<inTree*>& inputTrees, seqGenOptions *options, double max_path_length, double min_branch) 
{
	int num_discrete_steps = 1024;
	while (( (min_branch > 0.01) ? 0.01 : min_branch) < (max_path_length / num_discrete_steps)) 
		num_discrete_steps *= 2;
	double step_size = max_path_length / num_discrete_steps;
	if (step_size < 0.00001) step_size = 0.00001;	// To prevent underflow. For subst, 0=10e-06 //

	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
		if((*it)->my_tree->rooted) {
			(*it)->TimeScale_Tree(max_path_length, (*it)->my_tree->root);
		} else {
			options->SpoolWarnings("Cannot scale unrooted trees for epoch timepoints.");
		}
		(*it)->global_max_path_length = max_path_length;
	}

	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); it++) {
		(*it)->NudgeBranches(step_size,options->simulation_step_type);
		if (options->simulation_step_type == TIME_RELATIVE_STEPS)
			(*it)->my_tree->epoch_step_size = 1.0 / num_discrete_steps;
		else if (options->simulation_step_type == DISCRETE_EVOLUTIONARY_STEPS) {
			(*it)->my_tree->epoch_step_size = step_size;
			(*it)->my_tree->trs_example_step_size = 1.0 / num_discrete_steps;
		} else {		// GILLESPIE //
			//////////
			/// With Gillespie, epoch_step_size should be set as multiplier to dt?
			//////////
			(*it)->my_tree->epoch_step_size = 0.0;
		}
	}
}

void Leakage()
{
	long memoryLeaked = 0;

	//	cout << "@" << endl;
	//	cout << correct_motif_positions << " " << total_motif_positions << endl;
	//	cout << (double)correct_motif_positions/(double)total_motif_positions << endl;
	//	cout << num_template_violations << endl;
	//	cout << "*" << endl;

	memoryLeaked = 0
				 + print_memory("<TNode>:", 		   TNode::howMany(), 			sizeof(TNode))
				 + print_memory("<inClade>:", 		   inClade::howMany(), 			sizeof(inClade))
				 + print_memory("<inTree>:", 		   inTree::howMany(), 			sizeof(inTree))
				 + print_memory("<TTree>:", 		   TTree::howMany(), 			sizeof(TTree))
				 + print_memory("<Branch>:", 		   Branch::howMany(), 			sizeof(Branch))
				 + print_memory("<inMotif>:", 		   inMotif::howMany(), 			sizeof(inMotif))
				 + print_memory("<siteRegEx>:", 	   siteRegEx::howMany(), 		sizeof(siteRegEx))
				 + print_memory("<Sequence>:", 		   Sequence::howMany(), 		sizeof(Sequence))
				 + print_memory("<Site>:", 			   Site::howMany(), 			sizeof(Site))
				 + print_memory("<motifSite>:", 	   motifSite::howMany(), 		sizeof(motifSite))
				 + print_memory("<siteProperties>:",   siteProperties::howMany(),   sizeof(siteProperties))
				 + print_memory("<activeProperties>:", activeProperties::howMany(), sizeof(siteProperties))
				 + print_memory("<seqGenOptions>:",    seqGenOptions::howMany(), 	sizeof(seqGenOptions))
				 + print_memory("<Indel>:", 		   Indel::howMany(), 			sizeof(Indel))
				 + print_memory("<Substitution>:", 	   Substitution::howMany(), 	sizeof(Substitution))
				 + print_memory("<Insertion>:", 	   Insertion::howMany(), 		sizeof(Insertion))
				 + print_memory("<Deletion>:", 		   Deletion::howMany(), 		sizeof(Deletion))
				 + print_memory("<varSite>:", 		   varSite::howMany(), 			sizeof(varSite))
				 + print_memory("<eventTrack>:", 	   eventTrack::howMany(), 		sizeof(eventTrack))
				 ;

	cout << "Total memory leaked = " << memoryLeaked << " bytes." << endl << endl;
	cout << "Number of inserted positions: " << num_insert << " insertions for " << num_inserted_positions << endl;
	cout << "Number of deleted positions:  " << num_delete << " deletions for  " << num_deleted_positions << endl;
}

long print_memory(string message, size_t howMany, int object_size)
{
	cout << setw(20) << message << setw(10) << howMany << " x " << setw(5) << object_size
		 << " = " << setw(10) << howMany*object_size << endl;
	
	return howMany*object_size;
}

