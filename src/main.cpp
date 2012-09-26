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
#include "twister.h"
#include "random.h"
#include "trace.h"
#include "paleo.h"
#include "dependency.h"
#include "propose_path.h"
#include "stats.h"
#include "rate_type.h"

using namespace std;

bool	order_3_markov;
bool	Human_Data_simulation;

bool	empty_column_warning_spooled = false;
ofstream verbose_output;
extern 	ofstream verbose_output;
int 	num_copy_ctor = 0, num_dtor = 0, num_ctor = 0;
int 	num_insert = 0, num_delete = 0, len_insert = 0, len_delete = 0;
int 	actual_events = 0, virtual_events = 0;
int		eventNo; 
const string nucleotides="ACGT";
const string aminoAcids="ARNDCQEGHILKMFPSTWYV";
bool 	test_alternate_representations = true;
unsigned int		numStates_squared;
unsigned int		numStates_cubed;
int 	changed_site=0;
int		prev_state = 0;
bool	fast_simulation = true;
bool	optimize = true;
bool	rasmus_independent_proposals;
bool 	MCMC_sample_evenly = false;
bool	Qd=true, Pc=true, nij=false;
int 	print_stepwise_rates = 0;
bool	forward_simulation = false;
RateMatrix (*ptr2init)(TTree*, TNode*, TNode*, double, double, int);
void (*ptr2update)(RateMatrix*, TNode*, vector<Site>::iterator, double, double);


// MAIN prototypes
void setRates(list<inTree*>& inputTrees,seqGenOptions *options);
void Forward_Simulation(list<inTree*>& inputTrees, inClade *global_environment, int replicate, vector<ofstream*>& simulation_output_streams, seqGenOptions *options, list<eventTrack*> *events);
void Simulate(list<inTree*>& inputTrees, inClade *global_environment, vector<ofstream*>& simulation_output_streams, seqGenOptions *options);
void Path_Proposal(list<inTree*>& inputTrees, Statistics *stats, inClade *global_environment, vector<ofstream*>& simulation_output_streams, seqGenOptions *options, list<eventTrack*> *events);
void closeOutputStreams(list<inTree*>& inputTrees, vector<ofstream*>& simulation_output_streams, vector<bool> outfile_flags);
void openOutputStreams(vector<ofstream*>& simulation_output_streams, vector<bool> outfile_flags, string& outfile_name_root, bool simulation_type);
void Print_Root(list<inTree*> inputTrees, ostream& root_out, seqGenOptions *options);
void Print_Trace(list<inTree*> inputTrees, list<eventTrack*> events, ostream& trace_out, int step_type, int treeNo, int setNo);
void Print_Seq(list<inTree*> inputTrees, ostream& seq_out, seqGenOptions *options, int dataset_num);
void Print_MSA(list<inTree*> inputTrees, ostream& ma_out, seqGenOptions *options, int dataset_num, list<eventTrack*> events);
void Print_Motifs(list<inTree*>& inputTrees, ostream& motif_out);
void WriteAncestralMSA(vector<vector<char> >&, int *print_array_row, TNode *branch, int *n, TTree *tree, char **anc_names, seqGenOptions *options);
void WriteAncestralSequencesNode(vector<vector<char> >& print, int *print_array_row, TTree *tree, int *nodeNo, TNode *des, char **anc_names);
void PrintTitle();
void Print_NEXUS_Header(list<inTree*> inputTrees, seqGenOptions *options, ostream& out);
void Reset_Run(list<inTree*>& inTrees, seqGenOptions *options, bool reset);
void readInputTrees(list<inTree*>& inputTrees, 	seqGenOptions *options, inClade *global_environment, 
					vector<string>& OTU_names, int *numTrees, int numTaxa);
unsigned int rand_seed();
void setTRS(list<inTree*>& inputTrees, seqGenOptions *options, double max_path_length, double min_branch);
void CheckTrees(list<inTree*>& inputTrees, seqGenOptions *options);
size_t trimMSA(vector<vector<char> >& print, vector<bool>& empty_columns, seqGenOptions *options);
void fiddle_with_trees(list<inTree*>& inTrees, seqGenOptions *options);
void QuickTest();
void CheckMarkovCodonLikelihoods(Dependency *selective, Dependency *neutral);
int Nucleotide_Sequence_2_Index ( string sequence );
int Get_Nucleotide_Value ( char nucleotide );


int point_to_me (int return_val) { cerr << "pointed to me." << endl; return return_val; }
int point2me (int return_val) { cerr << "you chose me." << endl; return return_val; }

////////////////////////////////////////
////////////////////////////////////////

int main(int argc, char *argv[])
{	
	clock_t totalStart;
	double totalSecs;
	vector<ofstream*> simulation_output_streams;

/*	int (*ptr2func)(int) = NULL;
	
	ptr2func = &point_to_me;
	
	int check = (*ptr2func)(5);
	cerr << check << endl;

	ptr2func = &point2me;

	check = (*ptr2func)(10);
	cerr << check << endl;
	
	exit(0);
*/
	//QuickTest();
	//////////
	/// Major player variables.
	//////////
	inClade *global_environment;			// Environment
	list<inTree*> inputTrees;				// Tree + partition structure.
	totalStart = clock();

	// Options parsing, error-checking //
	seqGenOptions *options = new seqGenOptions();
	options->readOptions(argc, argv);

	if (options->gil_Seed.empty() ) {
		for (unsigned int i = 0; i < 4; i++) {
			options->gil_Seed.push_back(rand_seed());
		}
	}
	mt_srand(&options->gil_Seed[0], options->gil_Seed.size());
	SetSeed(options->randomSeed);
	global_environment = new inClade("Global Environment", options);	

	simulation_output_streams.resize(6);
	openOutputStreams(simulation_output_streams, options->output_file_flags, options->output_files, options->path_proposal);

	rasmus_independent_proposals = options->rasmus_independent_proposals;

	cerr << "Dependent Sites Model: " << options->dependence_model_counts << endl;
	cerr << "Neutral Model:         " << options->neutral_model_counts << endl;

	//////////
	/// Primary routine. Preprocesses data, then either calls Forward simulation or path proposals.
	//////////
	Simulate(
		     inputTrees, 
			 global_environment, 
			 simulation_output_streams, 
			 options
			);


	closeOutputStreams(inputTrees, simulation_output_streams, options->output_file_flags);
	options->SpoolWarnings("",true);

	//////////
	/// Cleanup all objects created in run.
	//////////
	delete global_environment;
	delete options;
	//Leakage();

	totalSecs = (double)(clock() - totalStart) / CLOCKS_PER_SEC;
	fprintf(stderr,"Time taken: %G seconds\n", totalSecs);
	
	//Leakage();
	return EXIT_SUCCESS;
}

////////////////////
//// FUNCTIONS /////
////////////////////
void QuickTest()
{
	/// This just allows me to test anything that I want to test that does not need to go through the
	/// Entire simulation routine.
	RateMatrix *rates;
	isNucModel = true;
	numStates = 4;
	numStates_squared = 16;
	numStates_cubed = 64;
	
	rates = new RateMatrix();
	
	rates->InitializeSubstitutionVectors();

//	for (int i = 0; i < numStates; i++) {
//		for (int j = 0; j < numStates; j++) {
//			if (i != j)
//				rates->Qij.at(i*numStates+j) = rates->Qij.at(i+j*numStates) = rndu();
//		}
//	}

//	double sum_neg;
//	for (int i = 0; i < numStates; i++) {
//		sum_neg = 0;
//		for (int j = 0; j < numStates; j++) {
//			if (i != j) {
//				sum_neg += rates->Qij.at(i*numStates+j);
//			} else rates->Qij.at(i*numStates+j) = 0;
//		}
//		rates->Qij.at(i*numStates+i) = -sum_neg;
//	}

//	rates->printQij();

	//////////
	// Testing Ziheng's exponentiation routine using A felsenstein representation...
	// Qij = S * pi_j, if i != j
	// Qii = -S * \sum_{j, j!=i} pi_j
	//
	// This leads to Pij(t) = e^{-st}+(1-e^{-st})*pi_j  if i=j
	//               Pij(t) = (1-e^{-st})*pi_j  if i != j
	//
	// S = constant, t = branch length.
	//////////
	rates->pi.at(0) = 0.1;
	rates->pi.at(1) = 0.2;
	rates->pi.at(2) = 0.3;
	rates->pi.at(3) = 0.4;

	double S = 0.1;
	double t = 10;
	double sum_neg;
	for (int i = 0; i < numStates; i++) 
		for (int j = 0; j < numStates; j++)
			if (i != j)
				rates->Qij.at(i*numStates+j) = S * rates->pi.at(j);

	for (int i = 0; i < numStates; i++) {
		sum_neg = 0;
		for (int j = 0; j < numStates; j++) 
			if (i != j) sum_neg += rates->Qij.at(i*numStates+j);
			else rates->Qij.at(i*numStates+j) = 0;
		rates->Qij.at(i*numStates+i) = -sum_neg;
	}

	rates->SetupMatrix(false);
	Site site;
	rates->setPij(site, t, NoRates);
	//rates->setPij(Site site, BL, rateHetero);
	rates->printQij();

	cerr << "    //{" << rates->pi.at(A) << "," << rates->pi.at(C) << "," << rates->pi.at(G) << "," << rates->pi.at(T) << "}" << endl;
	cerr << "    //S=" << S << " t=" << t << endl;
	rates->printPij();

	cerr << "    //P_{AA}=" << exp(-S*t)+(1-exp(-S*t))*rates->pi.at(A) << " ";
	cerr << "P_{AC}=" << (1-exp(-S*t))*rates->pi.at(C) << " ";
	cerr << "P_{AG}=" << (1-exp(-S*t))*rates->pi.at(G) << " ";
	cerr << "P_{AT}=" << (1-exp(-S*t))*rates->pi.at(T) << endl;

	cerr << "    //P_{CA}=" << (1-exp(-S*t))*rates->pi.at(A) << " ";
	cerr << "P_{CC}=" << exp(-S*t)+(1-exp(-S*t))*rates->pi.at(C) << " ";
	cerr << "P_{CG}=" << (1-exp(-S*t))*rates->pi.at(G) << " ";
	cerr << "P_{CT}=" << (1-exp(-S*t))*rates->pi.at(T) << endl;

	cerr << "    //P_{GA}=" << (1-exp(-S*t))*rates->pi.at(A) << " ";
	cerr << "P_{GC}=" << (1-exp(-S*t))*rates->pi.at(C) << " ";
	cerr << "P_{GG}=" << exp(-S*t)+(1-exp(-S*t))*rates->pi.at(G) << " ";
	cerr << "P_{GT}=" << (1-exp(-S*t))*rates->pi.at(T) << endl;

	cerr << "    //P_{TA}=" << (1-exp(-S*t))*rates->pi.at(A) << " ";
	cerr << "P_{TC}=" << exp(-S*t)+(1-exp(-S*t))*rates->pi.at(C) << " ";
	cerr << "P_{TG}=" << (1-exp(-S*t))*rates->pi.at(G) << " ";
	cerr << "P_{TT}=" << exp(-S*t)+(1-exp(-S*t))*rates->pi.at(T) << endl;

	exit(0);
}

void Simulate(
	   	      list<inTree*>& inputTrees, 
			  inClade *global_environment,
			  vector<ofstream*>& simulation_output_streams, 
			  seqGenOptions *options
			 )
{
	cerr << "Point-> Simulate() IN" << endl;

	paleobiology *paleontologicalProcess;	// Fossil deposition
	list<eventTrack*> events;				// Event tracking
	int numTaxa = 0;		// To check if all trees have the same number of taxa.
	vector<string> OTU_names;
	Statistics stats(options->number_of_datasets_per_tree);
	vector<SampleStatistics>::iterator stats_it = stats.sample_stats.begin();
 
	if (options->deposit_fossils) 
		paleontologicalProcess = new paleobiology(options->paleo_root_age, options->fossil_deposition_rate);


	//////////
	/// Gathering inTrees, checking sanity 
	//////////
	int numTrees = 0;
	readInputTrees(inputTrees,options,global_environment,OTU_names,&numTrees,numTaxa);
	CheckTrees(inputTrees, options);
	setRates(inputTrees, options);
	if(options->writeAncestors) {
		string node_names = "";
		inputTrees.front()->CheckPhylogeneticAncestralNodes(node_names, SET);
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
			if(!(*it)->CheckPhylogeneticAncestralNodes(node_names, CHECK)) {
				options->writeAncestors = false;
				options->SpoolWarnings("Cannot write ancestral sequences unless all partition tree topologies are equivalent.");
				break;
			}
		}
	}

	cerr << "Point-> Simulate() TREES READ" << endl;

	// Set the function pointers.
	//for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
	//	for (list<TNode*>::iterator jt = (*it)->my_tree->nodeList.begin(); jt != (*it)->my_tree->nodeList.end(); ++jt) {
	//		(*jt)->set_Qptr(Qd, Pc, nij);
	//	}
	//}

	if (Qd && Pc) ptr2init = &iQdPc;
	else if (Qd && !Pc && !nij) ptr2init = &iQdP; 
	else if (!Qd && Pc) ptr2init = &iQPc; 
	else if (Qd && nij) ptr2init = &iQN; 
	else ptr2init = &iQP;

	if (Qd && Pc) ptr2update = &uQdPc;
	else if (Qd && !Pc && !nij) ptr2update = &uQdP; 
	else if (!Qd && Pc) ptr2update = &uQPc; 
	else if (Qd && nij) ptr2update = &uQN; 
	else ptr2update = &uQP;

	cerr << "Point-> Simulate() FUNCTION PTR SET" << endl;

	if (!options->quiet && !options->path_proposal) {
		if (options->fileFormat == NEXUSFormat) 
			Print_NEXUS_Header(
							   inputTrees,
							   options,
							   (
							     (options->output_file_flags[seqGenOptions::ma])
							     ? *simulation_output_streams.at(seqGenOptions::ma)
							     : cout
							   )
							  );

		stringstream header;
		header << "Taxa bit representation:" << endl;
		inputTrees.front()->Print_Trace_Header(header, OTU_names, options);

		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
			header << "Clade List (<taxon members> = <clade name>): " << endl;
			for (list<TNode*>::iterator itn = (*it)->my_tree->nodeList.begin(); itn != (*it)->my_tree->nodeList.end(); itn++) {
				if ( !((*itn)->clade_label.empty()) ) {
					for (int j = 0; j < (*it)->my_tree->numTips; j++) header << (*itn)->bipartition.at(j);
					header << " = " << (*itn)->clade_label << endl;
				}
			}
		}

		if (!options->output_file_flags[seqGenOptions::trace]) 
			cout << header.str() << endl;
		else 
			*simulation_output_streams.at(seqGenOptions::trace) << header.str() << endl;
	}

	cerr << "Point-> Simulate() PREPROCESS SUCCESSFUL" << endl;

	//////////
	/// MAIN ROUTINE: All events that need to be done for each dataset are done inside this
	/// loop, and the Forward_Simulation or Propose_Path function is called.
	//////////
	for (int i = 0; i < options->number_of_datasets_per_tree; i++) {
		cerr << "---------- Replicate " << i+1 << "/" << options->number_of_datasets_per_tree << " ----------" << endl;
		double sum = 0;
		int numSites = 0;
		eventNo = 0;
		events.clear();

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
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
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
			(*it)->nameAncestralNodes();
		}
		//
		// (3) Finally, set TRS steps
		//
		setTRS(inputTrees, options, max_path_length, min_branch);
		//
		// Print out tree with perturbation (Perturbed tree will print out different trees for each
		// dataset) if the user wants the anc_trees.
		//
		if (!options->path_proposal) {
			if (options->output_file_flags[seqGenOptions::tree]) {
				if (options->simulation_step_type == TIME_RELATIVE_STEPS || options->simulation_step_type == UNIFORMIZATION) {
					for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it)
						(*it)->Print_Time_Rel_Tree(*simulation_output_streams.at(seqGenOptions::tree), options->writeAncestors);
				} else {
					for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it)
						(*it)->Print_Newick_Tree(*simulation_output_streams.at(seqGenOptions::tree), false);
				}
			}
		}
		//
		//////////

		//////////
		/// Setup of sequences for this run.
		//////////
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
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
						++rit;
					}
				}
			}
		}
		//
		//////////

		//////////
		/// One final check to see if all coding regions add up to 0 (mod 3)
		//////////
		if (isNucModel) {
			for(list<inTree*>::reverse_iterator rit = inputTrees.rbegin(); rit != inputTrees.rend(); ++rit) {
				if ( (*rit)->my_tree->root->nodeEnv->rateHetero == CodonRates ) {
					if ( ((*rit)->partitionLength + (*rit)->partitionLength) % 3 != 0 )
						options->SpoolWarnings("The lengths of all codon sequences do not add up to a factor of three.");
				}
			}
		}
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) 
			(*it)->partitionRate *= numSites / sum;

		//////////
		/// Need to calculate likelihoods under what conditions?
		//////////
		// i.   End-point conditioned runs
		// ii.  When specfile places probabilities of each state (Pfam-like?? Motif conservation)
		//////////
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
			//////////
			/// Check to see if there are dependencies specified.
			//////////
			if ( !options->dependency_file.empty() || order_3_markov || Human_Data_simulation ) {
				if (options->path_proposal) {
					(*it)->my_tree->dep.push_back(new Dependency(options->context_order, options->dependence_superscript));	// Current: hard-coded 3rd-order Markov.
					(*it)->my_tree->neutral_dep.push_back(new Dependency(options->context_order, global_environment));
				} else {
					cerr << "Not path proposal (Simulate())" << endl;
					if (order_3_markov) {
						cerr << "order_3_markov (Simulate())" << endl;
						(*it)->my_tree->neutral_dep.push_back(new Dependency(options->context_order, global_environment));
						(*it)->my_tree->dep.push_back(new Dependency(options->context_order, atof(options->dependence_superscript.c_str()), options->output_files));
					} else if (Human_Data_simulation) {
						cerr << "Human Data simulation (Simulate()) 1" << endl;
						(*it)->my_tree->dep.push_back(new Dependency(options->context_order, 3, options->dependence_model_counts));
						cerr << "Human Data simulation (Simulate()) 2" << endl;
						(*it)->my_tree->neutral_dep.push_back(new Dependency(options->context_order, 3, options->neutral_model_counts));
						cerr << "Human Data simulation (Simulate()) 3" << endl;
						//CheckMarkovCodonLikelihoods((*it)->my_tree->dep.front(), (*it)->my_tree->neutral_dep.front());
					} else {
						cerr << "Didn't enter a model? (-O <markov_sup> OR -2 <dep_counts> -3 <neutral_counts>) " << endl;
						exit(EXIT_FAILURE);
					}
					cerr << "Success." << endl;
				}
			}
		}

		cerr << "Point-> Simulate() entering Path_Proposal/Forward_Simulation" << endl;

		if (options->path_proposal) {
			Path_Proposal(
						  inputTrees, 
						  &stats,
						  global_environment, 
						  simulation_output_streams, 
						  options,
						  &events
						 );
		} else {
			cerr << "Point-> Simulate::Forward_Simulation() Beginning simulation." << endl;
			Forward_Simulation(
							   inputTrees, 
							   global_environment, 
							   i,
							   simulation_output_streams, 
							   options,
							   &events
							  );
		}

		//////////
		/// Remove objects from root sequence.
		//////////
		for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
			// The following line of code blows iSG up when a motif sequence is followed by a RANDOM,
			// sequence with the -1 option unset.
			if ( (*it)->rootSeqType != RANDOM || !(*it)->randomInvariableAssignment )
				(*it)->my_tree->root->Remove_Objects(true);
			//
			// Regardless of the motif situation, the general varSites still need to be removed. This is
			// done automatically if the if-stmt is true.
			else (*it)->my_tree->root->Remove_varSites();
			//
			//////////
			for (list<Dependency*>::iterator jt = (*it)->my_tree->dep.begin(); jt != (*it)->my_tree->dep.end(); ++jt)
				delete (*jt);
			(*it)->my_tree->dep.clear();
		}

		if (options->deposit_fossils) {
			paleontologicalProcess->max_path_length = max_path_length;
			cerr << "run " << i << ": Depositing fossils" << endl;
			paleontologicalProcess->doPaleontology(inputTrees.front()->my_tree, &events);
		}

		(*stats_it).calculateStatistics(&events);

		if (options->output_file_flags[seqGenOptions::verb]) verbose_output << "DATASET " << i << ":" << endl;
		if (options->tracking_defined()) {
			if (options->simulation_step_type == TIME_RELATIVE_STEPS || options->simulation_step_type == UNIFORMIZATION) {
				Print_Trace(
						    inputTrees, 
							events, 
							(
							  (options->output_file_flags[seqGenOptions::trace]) 
							  ? *simulation_output_streams.at(seqGenOptions::trace)
							  : cout
							), 
							options->simulation_step_type, 
							0, 
							i
						   );
			} else if (options->output_file_flags[seqGenOptions::trace]) 
				*simulation_output_streams.at(seqGenOptions::trace) << endl;
			else cout << endl;
		}

		if (!options->path_proposal) {
			Print_Root(
				 	   inputTrees, 
				 	   (
				 	     (options->output_file_flags[seqGenOptions::root]) 
				 	     ? *simulation_output_streams.at(seqGenOptions::root)
				 	     : cout
				 	   ), 
				 	   options 
				 	  );
			Print_Seq(
					  inputTrees, 
					  (
					    (options->output_file_flags[seqGenOptions::seq]) 
					    ? *simulation_output_streams.at(seqGenOptions::seq)
					    : cout
					  ), 
					  options, 
					  i
					 );
			Print_MSA(
					  inputTrees, 
					  (
					    (options->output_file_flags[seqGenOptions::ma]) 
					    ? *simulation_output_streams.at(seqGenOptions::ma)
					    : cout
					  ), 
					  options, 
					  i, 
					  events
					 );

		} else {
			//(*stats_it).calculateStatistics(&events);
			if (stats.returnM() < (*stats_it).m_i()) stats.setM((*stats_it).m_i());
			(*stats_it).set_numEvents(eventNo);
			++stats_it;
		}

		Reset_Run(inputTrees, options, (i+1 != options->number_of_datasets_per_tree));
	}

	//////////
	/// If this was a path proposal, calculate necessary statistics.
	//////////
	if (options->path_proposal) {
		stats.setWeights();
		stats.reportStatistics();
		//stats.calc_ESS();
	}

	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
		(*it)->my_tree->global_alignment->insert_sites.clear();
		for (list<TNode*>::iterator it2 = (*it)->my_tree->nodeList.begin(); it2 != (*it)->my_tree->nodeList.end(); ++it2) {
			delete (*it2)->branch->rates;
			delete (*it2)->branch;
			delete (*it2);
		}
		for (list<inClade*>::iterator it2 = (*it)->my_tree->treeEnv.begin(); it2 != (*it)->my_tree->treeEnv.end(); ++it2)
			delete (*it2);
		for (list<inMotif*>::iterator it2 = (*it)->motif_specs.begin(); it2 != (*it)->motif_specs.end(); ++it2)
			delete (*it2);
		delete (*it)->my_tree->global_alignment;
		delete (*it)->my_tree;
		delete (*it);
	}

	cerr << "Point-> Simulate() OUT" << endl;
	if (options->deposit_fossils) delete paleontologicalProcess;
}

/// This function is used to check the values of the lookup tables given various codon strings.
// INPUT: contextDependence data (lookup_table is private variable).
void CheckMarkovCodonLikelihoods(
								 Dependency *selective,
								 Dependency *neutral
								)
{
	ifstream is;
	string nucl_seq;
	int index;
	nucl_seq.assign(10, 'x');
	cerr << "Point-> CheckMarkovCodonLikelihood() IN" << endl;

	is.open("check_nucl_seqs.test");

	while(is.good()) {
		// Read in line of values. Element 1 is seed ktuplet (K) probability. Following are Pr(N|K).
		getline(is, nucl_seq);
		index = Nucleotide_Sequence_2_Index(nucl_seq);
		if (nucl_seq.size() == 6) {
			cerr << nucl_seq 
				 << "    " << index
				 << "    " << selective->context.return_lt_value(2, index)
				 << " (" << selective->context.return_lt_value(0, index) << ")"
				 << "    " << neutral->context.return_lt_value(2, index)
				 << " (" << neutral->context.return_lt_value(0, index) << ")"
				 << endl; 
		} else {
			cerr << nucl_seq 
				 << "    " << index
				 << "    " << selective->context.return_lt_value(1, index)
				 << "    " << neutral->context.return_lt_value(1, index)
				 << endl; 
		}
	}

	// int contextDependence::getOffset(
	//							 int environment,
	//							 int codon_position,
	//							 short i, 
	//							 short j
	//							)
	// CODON POSITION::: XXX123XXX
	// CODONS:::		 AAAAAAAAA
	// SEQ POSITION:::   012345678
	// MOD 3:::			 012012012
	//																	   A  C
	int AAAAAAAAA_index = Nucleotide_Sequence_2_Index("AAAAAAAAA");
	int AAAXAAAAA_index;
	int AAAAXAAAA_index;
	int AAAAAXAAA_index;
	int CATGATGTA_index = Nucleotide_Sequence_2_Index("CATGATGTA");
	int CATTATGTA_index = Nucleotide_Sequence_2_Index("CATTATGTA");
	cerr << "CATGATGTA->CATTATGTA: " << selective->context.return_lt_value(1, CATGATGTA_index)
		 << "   " << selective->context.return_lt_value(1, CATTATGTA_index)
		 << "   OFFSET_DIFF: " << selective->context.getOffset(1,0,stateCharacters.find("G"),stateCharacters.find("T"))
		 << "  VAL: " 
		 << selective->context.return_lt_value(1, CATGATGTA_index+selective->context.getOffset(1,0,stateCharacters.find("G"),stateCharacters.find("T"))) 
		 << endl << endl;

	cerr << "CATTATGTA->CATGATGTA: " << selective->context.return_lt_value(1, CATTATGTA_index)
		 << "   " << selective->context.return_lt_value(1, CATGATGTA_index)
		 << "   OFFSET_DIFF: " << selective->context.getOffset(1,0,stateCharacters.find("T"),stateCharacters.find("G"))
		 << "  VAL: " 
		 << selective->context.return_lt_value(1, CATTATGTA_index+selective->context.getOffset(1,0,stateCharacters.find("T"),stateCharacters.find("G"))) 
		 << endl << endl;

	cerr << "AAAAAAAAA->AAACAAAAA: " << selective->context.getOffset(1, 0, 0, 1) << endl;
	AAAXAAAAA_index = Nucleotide_Sequence_2_Index("AAACAAAAA");
	cerr << "     VALUES:  " << selective->context.return_lt_value(1, AAAAAAAAA_index)
		 << "   " << selective->context.return_lt_value(1, AAAXAAAAA_index)
		 << "   OFFSET_DIFF: " << selective->context.getOffset(1,0,0,1)
		 << "  VAL: " << selective->context.return_lt_value(1, AAAAAAAAA_index+selective->context.getOffset(1,0,0,1)) << endl;
	cerr << "AAAAAAAAA->AAAGAAAAA: " << selective->context.getOffset(1, 0, 0, 2) << endl;
	cerr << "AAAAAAAAA->AAATAAAAA: " << selective->context.getOffset(1, 0, 0, 3) << endl;
	cerr << "AAACAAAAA->AAAAAAAAA: " << selective->context.getOffset(1, 0, 1, 0) << endl;
	cerr << endl;
	cerr << "AAAAAAAAA->AAAACAAAA: " << selective->context.getOffset(1, 1, 0, 1) << endl;
	cerr << "AAAAAAAAA->AAAAGAAAA: " << selective->context.getOffset(1, 1, 0, 2) << endl;
	cerr << "AAAAAAAAA->AAAATAAAA: " << selective->context.getOffset(1, 1, 0, 3) << endl;
	cerr << "AAACAAAAA->AAAAAAAAA: " << selective->context.getOffset(1, 1, 1, 0) << endl;
	cerr << endl;
	cerr << "AAAAAAAAA->AAAAACAAA: " << selective->context.getOffset(1, 2, 0, 1) << endl;
	cerr << "AAAAAAAAA->AAAAAGAAA: " << selective->context.getOffset(1, 2, 0, 2) << endl;
	cerr << "AAAAAAAAA->AAAAATAAA: " << selective->context.getOffset(1, 2, 0, 3) << endl;
	cerr << "AAACAAAAA->AAAAAAAAA: " << selective->context.getOffset(1, 2, 1, 0) << endl;
	cerr << endl;

	exit(0);
}

int
Nucleotide_Sequence_2_Index ( string sequence )
{
	int char_val;
	int power_of_4 = 1;
	int index = 0;
	for (string::reverse_iterator it = sequence.rbegin(); it != sequence.rend(); ++it) {
		index += Get_Nucleotide_Value(*it) * power_of_4;
		power_of_4 *= 4;
	}
	return index;
}

int 
Get_Nucleotide_Value ( char nucleotide )
{
	return (int) stateCharacters.find(nucleotide);
}

void Path_Proposal(
	   	    	   list<inTree*>& inputTrees, 
	   	    	   Statistics *stats,
				   inClade *global_environment,
				   vector<ofstream*>& simulation_output_streams, 
				   seqGenOptions *options,
				   list<eventTrack*> *events
				  )
{
	vector<string> OTU_names;

	forward_simulation = false;
	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
		(*it)->path = new PathProposal();
		(*it)->path->parseEndpointInfile(options->end_point_condition_file);
		(*it)->path->setNodeSequences((*it)->my_tree);
		(*it)->path->setUnspecifiedNodeSequences((*it)->my_tree, (*it)->partitionLength);
		(*it)->my_tree->setTransitionProbabilities();

		//////////
		/// Propose paths. If we are emulating a forward run, then fill the event histories.
		//////////
		Create_Global_Arrays((*it)->my_tree, (*it)->my_tree->root->seq_evo.size());
		if (options->epc_emulate_forward_simulation) {
			(*it)->path->setEventHistory(events, options->forward_sim_event_file_to_emulate);
			(*it)->path->emulateForwardSimulation((*it)->my_tree, events);
			return;
		} else {
			// Set "current" path (for MCMC; otherwise, this is the only path).
			(*it)->path->Evolve((*it)->my_tree, events);
			for (list<eventTrack*>::iterator it2 = (*events).begin(); it2 != (*events).end(); ++it2) 
				(*it2)->Compute_MSA_Positions((*it)->my_tree, 0);
			// If MCMC, need to perform some steps to set up path and run MCMC.
			if (options->num_mcmc_steps) {
				(*it)->path->epc_events.assign((*events).begin(), (*events).end());
				stats->MCMC_run((*it)->my_tree, options->num_mcmc_steps, options->output_files, (*it)->path);
			}
		}

		for (list<eventTrack*>::iterator it2 = (*events).begin(); it2 != (*events).end(); ++it2) 
			(*it2)->Compute_MSA_Positions((*it)->my_tree, 0);

		delete (*it)->path;
	}
}

void Forward_Simulation(
	   	    	   		list<inTree*>& inputTrees, 
				   		inClade *global_environment,
				   		int replicate,
				   		vector<ofstream*>& simulation_output_streams, 
				   		seqGenOptions *options,
					 	list<eventTrack*> *events
					   )
{
	int size_prev_msa_positions = 0;

	forward_simulation = true;		// Declared at top of main. bad style, yes, but...
	cerr << "Made it to EvolveSequences!! But left lineage-specific rates untouched, so far" << endl;
	size_t x = 0;
	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it, x++) {
		(*it)->sim = new ForwardSimulation();
		(*it)->sim->EvolveSequences(*it, events, options);
		
		//(*it)->my_tree->global_alignment->Print();
		
		// Sort out the positions of the sequences.
		int msa_size = (*it)->my_tree->root_numSites;
		if (options->tracking_defined()) {
			for (list<eventTrack*>::iterator it2 = (*events).begin(); it2 != (*events).end(); ++it2) {
				msa_size = (*it2)->Compute_MSA_Positions((*it)->my_tree, size_prev_msa_positions);
			}
			size_prev_msa_positions += msa_size;			
		}
		if (options->simulation_step_type == DISCRETE_EVOLUTIONARY_STEPS && options->tracking_defined()) {
			Print_Trace(
						inputTrees, 
						*events, 
						( 
						  (options->output_file_flags[seqGenOptions::trace]) 
						  ? *simulation_output_streams.at(seqGenOptions::trace) 
						  : cout
						), 
						options->simulation_step_type, 
						(*it)->treeNum, 
						replicate
					   );
			(*events).clear();		// DES prints each partition separately, since they cannot be related.
		}
	}
}

//////////
/// This applies tree perturbations and treeScales to tree.
//////////
void fiddle_with_trees(
					   list<inTree*>& inTrees, 
					   seqGenOptions *options
					  ) 
{
	double all_paths_sum, mean, std_dev;

	//////////
	/// (1) Perturb and (2) scale tree
	//////////
	//
	// Set the perturbation values for the trees.
	for (list<inTree*>::iterator it = inTrees.begin(); it != inTrees.end(); ++it)
		(*it)->perturbTree(options->perturbTree);
	//
	//
	// If option is set, output tree stats:
	// * Average root-to-tip path lengths (+ std. dev).
	size_t x = 0;
	for (list<inTree*>::iterator it = inTrees.begin(); it != inTrees.end(); ++it, x++) {
		(*it)->calcDistancesFromRoot();
		if (options->treelength_check || options->treelength_scale || (*it)->scaleTree) {
			all_paths_sum = mean = std_dev = 0.0;
			for (vector<TNode*>::iterator jt = (*it)->my_tree->tips.begin(); jt != (*it)->my_tree->tips.end(); ++jt)
				all_paths_sum += (*jt)->DistanceFromRoot;
			mean = all_paths_sum / (*it)->my_tree->tips.size();
			for (vector<TNode*>::iterator jt = (*it)->my_tree->tips.begin(); jt != (*it)->my_tree->tips.end(); ++jt)
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
				for (vector<TNode*>::iterator jt = (*it)->my_tree->tips.begin(); jt != (*it)->my_tree->tips.end(); ++jt)
					all_paths_sum += (*jt)->DistanceFromRoot;
				mean = all_paths_sum / (*it)->my_tree->tips.size();
				for (vector<TNode*>::iterator jt = (*it)->my_tree->tips.begin(); jt != (*it)->my_tree->tips.end(); ++jt)
					std_dev += ((*jt)->DistanceFromRoot-mean)*((*jt)->DistanceFromRoot-mean);
				std_dev /= (*it)->my_tree->tips.size();
				cout << mean << " " << sqrt(std_dev) << endl;
			} else {
				cout << "If you wanted to scale average root-to-tip path length, input number must be negative." << endl; 
				cout << mean << " ... " << sqrt(std_dev) << endl;
				exit(EXIT_SUCCESS);
			}
		}
	}
	//
	//////////

	ofstream treescale_out;
	string treescale_outfile;
	treescale_outfile = options->output_files + ".scale_tree";
	treescale_out.open(treescale_outfile.c_str(), ios::trunc | ios::out);
	for (list<inTree*>::iterator it = inTrees.begin(); it != inTrees.end(); ++it) 
		(*it)->Print_Newick_Tree(treescale_out, false);

	return;
}

//////////
/// Remove single-replicate specific objects, which will be recreated in next run
//////////
void Reset_Run(
			   list<inTree*>& inTrees, 
			   seqGenOptions *options, 
			   bool reset
			  )
{
	//////////
	/// This small chunk of code applies ONLY to random sequences using the PROSITE library
	/// (option -1). 
	//////////
	for (list<inTree*>::iterator it = inTrees.begin(); it != inTrees.end(); ++it) {
		if ((*it)->rootSeqType == RANDOM && options->random_sequence_proportion_motif) {
			(*it)->motif_specs.clear();
			// Destroy inMotifs, will reset them for next run.
			for (list<inMotif*>::iterator jt = (*it)->my_tree->treeEnv.front()->my_motifs.begin(); 
										  jt != (*it)->my_tree->treeEnv.front()->my_motifs.end(); 
										  ++jt) 
			{
				delete *jt;
			}
			(*it)->my_tree->treeEnv.front()->my_motifs.clear();
		}

		for (list<TNode*>::iterator it2 = (*it)->my_tree->nodeList.begin(); it2 != (*it)->my_tree->nodeList.end(); ++it2) {
			(*it2)->evolvingSequence->evolutionaryAttributes.clear();
			delete (*it2)->evolvingSequence;
			// Variable region lists.
			if (reset) (*it2)->addGeneral_varSites();
		}
	}
	
	eventNo = 0;
}

//////////
/// Read the trees in the partition file, parse them, get statistics of trees, and perform
/// necessary bookkeeping.
//////////
void readInputTrees(
					list<inTree*>& inputTrees, 
					seqGenOptions *options, 
					inClade *global_environment, 
					vector<string>& OTU_names, 
					int *numTrees, 
					int numTaxa
				   )
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
			cerr << "Treefile error: Trailing characters after last tree. (Last tree may not end with ';\n')" << endl << treeInput << endl;
			exit(EXIT_FAILURE);
		}

		if(treeInput == "")
			continue;
		(*numTrees)++;
		inTree *tempTree = new inTree(options, *numTrees, global_environment);
		tempTree->perform_checks(treeInput, options);
		//////////
		/// Parsing lineages MUST come before parsing the tree. parseLineages collects motif data,
		/// all except the sitemap, which is done in parseTree-->Read_MA, which needs to know the
		/// motif marker in order to place the sitemap correctly.
		//////////
		if(options->lineageSpecificFile != "") tempTree->parseLineages(options);
		//////////
		/// Get the tree from the file
		//////////
		tempTree->parseTree(treeInput, options);
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
		tempTree->my_tree->my_iTree = tempTree;
		inputTrees.push_back(tempTree);
	}
}

void setRates(
			  list<inTree*>& inputTrees,
			  seqGenOptions *options
			 )
{
	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
		for (list<TNode*>::iterator jt = (*it)->my_tree->nodeList.begin(); jt != (*it)->my_tree->nodeList.end(); ++jt) {
			(*jt)->branch->rates = new RateMatrix();
			(*jt)->branch->rates->setModelConditions(
													 options->inputModel, 
													 options->rate_matrix_values, 
													 options->global_pi
													);
			(*jt)->branch->rates->SetModel((*jt)->nodeEnv);
			//////////
			/// Discrete category setup.
			//////////
			if (options->default_rateHetero == DiscreteGammaRates) {
				(*jt)->branch->rates->num_categories = options->num_discrete_gamma_categories;
				(*jt)->branch->rates->alphaGamma = options->alpha;
				DiscreteGamma(
							  (*jt)->branch->rates->freqRate, 
							  (*jt)->branch->rates->catRate, 
							  (*jt)->branch->rates->alphaGamma, 
							  (*jt)->branch->rates->alphaGamma, 
							  (*jt)->branch->rates->num_categories, 
							  false
							 );
				(*jt)->nodeEnv->rateHetero = DiscreteGammaRates;
			} else if (options->default_rateHetero == GammaRates) {
				cerr << "Setting gamma rates only (not discrete gamma) is left undealt with." << endl;
				(*jt)->branch->rates->alphaGamma = options->alpha;
				(*jt)->nodeEnv->rateHetero = GammaRates;
			}
		}
	}
}

//////////
/// This function checks to see if all lineages have marker matches, as well as other 
/// items of interest.
//////////
void CheckTrees(
			    list<inTree*>& inputTrees, 
			    seqGenOptions *options
			   )
{
	stringstream warning;
	bool set = false;
	// Now check each inMotif, find out if each marker defined is used or not. //
	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
		if ((*it)->rootSeqType != RANDOM) {
			for (list<inMotif*>::iterator jt = (*it)->motif_specs.begin(); jt != (*it)->motif_specs.end(); ++jt) {
				if ( (*jt)->sitemap.empty() ) {
					warning << "Motif " << (*jt)->marker << " is not used in any of the input root sequences." << endl;
					set = true;
				}
			}
		}
	}
	if (set) options->SpoolWarnings(warning.str());
}

void Print_Motifs(
				  list<inTree*>& inputTrees, 
				  ostream& motif_out
				 ) 
{
	size_t curr_tree = 0;

	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
		motif_out << curr_tree << ":" << endl;
		for (list<inMotif*>::iterator jt = (*it)->motif_specs.begin(); jt != (*it)->motif_specs.end(); ++jt) {
			motif_out << "Name: " << (*jt)->name << endl;
			motif_out << "Regular Expression: " << (*jt)->regex << endl;
			motif_out << "Placement on Root Sequence: " << (*jt)->sitemap << endl;
		}
	}
}

void Print_Trace(
				 list<inTree*> inputTrees, 
				 list<eventTrack*> events, 
				 ostream& trace_out, 
				 int step_type, 
				 int treeNo, 
				 int setNo
				)
{
	trace_out << ">Dataset_" << setNo;

	if(step_type == TIME_RELATIVE_STEPS || step_type == UNIFORMIZATION) {
		trace_out << endl;
		vector<eventTrack*> thisTip_event;
		thisTip_event.clear();
		int i = 0;
		for (list<eventTrack*>::iterator it = events.begin(); it != events.end(); ++it)
			thisTip_event.push_back((*it));

		if (thisTip_event.size() > 0) {
/*			double min_time;
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
*/
			string event_outstring;
			int i = 0;
			for (vector<eventTrack*>::iterator it = thisTip_event.begin(); it != thisTip_event.end(); ++it, i++) {
				event_outstring = (*it)->Print_Event();
				trace_out << event_outstring;
			}
			thisTip_event.clear();
		}
	} else {	// DISCRETE_EVOLUTIONARY_STEPS
		// Locate the current tree:
		string out_event;
		inTree *curr_tree = NULL;
		if(events.size() > 0) {
			for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
				if ((*it)->treeNum == treeNo) curr_tree = (*it);
			}
	
			if(curr_tree == NULL) {
				cerr << "Print_Trace: curr_tree has null value." << endl;
				exit(EXIT_FAILURE);
			}

			trace_out << "__partition_" << treeNo << endl;
			for (list<eventTrack*>::iterator it = events.begin(); it != events.end(); ++it) {
				out_event = (*it)->Print_Event();
				trace_out << out_event;
			}
		}  else {
			trace_out << "__partition_" << treeNo << endl << "-" << endl;
		}
	}
	
	// Free all of the event tracking.
	for (list<eventTrack*>::iterator it = events.begin(); it != events.end(); ++it) {
		delete (*it);
	}
}

void Print_Root(
				list<inTree*> inputTrees, 
				ostream& root_out, 
				seqGenOptions *options
			) 
{
	if (options->output_file_flags[seqGenOptions::verb]) verbose_output << "Root: ";

	for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
		for(int i = 0; i < (*it)->my_tree->root_numSites; i++) {
			root_out << stateCharacters[(*it)->my_tree->root->seq_evo.at(i).returnState()];
			if (options->output_file_flags[seqGenOptions::verb]) verbose_output << stateCharacters[(*it)->my_tree->root->seq_evo.at(i).returnState()];
		}
	}
	root_out << endl;
	if (options->output_file_flags[seqGenOptions::verb]) verbose_output << endl << endl;
}

void Print_NEXUS_Header(
						list<inTree*> inputTrees, 
						seqGenOptions *options, 
						ostream& out
					   ) 
{
	int root_length = 0;
	int partitions = 0;
	bool gene_flag = false;

	for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
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
	for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
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
			for (vector<double>::iterator pi_T = (*it)->my_tree->root->branch->rates->pi.begin();
										  pi_T!= (*it)->my_tree->root->branch->rates->pi.end();
										  pi_T++)
				out << setprecision(3) << setw(6) << (*pi_T) << " ";
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

	if (options->output_file_flags[seqGenOptions::verb]) verbose_output << "Sequences:" << endl;

	if(options->writeAncestors) {
		for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
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

			if (options->output_file_flags[seqGenOptions::verb]) {
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
				if (options->output_file_flags[seqGenOptions::verb]) 
					verbose_output << print.at(i).at(j);
			}
			seq_out << endl;
			if (options->output_file_flags[seqGenOptions::verb]) verbose_output << endl;
		}
		seq_out << endl;
		if (options->output_file_flags[seqGenOptions::verb]) verbose_output << endl;
	} else {
		// tips were correctly sorted in the beginning, so using the front tree only is OK.
		for (int i = 0; i < numTips; i++) {
			seq_out << ">" << inputTrees.front()->my_tree->names.at(i) << endl;
			if (options->output_file_flags[seqGenOptions::verb]) {
				int j;
				for(j = 0; j < 9 && j < inputTrees.front()->my_tree->names.at(i).size(); j++) {
					verbose_output << inputTrees.front()->my_tree->names.at(i).at(j);
				}
				while(j++ < 10) verbose_output << " ";
			}

			int aggregate = 0;
			for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
				for(int j = 0; j < (*it)->my_tree->tips.at(i)->seq_evo.size(); aggregate++, j++) {
					seq_out << stateCharacters[(*it)->my_tree->tips.at(i)->seq_evo.at(j).returnState()];
					if(options->fileFormat == FASTAFormat) {
						if((aggregate+1) % options->output_width == 0 && (aggregate+1) != (*it)->my_tree->tips.at(i)->seq_evo.size()) 
							seq_out << endl;
					}
					if (options->output_file_flags[seqGenOptions::verb]) verbose_output << stateCharacters[(*it)->my_tree->tips.at(i)->seq_evo.at(j).returnState()];
				}
			}
			seq_out << endl;
			if (options->output_file_flags[seqGenOptions::verb]) verbose_output << endl;
		}
		seq_out << endl;
		if (options->output_file_flags[seqGenOptions::verb]) verbose_output << endl;
	}
}      

void WriteAncestralSequencesNode(
								 vector<vector<char> >& print, 
								 int *print_array_row, 
								 TTree *tree, 
								 int *nodeNo, 
								 TNode *des, 
								 char **anc_names
								)
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

void Print_MSA(
			  list<inTree*> inputTrees, 
			  ostream& ma_out, 
			  seqGenOptions *options, 
			  int dataset_num,
			  list<eventTrack*> events
			 ) 
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
	//////////
	/// Filling the vector for printing.
	//////////
	for(list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
		if((*it)->my_tree->treeEnv.front()->rateHetero == CodonRates) gene_flag = true;
		if(options->writeAncestors) {
			print_array_row = taxa_ptr = 0;
			n = inputTrees.front()->my_tree->numTips+1;
			sprintf(anc_names[0], "%d", n);

			vector<insertSite>::iterator bt;
			list<siteModifier>::iterator ct;
			vector<Site>::iterator taxon_iterator = (*it)->my_tree->root->seq_evo.begin();
			bool isSite;
			for (
				 bt =  (*it)->my_tree->global_alignment->insert_sites.begin(); 
				 bt != (*it)->my_tree->global_alignment->insert_sites.end(); 
				 ++bt
				) 
			{
				isSite = true;
				//////////
				/// If: site is related to this insertion site, need to check if deleted
				/// Else: it is a gap character (a '+' if v is set)
				//////////
				if ( isAnc((*it)->my_tree->root, (*bt).fromAnc) ) {
					for (ct = (*bt).modifiers.begin(); ct != (*bt).modifiers.end(); ++ct) {
						if ( (*ct).action == 'd' )	// Site has potentially been removed.
							if ( isAnc((*it)->my_tree->root, (*ct).fromAnc) ) isSite = false;
					}
					//////////
					/// If it isSite is still true after all checks, then return character.
					//////////
					if (taxon_iterator == (*it)->my_tree->root->seq_evo.end()) 
						print.at(print_array_row).push_back('-');
					else 
						print.at(print_array_row).push_back ( 
											   (
											     (isSite) 
											     ? stateCharacters.at((*taxon_iterator).returnState())
											     : '-'
											   )
											  );
					if (isSite) taxon_iterator++;
				} else print.at(print_array_row).push_back( ((options->output_indel_gaps) ? '+' : '-') );
			}

			print_array_row++;
			if(!(*it)->my_tree->rooted)
				WriteAncestralMSA(print, &print_array_row, (*it)->my_tree->root->branch0, &n, (*it)->my_tree, anc_names, options);
			WriteAncestralMSA(print, &print_array_row, (*it)->my_tree->root->branch1, &n, (*it)->my_tree, anc_names, options);
			WriteAncestralMSA(print, &print_array_row, (*it)->my_tree->root->branch2, &n, (*it)->my_tree, anc_names, options);
		} else {
			vector<TNode*>::iterator at;
			vector<insertSite>::iterator bt;
			list<siteModifier>::iterator ct;
			vector<Site>::iterator taxon_iterator;
			bool isSite;
			int j = 0;	// This is a pointer to which print row we are filling. Follows tip ptr.
			for (at = (*it)->my_tree->tips.begin(); at != (*it)->my_tree->tips.end(); ++at, j++) {
				taxon_iterator = (*at)->seq_evo.begin();
				for (
					 bt = (*it)->my_tree->global_alignment->insert_sites.begin(); 
					 bt != (*it)->my_tree->global_alignment->insert_sites.end(); 
					 ++bt
					) 
				{
					isSite = true;
					//////////
					/// If: site is related to this insertion site, need to check if deleted
					/// Else: it is a gap character (a '+' if v is set)
					//////////
					if ( isAnc((*at), (*bt).fromAnc) ) {
						for (ct = (*bt).modifiers.begin(); ct != (*bt).modifiers.end(); ++ct) {
							if ( (*ct).action == 'd' )	// Site has potentially been removed.
								if ( isAnc((*at), (*ct).fromAnc) ) isSite = false;
						}
						//////////
						/// If it isSite is still true after all checks, then return character.
						//////////
						if (taxon_iterator == (*at)->seq_evo.end()) { 
							print.at(j).push_back('-');
							isSite = false;
						} else {
							print.at(j).push_back ( 
												   (
												     (isSite) 
												     ? stateCharacters.at((*taxon_iterator).returnState())
												     : '-'
												   )
												  );
						}
						if (isSite) taxon_iterator++;
					} else print.at(j).push_back( ((options->output_indel_gaps) ? '+' : '-') );
				}
			}
		}
	}

	//////////
	/// Detect empty columns
	//////////
	vector<bool> empty_columns(print.front().size(), false);
	int num_empty_col = 0;
	num_empty_col = trimMSA(print, empty_columns, options);
	num_empty_col = 0;
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
	/// MSA format footer data.
	//////////
	if(options->fileFormat == NEXUSFormat) ma_out << "\t;" << endl << "END;" << endl << endl;
	if(options->fileFormat == FASTAFormat) ma_out << endl;
}

void WriteAncestralMSA(
					  vector<vector<char> >& print, 
					  int *print_array_row, 
					  TNode *des, 
					  int *nodeNo, 
					  TTree *tree, 
					  char** anc_names, 
					  seqGenOptions *options
					 )
{
	if(des->tipNo == -1) {
		(*nodeNo)++;
		cerr << "des->ancestorNo: " << des->ancestorNo << "  nodeNo: " << *nodeNo << endl;
		sprintf(anc_names[*print_array_row], "%d", *nodeNo);
	} else {
		strcpy(anc_names[*print_array_row], (tree->names.at(des->tipNo).c_str()));
	}

	vector<insertSite>::iterator bt;
	list<siteModifier>::iterator ct;
	vector<Site>::iterator taxon_iterator = des->seq_evo.begin();
	bool isSite;

	for (
		 bt =  tree->global_alignment->insert_sites.begin(); 
		 bt != tree->global_alignment->insert_sites.end(); 
		 ++bt
		) 
	{
		isSite = true;
		//////////
		/// If: site is related to this insertion site, need to check if deleted
		/// Else: it is a gap character (a '+' if v is set)
		//////////
		if ( isAnc(des, (*bt).fromAnc) ) {
			for (ct = (*bt).modifiers.begin(); ct != (*bt).modifiers.end(); ++ct) {
				if ( (*ct).action == 'd' )	// Site has potentially been removed.
					if ( isAnc(des, (*ct).fromAnc) ) isSite = false;
			}
			//////////
			/// If it isSite is still true after all checks, then return character.
			//////////
			if (taxon_iterator == des->seq_evo.end()) {
				print.at(*print_array_row).push_back('-');
				isSite = false;
			} else 
				print.at(*print_array_row).push_back ( 
									   (
									     (isSite) 
									     ? stateCharacters.at((*taxon_iterator).returnState())
									     : '-'
									   )
									  );
			if (isSite) taxon_iterator++;
		} else print.at(*print_array_row).push_back( ((options->output_indel_gaps) ? '+' : '-') );
	}

    (*print_array_row)++;	// Go to the next row.

    if(des->tipNo==-1) {
		WriteAncestralMSA(print, print_array_row, des->branch1, nodeNo, tree, anc_names, options);
		WriteAncestralMSA(print, print_array_row, des->branch2, nodeNo, tree, anc_names, options);
    }
}

//////////
/// This function will detect the columns that are empty, and set those columns to true for use
/// in Print_MSA. If output_file_flags[trace] is unset, this function will not set empty_columns
/// to true in any place, since each site is important for the trace file. 
//////////
size_t trimMSA(
			   vector<vector<char> >& print, 
			   vector<bool>& empty_columns, 
			   seqGenOptions *options
			  )
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
		if (options->output_file_flags[seqGenOptions::trace]) {
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
	cerr << "(c) 2005--present Cory Strope" << endl;
}

//////////
/// Fast and dirty random seed algorithm; based off of NR's randqd1
//////////
unsigned int rand_seed()
{
	static unsigned int u = (unsigned int)(time(NULL)+3*getpid());
	return (u = u*1664525u + 1013904223u);
}

void setTRS(
			list<inTree*>& inputTrees, 
			seqGenOptions *options, 
			double max_path_length, 
			double min_branch
		   ) 
{
	//////////
	/// Calculate the max_path_length so that we can correctly place events relative to total sim time.
	//////////
	for (list<inTree*>::iterator it = inputTrees.begin(); it != inputTrees.end(); ++it) {
		if((*it)->my_tree->rooted) {
			(*it)->TimeScale_Tree(max_path_length, (*it)->my_tree->root);
		} else {
			options->SpoolWarnings("Cannot scale unrooted trees for epoch timepoints.");
		}
		(*it)->global_max_path_length = max_path_length;
	}
}

void openOutputStreams(
					   vector<ofstream*>& simulation_output_streams, 
					   vector<bool> outfile_flags, 
					   string& outfile_name_root,
					   bool epc
					  )
{
	string outfile_name;

	if (epc) outfile_name_root += ".epc";
	else outfile_name_root += ".sim";

	if (outfile_flags[seqGenOptions::tree]) {
		outfile_name = outfile_name_root + ".anc_tree";
		simulation_output_streams.at(seqGenOptions::tree) = new ofstream(outfile_name.c_str(), ios::trunc | ios::out);
	} else simulation_output_streams.at(seqGenOptions::tree) = NULL;
	if (outfile_flags[seqGenOptions::root]) {
		outfile_name = outfile_name_root + ".root";	
		simulation_output_streams.at(seqGenOptions::root) = new ofstream(outfile_name.c_str(), ios::trunc | ios::out);
	} else simulation_output_streams.at(seqGenOptions::root) = NULL;
	if (outfile_flags[seqGenOptions::seq]) {
		outfile_name = outfile_name_root + ".seq";
		simulation_output_streams.at(seqGenOptions::seq) = new ofstream(outfile_name.c_str(), ios::trunc | ios::out);
	} else simulation_output_streams.at(seqGenOptions::seq) = NULL;
	if (outfile_flags[seqGenOptions::verb]) {
		outfile_name = outfile_name_root + ".verb";
		simulation_output_streams.at(seqGenOptions::verb) = new ofstream(outfile_name.c_str(), ios::trunc | ios::out);
	} else simulation_output_streams.at(seqGenOptions::verb) = NULL;
	if (outfile_flags[seqGenOptions::ma]) {
		outfile_name = outfile_name_root + ".ma";
		simulation_output_streams.at(seqGenOptions::ma) = new ofstream(outfile_name.c_str(), ios::trunc | ios::out);
	} else simulation_output_streams.at(seqGenOptions::ma) = NULL;
	if (outfile_flags[seqGenOptions::trace]) {
		outfile_name = outfile_name_root + ".trace";
		simulation_output_streams.at(seqGenOptions::trace) = new ofstream(outfile_name.c_str(), ios::trunc | ios::out);
	} else simulation_output_streams.at(seqGenOptions::trace) = NULL;
}

void closeOutputStreams(
						list<inTree*>& inputTrees,
						vector<ofstream*>& simulation_output_streams, 
						vector<bool> outfile_flags
					   )
{
	size_t x = 1;

	if (outfile_flags[seqGenOptions::tree])
		delete simulation_output_streams.at(seqGenOptions::tree);
	if (outfile_flags[seqGenOptions::seq]) 
		delete simulation_output_streams.at(seqGenOptions::seq);
	if (outfile_flags[seqGenOptions::root]) 
		delete simulation_output_streams.at(seqGenOptions::root);
	if (outfile_flags[seqGenOptions::ma]) 
		delete simulation_output_streams.at(seqGenOptions::ma);
	if (outfile_flags[seqGenOptions::verb]) 
		delete simulation_output_streams.at(seqGenOptions::verb);
	if (outfile_flags[seqGenOptions::trace]) {
		//////////
		/// If there are no inputTrees, this is a path proposal. Close trace stream or events are underreported.
		//////////
		for (list<inTree*>::iterator xt = inputTrees.begin(); xt != inputTrees.end(); ++xt, x++)
			*simulation_output_streams.at(seqGenOptions::trace) << "GUIDE_TREE_PARTITION" << x << "=" << (*xt)->tree << ";" << endl;
		delete simulation_output_streams.at(seqGenOptions::trace);
	}
}
