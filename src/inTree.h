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

#ifndef INTREE_H_
#define INTREE_H_

//////////
/// Definitions
//////////
#define RANDOM 0
#define SINGLE_ROOT 1
#define MULTIPLE_ALIGNMENT_ROOT 2
#define DIST_CALC 0
#define TIME_SCALE 1
#define PREC_EPS 0.0001
#define ROOT_FLAG -2
#define SET 1
#define CHECK 2
#define PRECISION 7

#include <string>
#include <list>
#include <vector>
#include <exception>
#include <map>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <list>
#include <cmath>
#include "aamodels.h"
#include "evolve.h"
#include "inClade.h"
#include "model.h"
#include "motif.h"
#include "seqGenOptions.h"
#include "tree.h"
#include "treefile.h"

using namespace std;

// Forward definitions
class inClade;
class inMotif;

class inTree : private Counter<inTree>
{
public:
	using Counter<inTree>::howMany;
	// Variables
	bool randomInvariableAssignment;	
	bool scaleBranches;
	bool scaleTree;
	char ma_columnCollapseMethod;
	double branchScale;
	double global_max_path_length;
	double proportion_invariable;
	double partitionRate;
	int ma_numSeqsToUse;
	int rootSeqType;
	int partitionLength;
	int treeNum;
	list<inMotif*> motif_specs;
	short codon_offset;
	string filename;
	string invariable;
	string invariable_by_partition;
	string label;
	string ma_range;
	string rootseq_by_partition;
	string tree;
	TTree *my_tree;
	vector<string> input_MA;
	vector<string> input_MA_motifs;
	vector<string> motifs_by_partition;

	// Functions
	bool CheckPhylogeneticAncestralNodes(string& node_names, int action);
	double CalculatePathLengths();
	double SearchPath(TNode *anc, TNode *des, int flag);
	double SearchMinBranch(TNode *node, double curr_min);
	void Print_Trace_Header(stringstream& header, vector<string> otu_names, seqGenOptions* options);
	void Print_Branch(TNode *node, stringstream& header, vector<string>::iterator& taxon_names, int& num, seqGenOptions* options);
	double DoNudge(TNode *node, double step_size, int simulation_type);
	int  getRootSeq_and_Motif(int type, seqGenOptions *options);
	inMotif *FindMotif(char marker);
	inTree(const seqGenOptions *, int, inClade *);
	void Read_MA(string filename, seqGenOptions *options);
	string addFilePath(string suffix, string prefix);
	string tree_only(string input);
	void perform_checks(string, seqGenOptions *);
	void NudgeBranches(double step_size, int simulation_type);
	void ConvertRootSequence(TTree *tree);
	void TimeScale_Tree(double max_path, TNode *node);
	void BranchScale_Tree(TNode *anc, TNode *des, double scale);
	void calcDistancesFromRoot();
	void calcDist(TNode *anc, TNode *des);
	void perturbTree(double perturb_value);
	void setPerturb(TNode *node, double perturb_value);
	void Print_Time_Rel_Tree(ostream& out, bool writeAncestors);
	void PrintSubtree(TTree *tree, TNode *node, ostream& out, bool writeAncestors, int *nodeNo, bool time_rel, bool scale);
	void Print_Tree(ostream& out, bool scale);
	void Setup_Tree(seqGenOptions *options);
	void Define_Ancestors();
	void Define_Bipartitions(vector<string>& names);
	void Define_Clade_Bipartitions();
	void Push_Bipartition(TNode *node, int which_tip);
	void Define_Anc(TNode *anc, TNode *des, int *node_num);
	void CPAN(TNode *node, string& node_names, int action, bool *inval);
	void SortTaxonTips(vector<string>& otu_names);
	void parseTree(string, seqGenOptions *);
	void parseLineages(seqGenOptions *options);
	void parseClade(string& clade_string, seqGenOptions *options);
	void parseMotif(string& motif_string, seqGenOptions *options);
	void parseIndel(string& indel_info, inClade *clade, seqGenOptions *options);
	void parsePound(string& subtree_info, inClade *clade, bool local, seqGenOptions *options);
	void propagateTreePointers();
	void propagateCladePointers(TNode *node);
	void setupTreeStructs(int partitionLength, seqGenOptions *options);
	void Reset_inMotif();
	void setMultiStateCharacters();
	void apply_branchScales(double scalar);
	double apply2Subtree(TNode *node, double scalar);

	private:
};

list<string> split(string str, string delim, bool mixDelim = false, bool merge = false);
string& trim(string& str);
string& Remove_Whitespace(string& str);

#endif /*INTREE_H_*/
