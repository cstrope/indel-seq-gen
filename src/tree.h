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

#ifndef _TREE_H_
#define _TREE_H_

//////////
/// Definitions
//////////
//
#define MAX_NAME_LEN 256
#define MAX_INDEL_SIZE 100
#define INSERT 0
#define DELETE 1
#define traceEvents (options->output_file_flags_unset || options->output_file_flags[trace])

//////////

#include <string>
#include <cstring>
#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include <limits>
#include "gamma.h"

using namespace std;

// Forward definitions?
class TTree;
class TNode;
class inClade;
class seqGenOptions;
class motifSite;
class varSite;
class inTree;
class siteRegEx;
class inMotif;
class Branch;
class Sequence;
class Site;
class globalAlignment;

//////////
/// Variables used primarily to track data used in Strope09.
//////////
extern int total_motif_positions;
extern int correct_motif_positions;
extern int clone_iSGv1;
extern int num_template_violations;
extern int num_copy_ctor;
extern int num_ctor;
extern int num_dtor;

//////////
/// Enumerations
//////////
enum {
	codon11,
	codon31,
	protein,
	non_coding
};

//////////
/// rateHetero types.
//////////
enum {
	NoRates,
	CodonRates,
	GammaRates,
	DiscreteGammaRates
};

void crash(string message);

//////////
/// Counter for class instantiations (one for each class).
/// From: http://www.drdobbs.com/184403484;jsessionid=GEBN5SESG1APTQE1GHPCKHWATMY32JVN?pgno=1
//////////
template<typename T>
class Counter 
{
public:
	Counter() { ++count; num_ctor++; }
	Counter(const Counter&) { ++count; num_copy_ctor++; }
	~Counter() { --count; num_dtor++; }
	
	static size_t howMany() { return count; }
	
private:
	static size_t count;
};
template<typename T> 
size_t Counter<T>::count = 0;	// Set all counters to 0.


class TNode : private Counter<TNode> 
{
public:
	using Counter<TNode>::howMany;
	//////////
	/// Variables that obviously belong to the TNode class.
	//////////
	TNode *branch0, *branch1, *branch2, *anc;
	double trDistanceFromRoot;
	double DistanceFromRoot;
	int tipNo;
	int numEventsToSimulate;
	double BL_step_size;
	double atEpochTime;
	vector<bool> bipartition;
	string clade_label;
	inClade *nodeEnv;
	Branch *branch;
	int mytipNo;

	//////////
	/// Should these belong to the TNode... Only time (and the amount of effort needed to move 
	/// them) will tell.
	//////////
	varSite *one_site_varSite;
	varSite *unconstrained_varSite;
	list<varSite*> variable_region_list;

	//////////
	/// These variables more appropriately belong to a Sequence class.
	//////////
	// * Class for the sequence data
	Sequence *evolvingSequence;
	//////////

	//////////
	/// Functions:
	//////////
	// Constructor/Destructor
	TNode(TTree *tree);
	// Informational/Debugging
	void report_varSites();
	void Print_Active_Sites();
	void FULL_REPORT();
	void report(TTree *tree = NULL);
	void report_on_sites();
	// Other
	bool varSiteOverlap();
	bool isVarSite(varSite *chk_varSite);
	void constructUnconstrainedSequence(string &root_sequence);
	varSite *Generic_varSite();
	void calcMotifAccuracy();
	void setSitePointers();
	varSite *checkSitePointer(varSite *ptr, string message);
	void Inherit_varSites();
	void addGeneral_varSites();
	void InheritMotifSites();
	void clearMotifSites();
	void Site_postProcess(); 	// A routine to set all tracking pointers for sequence item.
	void dumpBetweenSitePointers(string root_sequence);
	void setInvariableArrayPos();
	void Remove_Objects(bool root_node = false);
	void Remove_varSites();
	void Print_Active_Properties();
	string output_sequence(string::iterator start, string::iterator end);
	size_t findFirstConstrained(size_t fromSite);
	size_t findForwardNumberOfConstrained(size_t fromSite, size_t numSites);
	size_t findBackwardNumberOfConstrained(size_t fromSite, size_t numSites);
	//////////
};

class Branch : private Counter<Branch>
{
public:
	using Counter<Branch>::howMany;
	double length0, length1, length2, param;
	double branch1_time_relative_length, branch2_time_relative_length, branch0_time_relative_length;
	double branch1_max_path, branch2_max_path;
	double perturbation;
	
	Branch() : length0(0.0), length1(0.0), length2(0.0), param(0.0),
			   branch1_time_relative_length(0.0), branch2_time_relative_length(0.0), branch0_time_relative_length(0.0),
			   branch1_max_path(0.0), branch2_max_path(0.0),
			   perturbation(1.0) 
    { }
	double anc_length() { return length0 * perturbation; }
	double des_length(int des) { return ( (des==1) ? length1 : length2 ) * perturbation; }

};

class TTree : private Counter<TTree>
{
public:
	using Counter<TTree>::howMany;
	// Variables
	int rooted, lengths;
	TNode *root; 
	list<TNode*> nodeList;
	int numTips, numNodes;
	double totalLength;
	vector<string> names;
	vector<TNode*> tips;
	int root_numSites;
	double epoch_step_size;
	double trs_example_step_size;
	list<inClade*> treeEnv;

	// Global arrays
	vector<globalAlignment> global_arrays;

	// Functions:
	TNode *ReadTip(string tree_str, int *pos, char ch, TTree *tree, bool first_partition);
	TNode *ReadNode(string tree_str, int *pos, TTree *tree, int detectPolytomies, bool first_partition);
	TNode *ReadBranch(string tree_str, int *pos, TTree *tree, bool first_partition);
	TNode *AddNode();
	inClade *AddClade(string& label);
	void constructRootMotifSites(list<inMotif*>& motif_specs, string &root_sequence, string &invariable);
	void ReadTree(string str, int *pos, TTree *tree, bool first_partition);
	void InitTree();
	void report_clades();
	void Place_Motif(list<siteRegEx*>& motif, char mark, list<inMotif*>& placed_motifs, inMotif* thisMotif, string &root_sequence, size_t start_site);
	void report_branches();

	list<inMotif*> DrawMotifs(string& root_sequence, double proportion_motifs);
	TTree();
};

class globalAlignment : private Counter<globalAlignment>
{
public:
	char	action;
	int		fromAnc;
	int		indelNo;

	globalAlignment(char a, int fA, int iN)
		: action(a),
		  fromAnc(fA),
		  indelNo(iN)
	{ }
};

#endif /* _TREE_H_ */

