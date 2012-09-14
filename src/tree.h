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
// Max size
#define MAX_NAME_LEN 256
#define MAX_INDEL_SIZE 100
#define RATE_AWAY_BUFFER 1.5
#define MAX_RATE_CATS 32
// Action related
#define NO_ACTION 0
#define INSERT 1
#define DELETE 2
#define SUBSTITUTION 4
#define FOSSIL 8
#define BRANCH_BEGIN 16
#define BRANCH_END 32
#define NO_EVENT 64
// Site constraint related
#define NO_CONSTRAINT 0
#define INVARIABLE 1
#define NO_INDEL 2
#define INVAR_AND_NOINDEL 3
#define EXACT 0
#define NODE_IS_ANC 1
#define NODE_IS_DES 2
#define NO_RELATION 3
// Function
#define MIN(a,b) ( (a<b) ? a : b )
#define MAX(a,b) ( (a<b) ? b : a )
#define from(x) (x*numStates)
#define to(x)   (x)
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
class globalArray;
class insertSite;
class siteModifier;
class RateMatrix;
class Likelihood;
class Dependency;
class eventTrack;
//class setRates;
class setRates;

//////////
/// Variables used to track data allocation.
//////////
extern int num_copy_ctor;
extern int num_ctor;
extern int num_dtor;
extern int numStates;
extern bool order_3_markov;
extern bool test_alternate_representations;
extern string stateCharacters;

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
	double atEpochTime;
	vector<bool> bipartition;
	string clade_label;
	string node_name;
	inClade *nodeEnv;
	Branch *branch;
	int mytipNo;
	int ancestorNo;
	setRates *Qptr;
//	RateMatrix (*ptr2init)(TTree*, TNode*, double, double, int);
//	void (*ptr2update)(RateMatrix*, vector<Site>::iterator, double, double);

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
	Sequence 		*evolvingSequence;
	double			rate_away_site_width;
	//////////

	//////////
	/// Functions:
	//////////
	// Constructor/Destructor
	TNode(TTree *tree);	/// TTree only so that it can be immediately added to the nodeList. Not copy constructor.
	TNode();
	// Informational/Debugging
	void report_varSites();
	void Print_Active_Sites();
	void FULL_REPORT();
	void report(TTree *tree = NULL);
	void report_on_sites();
	// Other
	bool varSiteOverlap();
	bool isVarSite(varSite *chk_varSite);
	void constructUnconstrainedSequence(string& root_sequence);
	varSite *Generic_varSite();
	void setSitePointers(string fromfunction);
	varSite *checkSitePointer(varSite *ptr, string message);
	void Inherit_varSites();
	void addGeneral_varSites();
	void InheritMotifSites();
	void clearMotifSites();
	void Site_postProcess(); 	// A routine to set all tracking pointers for sequence item.
	void dumpBetweenSitePointers(string& root_sequence);
	void setInvariableArrayPos();
	void Remove_Objects(bool root_node = false);
	void Remove_iz_Objects();
	void Remove_varSites();
	void Print_Active_Properties();
	void Print_Substitution_Properties();
	void setRateAway(int step_type);
	void inheritCategories(TNode *anc);
	void setStateLikelihood(size_t position, int category, double branch_length_scalar);
	vector<double> setCatLikelihood(size_t position, int category);
	void printGammaCategories();
	void setBranchTransitionProbabilities();
	double calculateForwardRateAwayFromSequence();
	double calculateForwardRateAwayFromSequence__order3Markov(TTree *tree, int event_site);
	double forwardRateAwayFromSequence();
	void updateSequence(TNode *update_node);
	bool EndpointCheck();
	bool isDescendant(TTree *tree, TNode *ancestor);
	string subtreePattern(int position);
	string printSequence(bool print_nucleotide_sequence = true);
	string printBipartition();
	string output_sequence(string::iterator start, string::iterator end);
	size_t findFirstConstrained(size_t fromSite);
	size_t findForwardNumberOfConstrained(size_t fromSite, size_t numSites);
	size_t findBackwardNumberOfConstrained(size_t fromSite, size_t numSites);
	int locateSite(double *rate_away_locator, int simulation_type);
	double calculateEndpointRateAwayFromSequence(TTree *tree, TNode *des, double T, double at_dt, int event_site);
	void calculateEndpointRateAwayFromSite(vector<Site>::iterator evolver_it, vector<Site>::iterator target_it, double *cumulative_site_sum_away, double *forward_site_sum_away, RateMatrix *rate_matrix);
	double TauIJ (TTree *tree);
	double TauIJ2 (TTree *tree, unsigned int position, unsigned int end_position);
	double TauIJ3 (TTree *tree, int env_index, int i_seq_index, short residue_i, short residue_j);
	void set_site_window(int order, int *start_site, unsigned int *end_site);
	void site_specific_Qmat(TTree *tree, int start_site, unsigned int end_site);
	RateMatrix initialize_rate_matrices(TTree *tree, TNode *des, double T, double at_dt, int event_site);
	void update_rate_matrices(RateMatrix *rates, vector<Site>::iterator site, double T,	double at_dt);
	double Rij (TTree *tree, unsigned int position, unsigned int end_position);
	void resetSequence(TNode *node);
	void printForwardRateAway();
	void write_sequence(ofstream& out);
	//void set_Qptr(bool Qd, bool Pc, bool Nij);
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
	RateMatrix *rates;		/// Each branch holds it's own rate matrix (inherits from ancestral branches)
	double S;		/// For EPC, value to scale independent transition probabilities.
	//////////
	/// n_{ij}(t) = # of sites where sequence at time t has nucleotide i and sequence at time T has nucleotide type j.
	/// n_{i.}(t) = \sum_{j={A,C,G,T}} n_{ij}(t)
	/// \tilde{P}_{ij}(T-t)=\frac{n_{ij}(t)}{n_{i.}(t)}
	///
	/// Pseudocounts may be necessary, since very short branches will have issues.
	/// \rho_{ij} = pseudocounts where sequence at time t has i and sequence at time T has j
	/// \rho_{i.} = \sum_j \rho_{ij}
	/// \tilde{\tilde{P}}_{ij}(t) = \frac{n_{ij}(t)+\rho_{ij}}{n_{i.}(t)+\rho_{i.}}
	//////////
	vector<double> nij;	/// For EPC again, transition probability vector approximations on a branch.
	vector<double> nij_pseudocounts;
	
	Branch() : 
		length0(0.0), 
		length1(0.0), 
		length2(0.0), 
		param(0.0),
		branch1_time_relative_length(0.0), 
		branch2_time_relative_length(0.0), 
		branch0_time_relative_length(0.0),
		branch1_max_path(0.0), 
		branch2_max_path(0.0),
		perturbation(1.0),
		rates ( NULL ),
		S (1.0)
    { }
	double anc_length() { return length0 * perturbation; }
	double des_length(int des) { return ( (des==1) ? length1 : length2 ) * perturbation; }
	void update_nij(short from_state, short to_state, short target_state);
	void print_nij(bool absolute = true);
	void report();
};

class TTree : private Counter<TTree>
{
public:
	using Counter<TTree>::howMany;

	// Variables
	int rooted, lengths;
	TNode *root; 
	inTree *my_iTree;
	list<TNode*> nodeList;
	int numTips, numNodes;
	double totalLength;
	vector<string> names;
	vector<TNode*> tips;
	int root_numSites;
	double epoch_step_size;
	double trs_example_step_size;
	list<inClade*> treeEnv;
	bool epc_root_is_set;
	list<Dependency*> dep;
	list<Dependency*> neutral_dep;
	int treeNum;		// The partition that this tree is assigned to.
	
	// Global arrays
	globalArray				*global_alignment;

	// Functions:
	TNode *ReadTip(string& tree_str, int *pos, char ch, TTree *tree, bool first_partition);
	TNode *ReadNode(string& tree_str, int *pos, TTree *tree, int detectPolytomies, bool first_partition);
	TNode *ReadBranch(string& tree_str, int *pos, TTree *tree, bool first_partition);
	TNode *AddNode();
	TNode *copyNode(TNode *node);
	inClade *AddClade(string& label);
	string sitePattern(int position);
	void setStateLikelihoods(double branch_length_scalar);
	void setCatLikelihoods(double branch_length_scalar);
	void constructRootMotifSites(list<inMotif*>& motif_specs, string& root_sequence, string& invariable);
	void ReadTree(string& str, int *pos, TTree *tree, bool first_partition);
	void InitTree();
	void report_clades();
	void Place_Motif(list<siteRegEx*>& motif, char mark, list<inMotif*>& placed_motifs, inMotif* thisMotif, string& root_sequence, size_t start_site);
	void report_branches();
	void report_sequences();
	void calculateJCLikelihoods();
	void setTransitionProbabilities();
	void sample_root_sequence();
	list<inMotif*> DrawMotifs(string& root_sequence, double proportion_motifs);
	TTree();
	void write_tree(ofstream& out, bool scale);
	void write_subtree(TNode *node, ofstream& out, bool writeAncestors, int *nodeNo, bool time_rel, bool scale);
};

class globalArray : private Counter<globalArray>
{
public:
	using Counter<globalArray>::howMany;
	vector<insertSite> insert_sites;

	void Print();
	vector<insertSite>::iterator locateEvent(TNode *des, size_t position);
	bool isDeleted(TNode *des, insertSite insert_site);
	void cleanArray( list<eventTrack*>& events );
};

class siteModifier : private Counter<siteModifier>
{
public:
	using 	Counter<siteModifier>::howMany;
	char	action;
	int		fromAnc;
	int		indelNo;
	bool	save;

	siteModifier(char a, int fA, int iN)
		: action(a),
		  fromAnc(fA),
		  indelNo(iN),
		  save(false)
	{ }
};

class insertSite : private Counter<insertSite>
{
public:
	using 	Counter<insertSite>::howMany;
	char	action;
	int		fromAnc;
	int		indelNo;
	
	list<siteModifier> modifiers;

	insertSite(char a, int fA, int iN)
		: action(a),
		  fromAnc(fA),
		  indelNo(iN)
	{ modifiers.clear(); }

};

#endif /* _TREE_H_ */

