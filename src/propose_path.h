#ifndef _PROPOSE_PATH_H_
#define _PROPOSE_PATH_H_

using namespace std;

#include<string>
#include<vector>
#include<list>
#include <assert.h>
#include "tree.h"
#include "model.h"
#include "inClade.h"
#include "evolve.h"
#include "motif.h"
#include "inTree.h"
#include "stats.h"

class RateMatrix;
class Sequence;
class eventTrack;
class Taxon;
class TTree;
class TNode;
class Branch;

class PathProbability : private Counter<PathProbability>
{
public:
	using Counter<PathProbability>::howMany;

	PathProbability() : forward (0), epc (0) { }
	double ForwardProbability( list<eventTrack*> sorted_events_by_time, bool do_last_event = true );
	double EPCProbability( list<eventTrack*> sorted_events_by_time, bool do_last_event = true );
	double EPC_Step(double Qij_k__T__, double Qidot_k__T__, int Qi2k, double t_z, double dt, double T);
	double EPC_NoEvent(double Qidot_k__T__,	double dt);
	double Forward_Step(double Qij, double Qidot, double dt);
	double Forward_NoEvent(double Qidot, double dt);
	double IndependentForwardProbability(list<eventTrack*> events, short start_state, RateMatrix *rates);

private:
	double forward;
	double epc;
};

class PathProposal : private Counter<PathProposal>
{
public:
	using Counter<PathProposal>::howMany;
	vector<Taxon*> taxa;
	list<eventTrack*> epc_events;
	PathProbability P_path;

	PathProposal() 
	{ }

	void Print_Path_Events();
	void write_path(ofstream& stream);
	void setNodeSequences(TTree *tree);
	void setUnspecifiedNodeSequences(TTree *tree, int sequence_length);
	void parseEndpointInfile(string& filename);
	void Evolve(TTree *tree, list<eventTrack*> *events);
	void EvolveToEndPoint(TTree *tree, TNode *anc, TNode *des, list<eventTrack*> *events);
	void htmlize(size_t location);
	void EvolveStep(TTree *tree, TNode *anc, TNode *des, double t_0, double T, list<eventTrack*> *events);
	void displayPaths(TTree *tree);
	vector<Taxon*> parseFastaFile(string& fasta_data);
	void emulateForwardSimulation(TTree *tree, list<eventTrack*> *events);
	void setEventHistory(list<eventTrack*> *events,	string event_history_file);
	void emulateToEndPoint(TTree *tree, TNode *anc, TNode *des, list<eventTrack*> *events);
	void EmulateStep(TTree *tree, TNode *i_z, TNode *k_0, double t_0, double T, list<eventTrack*> *events);
	list<eventTrack*>::iterator setSequenceAtTimePoint(TTree *tree, TNode *node, double to_time, list<eventTrack*>::iterator et, list<eventTrack*> *events);
	double select_next_dt(TNode *i_z, TNode *k_0, double current_dt, double lambda_T, double T, int *num_diff, bool isEndpoint);
	double rasmus_select_next_dt(double lambda_T, double T);

	void remove_site_events(int site_sampled);
	void evolve_independent_path(TTree *tree, TNode *i_z, TNode *k_0, double t_0, double T, list<eventTrack*> *events, vector<Site>::iterator i_z_site, vector<Site>::iterator k_0_site, int event_site);
	void EvolveIndependentStep(TTree *tree,	TNode *i_z,	TNode *k_0,	double t_0,	double T, list<eventTrack*> *events);
	list<eventTrack*> extract_site_events(int site);

private:
};

class Taxon : private Counter<Taxon>
{
public:
	using Counter<Taxon>::howMany;
	string taxon_name;
	string taxon_sequence;

	Taxon() 
		: taxon_name (""),
		  taxon_sequence ("")
	{  }
	

	
};

/// Function declaration
vector<bool> make_vector_bool (string str);
list<eventTrack*> getBranchEvents(list<eventTrack*>::iterator begin, list<eventTrack*>::iterator end);


#endif