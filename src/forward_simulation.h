#ifndef _FORWARD_SIMULATION_H_
#define _FORWARD_SIMULATION_H_

using namespace std;

#include <exception>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <iomanip>
#include <memory>
#include <string>
#include "tree.h"
#include "evolve.h"
#include "seqGenOptions.h"
#include "trace.h"
#include "inTree.h"

extern int num_subst;

//////////
/// Functions
//////////
template <class T> string to_string (const T& t);

class ForwardSimulation : private Counter<ForwardSimulation>
{
public:
	using Counter<ForwardSimulation>::howMany;

	ForwardSimulation() { }
	void EvolveNode(inTree *iTree, TNode *anc, TNode *des, int inNumSites, list<eventTrack*> *events, seqGenOptions *options);
	void EvolveSequences(inTree *iTree, list<eventTrack*> *events, seqGenOptions *options);
	double calcGillespieLambda(TTree *tree, TNode *des, double *I, double *D, double *S, int simulation_type, unsigned int event_site);
	void gillespie(inTree *iTree, TNode *des, double branch_len, list<eventTrack*> *events, seqGenOptions *options, int simulation_type, bool evolving_to_equilibrium);
	void evolve2Equilibrium(inTree *iTree, list<eventTrack*> *events, seqGenOptions *options);
private:

};

#endif