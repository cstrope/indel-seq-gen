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

#ifndef _INCLADE_H_
#define _INCLADE_H_

//////////
/// Definitions
//////////
#define NUM_CLADE_VARIABLES 32
#define NO_CHANGE 0
#define PSEUDOGENE 1
#define NEOFUNCTIONALIZATION 2

#include <string>
#include <cstring>
#include <list>
#include <vector>
#include <iostream>
#include "evolve.h"
#include "tree.h"
#include "inTree.h"
#include "twister.h"
#include "motif.h"

using namespace std;

class inMotif;
class varSite;

// Positions in set_in_clade;
enum {
	__P_ins___,
	__P_del___,
	__insert_lengthDistribution__,
	__delete_lengthDistribution__,
	__maxIndel__,
	__model__,
	__indelFlag__,
	__invariableSites__,
	__proportion_invariable__,
	__values2Export2Freq__,
	__rateHetero__,
	__constraintChange__
};

// Forward declarations.
class TTree;
class TNode;
class inTree;

class inClade : private Counter<inClade>
{
public:
	using Counter<inClade>::howMany;
	string environment_name;
	string clade_name;
	double P_ins_, P_del_;
	int maxIndel;
	vector<double> insert_lengthDistribution;
	vector<double> delete_lengthDistribution;
	int model;
	int rateHetero;
	bool indelFlag;
	bool processed;		// When enumerating clades in tree, all clades must be processed. If not, error.
	bool set_in_clade[NUM_CLADE_VARIABLES];
	bool invariableSites;
	double proportion_invariable;
	vector<double> values2Export2Freq;
	int numCats;
	double catRate[MAX_RATE_CATS];
	double gammaShape;
	short constraintChange;		// PSEUDOGENE, NEOFUNCTIONALIZATION, NO_CHANGE
	list<inMotif*> my_motifs;
	vector<bool> bipartition;

	// Functions
	inClade(string& label, TTree *tree);
	varSite *Generic_varSite(bool isTemplate);
	void rootEnvSetup(inClade *env);
	inClade(string label, const seqGenOptions *options);
	void Print_Environment();
	string report_rateHetero();
};

inClade *FindEnvironment(TTree *tree, string which_clade);
int CladeUnprocessed(TTree *tree);
string CladeName(TTree *tree, int which_clade_in_list);

#endif
