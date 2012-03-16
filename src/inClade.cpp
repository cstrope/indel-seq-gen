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

#include "inClade.h"

inClade::inClade(string& label, TTree *tree) 
{
	for (int i = 0; i < NUM_CLADE_VARIABLES; i++) 
		set_in_clade[i] = false;
	model = tree->treeEnv.front()->model;
	environment_name.assign(label);
	clade_name.assign(label);

	invariableSites=false;
	insert_lengthDistribution.clear();
	delete_lengthDistribution.clear();
	values2Export2Freq.clear();
	P_ins_=P_del_=maxIndel=0;
	indelFlag=false;
	proportion_invariable = 0;
	processed = false;
	rateHetero = NoRates;
	gammaShape = 1.0;
	numCats = 1;
	bipartition.clear();
	constraintChange = NO_CHANGE;
	for (int i = 0; i < MAX_RATE_CATS; i++)
		catRate[i] = 0;

	my_motifs.clear();
}

inClade::inClade(string label, const seqGenOptions *options) 
{
	clade_name.assign(label);
	environment_name.assign(label);
	P_ins_ = P_del_ = maxIndel = 0;
	indelFlag = false;
	model = options->global_model;
	bipartition.clear();
	values2Export2Freq.clear();
	for (int i = 0; i < numStates; i++) {
		if(isNucModel)
			values2Export2Freq.push_back(nucFreq[i]);
		else
			values2Export2Freq.push_back(aaFreq[i]);
	}
	insert_lengthDistribution.clear();
	delete_lengthDistribution.clear();
	proportion_invariable = options->default_proportion_invariable;	
	rateHetero=options->default_rateHetero;
	constraintChange = NO_CHANGE;
	for (int i = 0; i < MAX_RATE_CATS; i++)
		catRate[i] = options->default_catRate[i];
	for (int i = 0; i < numStates; i++) {
		values2Export2Freq.at(i) = ((numStates == 4) ? nucFreq[i] : aaFreq[i]);
	}
	if (rateHetero == DiscreteGammaRates) numCats = options->num_discrete_gamma_categories;
	if (rateHetero == GammaRates || rateHetero == DiscreteGammaRates) gammaShape = options->alpha;
	else gammaShape = 1.0;
	
	my_motifs.clear();
}

void inClade::rootEnvSetup(inClade *env) 
{
	processed = true;		// Global params for tree are always processed.
	environment_name.assign("tree_parameters");
	clade_name.assign("root_clade");
	P_ins_   = env->P_ins_;
	P_del_   = env->P_del_;
	maxIndel = env->maxIndel;
	model = env->model;
	insert_lengthDistribution.clear();
	delete_lengthDistribution.clear();
	values2Export2Freq.clear();
	values2Export2Freq = env->values2Export2Freq;
	proportion_invariable = env->proportion_invariable;
	rateHetero = env->rateHetero;
	gammaShape = env->gammaShape;
	if (rateHetero == DiscreteGammaRates || rateHetero == CodonRates)
		for (int i = 0; i < MAX_RATE_CATS; i++)
			catRate[i] = env->catRate[i];
	for (int i = 0; i < NUM_CLADE_VARIABLES; i++)
		set_in_clade[i] = true;
}

void inClade::Print_Environment() 
{
	if ( !environment_name.empty() ) cout << endl << "Environment name: " << environment_name << ", ";
	cout << "Clade: " << clade_name << "  Top affected node: ";
	for (vector<bool>::iterator it = bipartition.begin(); it != bipartition.end(); it++)
		cout << (*it);
	cout << endl;
	if (constraintChange) {
		cout << ((constraintChange == PSEUDOGENE)? " PSEUDOGENE" : " NEOFUNCTIONALIZATION");
		cout << endl;
	}
	cout << " INDEL PARAMS ";
	if (indelFlag) {
		cout << endl;
		cout << "   indelFlag = " << indelFlag << endl;
		cout << "   Pins = " << P_ins_ << endl;
		cout << "   Pdel = " << P_del_ << endl;
		cout << "   maxIndel = " << maxIndel << endl;
		if ( !(insert_lengthDistribution.empty()) ) {
			cout << "   inLD = " << endl << "    ";
			for (int i = 1; i <= maxIndel; i++) 
				cout << i << ":" << insert_lengthDistribution.at(i) << " ";
			cout << endl;
		}
		if ( !(delete_lengthDistribution.empty()) ) {
			cout << "   delLD = " << endl << "    ";
			for (int i = 1; i <= maxIndel; i++) 
				cout << i << ":" << delete_lengthDistribution.at(i) << " ";
			cout << endl;
		}
	} else cout << "[NONE]" << endl;
	cout << " SUBSTITUTION PARAMS " << endl;
	cout << "   model = " << modelNames[model] << endl;
	cout << "   Frequencies = ";
	if( !(values2Export2Freq.empty()) ) {
		for (int i = 0; i < numStates; i++)
			cout << stateCharacters[i] << ":" << values2Export2Freq.at(i) << " ";
	} else cout << " Unchanged from globals.";
	cout << endl;
	cout << "   Invariable = " << proportion_invariable << endl;
	cout << "   invariableSites = " << ((invariableSites) ? "true" : "false") << endl;
	cout << "   Rate Heterogeneity = ";
	switch (rateHetero) {
	case GammaRates:
		cout << "Gamma Rates, alpha = " << gammaShape;
		break;
	case DiscreteGammaRates:
		cout << "Discrete Gamma Rates, alpha = " << gammaShape << ", " << numCats << " categories.";
		break;
	case CodonRates:
		cout << "Codon Rates[pos] = ";
		for (int i = 0; i < 3; i++)
			cout << i+1 << ": " << catRate[i] << " ";
		break;
	case NoRates:
		cout << "Uniform Rates";
		break;
	}
	cout << endl;
		
	if ( !my_motifs.empty() ) {
		cerr << " MOTIFS: " << endl;
		for (list<inMotif*>::iterator it = my_motifs.begin(); it != my_motifs.end(); it++)
			(*it)->report();
	}
}

string inClade::report_rateHetero()
{
	switch (rateHetero) {
	case 0: return "Uniform Rates"; break;
	case 1: return "Codon Rates"; break;
	case 2: return "Gamma Rates"; break;
	case 3: return "Discrete Gamma Rates"; break;
	default: 
		cerr << "Invalid rateHetero assignment? How did you do this?" << endl;
		exit(EXIT_FAILURE);
		break;
	}	
}

string CladeName(TTree *tree, int which_clade_in_list) 
{
	string cladename = "";
	int clade_ctr = 0;

	for (list<inClade*>::iterator it = tree->treeEnv.begin(); it != tree->treeEnv.end(); it++, clade_ctr++)
		if (clade_ctr == which_clade_in_list) cladename.assign((*it)->clade_name);
		
	return cladename;
}

inClade *FindEnvironment(TTree *tree, string which_clade) 
{
	for (list<inClade*>::iterator it = tree->treeEnv.begin(); it != tree->treeEnv.end(); it++) {
		if ( which_clade.compare((*it)->clade_name) == 0) {
			return (*it);
		}
	}
	return NULL;
}

// Function to set clade params to be clade-specific or to point to tree parameters (if not)
int CladeUnprocessed(TTree *tree) 
{
	int unprocessed_clade = 0;
	int num_cycles = 0;

	for (list<inClade*>::iterator it = tree->treeEnv.begin(); it != tree->treeEnv.end(); it++) {
		if ((*it)->processed) unprocessed_clade++;
		else break;
		num_cycles++;
	}

	return ((num_cycles == unprocessed_clade) ? 0 : unprocessed_clade);
}


