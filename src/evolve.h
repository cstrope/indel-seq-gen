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

/*  Header file for evolve.c                               */
#ifndef _EVOLVE_H_
#define _EVOLVE_H_

//////////
/// Definitions
//////////
#define MAX_RATE_CATS 32
#define NO_ACTION 2
#define NO_CONSTRAINT 0
#define INVARIABLE 1
#define NO_INDEL 2
#define INVAR_AND_NOINDEL 3
#define EXACT 0
#define NODE_IS_ANC 1
#define NODE_IS_DES 2
#define NO_RELATION 3
#define MIN(a,b) ( (a<b) ? a : b )

//////////
/// WARNING: seq_evo is difficult to see, but it is a lot easier to make this than type such a
/// long chain of pointers. seq_evo is primarily (maybe only) used in Insert, Delete, and Mutate.
//////////
#define seq_evo evolvingSequence->evolutionaryAttributes

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
#include <string>
#include <iomanip>
#include <memory>

#include "aamodels.h"
#include "gamma.h"
#include "insert_freqs.h"
#include "model.h"
#include "motif.h"
#include "random.h"
#include "trace.h"
#include "tree.h"
#include "treefile.h"
#include "twister.h"

using namespace std;

// Forward Declarations
class inTree;
class seqGenOptions;

extern int writeAncestors;
extern int indel_fill_aas;
extern int indelNo;
extern bool output_trs_example;

/* prototypes */
void SetCategories(TNode *node, int inNumSites, seqGenOptions *options);
char SetState(vector<double>& P, string& caller);
void SetNucSequence(char *seq, char *source, int inFromSite, int inNumSites);
void SetAASequence(char *seq, char *source, int inFromSite, int inNumSites);

void CreateMatrix();
void DestroyMatrix();
void CreateRates(TNode *node, int inNumSites);
void Create_Global_Arrays(TTree *tree, int inNumSites);
void InheritCategories(TNode *node, int inNumSites);

void EvolveSequences(inTree *iTree, list<eventTrack*> *events, seqGenOptions *options);
void EvolveNode(inTree *iTree, TNode *anc, TNode *des, double scale, double invar_scale, int inNumSites, list<eventTrack*> *events, seqGenOptions *options);

bool Stop_Codon(int *codon);
short IsInvariable(double proportionInvariable);
void RandomSequence(char *seq, int inNumSites, seqGenOptions *options);
void MutateSequence(inTree *iTree, TNode *des, double len, double motif_correction, list<eventTrack*> *events);
void initializeCladeParameters(TTree *tree, TNode *node);
void initializeBranchRun(inTree *iTree, TNode *node, double *en_invar_scale, int inNumSites, seqGenOptions *options);
void initializeTreeScalingParameters(inTree *iTree, double *scale, double *invar_scale);
void initializeMotifPositions(TNode *des);
void InDel(TTree *tree, TNode *des, double len, list<eventTrack*> *events, int CnB_divisor, seqGenOptions *options);
void gillespieIndel(inTree *iTree, TNode *des, double indel_len, list<eventTrack*> *events, seqGenOptions *options);
int  Find_Indel_Size(vector<double>& dist,int max_indel);
int  Insert(TTree *tree,TNode *des, int indel_size);
int  Delete(TTree *tree,TNode *des, int *indel_size);
int  Find_Start_Location(TNode *des, int *indel_size, int final_aa, int action);
int  Find_Action(double len, double P_insert_, double P_delete_, int, int);
int  isAnc(TNode *thisNode, int anc);
int  isSpot(TTree *tree,int *pos, TNode *curr_taxa);
int  calcTrials(TNode *des, int *E_in_, int *E_del_);
string insertFillSequence(int inNumSites);
int testRelation(vector<bool>& n1, vector<bool>& n2);
void Print_Globals(TTree *tree);
void calculateMotifScale(inTree *iTree, TNode *des, seqGenOptions *options, double subst_len, double *motif_scale);
double CalculateMotifScale(TNode *des, int position, double subst_len, bitset<20>& motif_characters);
void Find_Motif_Positions(TTree *tree, TNode *des, motifSite *this_site, int indel_size, bool back, int start_aa, varSite **in_template_varSite, varSite **in_motif_varSite);
void Site_Find_Motif_Positions(TTree *tree, TNode *des, motifSite *this_site, bool back, int start_aa, varSite **in_template_varSite, varSite **in_motif_varSite);
double calcGillespieLambda(TNode *des, double *lambda_ins, int *ins_L, double *lambda_del, int *del_L);


#endif /* _EVOLVE_H_ */

