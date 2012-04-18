#ifndef _RATE_TYPE_H_
#define _RATE_TYPE_H_

using namespace std;

#include "tree.h"
#include "model.h"
#include "motif.h"

class setRates : private Counter<setRates>
{
public:
	using Counter<setRates>::howMany;

	void set_Qptr(bool Qd, bool Pc, bool Nij);

	// Pointers to functions for initializing and updating rates. //
	
	setRates() { }
};

extern RateMatrix (*ptr2init)(TTree*, TNode*, TNode*, double, double, int);
extern void (*ptr2update)(RateMatrix*, TNode*, vector<Site>::iterator, double, double);

RateMatrix iQdPc(TTree *tree, TNode *i_z, TNode *k_0, double T, double at_dt, int event_site);
RateMatrix iQdP (TTree *tree, TNode *i_z, TNode *k_0, double T, double at_dt, int event_site);
RateMatrix iQPc (TTree *tree, TNode *i_z, TNode *k_0, double T, double at_dt, int event_site);
RateMatrix iQP  (TTree *tree, TNode *i_z, TNode *k_0, double T, double at_dt, int event_site);
RateMatrix iQN  (TTree *tree, TNode *i_z, TNode *k_0, double T, double at_dt, int event_site);

void uQdPc(RateMatrix *rates, TNode *i_z, vector<Site>::iterator site, double T, double at_dt);
void uQdP (RateMatrix *rates, TNode *i_z, vector<Site>::iterator site, double T, double at_dt);
void uQPc (RateMatrix *rates, TNode *i_z, vector<Site>::iterator site, double T, double at_dt);
void uQP  (RateMatrix *rates, TNode *i_z, vector<Site>::iterator site, double T, double at_dt);
void uQN  (RateMatrix *rates, TNode *i_z, vector<Site>::iterator site, double T, double at_dt);

#endif
