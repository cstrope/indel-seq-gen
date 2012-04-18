#include "rate_type.h"
#include "dependency.h"

RateMatrix
iQN(
	TTree *tree,
    TNode *i_z,
	TNode *k_0,
	double T,
	double at_dt,
	int event_site
   )
{
	i_z->setRateAway(TIME_RELATIVE_STEPS);	// Resets site_sum_away, allocates e_QijDt.
	RateMatrix initial_rates(*k_0->branch->rates);

	//////////
	/// Independent sites: Copy the global model in for each site.
	//////////
	for (vector<Site>::iterator it = k_0->seq_evo.begin(); it != k_0->seq_evo.end(); ++it) 
		(*it).e_QijDt = k_0->branch->rates->Qij;

	vector<double> nij_row_sum (numStates, 0);
	int i = 0;
	vector<double>::iterator jt = nij_row_sum.begin();
	for (vector<double>::iterator it = i_z->branch->nij.begin(); it != i_z->branch->nij.end(); ++it, ++i) {
		if (i == numStates) { ++jt; i = 0; }
		(*jt) += (*it);
	}

	vector<double>::iterator pt = initial_rates.Pij.at(0).begin();	// Using nij, cannot use categories... yet. //
	jt = nij_row_sum.begin();
	i = 0;
	for (vector<double>::iterator it = i_z->branch->nij.begin(); it != i_z->branch->nij.end(); ++it, ++pt, ++i) {
		if (i == numStates) { ++jt; i = 0; }
		(*pt) = (*it) / (*jt);
	}

	i_z->branch->rates->setPij(i_z->seq_evo.front(), i_z->branch->S*(T-at_dt), i_z->nodeEnv->rateHetero);	// setPij dependes on the root Matrix

	pt = initial_rates.Pij.at(0).begin();
	cerr << endl << "Transition Probabilities using nij." << endl;
	for (i = 0; i < numStates; ++i) {
		for (int j = 0; j < numStates; ++j) {
			cerr << i_z->branch->rates->Pij.at(0).at(i*numStates+j) << " (" << initial_rates.Pij.at(0).at(i*numStates+j) << ") " << endl;
		}
		cerr << endl;
	}

	return initial_rates;
}

RateMatrix
iQPc(
	 TTree *tree,
     TNode *i_z,
	 TNode *k_0,
	 double T,
	 double at_dt,
	 int event_site
	)
{
	i_z->setRateAway(TIME_RELATIVE_STEPS);	// Resets site_sum_away, allocates e_QijDt.
	RateMatrix initial_rates(*k_0->branch->rates);
	unsigned int end_site;

	//////////
	/// Dependent sites. Recalculate Q matrices for positions at the beginning of the branch or
	/// affected by a change for subsequent calculation of the site-specific transition matrices.
	//////////
	i_z->set_site_window(tree->dep.front()->context.return_order(), &event_site, &end_site);
	i_z->site_specific_Qmat(tree, event_site, end_site);

	return initial_rates;
}

RateMatrix
iQP(
	TTree *tree,
    TNode *i_z,
	TNode *k_0,
	double T,
	double at_dt,
	int event_site
   )
{
	i_z->setRateAway(TIME_RELATIVE_STEPS);	// Resets site_sum_away, allocates e_QijDt.
	RateMatrix initial_rates(*k_0->branch->rates);
	//////////
	/// Independent sites: Copy the global model in for each site.
	//////////
	for (vector<Site>::iterator it = k_0->seq_evo.begin(); it != k_0->seq_evo.end(); ++it) 
		(*it).e_QijDt = k_0->branch->rates->Qij;

	// setPij uses the Root and Cijk matrices, already set up for branch->rates
	k_0->branch->rates->setPij(k_0->seq_evo.front(), k_0->branch->S*(T-at_dt), k_0->nodeEnv->rateHetero);	// setPij dependes on the root Matrix
	// Transfer independent rates into the initial_rates holding the dependent Qij's.
	vector<vector<double> >::iterator ppt = initial_rates.Pij.begin(); 
	for (vector<vector<double> >::iterator PPt = k_0->branch->rates->Pij.begin(); PPt != k_0->branch->rates->Pij.end(); ++PPt, ++ppt) {
		vector<double>::iterator pt = (*ppt).begin();
		for(vector<double>::iterator Pt = (*PPt).begin(); Pt != (*PPt).end(); ++Pt, ++pt) {
			(*pt) = (*Pt);
		}
	}

	return initial_rates;
}

RateMatrix
iQdPc(
	  TTree *tree,
	  TNode *i_z,
	  TNode *k_0,
	  double T,
	  double at_dt,
	  int event_site
	 )
{
	i_z->setRateAway(TIME_RELATIVE_STEPS);	// Resets site_sum_away, allocates e_QijDt.
	RateMatrix initial_rates(*k_0->branch->rates);
	unsigned int end_site;

	//////////
	/// Dependent sites. Recalculate Q matrices for positions at the beginning of the branch or
	/// affected by a change for subsequent calculation of the site-specific transition matrices.
	//////////
	i_z->set_site_window(tree->dep.front()->context.return_order(), &event_site, &end_site);
	i_z->site_specific_Qmat(tree, event_site, end_site);

	return initial_rates;
}

RateMatrix
iQdP(
	 TTree *tree,
	 TNode *i_z,
	 TNode *k_0,
	 double T,
	 double at_dt,
	 int event_site
	)
{
	i_z->setRateAway(TIME_RELATIVE_STEPS);	// Resets site_sum_away, allocates e_QijDt.
	RateMatrix initial_rates(*k_0->branch->rates);
	unsigned int end_site;

	//////////
	/// Dependent sites. Recalculate Q matrices for positions at the beginning of the branch or
	/// affected by a change for subsequent calculation of the site-specific transition matrices.
	//////////
	i_z->set_site_window(tree->dep.front()->context.return_order(), &event_site, &end_site);
	i_z->site_specific_Qmat(tree, event_site, end_site);

	// setPij uses the Root and Cijk matrices, already set up for branch->rates
	k_0->branch->rates->setPij(k_0->seq_evo.front(), k_0->branch->S*(T-at_dt), k_0->nodeEnv->rateHetero);	// setPij dependes on the root Matrix
	// Transfer independent rates into the initial_rates holding the dependent Qij's.
	vector<vector<double> >::iterator ppt = initial_rates.Pij.begin(); 
	for (vector<vector<double> >::iterator PPt = k_0->branch->rates->Pij.begin(); PPt != k_0->branch->rates->Pij.end(); ++PPt, ++ppt) {
		vector<double>::iterator pt = (*ppt).begin();
		for(vector<double>::iterator Pt = (*PPt).begin(); Pt != (*PPt).end(); ++Pt, ++pt) {
			(*pt) = (*Pt);
		}
	}

	return initial_rates;
}

void
uQN(
	RateMatrix *rates, 
	TNode *i_z, 
	vector<Site>::iterator site, 
	double T,
	double at_dt
   )
{
	/// Passing global probabilities. This route is a great deal faster.
}

void
uQdPc(
	  RateMatrix *rates, 
	  TNode *i_z, 
	  vector<Site>::iterator site, 
	  double T,
	  double at_dt
     )
{
	/// Calculating and passing the site-specific transition probabilities.
	(*rates).Qij = (*site).e_QijDt;							// Entire Qij matrix. For exponentiating.
	(*rates).SetupMatrix(false);							// Set Cijk matrix, for exponentiating.
	(*rates).setPij((*site), (T-at_dt), i_z->nodeEnv->rateHetero); // Exponentiation.
}

void
uQdP(
	 RateMatrix *rates, 
	 TNode *i_z, 
	 vector<Site>::iterator site, 
	 double T,
	 double at_dt
    )
{
	// Fast option
	(*rates).Qij = (*site).e_QijDt;
}

void
uQPc(
	 RateMatrix *rates, 
	 TNode *i_z, 
	 vector<Site>::iterator site, 
	 double T,
	 double at_dt
    )
{
	(*rates).Qij = (*site).e_QijDt;							// Entire Qij matrix. For exponentiating.
	(*rates).SetupMatrix(false);							// Set Cijk matrix, for exponentiating.
	(*rates).setPij((*site), (T-at_dt), i_z->nodeEnv->rateHetero); // Exponentiation.
	(*rates).Qij = i_z->branch->rates->Qij;					// Reset dependent rates to independent rates.
}

void
uQP(
	RateMatrix *rates, 
	TNode *i_z, 
	vector<Site>::iterator site, 
	double T,
	double at_dt
   )
{
	/// Passing global probabilities. This route is a great deal faster.
}

