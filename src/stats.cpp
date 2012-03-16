#include "stats.h"
#include "trace.h"
#include "dependency.h"
#include "paleo.h"

extern bool forward_simulation;
extern bool rasmus_independent_proposals;
extern bool MCMC_sample_evenly;

vector<double> t_path_start (100001, 0);
vector<double> t_path_end (100001, 0);
stringstream accepted_subpath_times;	// To compare function vs. current representations. // 

SampleStatistics::SampleStatistics() : forward_probability(0), epc_probability(0), w_i(0), numEvents(0), P_path (new PathProbability()) { }

void 
SampleStatistics::calculateStatistics (
									  		list<eventTrack*> *events
									  	   )
{
	vector<eventTrack*> sorted_events_by_time;

	cerr << "--------------calculateStatistics" << endl;

	//////////
	/// EPC rate vectors: Q.at(idot_k__T__) and Q.at(ij_k__T__)
	//////////
	cerr << "num events:: " << sorted_events_by_time.size() << endl;
	cerr << (*events).front()->Q.idot_k__T__ << endl;

	if ( !(*events).empty() ) {
		if ( forward_simulation ) {
			forward_probability = P_path->ForwardProbability(*events);
			cerr << "FORWARD_PROBABILITY: " << forward_probability << endl;
		} else {
			forward_probability = P_path->ForwardProbability(*events);
			epc_probability = P_path->EPCProbability(*events);
		}
	} else {
		cerr << "NO EVENTS IN SampleStatistics::calculateStatistics()." << endl;
		exit(EXIT_FAILURE);
	}
}		

void 
SampleStatistics::reportSampleStatistics()
{
	cerr << "----------" << endl;
	cerr << "Number of Substitutions:                 " << numEvents << endl;
	cerr << "Forward Probability        [Qi. ]:       " << forward_probability << endl;
	cerr << "EPC probability            [Qi.|k(T-t)]: " << epc_probability << endl;
	cerr << "m_i:                                   : " << m_i() << endl;
	cerr << "w_i:                                   : " << return_w_i() << endl;
}

void 
Statistics::reportStatistics()
{
	cerr << "Relative weights: " << endl;
	for (vector<SampleStatistics>::iterator stats_it = sample_stats.begin(); stats_it != sample_stats.end(); ++stats_it) 
		(*stats_it).reportSampleStatistics();

	cerr << "Sample set statistics: " << endl;
	cerr << "  WEIGHTS" << endl;
	cerr << "    mean:     " << importance_sampling.return_mu() << endl;
	cerr << "    stdev:    " << importance_sampling.return_sigma() << endl;
	cerr << "    median:   " << importance_sampling.return_median() << endl;
	cerr << "    range:    " << importance_sampling.range.return_min() << ":" << importance_sampling.range.return_max() << endl;

}

void 
Statistics::setWeights()
{
	vector<double> mod_m_i (sample_stats.size(), 0);
	vector<double>::iterator jt = mod_m_i.begin();

	//////////
	/// Still in log space. This loop reverts from log space. In order to avoid underflow, we subtract
	/// the maximum log-space value from each m_i, making the values relative.
	//////////
	for (vector<SampleStatistics>::iterator it = sample_stats.begin(); it != sample_stats.end(); ++it, ++jt) {
		(*jt) = exp((*it).m_i() - returnM());	/// Exponentiating the difference of m_i and maximum m_i.
		m_i_sum += (*jt);
		cerr << (*jt) << endl;
	}

	jt = mod_m_i.begin();
	for (vector<SampleStatistics>::iterator it = sample_stats.begin(); it != sample_stats.end(); ++it, ++jt) {
		(*it).set_w_i( (*jt) / m_i_sum );
		(*jt) = (*it).return_w_i();
	}

	importance_sampling.mu_sigma(mod_m_i);

}

void 
Statistics::ImportanceSampling::mu_sigma( vector<double> X )
{
	int N = X.size();
	double mu = 0;
	double sigma = 0;
	
//	std::sort(X.begin(), X.end());
	med = X.at(X.size()/2);
	range.set_min(X.front());
	range.set_max(X.back());

	/// Calculate average
	for (vector<double>::iterator it = X.begin(); it != X.end(); ++it) {
		mu += (*it);
	}

	mu /= X.size();
	for (vector<double>::iterator it = X.begin(); it != X.end(); ++it) sigma += (mu - (*it)) * (mu - (*it));
	sigma /= X.size();

	u = mu;
	r = sqrt(sigma);
}

void 
Statistics::calc_ESS()
{
	int N = sample_stats.size();
	effective_sample_size = N / ( 1 + ESS() );

	cerr << N << " / ( " << 1 + ESS() << " )" << endl;
	
	cerr << "ESS = " << effective_sample_size << endl;
}

double 
Statistics::ESS()
{
	int N = sample_stats.size();
	double denominator;
	double mu_m_i = 1.0/N * m_i_sum;
	
	double numerator = 0, squareme;
	for (vector<SampleStatistics>::iterator it = sample_stats.begin(); it != sample_stats.end(); ++it) {
		squareme = ( exp((*it).m_i()-returnM() ) - mu_m_i );
		numerator += squareme*squareme;
		cerr << " -> exp(" << (*it).m_i() << "-" << returnM() << ") - " << mu_m_i << endl;
	}
	cerr << "1.0/" << N << " * " << m_i_sum << ": " << mu_m_i << endl;
	cerr << numerator << " / " << (N-1.0) * mu_m_i * mu_m_i << endl;
	denominator = 1.0 / ((N-1.0) * mu_m_i * mu_m_i);

	return numerator * denominator;
}

void 
Statistics::MCMC_run(
					 TTree *tree,
					 int number_of_steps, 
					 PathProposal *path
					)
{
	double t_0, T;	/// The begin and end of the branch. Need to change once we do a tree.
	double t_B, t_E;	/// The two times between which we want to simulate.
	TNode *i_B, *i_E, *work;	/// Virtual nodes within the branch that we are simulating, one at t_B and one at t_E.
	list<eventTrack*>::iterator e_B, e_E;	/// Iterators into the event history marking the begin and end times.
	list<eventTrack*> proposed_subpath_events, current_subpath_events;
	double cP, cJ, pP, pJ;
	stringstream report;
	vector<double> mcmc_chain_forward;
	int accepted = 0;
	double post_subpath_event_time;
	int cycle_sample_length = 1;

	//////////
	/// Set originally simulated path to the current path (along with events).
	//////////
	mcmc.current = setPath(path);
	mcmc.current->epc_events = getBranchEvents(path->epc_events.begin(), path->epc_events.end());
	t_0 = mcmc.current->epc_events.front()->eventTime;
	T = mcmc.current->epc_events.back()->eventTime;

	TNode *i_0 = tree->root;	// tree->root is temporary. Use this while we are working only on the single branch.
	TNode *k_0 = tree->root->branch1;
	work = tree->copyNode(i_0);
	delete (*mcmc.current->epc_events.begin()); delete (*mcmc.current->epc_events.rbegin());
	mcmc.current->epc_events.pop_front(); mcmc.current->epc_events.pop_back();
	mcmc.current->EmulateStep(tree, work, k_0, t_0, T, &mcmc.current->epc_events);
	cerr << "--------------------------------------------------------" << endl;
	cerr << " Calculation of original path" << endl;
	cerr << "--------------------------------------------------------" << endl;
	cerr << "Original Path Events: " << endl;
	cerr << "i_0: " << i_0->printSequence(true) << endl;
	cerr << "k_0: " << k_0->printSequence(true) << endl;
	//mcmc.current->Print_Path_Events();
	mcmc_chain_forward.push_back(mcmc.current->P_path.ForwardProbability(mcmc.current->epc_events));

	cerr << "--------------------------------------------------------" << endl;
	Leakage();
	cerr << "+++++++++===============++++++++++++++++++++++++++++++" << endl;

	//cout << "i_0: " << i_0->printSequence(true) << endl;
	//cout << "k_0: " << k_0->printSequence(true) << endl;

	for (int i = 0; i < number_of_steps; i++) {
		cerr << "***********************************************************************************" << endl;
		cerr << "****** STEP " << i << " *** STEP " << i << " *** STEP " << i << "************************************************" << endl;
		cerr << "***********************************************************************************" << endl;

		cout << i << "  ";
		if (rasmus_independent_proposals) 
			mcmc_chain_forward.push_back( mcmc.rasmus_resample(tree, i_0, k_0, t_0, T, mcmc_chain_forward.back()) );
		else {
			mcmc_chain_forward.push_back( mcmc.resample_subpath(tree, i_0, k_0, t_0, T, mcmc_chain_forward.back()) );
		}
		cout << "  " << mcmc_chain_forward.back() << endl;

		if (i % cycle_sample_length == 0) tree->global_alignment->cleanArray(mcmc.current->epc_events);	
		cerr << "DONE." << endl;
	}

	work->Remove_Objects();
	delete work->evolvingSequence;
	delete work;
	i_0->Remove_Objects();
	delete i_0->evolvingSequence;
	delete i_0;
	delete k_0->evolvingSequence;
	k_0->Remove_varSites();
	delete k_0;

	cerr << "+++++++++===============++++++++++++++++++++++++++++++" << endl;
	Leakage();

	cerr << "Times on branch of accepted proposals: " << endl;
	cerr << accepted_subpath_times.str() << endl;

	exit(0);
}

double
Statistics::MCMC::rasmus_resample(
								  TTree *tree,
								  TNode *i_0,
								  TNode *k_0,
								  double t_0,
								  double T,
								  double current_cycle_probability
								 )
{
	int site_sampled;
	TNode *work = tree->copyNode(i_0);
	bool accepted = false;
	double return_probability;

	site_sampled = (int)(rndu()*i_0->seq_evo.size());
	cerr << "Sampled site " << site_sampled << endl;

	vector<Site>::iterator site_work = work->seq_evo.begin()+site_sampled;
	vector<Site>::iterator site_i_0 = i_0->seq_evo.begin()+site_sampled;
	vector<Site>::iterator site_k_0 = k_0->seq_evo.begin()+site_sampled;

	proposed = setPath(current);
	proposed->remove_site_events(site_sampled);
	order_3_markov = false;

	// This is the black box routine that evolves the path until a time is chosen that exceeds T. Rejects if
	// the state doesn't match, otherwise accepts path.
	proposed->evolve_independent_path(tree, work, k_0, t_0, T, &proposed->epc_events, site_work, site_k_0, site_sampled);
	order_3_markov = true;
	delete (*proposed->epc_events.begin()); delete (*proposed->epc_events.rbegin());
	proposed->epc_events.pop_front(); proposed->epc_events.pop_back();

	work->resetSequence(i_0);
	proposed->EmulateStep(tree, work, k_0, t_0, T, &proposed->epc_events);
	double pForwardDep, pForwardInd, cForwardDep, cForwardInd;
	cForwardDep = current->P_path.ForwardProbability(current->epc_events, true);
	pForwardDep = proposed->P_path.ForwardProbability(proposed->epc_events, true);
	list<eventTrack*> current_site_path = current->extract_site_events(site_sampled);
	list<eventTrack*> proposed_site_path = proposed->extract_site_events(site_sampled);

	//cerr << "Current site path for site " << site_sampled << endl;
	//for (list<eventTrack*>::iterator it = current_site_path.begin(); it != current_site_path.end(); ++it)
	//	cerr << (*it)->Print_Event();
	//cerr << "Proposed site path for site " << site_sampled << endl;
	//for (list<eventTrack*>::iterator it = proposed_site_path.begin(); it != proposed_site_path.end(); ++it)
	//	cerr << (*it)->Print_Event();

	cForwardInd = current->P_path.IndependentForwardProbability(current_site_path, (*site_i_0).returnState(), k_0->branch->rates);
	pForwardInd = proposed->P_path.IndependentForwardProbability(proposed_site_path, (*site_i_0).returnState(), k_0->branch->rates);

	//cerr << endl << "cForwardDep: "  << cForwardDep << endl << "pForwardDep: " << pForwardDep << endl;
	//cerr << "cForwardInd: "  << cForwardInd << endl << "pForwardInd: " << pForwardInd << endl;

	if (accept_proposal(cForwardDep, cForwardInd, pForwardDep, pForwardInd)) {
		for (list<eventTrack*>::iterator iiit = current->epc_events.begin(); iiit != current->epc_events.end(); ++iiit)
			delete (*iiit);
		delete current;
		current = setPath(proposed);
		tree->global_alignment->cleanArray(current->epc_events);
		return_probability = current_cycle_probability - cForwardInd + pForwardInd;
	} else {
		return_probability = current_cycle_probability;
	}
	
	for (list<eventTrack*>::iterator jjjt = proposed->epc_events.begin(); jjjt != proposed->epc_events.end(); ++jjjt) 
		delete (*jjjt);
	delete proposed;
	work->Remove_Objects();
	delete work->evolvingSequence;
	delete work;

	cout << "  " << current->epc_events.size()-2;

	return return_probability;
}

double
Statistics::MCMC::resample_subpath(
								   TTree *tree,
								   TNode *i_0,
								   TNode *k_0,
								   double t_0,
								   double T,
								   double current_cycle_probability
								  )
{
	double t_B, t_E;	/// The two times between which we want to simulate.
	TNode *i_B, *i_E, *work;	/// Virtual nodes within the branch that we are simulating, one at t_B and one at t_E.
	list<eventTrack*>::iterator e_B, e_E;	/// Iterators into the event history marking the begin and end times.
	list<eventTrack*> proposed_subpath_events, current_subpath_events;
	double cP, cJ, pP, pJ;
	double return_probability;
	double post_subpath_event_time;
	bool error_debug = false;

	cP = cJ = pP = pJ = 0;	// This resets P(\rho) and J(\rho | \rho'), P(\rho') and J(\rho | \rho') to zero for this proposal step.
	// Create a PathProposal object for this round of proposals.
	proposed = setPath(current);

	//////////
	/// Setup all necessary sequences for a resimulation of subpath.
	//////////
	// 4TREE: Select branch.
	// Collect events that occurred on branch.
	t_B = t_0 + rndu() * (T-t_0) + t_0;
	t_E = t_0 + rndu() * (T-t_0) + t_0;
	if (t_B > t_E) swap(t_B, t_E);	// STL swap routine.

	//////////
	/// Subpath sampling test:
	/// 1/3 probability of t_0->t_B
	/// 1/3 probability of t_B->t_E
	/// 1/3 probability of t_E->T
	//////////
	if (MCMC_sample_evenly) {
		double RN = rndu();
//		if (RN < 0.5) { t_E = t_B; t_B = t_0; }
//		else { t_B = t_E; t_E = T; }

//		t_E = t_B; t_B = t_0;	// First interval
//		t_B = t_E; t_E = T;		// third interval
//		;						// second interval
		
		// Interval 123 combined
//		if (RN < 1.0/3.0) { t_E = t_B; t_B = t_0; }
//		else if (RN > 2.0/3.0) { t_B = t_E; t_E = T; }
	}

	list<eventTrack*>::iterator xf = proposed->epc_events.begin(); xf++;
	list<eventTrack*>::reverse_iterator xb = proposed->epc_events.rbegin(); xb++;
	cout << t_B << "  " << t_E << "  " << "  " << (*xf)->eventTime << "  " << (*xb)->eventTime;

	// Set up nodes
	i_B = tree->copyNode(i_0);
	i_E = tree->copyNode(i_0);
	work = tree->copyNode(i_0);

	//////////
	/// Set path pointers in sequence time.
	//////////
	//// 3-8-2012: SPEED THIS UP!!!! e_E CAN BE SET starting FROM E_B.
	e_B = proposed->setSequenceAtTimePoint(tree, i_B, t_B, &(proposed->epc_events));
	e_E = proposed->setSequenceAtTimePoint(tree, i_E, t_E, &(proposed->epc_events));
	if ( e_E == proposed->epc_events.end() ) post_subpath_event_time = T;
	else { e_E++; post_subpath_event_time = (*e_E)->eventTime; e_E--; }

	/// Current setup
	/// eventID:     0     1   2         3  45      6         7       8
	/// Branch:  t_0 |-----|---|--*------|--||------|--*------|-------|T
	///                          t_B    e_B                  e_E
	/// Set to the event AFTER t_E (ID=7), outside of subpath. e_B set to event after t_B (ID=3), which is correct.
	/// One consideration: Is this event at the end of the path?
	work->resetSequence(i_B);	// Set scratch sequence to begin state.

	//////////
	/// p = path; p'= proposed path; p_L, p_M, p_R = probabilities of left, middle, and right subpaths;
	///
	/// After randomly choosing the subpath interval, we now simulate a new subpath, in the end getting
	/// current_path = { p_L, p_M, p_R }
	/// proposed_path = { p_L, p'_M, p_R }
	/// where
	///	 |  p_L     |  p_M   |        p_R         |
	///  |          | subpath|					  |
	///  |----------|--------|--------------------|
	///  |         t_B      t_E					  |
	///  |  p_L     |  p'_M  |        p_R		  |
	///
	/// We use the probabilities of p_M and p'_M as the indicator of accepting the proposed path,
	/// since p_L and p_R cancel out.
	///
	///  P(p)  = P(p_L) x P(p_M|p_L)  x P(p_R|p_M)
	///  P(p') = P(p_L) x P(p'_M|p_L) x P(p_R|p'_M)
	//////////

	cerr << "(SUB)PATH SETUP: " << t_0 << " --> " << t_B << "  " << t_E << "<- " << T << endl;
	//cerr << "i_0: " << i_0->printSequence(true) << endl;
	//cerr << "i_B: " << i_B->printSequence(true) << endl;
	//cerr << "  Rate away from sequence: " << i_B->evolvingSequence->Qidot << endl;
	//cerr << "i_E: " << i_E->printSequence(true) << endl;
	//cerr << "  Rate away from sequence: " << i_E->evolvingSequence->Qidot << endl;
	//cerr << "k_0: " << k_0->printSequence(true) << endl;

	list<eventTrack*>::iterator pt, ct;
	//////////
	/// PROPOSED SUBPATH
	//////////
	// Set the site-state likelihoods for the subpath. 
	for (vector<Site>::iterator qt = i_E->seq_evo.begin(); qt != i_E->seq_evo.end(); ++qt)
		(*qt).L_i.Li_xi_.at( (*qt).returnState() ) = 1;
	// Propose events for the subpath. Do not need to emulate path, since the proposal step will calculate
	// the probabilities of the events (Qij, Qi., Qij|k(t), Qi.|k(t)).
	proposed_subpath_events.clear();
	work->resetSequence(i_B);
	tree->dep.front()->context.set_sequence_indices(work);
	
	cout << "  " << work->evolvingSequence->compare_sequence(i_E->evolvingSequence);	// Gets the number of differences. If times are really short, and no differences are observed, this may explain some problems that the method is having.
	
	proposed->EvolveStep(tree, work, i_E, t_B, t_E, &proposed_subpath_events);
	for (pt = proposed_subpath_events.begin(); pt != proposed_subpath_events.end(); ++pt)
		(*pt)->Compute_MSA_Positions(tree, 0);

	//////////
	/// CURRENT SUBPATH
	//////////
	// Excise the events that were proposed for the previous round in the subpath.
	current_subpath_events.clear();
	if (e_B != e_E) current_subpath_events.assign(e_B, e_E);
	// i_B holds the state at the beginning of the subpath, work is a scratch sequence. Reset
	// work to appropriate sequence for emulations.
	work->resetSequence(i_B);
	PathProposal *current_subpath = new PathProposal(*current);
	// Both paths emulated. Now we should be able to calculate both the forward and epc probabilities. Check full events, first.
	current_subpath->EmulateStep(tree, work, i_E, t_B, t_E, &current_subpath_events);

	cout << "  " << current_subpath_events.size() << "," << proposed_subpath_events.size();

	//if (current_subpath_events.size() == 2 && proposed_subpath_events.size() == 4) {
	//	for (pt = proposed_subpath_events.begin(); pt != proposed_subpath_events.end(); ++pt)
	//		cout << "pt: " << (*pt)->Print_Event();
	//	cout << "---------------------------------------------------------" << endl;
	//	for (ct = current_subpath_events.begin(); ct != current_subpath_events.end(); ++ct)
	//		cout << "ct: " << (*ct)->Print_Event();
	//	cout << "---------------------------------------------------------" << endl;
	//	error_debug = true;
	//}

	// Calculate probabilities.
	pP = proposed->P_path.ForwardProbability(proposed_subpath_events);
	pJ = proposed->P_path.EPCProbability(proposed_subpath_events);
	cP = current_subpath->P_path.ForwardProbability(current_subpath_events);
	cJ = current_subpath->P_path.EPCProbability(current_subpath_events);
	// current_subpath_events has eventTracks from mcmc.current sandwiched by independent front and back eventTrack pointers from emulateStep. //

	// Clear up first and last events (non-subst events)
	delete (*current_subpath_events.begin()); delete (*current_subpath_events.rbegin());

	// Emulate the two subpaths to get their forward probability (P(p)) and their EPC probabilities
	// (J(p|p')) necessary to calculate r = min{1, (P(p')/J(p'|p)) / (P(p)/J(p|p'))
//	double logRN = log(rndu());
//	double log_r = logr(cP, cJ, pP, pJ);
//	bool accept_new_path = ( (logRN < log_r ) ? true : false );
	bool accept_new_path = accept_proposal(cP, cJ, pP, pJ);

	// Now that events have been emulated, probabilities have been calculated, and path has been accepted,
	// modify the times of the proposed path to add in the beginning of the subpath interval.
	if ( accept_new_path ) {
		list<eventTrack*>::iterator current_start, current_end, proposed_start, proposed_end;

		//cout << "Time resimulated: " << t_B << "..." << t_E << endl;
		//cout << "i_B: " << i_B->printSequence() << endl;
		//cout << "i_E: " << i_E->printSequence() << endl;

		cout << " A";

		// Get the beginning and end events for both the current and proposed paths.
		work->resetSequence(i_B);
		current_start = current->setSequenceAtTimePoint(tree, i_B, t_B, &(current->epc_events));
		work->resetSequence(i_E);
		current_end = current->setSequenceAtTimePoint(tree, i_E, t_E, &(current->epc_events));
		work->resetSequence(i_B);
		proposed_start = proposed->setSequenceAtTimePoint(tree, i_B, t_B, &(proposed->epc_events));
		work->resetSequence(i_E);
		proposed_end = proposed->setSequenceAtTimePoint(tree, i_E, t_E, &(proposed->epc_events));

		//cerr << "Time resimulated: " << t_B << "..." << t_E << endl;
		//cerr << "current_event start iterator eventTime:  " << (*current_start)->eventTime << endl;
		//cerr << "current_event end iterator eventTime:    " << (*current_end)->eventTime << endl;
		//cerr << "proposed_event start iterator eventTime: " << (*proposed_start)->eventTime << endl;
		//cerr << "proposed_event end iterator eventTime:   " << (*proposed_end)->eventTime << endl;
		//cerr   << "*P(p):    " << cP << endl << "*J(p|p'): " << cJ << endl << "*P(p'):   " << pP << endl << "*J(p'|p): " << pJ << endl;
		for (list<eventTrack*>::iterator erase_it = proposed_start; erase_it != proposed_end; ++erase_it)
			delete (*erase_it);
		proposed->epc_events.erase(proposed_start, proposed_end);
		delete (*proposed_subpath_events.begin()); delete (*proposed_subpath_events.rbegin());
		proposed_subpath_events.pop_front(); proposed_subpath_events.pop_back();
		proposed->epc_events.insert(proposed_end, proposed_subpath_events.begin(), proposed_subpath_events.end());

		for (list<eventTrack*>::iterator iiit = current->epc_events.begin(); iiit != current->epc_events.end(); ++iiit)
			delete (*iiit);
		delete current;
		current = setPath(proposed);

		//if (error_debug) {
		//	cout << "ERROR DEBUG: pP = " << pP << " cP = " << cP << endl;
		//
		//	cout << "// PROPOSED EVENTS TO INSERT INTO SUBPATH //" << endl;
		//	for (list<eventTrack*>::iterator yt = proposed_subpath_events.begin(); yt != proposed_subpath_events.end(); ++yt)
		//		cout << (*yt)->Print_Event();
		//	cout << "// PROPOSED SUBPATH //" << endl;
		//	proposed->Print_Path_Events();
		//
		//	// Don't do this. If error_debug is set, 
		//	//cout << "// CURRENT EVENTS IN SUBPATH //" << endl;
		//	//for (list<eventTrack*>::iterator yt = current_subpath_events.begin(); yt != current_subpath_events.end(); ++yt)
		//	//	cout << (*yt)->Print_Event();
		//	cout << "// CURRENT SUBPATH //" << endl;
		//	current->Print_Path_Events();
		//}

		return_probability = current_cycle_probability + pP - cP;
		//cerr << "Forward fullpath difference: " << forward_prob_fullpath << "   subpath: " << cP - pP << endl;
		//cerr << "Forward probability: " << current_cycle_probability << " + " << pP << " - " << cP << " = " << return_probability << endl;

		//current->Print_Path_Events();

		/// Checking routine, to see if the values from the subpath insertion are correct.
		//double fwd_curr = current->P_path.ForwardProbability(current->epc_events);
		//double fwd_prop = proposed->P_path.ForwardProbability(proposed->epc_events);
		//double forward_prob_fullpath = fwd_curr - fwd_prop;
		//if (forward_prob != fwd_prop) {
		//	cerr << "Not equal, step " << i << endl;
		//	cerr << "Forward: " << return_probability << "    " << fwd_prop << endl;
		//	cerr << "     Differences: " << forward_prob_fullpath << " vs: " << cP-pP << endl;
	
		//	exit(0);
		//}

		accepted_subpath_times << " " << t_B << " " << t_E << endl;
	} else {
		cout << " R";
		for (list<eventTrack*>::iterator iit = proposed_subpath_events.begin(); iit != proposed_subpath_events.end(); ++iit)
			delete (*iit);
		return_probability = current_cycle_probability;
	}

	delete current_subpath;
	for (list<eventTrack*>::iterator it = proposed->epc_events.begin(); it != proposed->epc_events.end(); ++it)
		delete (*it);
	delete proposed;	// Speed up: don't need to delete this on an accepted path.
	i_B->Remove_Objects();
	delete i_B->evolvingSequence;
	delete i_B;
	i_E->Remove_Objects();
	delete i_E->evolvingSequence;
	delete i_E;
	work->Remove_Objects();
	delete work->evolvingSequence;
	delete work;

	cout << "  " << current->epc_events.size()-2;

	cerr << "DONE." << endl;
	
	return return_probability;
}

PathProposal *setPath(PathProposal *path)
{
	eventTrack *new_event;
	PathProposal *new_path = new PathProposal(*path);
	new_path->epc_events.clear();
	for (list<eventTrack*>::iterator it = path->epc_events.begin(); it != path->epc_events.end(); ++it) {
		new_event = new eventTrack(*(*it));
		new_path->epc_events.push_back(new_event);
	}
	return new_path;
}

double Statistics::MCMC::logr(
							  double P,
							  double J,
							  double Pstar,
							  double Jstar
							 )
{
	double r = (Pstar-Jstar)-(P-J);
	cerr << "Statistics::MCMC::rlog: (" << Pstar << "-" << Jstar << ") - (" << P << "-" << J << ") = " << r << endl;
	return min(0.0,r);
}

bool
Statistics::MCMC::accept_proposal (
								   double P,
								   double J,
								   double Pstar,
								   double Jstar
								  )
{
	// (J(p|p')) necessary to calculate r = min{1, (P(p')/J(p'|p)) / (P(p)/J(p|p'))
	double logRN = log(rndu());
	double log_r = logr(P, J, Pstar, Jstar);
	bool accept_new_path = ( (logRN < log_r ) ? true : false );

	cerr   << "*P(p):    " << P << endl << "*J(p|p'): " << J << endl 
		   << "*P(p'):   " << Pstar << endl << "*J(p'|p): " << Jstar << endl 
		   << "*log_r = min{0, -->}" << "Acceptance test: (logRN < log(r)) = " << logRN 
		   << " > " << log_r << " = subpath " << ( (accept_new_path) ? "accepted" : "rejected" ) 
		   << "." << endl;

	return accept_new_path;
}

void Leakage()
{
	long memoryLeaked = 0;

	memoryLeaked = 0
				 + print_memory("<TNode>:", 		   			TNode::howMany(), 					sizeof(TNode))
				 + print_memory("<inClade>:", 		   			inClade::howMany(), 				sizeof(inClade))
				 + print_memory("<inTree>:", 		   			inTree::howMany(), 					sizeof(inTree))
				 + print_memory("<TTree>:", 		   			TTree::howMany(), 					sizeof(TTree))
				 + print_memory("<Branch>:", 		   			Branch::howMany(), 					sizeof(Branch))
				 + print_memory("<inMotif>:", 		   			inMotif::howMany(), 				sizeof(inMotif))
				 + print_memory("<siteRegEx>:", 	   			siteRegEx::howMany(), 				sizeof(siteRegEx))
				 + print_memory("<Sequence>:", 		   			Sequence::howMany(), 				sizeof(Sequence))
				 + print_memory("<Site>:", 			   			Site::howMany(), 					sizeof(Site))
				 + print_memory("<motifSite>:", 	   			motifSite::howMany(), 				sizeof(motifSite))
				 + print_memory("<siteProperties>:",   			siteProperties::howMany(),   		sizeof(siteProperties))
				 + print_memory("<activeProperties>:", 			activeProperties::howMany(), 		sizeof(siteProperties))
				 + print_memory("<seqGenOptions>:",    			seqGenOptions::howMany(), 			sizeof(seqGenOptions))
				 + print_memory("<Indel>:", 		   			Indel::howMany(), 					sizeof(Indel))
				 + print_memory("<Substitution>:", 	   			Substitution::howMany(), 			sizeof(Substitution))
				 + print_memory("<Insertion>:", 	   			Insertion::howMany(), 				sizeof(Insertion))
				 + print_memory("<Deletion>:", 		   			Deletion::howMany(), 				sizeof(Deletion))
				 + print_memory("<varSite>:", 		   			varSite::howMany(), 				sizeof(varSite))
				 + print_memory("<eventTrack>:", 	   			eventTrack::howMany(), 				sizeof(eventTrack))
				 + print_memory("<eventTrack::ratesAway>:",		eventTrack::ratesAway::howMany(),	sizeof(eventTrack::ratesAway))
				 + print_memory("<globalArray>:", 	   			globalArray::howMany(), 			sizeof(globalArray))
				 + print_memory("<insertSite>:", 	   			insertSite::howMany(), 				sizeof(insertSite))
				 + print_memory("<siteModifier>:",	   			siteModifier::howMany(),			sizeof(siteModifier))
				 + print_memory("<RateMatrix>:",	   			RateMatrix::howMany(),       		sizeof(RateMatrix))
				 + print_memory("<PathProposal>:",				PathProposal::howMany(),			sizeof(PathProposal))
				 + print_memory("<PathProbability>:",			PathProbability::howMany(),			sizeof(PathProbability))
				 + print_memory("<Taxon>:",						Taxon::howMany(),					sizeof(Taxon))
				 + print_memory("<paleobiology>:",				paleobiology::howMany(),			sizeof(paleobiology))
				 + print_memory("<Statistics>:",       			Statistics::howMany(),             	sizeof(Statistics))
				 + print_memory("<Dependency>:",       			Dependency::howMany(),       		sizeof(Dependency))
				 + print_memory("<contextDependence>:",			contextDependence::howMany(),		sizeof(contextDependence))
				 + print_memory("<LookUp>:",                    LookUp::howMany(),					sizeof(LookUp))
				 ;

	cerr << "Total memory leaked = " << memoryLeaked << " bytes." << endl << endl;
}

long print_memory(
				  string message, 
				  size_t howMany, 
				  int object_size
				 )
{
	cerr << setw(25) << message << setw(10) << howMany << " x " << setw(5) << object_size
		 << " = " << setw(10) << howMany*object_size << endl;
	
	return howMany*object_size;
}
