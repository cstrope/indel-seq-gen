#include "propose_path.h"
#include "dependency.h"

extern vector<short> acgt;
extern vector<int> acagatcgctgt;
extern int print_stepwise_rates;

#define PRE 0
#define IN 1
#define POST 2

extern bool round1;
extern int changed_site;
extern bool rasmus_independent_proposals;
extern bool Qd, Pc, nij;
double sum_Pr_ratio = 0;
double sum_QPr = 0;
double sum_Pr_away = 0;
double sum_Pr_away2 = 0;
int num_changes = 0;
double rate_away_same = 0, rate_away_diff = 0;
double sum_Pij_same = 0, sum_Pij_diff = 0;
double _epc_S_average_ = 0;

// For a specific time point, this array collects the number of changes that occur to each sequence
// site after a time point X={0.5}
vector<double> exp_mean_draws;

////////////////////
////// Sets the sequences for all nodes that are specified in the fasta-format input file.
////////////////////
void 
PathProposal::setNodeSequences( 
									TTree *tree
							   	   )
{
	vector<Taxon*>::iterator jt;
	bool found_taxon;
	bool leaf_node;
	size_t jt_count = 0;
	TNode *found;

	for (vector<Taxon*>::iterator jt = taxa.begin(); jt != taxa.end(); ++jt) {
		found = NULL;
		for (list<TNode*>::iterator it = tree->nodeList.begin(); it != tree->nodeList.end(); ++it) {
			// Leaf taxon found.
			if ( (*jt)->taxon_name.compare((*it)->node_name.c_str()) == 0 ) {
				found = (*it);
				break;
			} else if ( (*it)->ancestorNo == atoi((*jt)->taxon_name.c_str()) ) {
				/// Internal node found.
				if (it == tree->nodeList.begin()) tree->epc_root_is_set = true;
				found = (*it);
			} 
		}

		if (found) {
			TNode *tmp_anc = found->anc;
			found->anc = found;
			found->constructUnconstrainedSequence((*jt)->taxon_sequence);
			found->evolvingSequence->ConvertRootSequence();
			found->anc = tmp_anc;
		}
	}
}

void 
PathProposal::setUnspecifiedNodeSequences(
											   TTree *tree,
											   int sequence_length
											  )
{
	string sequence;
	sequence.assign(sequence_length, 0);	/// All A's

	if (!tree->epc_root_is_set) {
		for (vector<Site>::iterator jt = tree->root->seq_evo.begin(); jt != tree->root->seq_evo.end(); ++jt) 
			(*jt).setState(-1);
	}

	for (list<TNode*>::iterator it = tree->nodeList.begin(); it != tree->nodeList.end(); ++it) {
		//////////
		/// if there is not evolvingSequence, then this sequence has not yet been set.
		/// Should be okay to do this for any type of node (even leaf), for, e.g., forward sim.
		//////////
		if ((*it)->evolvingSequence == NULL) {
			(*it)->constructUnconstrainedSequence(sequence);
			for (vector<Site>::iterator jt = (*it)->seq_evo.begin(); jt != (*it)->seq_evo.end(); ++jt) {
				(*jt).setState(-1);
			}
		}
	}
}

////////////////////
////// Lists all of the paths that will be simulated
////////////////////
void 
PathProposal::displayPaths(
								TTree *tree
							   )
{
	for (vector<TNode*>::iterator it = tree->tips.begin(); it != tree->tips.end(); ++it) {
		cerr << "From sequence: " 
			 << (*it)->anc->printSequence()
			 << endl << "  to sequence: "
			 << (*it)->printSequence()
			 << endl << "  with BL = " << (*it)->branch->length0 << endl << endl;
	}
}

void 
PathProposal::parseEndpointInfile(
									   string& filename
									  )
{
	ifstream infile;
	string in_buffer = "";
	infile.open(filename.c_str());

	//////////
	/// Use to read in endpoints (or more sequences, later).
	//////////
	if (infile) {
		in_buffer += " ";
		while (infile.good()) {
			in_buffer += infile.get();
		}
		infile.close();
		in_buffer += "\n>";
		//////////
		/// This will be where a more advanced version of parsing will take place, eventually.
		//////////
		trim(in_buffer);
		taxa = parseFastaFile(in_buffer);
	} else {
		cerr << "Unable to open file \"" << filename << "\"." << endl;
		exit(EXIT_FAILURE);
	}
}

vector<Taxon*> 
PathProposal::parseFastaFile(
											string& fasta_data
										   )
{
	vector<Taxon*> return_taxa;
	string sequence_build;
	size_t found;
	list<string> each_taxon = split (fasta_data, ">");
	each_taxon.pop_front();		// Get rid of leading space.
	each_taxon.pop_back();		// Added extra ">". See what happens...

	list<string> each_line;
	Taxon *current_taxon;
	for (list<string>::iterator it = each_taxon.begin(); it != each_taxon.end(); ++it) {
		current_taxon = new Taxon();
		each_line = split((*it),"\n");
		trim(each_line.front());
		current_taxon->taxon_name = each_line.front();
		each_line.pop_front();
		for (list<string>::iterator jt = each_line.begin(); jt != each_line.end(); ++jt) {
			sequence_build.clear();
			//////////
			/// Getting odd '?' characters at the begin and end of my file. This will blast them.
			//////////
			trim(*jt);
			Remove_Whitespace(*jt);
			for (string::iterator kt = (*jt).begin(); kt != (*jt).end(); ++kt)
				if ( stateCharacters.find(*kt) != string::npos ) sequence_build += (*kt);
			//
			//////////
			current_taxon->taxon_sequence += sequence_build;
		}
		return_taxa.push_back(current_taxon);
	}

	return return_taxa;
}

void 
PathProposal::Evolve(
						  TTree *tree,
						  list<eventTrack*> *events
						 )
{
	if (tree->root->nodeEnv->rateHetero == DiscreteGammaRates) {
		//tree->calculateJCLikelihoods();
		tree->setCatLikelihoods(1);
	} else if (tree->root->nodeEnv->rateHetero == GammaRates) {
		cerr << "At this time, path proposals do not work with non-discrete gamma rates." << endl;
		exit(EXIT_FAILURE);
	}

	//////////
	/// Setting the state likelihoods depends on the gamma categories. Two cases reach this point.
	///   (1) No gamma categories, so need to set likelihoods.
	///   (2) Gamma cats set above, so need to recalculate conditional probs for the given cat
	///       (this is extra work, so may re-do it in the setting of cat likelihoods).
	//////////
	tree->setStateLikelihoods(1);

	//////////
	/// Now have conditional likelihoods based on categories. Choose root sequence. If this is
	/// a single-branch tree, the likelihoods for each site will be 1.0 for the correct state, so
	/// the sequence will simply be copied.
	//////////
	tree->sample_root_sequence();
	//exit(0);

	//////////
	/// Evolve both branches.
	//////////
	this->EvolveToEndPoint(tree, tree->root, tree->root->branch1, events);
	this->EvolveToEndPoint(tree, tree->root, tree->root->branch2, events);

	for (vector<Taxon*>::iterator it = taxa.begin(); it != taxa.end(); ++it) 
		delete (*it);

//	cout << "EVENTS:" << eventNo << endl;	//XOUT
	cerr << "EVENTS:" << eventNo << endl;	//XOUT
}

void 
PathProposal::EvolveToEndPoint(
							   TTree *tree,
							   TNode *anc,
							   TNode *des,
							   list<eventTrack*> *events
							  )
{
	if (des->branch->length0 == 0) return;

	TNode *work = new TNode();		// Scratch sequence, used as substitutions accumulate down branch.
	work->evolvingSequence = new Sequence(anc->evolvingSequence, des);
	work->nodeEnv = anc->nodeEnv;
	work->bipartition = des->bipartition;
	work->branch = des->branch;

	des->branch->rates->Qij = anc->branch->rates->Qij;
	work->branch->rates->Qij = anc->branch->rates->Qij;

	//////////
	/// Des will need to inherit the gamma categories.
	//////////
	des->inheritCategories(anc);
	if (rasmus_independent_proposals) {
		EvolveIndependentStep(tree, work, des, 0, des->branch->length0, events); // NEED TO WATCH... OPTIONS MAY CONFLICT WITH DEPENDENT SITES STUFF
	} else {
		EvolveStep(tree, work, des, 0, des->branch->length0, events);
	}
	des->updateSequence(work);
    if (des->tipNo==-1) { 
	    EvolveToEndPoint(tree, des, des->branch1, events);
	   	EvolveToEndPoint(tree, des, des->branch2, events);
    }

	des->Remove_Objects();
	work->Remove_iz_Objects();
	work->evolvingSequence->evolutionaryAttributes.clear();
	delete work->evolvingSequence;
	delete work;
}

void 
PathProposal::EvolveIndependentStep(
						 			TTree *tree,
						 			TNode *i_z,
						 			TNode *k_0,
						 			double t_0,
						 			double T,
						 			list<eventTrack*> *events
								   )
{
	cerr << "PathProposal::EvolveIndependentStep: " << endl;
	order_3_markov = false;

	bool isEndpoint = false;
	eventTrack *begin_event, *end_event;

	isEndpoint = false;		// Set isEndpoint to false. This allows for non-acceptable proposals in the loop below.
	//////////
	/// Quick calculation of evolvingSequence->Qidot (forward rate away);
	//////////
	(*(i_z->seq_evo.begin())).forward_rate_away_from_site(i_z->branch);
	i_z->calculateForwardRateAwayFromSequence__order3Markov(tree, 0);

	// Independent. Rate away from sequence is based on forward rates, and there are no EPC rates.
	begin_event = branch_terminal_event(k_0->mytipNo, BRANCH_BEGIN, k_0->DistanceFromRoot, k_0->bipartition, t_0, T); 
	begin_event->assign_Q(i_z, i_z->seq_evo.at(0), SUBSTITUTION, 1);
	begin_event->Q.idot_k__T__ = 0;
	(*events).push_back(begin_event);

	//////////
	// Indel routine	
	//  * Should codonRates affect the probability of an indel occurring?
	//////////
	cerr << i_z->printBipartition() << "(0)  lambda_T(epc) SEQUENCE: " << i_z->evolvingSequence->returnRij() << endl; //XOUT
	//i_z->branch->rates->printCategories();

	cerr << "i_z = " << i_z->printSequence(true) << endl;
	cerr << "k_0 = " << k_0->printSequence(true) << endl;
	// With independence, we simulate each site, rather than whole sequence.
	vector<Site>::iterator kt = k_0->seq_evo.begin();
	int site = 0;
	for (vector<Site>::iterator it = i_z->seq_evo.begin(); it != i_z->seq_evo.end(); ++it, ++kt, ++site) {
		(*it).forward_rate_away_from_site(i_z->branch);
		evolve_independent_path(tree, i_z, k_0, t_0, T, events, it, kt, site);
	}	// NEXT SITE //

	k_0->calculateForwardRateAwayFromSequence__order3Markov(tree, 0);
	end_event = branch_terminal_event(k_0->mytipNo, BRANCH_END, k_0->DistanceFromRoot, k_0->bipartition, t_0, T); 
	end_event->assign_Q(i_z, i_z->seq_evo.at(0), SUBSTITUTION, 1);
	(*events).push_back(end_event);

	tree->global_alignment->cleanArray(*events);

	sortEventsByTime(events);
	cerr << endl << "Rasmus Nielsen path:" << endl;
	for (list<eventTrack*>::iterator it = (*events).begin(); it != (*events).end(); ++it)
		cerr << (*it)->Print_Event();

	order_3_markov = true;
}

void 
PathProposal::evolve_independent_path(
							  		  TTree *tree,
							  		  TNode *i_z,
							  		  TNode *k_0,
							  		  double t_0,
							  		  double T,
							  		  list<eventTrack*> *events,
							  		  vector<Site>::iterator i_z_site,
							  		  vector<Site>::iterator k_0_site,
							  		  int event_site
									 )
{
	order_3_markov = false;

	list<eventTrack*> proposal_events;
	bool success;
	double next_dt, lambda_T;
	string event;
	int action = INSERT;
	int devnull = 1;
	int proposalNo = 0;
	short state_x0, state_xT;
	state_x0 = (*i_z_site).returnState();
	state_xT = (*k_0_site).returnState();
	double t0_pr = t_0;

	cerr << "Proposal from site " << event_site << ": " << stateCharacters.at((*i_z_site).returnState()) << " to " << stateCharacters.at((*k_0_site).returnState()) << endl;

	do {
		if (proposalNo) cerr << "proposal " << proposalNo << " XXXXXXXXX" << endl;

		// If there is a rejection by the end of the branch, reset the state to it's original.
		(*i_z_site).setState(state_x0);
		t0_pr = t_0;

		// Evolve path using Nielsen's method.
		// If there is a difference between the start and end states, need to use special case to make sure
		// there is a change
		if (state_x0 != state_xT) {
			// To determine the first change, we calculate the probability of a change occurring at time 
			// t1, normalized by the probability of that a change will occur in time less than T:
			// f(t1|t1<T) = \frac{Q_{X(0).}\exp{-Q_{X(0)t_1}} {1-\exp{-Q_{X(0).}T}}
			lambda_T = i_z_site->forward_rate_away_from_site(i_z->branch);
			t0_pr = rasmus_select_next_dt(lambda_T, T);

			cerr << "t0_pr = " << t0_pr << " from branch length " << T << endl;
		}

		for (list<eventTrack*>::iterator it = proposal_events.begin(); it != proposal_events.end(); ++it)
			delete (*it);
		proposal_events.clear();
		lambda_T = i_z_site->forward_rate_away_from_site(i_z->branch);
		next_dt = select_next_dt(i_z, k_0, t0_pr, lambda_T, T, &devnull, false);
		devnull++;
		for (double dt = next_dt; dt <= T; dt = next_dt) {
			success = false;
			//////////
			/// If we sent rate_away matrix, then we simply need to select a substitution from
			/// the probabilities that are calculated.
			//////////
			event.clear();
			// Passing in rate away from sequence as lambda_T+event_site. lambda_T \in {0, 1], which is the independent
			// rate away (e.g., number of substitutions per site), and event_site is the current site. This is necessary
			// because the substitute routine assumes that the state is a sequence, whereas in the independent model, the
			// site is the current state. Substitute chooses the site to change based on the rate away from the sequence.
			// Adding event_site to lambda_T will cause the substitute routine to make changes based on the site under
			// consideration.
			success = substitute(tree, i_z, &event_site, INDEPENDENT_ENDPOINT_CONDITIONED, event, true, lambda_T);
			action = SUBSTITUTION;
			if (success) {
				eventTrack *new_event;
				//////////
				/// Calculate the relative time that the event occurs at.
				//////////
				// * i_z->trDistanceFromRoot: The scaled distance from the root
				// * (dt/indel_len): Percentage of the branch that has been simulated
				// * k_0->branch0_time_relative_length: Time rel length of the branch being simulated.
				// * iTree->global_max_path_length: Scalar for the sum of above terms to set time of occurrence between 0 and 1.
				//////////
				k_0->atEpochTime = i_z->DistanceFromRoot + dt;
				new_event = new eventTrack(
										   eventNo,
										   action,
										   k_0->atEpochTime,
										   k_0->bipartition,
										   event,
										   T
										  );
				// Since only one change occurs each iteration, we do not need to compare the sequences
				// until at least num_diff changes have occurred.
				new_event->assign_Q(i_z, i_z->seq_evo.at(event_site), SUBSTITUTION, 1);
				lambda_T = i_z_site->forward_rate_away_from_site(i_z->branch);
				proposal_events.push_back(new_event);
				cerr << "IND-EVO subst pos: " << event_site << " " << event << "->" << stateCharacters.at((*k_0_site).returnState()) << ", proposal " << proposalNo << ":  "; //XOUT
				eventNo++;
			}
				cerr << "---------- Time left: " << T-dt << "/" << T << "----------"; //XOUT
	//		cerr << "EVO: " << i_z->printSequence(true) << "  Qi. = " << i_z->evolvingSequence->Qidot << endl;
			cerr << "  Qi. = " << lambda_T << endl;
	//		cout << "DT:" << dt << endl;			//XOUT
	//		cout << "LAMBDA_T:" << lambda_T << endl; //XOUT
	//		cout << "Qi.|k:" << i_z->evolvingSequence->Qidot_k__T__ << endl; //XOUT
	//		cout << "JC_Rij:" << JC_Rij << endl;	//XOUT
	//		cout << "Qi.:" << i_z->evolvingSequence->Qidot << endl;	//XOUT
	//		cout << "SITE_RATE_AWAY:" << i_z->seq_evo.at(changed_site).site_rate_away.back() << endl; //XOUT
	//		cout << endl;		//XOUT

			// isEndpoint = false, doesn't resample waiting time from exp for failed paths.
			next_dt = select_next_dt(i_z, k_0, dt, lambda_T, T, &devnull, false); 	
			devnull++;
		}
		proposalNo++;
	} while ( (*i_z_site).returnState() != (*k_0_site).returnState() );

	for (list<eventTrack*>::iterator it = proposal_events.begin(); it != proposal_events.end(); ++it)
		(*it)->Compute_MSA_Positions(tree, 0);

	//////////
	/// Accepted proposal. Push events into the event array.
	//////////
	for (list<eventTrack*>::iterator it = proposal_events.begin(); it != proposal_events.end(); ++it)
		(*events).push_back(*it);
	sortEventsByTime(events);
	cerr << "  " << proposalNo << " proposals before success." << endl;
}

double
PathProposal::rasmus_select_next_dt(
									double lambda_T,
									double T
								   )
{
	double next_dt, exp_mean;

	// Simulation of the first substitution for the Nielsen method when A != E.
	// t_1 = -\log(1-U(1-\exp{-q_AT})/q_A
	// where q_A is the rate away from the ancestral state A, U~Uniform(0,1), and T is the branch length.
	
	next_dt = -log(1-rndu()*(1-exp(-lambda_T*T)))/lambda_T;

	return next_dt;
}

void 
PathProposal::EvolveStep(
						 TTree *tree,
						 TNode *i_z,
						 TNode *k_0,
						 double t_0,
						 double T,
						 list<eventTrack*> *events
						)
{
	double lambda_T, next_dt;
	bool success, isEndpoint = false;
	int action = INSERT, event_site, num_diff;
	string event;
	eventTrack *begin_event, *end_event;
	double i_z_rate_away;

	isEndpoint = k_0->EndpointCheck();

	//////////
	/// Set up indices of dependency for the sequence on this branch.
	//////////
	if (order_3_markov) tree->dep.front()->context.set_sequence_indices(i_z);
	//////////
	/// Quick calculation of evolvingSequence->Qidot (forward rate away);
	//////////
	if (order_3_markov) i_z_rate_away = i_z->calculateForwardRateAwayFromSequence__order3Markov(tree, -1);
	else i_z->evolvingSequence->Qidot = k_0->seq_evo.size();
	num_diff = k_0->evolvingSequence->compare_sequence(i_z->evolvingSequence);

	// For independent sites, branch needs to be scaled from substitutions per site to an approximation
	// of the dependent number of substitutions per site. S is the value that we use.
	if (!Pc && !nij) i_z->branch->S = i_z_rate_away / i_z->seq_evo.size();
	else i_z->branch->S = 1;

	// Initial calculation of nij's
	if (nij) {
		i_z->branch->nij.assign(numStates*numStates, 0);
		i_z->branch->nij_pseudocounts.assign(numStates*numStates,0.01);
		vector<Site>::iterator kt = k_0->seq_evo.begin();
		for (vector<Site>::iterator it = i_z->seq_evo.begin(); it != i_z->seq_evo.end(); ++it, ++kt) {
			i_z->branch->nij.at( (*it).returnState()*numStates + (*kt).returnState() )++;
		}
		vector<double>::iterator pk = i_z->branch->nij_pseudocounts.begin();
		for (vector<double>::iterator it = i_z->branch->nij.begin(); it != i_z->branch->nij.end(); ++it, ++pk) {
			(*it)+=(*pk);
		}
	}

	lambda_T = i_z->calculateEndpointRateAwayFromSequence(tree, k_0, T, t_0, -1);
	begin_event = branch_terminal_event(k_0->mytipNo, BRANCH_BEGIN, k_0->DistanceFromRoot, k_0->bipartition, t_0, T); 
	begin_event->assign_Q(i_z, i_z->seq_evo.at(0), SUBSTITUTION, num_diff);
	begin_event->Q.idot_k__T__ = lambda_T;
	(*events).push_back(begin_event);

	//////////
	// Indel routine	
	//  * Should codonRates affect the probability of an indel occurring?
	//////////
	cout << lambda_T << endl;	///TESTING_ONLY: FOR HEAT MAP.pl
	cerr << i_z->printBipartition() << "(0)  lambda_T(epc): " << lambda_T << endl; //XOUT
	//i_z->branch->rates->printCategories();

	next_dt = select_next_dt(i_z, k_0, t_0, lambda_T, T, &num_diff, isEndpoint);

//	double JC_Rij; int ID2DIFF = 0, DIFF2DIFF = 0, DIFF2ID = 0, ID2ID = 0;
	for (double dt = next_dt; dt <= T; dt = next_dt) {
		success = false;
		//////////
		/// If we sent rate_away matrix, then we simply need to select a substitution from
		/// the probabilities that are calculated.
		//////////
		event.clear();
		success = substitute(tree, i_z, &event_site, ENDPOINT_CONDITIONED, event, true, lambda_T);
		action = SUBSTITUTION;
		if (success) {
			eventTrack *new_event;
			if (nij) i_z->branch->update_nij(stateCharacters.find_first_of(event.at(0)), stateCharacters.find_first_of(event.at(1)), k_0->seq_evo.at(event_site).returnState());		// Update the nij values based on the change.
			//////////
			/// Calculate the relative time that the event occurs at.
			//////////
			// * i_z->trDistanceFromRoot: The scaled distance from the root
			// * (dt/indel_len): Percentage of the branch that has been simulated
			// * k_0->branch0_time_relative_length: Time rel length of the branch being simulated.
			// * iTree->global_max_path_length: Scalar for the sum of above terms to set time of occurrence between 0 and 1.
			//////////
			k_0->atEpochTime = i_z->DistanceFromRoot + dt;
			new_event = new eventTrack(
									   eventNo,
									   action,
									   k_0->atEpochTime,
									   k_0->bipartition,
									   event,
									   T
									  );
			// Since only one change occurs each iteration, we do not need to compare the sequences
			// until at least num_diff changes have occurred.
			if (--num_diff <= 0)
				num_diff = k_0->evolvingSequence->compare_sequence(i_z->evolvingSequence);
			new_event->assign_Q(i_z, i_z->seq_evo.at(event_site), SUBSTITUTION, num_diff);
			tree->dep.front()->context.reset_sequence_indices(i_z, event_site, event);
			lambda_T = i_z->calculateEndpointRateAwayFromSequence(tree, k_0, T, dt, event_site);
			new_event->Q.idot = i_z->evolvingSequence->Qidot;
			new_event->Q.idot_k__T__ = lambda_T;
			//cerr << "New event: "  << new_event->Print_Event() << endl;
			(*events).push_back(new_event);
			cerr << "EVO subst pos: " << event_site << " " << event << "|" << stateCharacters.at(k_0->seq_evo.at(event_site).returnState()) << ":  "; //XOUT
			if (nij) {
				i_z->branch->print_nij();
				i_z->branch->print_nij(false);
				
				vector<double>::iterator pt = i_z->branch->rates->Pij.at(0).begin();
				int i = 0;
				cerr << endl << "Transition Probabilities (Jukes Cantor).";
				for (; pt != i_z->branch->rates->Pij.at(0).end(); ++pt, ++i) {
					if (i % numStates == 0) cerr << endl;
					cerr << (*pt) << " ";	
				}

			}
			eventNo++;
		}

		//JC_Rij 
		//= Ns * (0.75 * ((1-exp(-(T-dt))) / (1+3*exp(-(T-dt)))) )
		//+ Nd * (0.25 * ( (exp(-(T-dt)) + 3.0) / (1-exp(-(T-dt))) ) );

		cerr << "---------- Time left: " << T-dt << "/" << T << "----------"; //XOUT
//		cerr << "EVO: " << i_z->printSequence(true) << "  Qi. = " << i_z->evolvingSequence->Qidot << endl;
//		cerr << "  Qi. = " << i_z->evolvingSequence->Qidot << endl;
//		cout << "DT:" << dt << endl;			//XOUT
//		cout << "LAMBDA_T:" << lambda_T << endl; //XOUT
//		cout << "Qi.|k:" << i_z->evolvingSequence->Qidot_k__T__ << endl; //XOUT
//		cout << "JC_Rij:" << JC_Rij << endl;	//XOUT
//		cout << "Qi.:" << i_z->evolvingSequence->Qidot << endl;	//XOUT
//		cout << "SITE_RATE_AWAY:" << i_z->seq_evo.at(changed_site).site_rate_away.back() << endl; //XOUT
//		cout << endl;		//XOUT

		next_dt = select_next_dt(i_z, k_0, dt, lambda_T, T, &num_diff, isEndpoint);
	}

	num_diff = k_0->evolvingSequence->compare_sequence(i_z->evolvingSequence);
	if (num_diff > 0) {	cerr << "num_diff = " << num_diff << endl; exit(0); }
	
	end_event = branch_terminal_event(k_0->mytipNo, BRANCH_END, k_0->DistanceFromRoot, k_0->bipartition, t_0, T); 
	end_event->assign_Q(i_z, i_z->seq_evo.at(0), SUBSTITUTION, num_diff);
	(*events).push_back(end_event);
}

double 
PathProposal::select_next_dt(
							 TNode *i_z,
							 TNode *k_0,
							 double current_dt, 
							 double lambda_T,
							 double T,
							 int *num_diff,
							 bool isEndpoint		/// If no, then we simulate towards conditional Likelihoods, and will not reject a change.
							)
{
	double next_dt, exp_mean;

	exp_mean = 1.0/lambda_T;
	next_dt = current_dt + rand_exp(exp_mean);
	if (next_dt > T && isEndpoint) {
		if ( *num_diff > 0 ) {
			int num_trials = 0;
			cerr << *num_diff << "->";
			*num_diff = k_0->evolvingSequence->compare_sequence(i_z->evolvingSequence);
			cerr << *num_diff << "  ";
			cerr << "Ooops... still " << *num_diff << " change(s) to do: " << current_dt << " " << next_dt << "->";
			do {
				next_dt = current_dt + rand_exp(exp_mean); 
				if (num_trials++ > 1000) {
					cerr << "Proposed parameters for end-point conditioning may be unsuitable."
						 << endl;
					return current_dt + 0.0000000000000000000000000001;	// For testing purposes, need to choose a really small time slot.
					//exit(EXIT_FAILURE); // Normal procedure. Bad EPC models should normally do this.
				}
			} while (next_dt > T);
			cerr << next_dt << " " << T << endl;
		}
	}
	
	return next_dt;
}

void 
PathProposal::emulateForwardSimulation (
						  				TTree *tree,
										list<eventTrack*> *events 
									   )
{
	tree->setStateLikelihoods(1);
	tree->sample_root_sequence();

	//////////
	/// Evolve both branches.
	//////////
	this->emulateToEndPoint(tree, tree->root, tree->root->branch1, events);
//	this->emulateToEndPoint(tree, tree->root, tree->root->branch2, events);

	
}

list<eventTrack*>::iterator 
PathProposal::setSequenceAtTimePoint(
									 TTree *tree,
									 TNode *node,
									 double to_time,
									 list<eventTrack*>::iterator et,
									 list<eventTrack*> *events
									)
{
	while ( (*et)->ID == -1 ) ++et;
	while ( et != (*events).end() && ((*et)->eventTime < to_time)) {	// Check if et is at end, short-circuit and avoid seg fault.
		if ( (*et)->ID == -1 ) { 
			et = (*events).end(); 	// Sets et 1 event past the last event.
			et--; 					// Set to last event.
			break;
		} else {
			// Make replacement of sequence position.
			node->seq_evo.at( (*et)->MSA_positions.at(0) ).setState(stateCharacters.find_first_of((*et)->size.at(1)) );
			/// Starting at i_0 and setting at i_z; need to also set the Qij matrices at each position, or else strange
			/// values will be returned.
			++et;
		}
	}

	return et;
}

void 
PathProposal::emulateToEndPoint(
								TTree *tree,
								TNode *anc,
								TNode *des,
								list<eventTrack*> *events
							   )
{
	TNode *work = tree->copyNode(des);	// Scratch sequence, used as substitutions accumulate down branch.
//	TNode *work = new TNode();		// Scratch sequence, used as substitutions accumulate down branch.
//	work->evolvingSequence = new Sequence(anc->evolvingSequence, des);
//	work->nodeEnv = anc->nodeEnv;
//	work->bipartition = des->bipartition;
//	work->branch = des->branch;

	des->branch->rates->Qij = anc->branch->rates->Qij;
	work->branch->rates->Qij = anc->branch->rates->Qij;

	EmulateStep(tree, work, des, 0, des->branch->length0, events);
//  For now, I am happy simple emulating a single branch as testcase. 
//	des->updateSequence(work);
//    if (des->tipNo==-1) { 
//	    emulateToEndPoint(tree, des, des->branch1, events);
//	   	emulateToEndPoint(tree, des, des->branch2, events);
//	  }

}

void 
PathProposal::EmulateStep(
						  TTree *tree,
						  TNode *i_z,
						  TNode *k_0,
						  double t_0,
						  double T,						// True branch end, not subpath branch end in subpath simulation.
						  list<eventTrack*> *events
						 )
{
	double lambda_T;
	eventTrack *begin_event, *end_event;
	int num_diff = 0;
	double i_z_rate_away;
	
	assert(T-t_0>0);
	if (!(*events).empty())
	assert((*events).front()->ID != -1);

	if (order_3_markov) tree->dep.front()->context.set_sequence_indices(i_z);

	//for (list<eventTrack*>::iterator it = (*events).begin(); it != (*events).end(); ++it) {
	//	cerr << (*it)->Print_Event();
	//}

	for (vector<Site>::iterator it = i_z->seq_evo.begin(); it != i_z->seq_evo.end(); ++it) (*it).calcTauIJ = true;
	//////////
	/// Quick calculation of evolvingSequence->Qidot (forward rate away);
	//////////
	if (order_3_markov) i_z_rate_away = i_z->calculateForwardRateAwayFromSequence__order3Markov(tree, -1);
	else i_z->evolvingSequence->Qidot = k_0->seq_evo.size();

	// For independent sites, branch needs to be scaled from substitutions per site to an approximation
	// of the dependent number of substitutions per site. S is the value that we use.
	if (!Pc && !nij) i_z->branch->S = i_z_rate_away / i_z->seq_evo.size();
	else i_z->branch->S = 1;

	// Initial calculation of nij's
	if (nij) {
		i_z->branch->nij.assign(numStates*numStates, 0);
		i_z->branch->nij_pseudocounts.assign(numStates*numStates,0.01);
		vector<Site>::iterator kt = k_0->seq_evo.begin();
		for (vector<Site>::iterator it = i_z->seq_evo.begin(); it != i_z->seq_evo.end(); ++it, ++kt) {
			i_z->branch->nij.at( (*it).returnState()*numStates + (*kt).returnState() )++;
		}
		vector<double>::iterator pk = i_z->branch->nij_pseudocounts.begin();
		for (vector<double>::iterator it = i_z->branch->nij.begin(); it != i_z->branch->nij.end(); ++it, ++pk) {
			(*it)+=(*pk);
		}
	}

	// BRANCH BEGIN:
	lambda_T = i_z->calculateEndpointRateAwayFromSequence(tree, k_0, T, t_0, -1);
	begin_event = branch_terminal_event(k_0->mytipNo, BRANCH_BEGIN, k_0->DistanceFromRoot, k_0->bipartition, t_0, T); 
	begin_event->assign_Q(i_z, i_z->seq_evo.at(0), SUBSTITUTION, num_diff);

//	cerr << "Pre-looping thru events:" << endl;
	cerr << "  Qi. = " << i_z->evolvingSequence->Qidot << endl;
	if (print_stepwise_rates) cout << print_stepwise_rates << " " << i_z->evolvingSequence->Qidot << endl;
//	cerr << "  Qi.|k(" << T-t_0 << ") = " << lambda_T << endl;

	//////////
	/// EVENTS: For collecting stats, it is important to maintain order of operations. See EvolveStep.
	//////////
	string event = "XX";
	for (list<eventTrack*>::iterator et = (*events).begin(); et != (*events).end(); ++et) {
		// Emulate substitution: 
		// * MSA_position where this substitution occurred w.r.t. MSA.
		// * size.at(1), for a substitution, is the nature of the change, where 0th position is old state, first is new state.
		event.at(0) = stateCharacters.at(i_z->seq_evo.at( (*et)->MSA_positions.at(0)).returnState());
		i_z->seq_evo.at( (*et)->MSA_positions.at(0) ).setState(stateCharacters.find_first_of( (*et)->size.at(1) ) );
		event.at(1) = stateCharacters.at(i_z->seq_evo.at( (*et)->MSA_positions.at(0)).returnState());
		cerr << "EMU pos " << (*et)->MSA_positions.at(0) << " " << event << "---------- Time (" << (*et)->eventTime << "), left: " << T-(*et)->eventTime << "/" << T-t_0 << "----------"; //XOUT
		num_diff = k_0->evolvingSequence->compare_sequence(i_z->evolvingSequence);
		(*et)->assign_Q(i_z, i_z->seq_evo.at((*et)->MSA_positions.at(0)), SUBSTITUTION, num_diff);
		(*et)->event_occurrence_branch_length = T-t_0;
		if (nij) i_z->branch->update_nij(stateCharacters.find_first_of(event.at(0)), stateCharacters.find_first_of(event.at(1)), k_0->seq_evo.at((*et)->MSA_positions.at(0)).returnState());		// Update the nij values based on the change.

		tree->dep.front()->context.reset_sequence_indices(i_z, (*et)->MSA_positions.at(0), event);
		lambda_T = i_z->calculateEndpointRateAwayFromSequence(tree, k_0, T, (*et)->eventTime, (*et)->MSA_positions.at(0) );
		(*et)->Q.idot = i_z->evolvingSequence->Qidot;
		(*et)->Q.idot_k__T__ = lambda_T;
		if (print_stepwise_rates) cout << print_stepwise_rates << " " << (*et)->eventTime << " " << i_z->evolvingSequence->Qidot << endl;
//		cerr << "EMU: " << i_z->printSequence(true) << "  Qi. = " << i_z->evolvingSequence->Qidot << endl;
//		cerr << "  Qi.|k(" << T-(*et)->eventTime << ") = " << lambda_T << endl;
//		cout << "DT:" << (*et)->eventTime << endl;			//XOUT
//		cout << "LAMBDA_T:" << lambda_T << endl; //XOUT
//		cout << "Qi.|k:" << i_z->evolvingSequence->Qidot_k__T__ << endl; //XOUT
//		cout << "Qi.:" << i_z->evolvingSequence->Qidot << endl;	//XOUT
//		cout << endl;		//XOUT

		eventNo++;
	}	

	// BRANCH END
	end_event = branch_terminal_event(k_0->mytipNo, BRANCH_END, k_0->DistanceFromRoot, k_0->bipartition, t_0, T); 
	end_event->assign_Q(i_z, i_z->seq_evo.at(0), SUBSTITUTION, num_diff);
	(*events).push_front(begin_event);
	(*events).push_back(end_event);
}

void PathProposal::setEventHistory (
									list<eventTrack*> *events,
									string event_history_file
								   )
{
	ifstream infile;
	string in_buffer;
	infile.open(event_history_file.c_str());
	list<string> event_data;
	list<string>::iterator event_it;

	//////////
	/// Use to read in endpoints (or more sequences, later).
	//////////
	if (infile) {
		while (infile.good()) {
			getline(infile, in_buffer);
			if ( isEvent(in_buffer) ) {
				event_data.clear();
				event_data = split (in_buffer, ",");	// Split all data items.
				event_it = event_data.begin();

				eventTrack *new_event;
				new_event = new eventTrack(
										   atoi((*(event_it++)).c_str()),
										   (*(event_it++)).at(0),		// A single character.
										   strtod((*(event_it++)).c_str(), NULL),
										   make_vector_bool(*(event_it++)),
										   (*(event_it++)),
										   atoi((*(event_it++)).c_str()),
										   strtod((*(event_it++)).c_str(), NULL),
										   strtod((*(event_it++)).c_str(), NULL),
										   strtod((*(event_it++)).c_str(), NULL),
										   strtod((*(event_it++)).c_str(), NULL),
										   strtod((*(event_it++)).c_str(), NULL),
										   strtod((*(event_it++)).c_str(), NULL),
										   strtod((*(event_it++)).c_str(), NULL)
										  );
				(*events).push_back(new_event);
			}
		}
		infile.close();
		//////////
		/// This will be where a more advanced version of parsing will take place, eventually.
		//////////
	} else {
		cerr << "Unable to open file \"" << event_history_file << "\"." << endl;
		exit(EXIT_FAILURE);
	}	

	cerr << "Number of events: " << (*events).size() << endl;
//	for (list<eventTrack*>::iterator it = (*events).begin(); it != (*events).end(); ++it)
//		cerr << (*it)->Print_Event();
}

vector<bool> 
make_vector_bool (string str)
{
	string::iterator it;
	vector<bool> bool_vector;
	
	for (it = str.begin(); it != str.end(); ++it)
		if ((*it) == '1')
			bool_vector.push_back(true);
		else if ((*it) == '0')
			bool_vector.push_back(false);
		else {
			cerr << "String value " << (*it) << " is not a boolean character." << endl << "string = \"" << str << "\"" << endl;
			exit(EXIT_FAILURE);
		}

	return bool_vector;
}

double 
PathProbability::EPCProbability (
						     	 list<eventTrack*> events,
						     	 bool do_last_event
						    	)
{
	long double log_sum_epc = 0;
	double exp_value_dt;
	double dt;
	bool profile=false;							//XOUT
	list<eventTrack*>::iterator it, prev2it, end_event;
	list<eventTrack*> branch_events;

	double T, t_z, i_z_Qidot_k__T__;
	
	// Get single branch events, for now. NEED TO DO FOR ALL BRANCHES EVENTUALLY.
	it = events.begin();
	while (it != events.end() && (*it)->eventType != BRANCH_END) {
		branch_events.push_back(*it);
		++it;
	}
	branch_events.push_back(*it);
	// This is one branch. Data I have:
	// front() -> [-1,B,t,...]	-1 flags as non-event, B = beginning of branch, t = time at start, Qi. also held.
	// event1
	// event2
	// ...
	// back()  -> [-1,E,t,...]  As above.
	// This is all that is necessary to calculate the forwardProbability.
	T = branch_events.back()->eventTime;

	if (branch_events.size() > 2) {
		it = branch_events.begin();
		prev2it = branch_events.begin();
		end_event = branch_events.end(); end_event--;
		it++;
		int xx = 0;
		for (; it != end_event; ++it, xx++, ++prev2it) {
			dt = (*it)->eventTime - (*prev2it)->eventTime;
			if (dt == 0) dt = numeric_limits<double>::min();
			log_sum_epc	+= EPC_Step((*it)->Q.ij_k__T__, (*prev2it)->Q.idot_k__T__, (*it)->Q.i2k, (*prev2it)->eventTime, dt,	T);
			//cerr << " = " << log_sum_epc << endl;
		}
		t_z = (*prev2it)->eventTime;
		i_z_Qidot_k__T__ = (*prev2it)->Q.idot_k__T__;
	} else {
		t_z = events.front()->eventTime;
		i_z_Qidot_k__T__ = branch_events.front()->Q.idot_k__T__;
	}

	// it and prev2it are exactly where we want them to be for this calculation.
	if (do_last_event) { 
		dt = T - t_z;
		log_sum_epc += EPC_NoEvent(i_z_Qidot_k__T__, dt);
	}

	cerr << "LOG_EPC: " << log_sum_epc << endl;
	return log_sum_epc;
}

double 
PathProbability::EPC_NoEvent(
									double Qidot_k__T__,
									double dt
								   )
{
	if (dt == 0) dt = numeric_limits<double>::min();
	//cerr << -Qidot_k__T__ << " * " << dt;
	return -Qidot_k__T__ * dt;
}

double 
PathProbability::EPC_Step(
								double Qij_k__T__,
								double Qidot_k__T__,
								int Qi2k,
								double t_z,
								double dt,
								double T
							   )
{
	// Calculation of the overstep penalty (t_z+dt > T) when i!=k,
	// i.e., log( 1 - ( ((*it)->Q.i2k != 0) ? exp_value_dt : 0 )) 
	// NOTE: if exp_value_dt is effectively 0, then this function is log(1) = 0, thus overstep_penalty
	// is originally set to 0. However, when exp_value is significant enough to care about calculating it,
	// (to XX decimal places of accuracy), then we need to calculate both the log and the exp. Otherwise,
	// skip it to save computation time.

	double exp_value = -Qidot_k__T__ * (T-t_z);
	double overstep_penalty = 0;
	if (exp_value > -115) // -115: 50 decimal places of accuracy; -230: 100 decimal places.
		if (Qi2k != 0) 
			overstep_penalty = log( 1 - exp( exp_value ) );
	//cerr <<	"log(" << Qij_k__T__  << ") - (" << Qidot_k__T__ << " * " << dt << ") - " << overstep_penalty;
	return 	log( Qij_k__T__ ) - (Qidot_k__T__ * dt) - overstep_penalty;
}

double 
PathProbability::Forward_NoEvent(
								 double Qidot,
								 double dt
								)
{
	//cerr << -Qidot << " * " << dt;
	return -Qidot * dt;
}

double 
PathProbability::Forward_Step(
							  double Qij,
							  double Qidot,
							  double dt
							 )
{
	assert(Qij != 0);
	
	//cerr << log(Qij) << " + ";
	return log(Qij) + Forward_NoEvent(Qidot, dt);
}

double 
PathProbability::ForwardProbability (
						   	   	 	 list<eventTrack*> events,
						   	   	 	 bool do_last_event
							  		)
{
	double log_sum_forward = 0;
	double dt;
	list<eventTrack*> branch_events;
	list<eventTrack*>::iterator prev2it, end_event, it;
	double T, t_z;
	double i_z_Qidot;
	

	// Get single branch events, for now. NEED TO DO FOR ALL BRANCHES EVENTUALLY.
	branch_events = getBranchEvents(events.begin(), events.end()/*?*/);

	// This is one branch. Data I have:
	// front() -> [-1,B,t,...]	-1 flags as non-event, B = beginning of branch, t = time at start, Qi. also held.
	// event1
	// event2
	// ...
	// back()  -> [-1,E,t,...]  As above.
	// This is all that is necessary to calculate the forwardProbability.
	T = branch_events.back()->eventTime;

	if (branch_events.size() > 2) {		// 2: 1 branch_begin + 1 branch_end + 0 events.
		//////////
		/// Forward rate vectors: Q.idot Q.ij
		//////////
		// Start at first event occurrence (first reported event is the beginning of the branch), and
		// last event occurrence (end of the branch). Will need to test for -1 as the event ID eventually.
		//
		// it: holds the Qij of the current change, time of current change, rate away of sequence AFTER change
		// prev2it: holds the time of previous change, rate away of sequence BEFORE change -> Qi. for calculating current probability 
		//it = branch_events.begin();
		//prev2it = branch_events.begin();
		end_event = branch_events.end(); end_event--;
		//++it;
		int xx = 0;
		//for (; it != end_event; ++it, xx++, ++prev2it)
		//	cerr << (*it)->Print_Event();
		it = branch_events.begin();
		prev2it = branch_events.begin();
		++it;
		for (; it != end_event; ++it, xx++, ++prev2it) {
			//cerr << (*it)->Print_Event();
			dt = (*it)->eventTime - (*prev2it)->eventTime;
			//cerr << "Qij: " << (*it)->Q.ij << "  Qi.: " << (*prev2it)->Q.idot << "  dt: " << dt << " ---> ";
			log_sum_forward += Forward_Step((*it)->Q.ij, (*prev2it)->Q.idot, dt);
			//cerr << " = " << log_sum_forward << endl;
		}
		t_z = (*prev2it)->eventTime;
		i_z_Qidot = (*prev2it)->Q.idot;
	} else {
		t_z = events.front()->eventTime;
		i_z_Qidot = branch_events.front()->Q.idot;
	}

	//////////
	/// Poisson(N=0|T-t_z), where N is the number of events, T is end of branch, t_z is time at state z.
	//////////
	// Probability of no events occurring until the end of the branch
	if (do_last_event) {
		//////////
		/// Finally, the probability of no events happening from the last event to the end of the
		/// branch.
		//////////
		dt = T - t_z;		// For full path, T. For subpath, next event (possibly T).
		log_sum_forward += Forward_NoEvent(i_z_Qidot, dt);
	}
	
	// assert(log_sum_forward < 0);	// For high dependency models, this may not be true, e.g., BL=1, depsup=2.
	cerr << "LOG_SUM_FORWARD:  " << log_sum_forward << endl;	
	return log_sum_forward;
}

double
PathProbability::IndependentForwardProbability (
											    list<eventTrack*> events,
											    short start_state,		// Only necessary for no subst events.
											    RateMatrix *rates
											   )
{
	double log_sum_forward = 0;
	double dt;
	double t_0, T, t_z;
	list<eventTrack*>::iterator prev2it, it, end_event;
	double i_z_Qidot;
	double Qij, Qidot;
	size_t from, to;

	//cerr << "Events coming into IFP: " << endl;
	//for (list<eventTrack*>::iterator it = events.begin(); it != events.end(); ++it)
	//	cerr << (*it)->Print_Event();

	t_0 = events.front()->eventTime;
	T = events.back()->eventTime;
	
	if (events.size() > 2) {
		it = events.begin();
		end_event = events.end(); end_event--;
		prev2it = events.begin();
		++it;
		t_z = t_0;

		for (; it != end_event; ++it, ++prev2it) {
			from = stateCharacters.find_first_of((*it)->size.at(0));	// NOTE: size, for a substitution, is the nature of the event.
			to = stateCharacters.find_first_of((*it)->size.at(1));	
			Qij = rates->Qij.at(from*numStates + to);
			Qidot = -rates->Qij.at(from*numStates + from);	// -Qii
			dt = (*it)->eventTime - (*prev2it)->eventTime;
			log_sum_forward += Forward_Step(Qij, Qidot, dt);
		}
		t_z = (*prev2it)->eventTime;
		from = stateCharacters.find_first_of((*prev2it)->size.at(1));
		i_z_Qidot = rates->Qij.at(from*numStates+from);
	} else {
		t_z = 0;
		i_z_Qidot = -rates->Qij.at(start_state*numStates+start_state);
	}
	
	//////////
	/// Finally, the probability of no events happening from the last event to the end of the
	/// branch.
	//////////
	dt = T - t_z;		// For full path, T. For subpath, next event (possibly T).
	log_sum_forward += Forward_NoEvent(i_z_Qidot, dt);

	return log_sum_forward;
}

void
PathProposal::remove_site_events(
								 int site_sampled
								)
{
	int i = 0;
	for (list<eventTrack*>::iterator it = epc_events.begin(); it != epc_events.end(); ++it, ++i) {
		if ((*it)->ID != -1)
			if ((*it)->MSA_positions.at(0) == site_sampled) {
				epc_events.erase(it);
			}
	}
}

list<eventTrack*> 
PathProposal::extract_site_events (int site)
{
	list<eventTrack*> site_events;
	for (list<eventTrack*>::iterator it = epc_events.begin(); it != epc_events.end(); ++it) {
		if ((*it)->ID == -1 || (*it)->MSA_positions.at(0) == site) {
			site_events.push_back(*it);
		}
	}
	return site_events;
}

list<eventTrack*>  
getBranchEvents(
				list<eventTrack*>::iterator begin,
				list<eventTrack*>::iterator end
			   )
{
	// Get single branch events, for now. NEED TO DO FOR ALL BRANCHES EVENTUALLY.

	list<eventTrack*>::iterator it = begin;
	list<eventTrack*> branch_events;
	while (it != end && (*it)->eventType != BRANCH_END) {
		branch_events.push_back(*it);
		++it;
	}
	branch_events.push_back(*it);
	return branch_events;
}

void 
PathProposal::write_path(ofstream& stream, int pathID)
{
	stream << "#" << pathID << endl;
	list<eventTrack*>::iterator it = epc_events.begin(); it++;
	for (; it != epc_events.end(); ++it) {
		stream << (*it)->write_path_event();
	}
}

void 
PathProposal::Print_Path_Events()
{
	for (list<eventTrack*>::iterator it = epc_events.begin(); it != epc_events.end(); ++it) {
		cerr << (*it)->Print_Event();
		cout << (*it)->Print_Event();
	}
}
