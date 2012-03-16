#include "forward_simulation.h"
#include "dependency.h"

void ForwardSimulation::EvolveNode(
								   inTree *iTree, 
								   TNode *anc, 	// Currently set in this function.
								   TNode *des, 
								   int inNumSites, 
								   list<eventTrack*> *events, 
								   seqGenOptions *options
			   					  )
{ 
	bool chang_and_benner = false;
	ofstream step_out;
	string trs_outfile;

    des->anc=anc;
	des->atEpochTime = anc->atEpochTime;

	initializeMotifPositions(des);
	initializeCladeParameters(iTree->my_tree, des);
	initializeBranchRun(iTree, des, anc->seq_evo.size(), options);

	//des->nodeEnv->Print_Environment();
	des->branch->rates->Qij = anc->branch->rates->Qij;

	//////////
	/// 050311: Breaking chang and benner completely!!!!!!!
	/// The Chang and Benner model does not reconcile with the Gillespie model. Run
	/// the C&B model as is done in iSGv1.
	//////////
	chang_and_benner = (des->nodeEnv->P_ins_ == 0 && des->nodeEnv->P_del_ == 0 && des->nodeEnv->indelFlag);
	if (options->simulation_step_type == UNIFORMIZATION && !chang_and_benner) {
		des->setRateAway(UNIFORMIZATION);
		gillespie(iTree, des, des->branch->anc_length(), events, options, UNIFORMIZATION, FORWARD_EVOLUTION);
	} else {
		//des->setInteractions();
		gillespie(iTree, des, des->branch->anc_length(), events, options, TIME_RELATIVE_STEPS, FORWARD_EVOLUTION);
	}

	//////////
	/// If not a leaf, proceed down both subtrees.
	//////////
    if (des->tipNo==-1) { 
    	//des->setSitePointers("EvolveNode");
        EvolveNode(iTree, des, des->branch1, inNumSites, events, options);
    	EvolveNode(iTree, des, des->branch2, inNumSites, events, options);
    }
	des->Remove_Objects();
} 

void ForwardSimulation::EvolveSequences(
					 					inTree *iTree, 
					 					list<eventTrack*> *events, 
					 					seqGenOptions *options
									   )
{
	num_subst = 0;

	//////////
	/// FORWARD SIMULATION SPECIFIC SETUP
	//////////
	cerr << "Forward Category breakdown: " << endl;	//XOUT
	vector<short> cats (4,0);
	vector<Site>::iterator root_it;
	for (root_it = iTree->my_tree->root->seq_evo.begin(); root_it != iTree->my_tree->root->seq_evo.end(); ++root_it) 
		cats.at((*root_it).returnCategory())++;

	double overall_rate_multiplier = 0;
	for (int i = 0; i < 4; i++) {
		overall_rate_multiplier += cats.at(i) * iTree->my_tree->root->nodeEnv->catRate[i];
		cerr << "cat" << i+1 << " occupancy: " << cats.at(i) << endl;//XOUT
	}

	overall_rate_multiplier /= iTree->my_tree->root->seq_evo.size();
	cerr << "Overall multiplier: " << overall_rate_multiplier << endl;	//XOUT
	//cout << overall_rate_multiplier;		//XOUT

	if ( order_3_markov ) {
		iTree->my_tree->dep.front()->context.set_sequence_indices(iTree->my_tree->root);
		//iTree->my_tree->dep.front()->context.set_Qmat(iTree->my_tree->root);
		cerr << "EVOLVING TO EQUILIBRIUM>>>>>" << endl;
		evolve2Equilibrium(iTree, events, options);
		cerr << "EVOLVED TO EQUILIBRIUM>>>>>" << endl;
		//exit(0);
	}
	
	EvolveNode(iTree, iTree->my_tree->root, iTree->my_tree->root->branch1, iTree->partitionLength, events, options);
    EvolveNode(iTree, iTree->my_tree->root, iTree->my_tree->root->branch2, iTree->partitionLength, events, options);
    if (!iTree->my_tree->rooted) {
		cerr << "Unrooted tree??" << endl; 
		exit(EXIT_FAILURE);
	}

	//////////
	/// Remove objects from root sequence.
	//////////
	// The following line of code blows iSG up when a motif sequence is followed by a RANDOM,
	// sequence with the -1 option unset.
//	if ( iTree->rootSeqType != RANDOM || !iTree->randomInvariableAssignment )
//		iTree->my_tree->root->Remove_Objects(true);
	//
	// Regardless of the motif situation, the general varSites still need to be removed. This is
	// done automatically if the if-stmt is true.
//	else iTree->my_tree->root->Remove_varSites();
	//
	//////////
	
	cerr << "EVENTS:" << eventNo << endl;
	PathProposal *path = new PathProposal();
	cout << path->P_path.ForwardProbability(*events, true) << endl;

	cout << "EVENTS:" << eventNo << endl;
} 

double ForwardSimulation::calcGillespieLambda(
						   					  TTree *tree,
						   					  TNode *des, 
						   					  double *I, 
						   					  double *D, 
						   					  double *S,
						   					  int simulation_type,
						   					  unsigned int event_site
						  					 )
{
	double lambda = 0;
	double del_freq_by_size=0, del_freq=0, u_del=0;
	int ins_L, del_L;
	int i = 0;

	if (simulation_type == UNIFORMIZATION) {
		for (vector<double>::iterator it = des->nodeEnv->delete_lengthDistribution.begin(); it != des->nodeEnv->delete_lengthDistribution.end(); ++it, i++) {
			del_freq += (*it);
			del_freq_by_size += (*it)*i;
		}

		if (del_freq == 0) u_del = 0;
		else u_del = del_freq_by_size / del_freq;

		calcTrials(des, &ins_L, &del_L);
		*I = des->nodeEnv->P_ins_ * (ins_L + 1);
		*D = des->nodeEnv->P_del_ * (u_del - 1 + del_L);
		*S = des->seq_evo.size() * des->rate_away_site_width;
	} else {
		if (order_3_markov) *S = des->calculateForwardRateAwayFromSequence__order3Markov(tree, event_site);
		else {
			for (vector<double>::iterator it = des->nodeEnv->delete_lengthDistribution.begin(); it != des->nodeEnv->delete_lengthDistribution.end(); ++it, i++) {
				del_freq += (*it);
				del_freq_by_size += (*it)*i;
			}

			if (del_freq == 0) u_del = 0;
			else u_del = del_freq_by_size / del_freq;

			calcTrials(des, &ins_L, &del_L);
			*I = des->nodeEnv->P_ins_ * (ins_L + 1);
			*D = des->nodeEnv->P_del_ * (u_del - 1 + del_L);
			*S = des->calculateForwardRateAwayFromSequence();
		}
	}

//	cerr << endl;
//	cerr << "*I = " << des->nodeEnv->P_ins_ << " * (" << ins_L << " + 1) = " << *I << endl;
//	cerr << "*D = " << des->nodeEnv->P_del_ << " * (" << u_del << "- 1 + " << del_L << ") = " << *D << endl;
//	cerr << "*S = " << des->seq_evo.size() * des->rate_away_site_width << " = " << *S << endl;

	lambda = *I + *D + *S;
	return lambda;
}
               
void ForwardSimulation::gillespie(
			   					  inTree *iTree, 
			   					  TNode *des,
								  double branch_len,
								  list<eventTrack*> *events,
								  seqGenOptions *options,
			   					  int simulation_type,			/// may be options->simulation_step_type, but who cares.
			 				 	  bool evolving_to_equilibrium
			  )
{
	double lambda_T, exp_mean, I=0, D=0, S=0;
	bool success;
	int event_size = 0, action = INSERT;
	double RN;
	string event;
	bool trackEvent = false;
	double max_path_length_inverse = 1.0 / iTree->global_max_path_length;
	double branch_len_inverse = 1.0 / branch_len;
	int event_site;		/// Placeholder for substitution routine, at the moment. ///
	bool profile=false;
		eventTrack *begin_event, *end_event;

	//////////
	// Indel routine	
	//  * Should codonRates affect the probability of an indel occurring?
	//////////
	lambda_T = calcGillespieLambda(iTree->my_tree, des, &I, &D, &S, simulation_type, -1);
	cerr << "lambda_T(sim): " << lambda_T << endl;		//XOUT
	exp_mean = 1.0/lambda_T;
	if (!evolving_to_equilibrium) {
		// Forward simulation, tree & branches fully made, thus, begin and end events do not need to know the branch lengths (last variable = 0).
		// des->anc->DistanceFromRoot: beginning of branch starts at ancestor's DistanceFromRoot.
		begin_event = branch_terminal_event(des->mytipNo, BRANCH_BEGIN, des->anc->DistanceFromRoot, des->bipartition, 0, branch_len); 
		begin_event->assign_Q(des, des->seq_evo.at(0), SUBSTITUTION, 1);
		(*events).push_back(begin_event);
	}
	for (double dt = rand_exp(exp_mean); dt <= branch_len; dt += rand_exp(exp_mean)) {
		trackEvent = false;
		success = false;
		RN = (double)rndu();
		if ( RN < (I / lambda_T) ) {	// Insertion.
        	event_size=Find_Indel_Size(des->nodeEnv->insert_lengthDistribution,des->nodeEnv->maxIndel);	
			action = INSERT;
			if ( (success=Insert(iTree->my_tree,des,event_size)) ) {
				event = to_string(event_size);
				if (options->events_to_track[seqGenOptions::track_insertion] && !evolving_to_equilibrium) trackEvent = true;
			}
		} else if (RN <= (I + S) / lambda_T ) { // Substitution
			//////////
			/// If we sent rate_away matrix, then we simply need to select a substitution from
			/// the probabilities that are calculated.
			//////////
			event.clear();
			success = substitute(
								 iTree->my_tree, 
								 des, 
								 &event_site, 
								 simulation_type, 
								 event, 
								 options->events_to_track[seqGenOptions::track_substitution], 
								 S
								);
			if (options->events_to_track[seqGenOptions::track_substitution]) trackEvent = true;
			action = SUBSTITUTION;
			num_subst++;
		} else { // Deletion.
        	event_size=Find_Indel_Size(des->nodeEnv->delete_lengthDistribution,des->nodeEnv->maxIndel);
			action = DELETE;
			if ( (des->seq_evo.size()-event_size) > 0) {
				//des->setSitePointers("gillespie3");
	    		if ( (success=Delete(iTree->my_tree,des,&event_size)) ) {
					event = to_string(event_size);
					if (options->events_to_track[seqGenOptions::track_deletion]) trackEvent = true;
				}
			}
		}
		if (success /*&& trackEvent*/) {
			eventTrack *new_event;
			cerr << "FWD subst pos: " << event_site << ":  "; //XOUT
			if (!evolving_to_equilibrium) {
				//////////
				/// Calculate the relative time that the event occurs at.
				//////////
				// * anc->trDistanceFromRoot: The scaled distance from the root
				// * (dt/indel_len): Percentage of the branch that has been simulated
				// * des->branch0_time_relative_length: Time rel length of the branch being simulated.
				// * iTree->global_max_path_length: Scalar for the sum of above terms to set time of occurrence between 0 and 1.
				//////////
				//if (order_3_markov) {
				des->atEpochTime = des->anc->DistanceFromRoot + dt;
				//} else { des->atEpochTime = (des->anc->trDistanceFromRoot+(dt*branch_len_inverse)*des->branch->branch0_time_relative_length)*max_path_length_inverse;}
				new_event = new eventTrack(
										   eventNo,
										   action,
										   des->atEpochTime,
										   des->bipartition,
										   event,
										   branch_len
										  );
				//////////
				/// In an indel, the event_site is inconsequential (Qij not calculated). Also, in deletion
				/// at the end of the sequence, site no longer exists.
				//////////
				if (action == DELETE || action == INSERT) event_site = 0;
				new_event->assign_Q(des, des->seq_evo.at(event_site), action, 1);
			}
			if (action == DELETE || action == INSERT) event_site = -1;		// Temp. Flags to recalculate TauIJ for entire sequence.
			iTree->my_tree->dep.front()->context.reset_sequence_indices(des, event_site, event);
			lambda_T = calcGillespieLambda(iTree->my_tree, des, &I, &D, &S, simulation_type, event_site);
			if (!evolving_to_equilibrium) {
				new_event->Q.idot = lambda_T;
				(*events).push_back(new_event);
				eventNo++;
				if (profile) cerr << "gillespie::post_change_sequence: " << des->printSequence() << endl;	//XOUT
			}
		}
		if (!evolving_to_equilibrium) {	//XOUT
			cerr << "---------- Time left: " << branch_len-dt << "/" << branch_len << "----------"; //XOUT
			cerr << "Qi.: " << des->evolvingSequence->Qidot; //XOUT
//			cout << "DT:" << dt << endl;	//XOUT
//			cout << "LAMBDA_T:" << lambda_T << endl;			//XOUT
//			cout << "Qi.:" << des->evolvingSequence->Qidot << endl;	//XOUT
//			cout << "Qi.|k:" << des->evolvingSequence->Qidot_k__T__ << endl;	//XOUT
//			cout << "SITE_RATE_AWAY:" << des->seq_evo.at(changed_site).site_rate_away.back() << endl; //XOUT
		} else {	//XOUT
			cerr << "---------- E2E->Time left: " << branch_len-dt << "/" << branch_len << "----------"; //XOUT
			cerr << "Qi.: " << des->evolvingSequence->Qidot; //XOUT
			cerr << endl;
		}

		exp_mean = 1.0/lambda_T;

		if (!evolving_to_equilibrium) {
			//cout << endl;	//XOUT
			cerr << endl;
		}
	}

	if (!evolving_to_equilibrium) {
		// Forward simulation, tree & branches fully made, thus, begin and end events do not need to know the branch lengths (last variable = 0).
		end_event = branch_terminal_event(des->mytipNo, BRANCH_END, des->DistanceFromRoot, des->bipartition, 0, branch_len); 
		end_event->assign_Q(des, des->seq_evo.at(0), SUBSTITUTION, 1);
		(*events).push_back(end_event);
	}
}

void ForwardSimulation::evolve2Equilibrium(
						inTree *iTree,
					 	list<eventTrack*> *events, 
					 	seqGenOptions *options
					   )
{
	double lead_in_branch_length = 10;

	iTree->my_tree->root->evolvingSequence->print_sequence();
	options->events_to_track[seqGenOptions::track_substitution] = false;
	gillespie(iTree, iTree->my_tree->root, lead_in_branch_length, events, options, TIME_RELATIVE_STEPS, EVOLVING_TO_EQUILIBRIUM);
	options->events_to_track[seqGenOptions::track_substitution] = true;
	iTree->my_tree->root->evolvingSequence->print_sequence();
	num_subst = 0;
	for (vector<Site>::iterator it = iTree->my_tree->root->seq_evo.begin(); it != iTree->my_tree->root->seq_evo.end(); ++it) {
		(*it).interactions.clear();
	}

	cerr << "evolve2Equilibrium::Number of substitutions: " << num_subst << endl;
	cerr << "                    Number of recorded events: " << (*events).size() << " " << eventNo << endl;
	cerr << "evolve2Equilibrium::sequence: " << iTree->my_tree->root->printSequence() << endl;
}

template <class T> string to_string (
								     const T& t
								    )
{
	stringstream ss;
	ss << t;
	return ss.str();
}
 
int quick_lookup(const std::vector<short>& letters)
{
	int temp=1;
	for(int i=0;i<letters.size();i++)
		temp *= letters[i];
	return temp;
}

template <typename T, typename U>
class quick_map
{
	std::vector<U> data;
public:
// Index by T
	const U& operator[](const T& letters) const {
		int index = quick_lookup(letters);
		return data[index];
	}

	U& operator[](const T& letters)       {
		int index = quick_lookup(letters);
		return data[index];
	}

	void set_value(const T& letters)
	{
		int index = quick_lookup(letters);
		if (index < data.size())
			data.resize(index+1);
		return data[index];
	}

// Index by int
	const U& operator[](int index) const {
		return data[index];
	}


	U& operator[](int index)
	{
		return data[index];
	}

// Handle sizing
	void resize(unsigned s) {data.resize(s);}

	unsigned size() const {return data.size();}

//Initialize
	quick_map() {}
	quick_map(int size):data(size) { }
};

int f()
{
	quick_map< std::vector<short>, int > M(100);
	std::vector<short> v;
	v.push_back(0);
	v.push_back(1);

//	quick_map< std::vector<int>, int > M2(100);


	int out = M[v];
	return out;
}