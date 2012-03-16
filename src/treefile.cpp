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

#include "limits.h"
#include "treefile.h"
#include "prosite_motif.h"
#include "dependency.h"
#include "stats.h"

char treeErrorMsg[256];
int treeError;

//vector<double> save_Pr (64,0);

extern bool round1;
extern double rate_away_same, rate_away_diff;
extern double sum_Pij_same, sum_Pij_diff;
extern int changed_site;
extern int prev_state;
extern int changed_site;
extern bool optimize;
extern bool Pc, Qd, nij;

int TauIJ_calls = 0;
int TauIJ2_calls = 0;

/* functions */

TNode::TNode(
			 TTree *tree
			) 
	  : evolvingSequence(NULL),
	    branch0(NULL),
	    branch1(NULL),
	    branch2(NULL),
	    anc(NULL),
	    trDistanceFromRoot(0.0),
	    DistanceFromRoot(0.0),
	    tipNo(-1),
	    atEpochTime(0.0),
	    nodeEnv(NULL),
		node_name(""),
		branch (new Branch()),
	    mytipNo(-1)
{
	clade_label.clear();
	tree->numNodes++;
	bipartition.clear();
	addGeneral_varSites();
} 

TNode::TNode() 
	  : evolvingSequence(NULL),
	    branch0(NULL),
	    branch1(NULL),
	    branch2(NULL),
	    anc(NULL),
	    trDistanceFromRoot(0.0),
	    DistanceFromRoot(0.0),
	    tipNo(-1),
	    atEpochTime(0.0),
		node_name(""),
	    nodeEnv(NULL),
		branch (new Branch()),
	    mytipNo(-1)
{
	clade_label.clear();
	bipartition.clear();
	addGeneral_varSites();
} 

void TTree::setTransitionProbabilities()
{
	int position = 0;
	int branch_length_scalar = 1;

//	cerr << "At root: BLs: " << root->branch->length1 << " " << root->branch->length2 << endl;

	root->branch1->branch->rates->setPij(
	 							root->seq_evo.at(position), 
								root->branch->length1 * branch_length_scalar, 
								root->nodeEnv->rateHetero
						 	   );

//	cerr << "root->branch1->branch->rates->setPij(," << root->branch->length1*branch_length_scalar << ",)" << endl;

	root->branch1->setBranchTransitionProbabilities();
	root->branch2->branch->rates->setPij(
	 							root->seq_evo.at(position), 
								root->branch->length2 * branch_length_scalar, 
								root->nodeEnv->rateHetero
						 	   );
	root->branch2->setBranchTransitionProbabilities();
}

void TNode::setBranchTransitionProbabilities()
{
	int position = 0;
	int branch_length_scalar = 1;

	if (tipNo == -1) {
		branch1->branch->rates->setPij(
		 					  seq_evo.at(position), 
							  branch->length1 * branch_length_scalar, 
							  nodeEnv->rateHetero
							 );
		branch1->setBranchTransitionProbabilities();
		branch2->branch->rates->setPij(
		 					  seq_evo.at(position), 
							  branch->length2 * branch_length_scalar, 
							  nodeEnv->rateHetero
							 );
		branch2->setBranchTransitionProbabilities();
	}
}

void TTree::setStateLikelihoods( 
								double branch_length_scalar
							   )
{
	size_t position = 0;
	int cat;

	//root->branch->rates->printPij();

	for (vector<Site>::iterator it = root->seq_evo.begin(); it != root->seq_evo.end(); ++it, position++) {
		cat = (*it).returnCategory();
		root->branch1->setStateLikelihood(position, cat, branch_length_scalar);
		root->branch2->setStateLikelihood(position, cat, branch_length_scalar);

		root->seq_evo.at(position).L_i.calculateStateLikelihood(
																root, 
																cat,
																position,
																branch_length_scalar
															   );
		//
		size_t pi_position = 0;
		for (vector<double>::iterator it = root->seq_evo.at(position).L_i.Li_xi_.begin(); it != root->seq_evo.at(position).L_i.Li_xi_.end(); ++it, pi_position++) {
			(*it) *= root->branch->rates->pi.at(pi_position);
		//	cerr << (*it) << " "; 
		}
		//cerr << endl;
	}
}

void TNode::setStateLikelihood( 
						       size_t position,
						       int category,
						       double branch_length_scalar
						      )
{

	//////////
	/// Leaf node likelihood setting
	//////////
	// Below: If a sequence position is unset, it will have the state -1. If the position is set,
	// then there is a 100% probability that that position will be the value of it's current state.
	if ( seq_evo.at(position).isSet() == true ) {
		bool a_one_set = false;
		for (int i = 0; i < numStates; i++)
			if (i == seq_evo.at(position).returnState()) 
				seq_evo.at(position).L_i.Li_xi_.at(i) = 1;
			else 
				seq_evo.at(position).L_i.Li_xi_.at(i) = 0;
	}
	//
	//////////

	//cerr << "Node: " << mytipNo << endl;
	//for (vector<double>::iterator it = seq_evo.at(position).L_i.Li_xi_.begin(); it != seq_evo.at(position).L_i.Li_xi_.end(); ++it) {
	//	cerr << (*it) << " "; 
	//}
	//cerr << endl;

	// Stopping condition
	if (tipNo == -1) {
		//////////
		/// branch1
		//////////
		/// Make sure that both branches nodes have their likelihood vectors set:
		branch1->setStateLikelihood(position, category, branch_length_scalar);
		branch2->setStateLikelihood(position, category, branch_length_scalar);

		if ( !seq_evo.at(position).isSet() )
			seq_evo.at(position).L_i.calculateStateLikelihood(
															  this, 
															  category, 
															  position, 
															  branch_length_scalar
															 );
	} 
}

//////////
/// Calculate the probability of each category given the observed data and the model.
///
/// Calculate P( C=j | X_1,X_2,X_3,\theta ) and set categories.
///              P( C=j,X_1,X_2,X_3 | \theta )
///			  = -------------------------------
///				   P( X_1,X_2,X_3 | \theta )
///
///              P( X_1,X_2,X_3 | C=j,\theta ) x P( C=j | \theta )
///			  = ---------------------------------------------------
///				   P( X_1,X_2,X_3 | \theta )
///
///              P( X_1,X_2,X_3 | C=j,\theta ) x 1 / #_cats
///			  = --------------------------------------------
///				   P( X_1,X_2,X_3 | \theta )
///
/// \theta = { \tau, branch_lengths, substitution model}
//////////
void TTree::setCatLikelihoods(
							  double branch_length_scalar
							 )
{
	vector<vector<double> > Li_per_cat;
	size_t position = 0;
	int num_cats = root->branch->rates->num_categories;
	vector<double> b1_Li, b1_Pij;
	vector<double> b2_Li, b2_Pij;
	vector<double> Pr (root->seq_evo.size(), 0);
	vector<double> my_Li (4,0);
	double total_Li = 0;
	vector<double> cat_likelihood;
	vector<double>::iterator pi_it;
	double sum_cats;

	cerr << "DISCRETE GAMMA RATES CALCULATION::::" << endl;

	//////////
	/// Start from root, calculate likelihoods for the left node and right node recursively
	/// for each site
	//////////
	for (vector<Site>::iterator it = root->seq_evo.begin(); it != root->seq_evo.end(); ++it, position++) {
		//cerr << endl << "SITE " << position << ":" << endl;
		//////////
		/// Calculate P( X_1,X_2,X_3 | \theta )
		/// The site likelihoods; the probability of the observed data given the model.
		//////////
		sum_cats = 0;
		cat_likelihood.assign(num_cats, 0);
		for (int cat = 0; cat < num_cats; cat++) {
			//cerr << "Cat " << cat << ": " << endl;
			//////////
			/// Calculate conditional likelihoods of left and right descendant nodes.
			//////////
			//
			// Call recursive routines to set left node conditional likelihoods
			b1_Li = root->branch1->setCatLikelihood(position, cat);
			// and right node conditional likelihoods
			b2_Li = root->branch2->setCatLikelihood(position, cat);
			//cerr << " ROOT_CALC: " << endl;
			// and utilize these to calculate the likelihoods at the current node.
			// ** P( X_1,X_2,X_3 | C=j,\theta ) in the numerator.
		    my_Li =(*it).L_i.calculateCatLikelihood(
											 		root, 
													cat, 
													position,
													b1_Li,
													b2_Li
												   );

			//////////
			/// Test code to calculate the proper likelihoods.
			//////////
			//cerr << "b1_node_likelihood, cat " << cat
			//	 << "  A: " << my_Li.at(0) 
			//	 << "  C: " << my_Li.at(1) 
			//	 << "  G: " << my_Li.at(2) 
			//	 << "  T: " << my_Li.at(3) 
			//	 << endl;

			//////////
			/// Multiply by stationary frequencies and normalize by number of categories.
			//////////
			pi_it = root->branch->rates->pi.begin();
			double this_cat_sum = 0;
			for (vector<double>::iterator jt = my_Li.begin(); jt != my_Li.end(); ++jt, ++pi_it) {
				//////////
				///      P( X_1,X_2,X_3 | C=j,\theta ) x 1 / #_cats
				///	  = --------------------------------------------
				///		   		P( X_1,X_2,X_3 | \theta )
				//////////
				// (*jt) is the root state likelihood vector:
				//   * P( X_1,X_2,X_3 | X_R=i,C=j,\theta )
				//
				// Multiply by ¹_i
				(*jt) *= (*pi_it);
				// and 1 / #_cats
				(*jt) *= (1.0 / root->branch->rates->num_categories);
				//
				// and add this to the current category likelihood.
				cat_likelihood.at(cat) += (*jt);
				// sum_cats is the probability of observed data given the model (X_R and C are summed
				// away).
				sum_cats += (*jt);
			}
		}

		double RND = rndu() * sum_cats;
		vector<double>::iterator jt = cat_likelihood.begin();
		double running_sum = (*jt);
		++jt;
		size_t category = 0;
		while (RND > running_sum) {
			running_sum += (*jt);
			++jt;
			category++;
		}

//		cerr << sitePattern(position) << ":   "
//			 << ((category == 0) ? "*" : "") << cat_likelihood.at(0) << " " 
//			 << ((category == 1) ? "*" : "") << cat_likelihood.at(1) << " " 
//			 << ((category == 2) ? "*" : "") << cat_likelihood.at(2) << " " 
//			 << ((category == 3) ? "*" : "") << cat_likelihood.at(3) 
//			 << endl;

		(*it).setCategory(category);
	}
	
	cerr << "Path Category breakdown: " << endl;
	vector<short> cats (4,0);
	vector<Site>::iterator anc_it;
	for (anc_it = root->seq_evo.begin(); anc_it != root->seq_evo.end(); ++anc_it) 
		cats.at((*anc_it).returnCategory())++;

	double overall_rate_multiplier = 0;
	for (int i = 0; i < root->branch->rates->num_categories; i++) {
		overall_rate_multiplier += cats.at(i) * root->branch->rates->catRate.at(i);
		cerr << "cat" << i+1 << " occupancy: " << cats.at(i) << endl;
	}

	overall_rate_multiplier /= root->seq_evo.size();
	cerr << "Overall multiplier: " << overall_rate_multiplier << endl;
	cout << overall_rate_multiplier << endl;
	//exit(0);
}

vector<double> TNode::setCatLikelihood(
									   size_t position,
									   int category
									  )
{
	vector<double> my_Li (numStates, 0);	/// This is the returned likelihood from me.
	vector<double> b1_Li, b2_Li;			/// Keep track of my_Li of each child.

	//cerr << "starting node " << node_name << "(" << mytipNo << ")" << endl;

	//////////
	/// If this site is set, the likelihood of the state at the site is set to 1, other states = 0.
	//////////
	if (seq_evo.at(position).isSet() == true) 
		my_Li.at(seq_evo.at(position).returnState()) = 1;
	if (tipNo == -1) {
		//////////
		/// Internal node, calculate left and right branch likelihoods.
		//////////
		b1_Li = branch1->setCatLikelihood(position, category);
		b2_Li = branch2->setCatLikelihood(position, category);
		//
		// If this site is unset previously, then set the conditional likelihoods.
		if ( !seq_evo.at(position).isSet() ) {
			my_Li 
			= seq_evo.at(position).L_i.calculateCatLikelihood(
														      this, 
														      category, 
														      position,
														      b1_Li,
														      b2_Li
														     );
		}
	}

	//cerr << "b1_node_likelihood, cat " << category
	//	 << "  A: " << my_Li.at(0) 
	//	 << "  C: " << my_Li.at(1) 
	//	 << "  G: " << my_Li.at(2) 
	//	 << "  T: " << my_Li.at(3) 
	//	 << endl;

	//cerr << "Node: " << mytipNo << endl;
	//for (vector<double>::iterator it = my_Li.begin(); it != my_Li.end(); ++it) {
	//	cerr << (*it) << " "; 
	//}
	//cerr << endl;

	return my_Li;
}

string 
TTree::sitePattern(
				   int position
				  )
{
	string pattern = "";

	pattern += root->branch1->subtreePattern(position);
	pattern += root->branch2->subtreePattern(position);

	return pattern;
}

string 
TNode::subtreePattern(
					   		 int position
					  		) 
{
	string pattern = "";
	
	if (tipNo == -1) {
		pattern += branch1->subtreePattern(position);
		pattern += branch2->subtreePattern(position);
	} else pattern = stateCharacters.at(seq_evo.at(position).returnState());

	return pattern;
}					  

void 
TTree::calculateJCLikelihoods()
{
	double X, Y;
	double branch_length = 0.5;
	double gamma_rate;
	int num_same, ptr;
	vector<double> Pr (numStates_cubed, 0);
	double total_Pr = 0;
	vector<double> cat_Pr (TTree::root->TNode::branch->Branch::rates->RateMatrix::num_categories, 0);
	vector<double> site_Li (numStates_cubed, 0);
	unsigned int squared[numStates], power_1[numStates];

	for (unsigned int i = 0; i < numStates; i++) {
		squared[i] = i * numStates_squared;
		power_1[i] = i * numStates;
	}
	//root->branch->rates->setPij(root->seq_evo.front(), branch_length, DiscreteGammaRates);	
	for (int x = 0; x < TTree::root->TNode::branch->Branch::rates->RateMatrix::num_categories; x++) {
		gamma_rate = root->branch->rates->catRate.at(x);
		//cerr << "Category " << x << " rate: " << gamma_rate << endl;
		X = 0.25 + 0.75 * exp( -4.0/3.0 * gamma_rate * branch_length);
		Y = 0.25 - 0.25 * exp( -4.0/3.0 * gamma_rate * branch_length);

		//////////
		/// Calculate for all possible site patterns
		//////////
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++) {
				for (int k = 0; k < numStates; k++) {
					num_same = ( (i==j) ? 1 : 0 ) + ( (i==k) ? 1 : 0 ) + ( (j==k) ? 1 : 0 );
					ptr = squared[i]+power_1[j]+k;
					switch(num_same)
					{
						case 0:
							//////////
							/// All different: XYZ
							//////////
							Pr.at(ptr) = 0.75*(X*Y*Y) + 0.25*(Y*Y*Y);
							break;
						case 1:
							//////////
							/// Two are the same: XXY, XYX, YXX
							//////////
							Pr.at(ptr) = 0.25*(X*X*Y) + 0.25*(X*Y*Y) + 0.5*(Y*Y*Y);
							break;
						case 3:
							//////////
							/// All same: XXX
							//////////
							Pr.at(ptr) = 0.25*(X*X*X) + 0.75*(Y*Y*Y);
							break;
						default:
							cerr << "Impossible to have more than 3; num_same = " << num_same << endl;
							exit(EXIT_FAILURE);
					}


					cerr << "cat(" << x << ") " 
						 << stateCharacters.at(i) 
						 << stateCharacters.at(j) 
						 << stateCharacters.at(k) 
						 << "  num_same: " << num_same
						 << "  Li: " << Pr.at(ptr)
					;

					cerr << endl;
				}
			}
		}

		vector<double>::iterator st = site_Li.begin();
		for (vector<double>::iterator it = Pr.begin(); it != Pr.end(); ++it, ++st) {
			(*st) += (*it);
			cat_Pr.at(x) += (*it);
			total_Pr += (*it);
		}
	}

	int x = 0;
	cerr << "Site likelihoods: " << endl;
	vector<double>::iterator st = site_Li.begin();
	int sites;
	int i, j, k;
	for (st = site_Li.begin(); st != site_Li.end(); ++st, x++) {
		i = x / 16;
		sites = x - i*16;
		j = sites / 4;
		sites -= j*4;
		k = sites;
		
		(*st) /= 4.0;

		cerr << stateCharacters.at(i) 
			 << stateCharacters.at(j) 
			 << stateCharacters.at(k) 
			 << " " << (*st) << endl;
	}
	

	cerr << "Sanity checks: " << endl;
	for (int x = 0; x < TTree::root->TNode::branch->Branch::rates->RateMatrix::num_categories; x++) {
		cat_Pr.at(x) *= 0.25;
		cerr << "  cat " << x << " Likelihood: " << cat_Pr.at(x) << endl;
	}
	total_Pr *= 0.25;
	cerr << "  total Likelihood = " << total_Pr << endl;
}

void 
TTree::sample_root_sequence()
{
	double sum;
	double RND;
	vector<double>::iterator jt;
	
	for (vector<Site>::iterator it = root->seq_evo.begin(); it != root->seq_evo.end(); ++it) {
		sum = 0;
		for (jt = (*it).L_i.Li_xi_.begin(); jt != (*it).L_i.Li_xi_.end(); ++jt)
			sum += (*jt);

		RND = (double)rndu() * sum;
		jt = (*it).L_i.Li_xi_.begin();
		size_t state;
		sum = 0;
		for (state = 0; state < numStates; state++, ++jt) {
			sum += (*jt);
			if (RND < sum) break;
		}

		if (state == 4) {
			cerr << "TTree::sample_root_sequence(): probably did not use option \"-o f\" in forward simulation, or input file is not in fasta format." << endl;
			exit(EXIT_FAILURE);
		} else {
			(*it).setState(state);
		}

		//sum = 0;
		//for (jt = (*it).L_i.Li_xi_.begin(); jt != (*it).L_i.Li_xi_.end(); ++jt, sum++) {
		//	if (sum == state) cerr << "*";
		//	cerr << (*jt) << "  ";
		//} cerr << endl;
		
	}
}

void 
TNode::setRateAway(
				   int step_type
				  )
{
	double max_Qii = 0, max_gamma = 0;
	//////////
	/// First discover the maximum values of Qii and gamma.
	//////////
	if (step_type == UNIFORMIZATION) {
		for (int i = 0; i < numStates; i++) 
			if (-branch->rates->Qij[i*numStates+i] > max_Qii) 
				max_Qii = -branch->rates->Qij[i*numStates+i];
		for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it) {
			if ( (*it).returnGamma() == 0 ) {
				max_gamma = nodeEnv->catRate[nodeEnv->numCats-1];
				it = seq_evo.end()-1;
			} else {
				if ( (*it).returnGamma() > max_gamma ) 
					max_gamma = (*it).returnGamma();
			}
		}
		rate_away_site_width = max_Qii * max_gamma;
		if ( (*seq_evo.begin()).returnGamma() != 0 ) rate_away_site_width *= RATE_AWAY_BUFFER;
		for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it)
			(*it).setSiteRateAway
				(
				 rate_away_site_width, 
				 (
				   (nodeEnv->rateHetero == DiscreteGammaRates)
				   ? nodeEnv->catRate[(*it).returnCategory()]
				   : (*it).returnGamma()
				 ),
				 branch->rates
				);
	} else {
		//////////
		/// Calculate the rate away, instead of des->seq_evo.size(), e.g., if A:
		///  Rate away = Q_{AC} + Q_{AG} + Q_{AT}
		///  Presumably, this means that the transition probabilities need not be included in calculation
		//////////
		double prev = 0;
		for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it)
			prev = (*it).setSiteRateAway(branch->rates->Qij, branch->rates);
	}
}

double 
TNode::calculateForwardRateAwayFromSequence( void )
{
	double R_D = 0;
	double rate_away;
	bool is_dependent;
	int times_thru_D = 0;
	double forward_sum_away = 0;
	
	for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it) {
		is_dependent = false;
		(*it).forward_rate_away.assign(numStates, 0);
		for (int i = 0; i < numStates; i++) {
			(*it).forward_rate_away.at(i) = branch->rates->pi.at(i);
			forward_sum_away += (*it).forward_rate_away.at(i);
		}
		if ( !(*it).isSet() ) {
			cerr << "TNode::calculateForwardRateAwayFromSequence:" << endl;
			cerr << "  In order to do calculate the rate away from a sequence, all sites in the sequence "
				 << "must be set." << endl;
			exit(EXIT_FAILURE);
		}
		for (list<siteDependencies*>::iterator jt = (*it).interactions.begin(); jt != (*it).interactions.end(); ++jt) {
			if ((*jt)->isDependentInteraction()) {
				rate_away = (*it).setSiteRateAway(branch->rates->Qij, branch->rates);	
				is_dependent = true;
				times_thru_D++;
			}
		}
		if (!is_dependent) rate_away = (*it).setSiteRateAway(branch->rates->Qij, branch->rates);

		R_D += rate_away;
	}

	evolvingSequence->Qidot_k__T__ = -1;	// Forward simulation doesn't do this.
	evolvingSequence->Qidot = forward_sum_away;

	return R_D;
}

void 
TNode::set_site_window(
					   int order, 
					   int *start_site, 
					   unsigned int *end_site
					  )
{
	if (*start_site == -1) { 
		*start_site = 0; 
		*end_site = seq_evo.size(); 
	} else if ( *start_site > seq_evo.size() - (order+1) ) {
		*end_site = seq_evo.size();
		*start_site -= order;
	} else {
		*end_site = *start_site + order+1;
		*start_site -= order;
		if (*start_site < 0) *start_site = 0;
	}
}

double 
TNode::calculateForwardRateAwayFromSequence__order3Markov( 
														  TTree *tree, 
														  int event_site 
														 )
{
	double R_D = 0;
	unsigned int end_site;

	if (order_3_markov) {
		set_site_window(tree->dep.front()->context.return_order(), &event_site, &end_site);
		if (!Pc) Qd = true;		// Dependent sites is definitely true for forward simulation. (!Pc flags EPC.)
	} else end_site = event_site+1;

	//////////
	/// Calculating the rates away from sequence i per 
	/// S.-C. Choi, B.D. Redelings, and J.L. Thorne. 2008. "Basing population genetic inferences 
	/// 	and models of molecular evolution upoon desired stationary distributions of DNA or 
	/// 	protein sequences". Phil. Trans. R. Soc. B 363:3931-3939. doi:10.1098/rstb.2008.0167.
	///
	///	This representation is different since we are (currently) not calculating the VLMM on
	/// protein sequences, but are rather leaving everything at the nucleotide level.
	//////////
	R_D = Rij(tree, event_site, end_site);

	//////////
	/// Make CDF of site_rate_away, but only for the sites that have changed!
	//////////
	evolvingSequence->forward_rate_away_from_sequence(branch, event_site, end_site);

	return R_D;
}

double 
TNode::Rij (
			TTree *tree, 
			unsigned int start_position, 
			unsigned int end_position
		   )
{
	vector<Site>::iterator start = seq_evo.begin()+start_position;
	vector<Site>::iterator end = seq_evo.begin()+end_position;
	double tau_ij;
	double site_rij;
	int position = start_position;
	double sum_rate_away = 0;
	bool initial_rates_calculation = false;
	double rate_ij;

	//////////
	/// Each site carries its own rate away, but on the first pass, those rates away need to be
	/// set. Thus, if the initial_rates_calculation flag is set, the sequence rate away is to be
	/// recalculated from 0. If the flag is not set, then the site rates away are pre-set (from the
	/// previous round), and we simply need to subtract that previous rate away and add the new
	/// rate away for the site.
	//////////
	if (start == seq_evo.begin() && end == seq_evo.end()) {
		initial_rates_calculation = true;
		sum_rate_away = 0;
	} else sum_rate_away = evolvingSequence->returnRij(); 

	position = start_position;
	for (vector<Site>::iterator i = seq_evo.begin()+start_position; i != seq_evo.begin()+end_position; ++i, ++position) {
		(*i).site_rate_away.clear();
		site_rij = 0;
		for (int j = 0; j < numStates; ++j) {
			if (j != (*i).returnState()) {
				if (Qd) {
					//////////
					/// Calculating the rates away from sequence i per 
					/// S.-C. Choi, B.D. Redelings, and J.L. Thorne. 2008. "Basing population genetic inferences 
					/// 	and models of molecular evolution upoon desired stationary distributions of DNA or 
					/// 	protein sequences". Phil. Trans. R. Soc. B 363:3931-3939. doi:10.1098/rstb.2008.0167.
					//////////
					tau_ij = TauIJ3(
								    tree, 
								    (*i).return_lookup_table_environment_index(), 
								    (*i).return_lookup_table_sequence_index(), 
								    (*i).returnState(), 
								    j
								   );
					//cerr << "Rij:: tau_ij value: " << tau_ij << endl << endl;
					///
					/// Equation 1.9 from Choi et al, utilizing the equation from Yang & Nielsen:
					/// 2N*Pr(Z_ij) ~ log(tau_ij)(1-1/tau_ij).
					if (tau_ij != 1) {
						rate_ij 
						= branch->rates->pi.at(j)
						  *
						  ( log(tau_ij)/(1 - 1.0/tau_ij) );
					} else rate_ij = tree->root->branch->rates->pi.at(j); 
				} else rate_ij = tree->root->branch->rates->pi.at(j);
			//////////
			/// i == j, therefore, is not a rate away.
			//////////
			} else rate_ij = 0;
			site_rij += rate_ij;
			(*i).site_rate_away.push_back(rate_ij);
			if (initial_rates_calculation) sum_rate_away += rate_ij;
		}
		if (!initial_rates_calculation) sum_rate_away +=  site_rij - (*i).returnrij();
		(*i).setrij(site_rij);
	}
	evolvingSequence->setRij(sum_rate_away);

	return sum_rate_away;
}

double 
TNode::TauIJ3 (
			  TTree *tree,
			  int env_index,
			  int i_seq_index,
			  short residue_i,
			  short residue_j
			 )
{
	double diffPji, diffP0ji;
	register double tau_ij; // The tau_ij parameter, equation 1.7 Choi et al.
	int j_seq_index = i_seq_index + tree->dep.front()->context.getOffset(env_index, residue_i, residue_j);

//	cerr << "env: " << env_index << "  ";
//	cerr << "j_seq_index = " << i_seq_index << " + " << tree->dep.front()->context.getOffset(env_index, residue_i, residue_j) << endl;

	//////////
	///	This representation is different since we are (currently) not calculating the VLMM on
	/// protein sequences, but are rather leaving everything at the nucleotide level.
	///
	/// Equation (1.1) Independent sites probability of sequence i.
	/// * Unnecessary to calculate: only need to find the difference.
	/// * Calculate P(i) given the 3rd order Markov model. This has no equation in Choi et al, but
	///   it plays an integral part in Equation (1.7).
	/// * First triplet stationary probability for the state of the sequence.
	/// *Calculating the probabilities for sequence j. This is different than sequence i, since
	///  j will be different in one position. This difference is what is difficult to do, however,
	///  since we will need to calculate the probabilities of all sequences that have that position.
	///  The first thing to do, in any case, is to keep the P_i_ calculation, and re-calculating
	///  only the necessary positions.
	///    * Inheriting 2 values: the probability of the triplet for i, and
	///    * The multiplicand for each site in j that is unaffected by the changed site.
	///
	/// Position 1:
	/// P(i)P(j_1,j_2,j_3)P(j_4|j_1,j_2,j_3)
	/// ------------------------------------
	/// P(i_1,i_2,i_3)P(i_4|i_1,i_2,i_3)
	///
	/// Rest of sequences cancel out, since none depend on the i_1->j_1 change. Upon close 
	/// inspection, one can see that this formula is simply tau_{ij} from Choi et al., solved
	/// for P(J). Thus, P(J) is P(I) times a ratio of P(J)'/P(I)', where the ' indicates the
	/// probabilities changed by the i->j change.
	//////////
	diffPji 
	= tree->dep.front()->context.lt_markov_ratio(
											     env_index,			// env.
											     i_seq_index,		// i
											     j_seq_index		// j
											    );
	//////////
	/// Final value needed to complete TauIJ for this particular position is
	/// the calculation of the stationary frequency of P0_j_ under the independent
	/// model.
	//////////
	diffP0ji = branch->rates->pi.at(residue_j) / branch->rates->pi.at(residue_i);
	
	//////////
	/// Equation 1.7 from Choi et al, tau_ij.
	//////////
	double diffP0ji_inv = 1.0 / diffP0ji;
	tau_ij = diffPji;
    tau_ij *= diffP0ji_inv;
	//////////

	return tau_ij;
}

void 
TNode::printForwardRateAway()
{
	double rate_away = 0;
	for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it) {
		rate_away += (*it).site_rate_away.at((*it).returnState());
	}
}

void 
TNode::updateSequence(
					  TNode *update_node
					 )
{
	vector<Site>::iterator it, jt;
	jt = update_node->seq_evo.begin();
	for (it = seq_evo.begin(); it != seq_evo.end(); ++it, ++jt) {
		(*it).setState((*jt).returnState());
	}
}

string 
TNode::printBipartition()
{
	string out = "";
	for (vector<bool>::iterator it = bipartition.begin(); it != bipartition.end(); ++it)
		out += ( (*it) ? "1" : "0" );
	return out;
}

////////////////////
/// Calculates rate away from site for positions that have specific states at their endpoints.
///  * Currently, only applicable for EPC runs,
///  * Function designed so that sequences can be guided to motifs more realistically in FWD sim.
////////////////////
void 
TNode::calculateEndpointRateAwayFromSite(
										 vector<Site>::iterator evolver_it,
										 vector<Site>::iterator target_it,
										 double *cumulative_site_sum_away,
										 double *forward_site_sum_away,
										 RateMatrix *rate_matrix
										)
{
	register double Pik_inverse, Pik, Qij, Pjk, Lk;
	double sum_away = 0;
	vector<double>::iterator seq_away_it;
	vector<double>::iterator forward_away_it;
	unsigned int site_no = 0;
	bool profile=false;			//XOUT

	//////////
	/// Qij * Pjk(T)/Pik(T)
	/// i = begin state, j = middle state, k = target state
	///      i   j                               k
	///      |---+-------------------------------|
	/// time 0   dt                              T
	//////////
	// Sum all Pij where i!=j
	// Pik will be the divisor for each other transition.
	// k is the target, to which all of these comparisons are conditioned.
	//////////
	if (!order_3_markov) (*evolver_it).site_rate_away.assign(numStates,0);
	(*evolver_it).Pjk0.assign(numStates, 0);
	(*evolver_it).forward_rate_away.assign(numStates,0);
	//////////
	/// Calculate sum_away
	//////////
	//
	// Current state (observed)
	int i = (*evolver_it).returnState();
	//
	//////////
	/// Main loop, taking into account conditional likelihoods for unset internal nodes.
	//////////
	// calculate the "bin size" for each particular value of j:
	//
	// sum_away_j_ = Qij * Lk
	//
	// Where 
	// (i) sum_away_j_ is the ?probability? that the next substitution step will occur on
	// the current site TO the value of j under consideration. To get the true probability, this
	// value is divided by the total sum away for the sequence, and 
	// (ii) Lk is the likelihoods of the target state k, calculated by pruning.
	//////////
	seq_away_it = (*evolver_it).site_rate_away.begin();
	forward_away_it = (*evolver_it).forward_rate_away.begin();

	unsigned int j, k;
	if ( (*target_it).isDefinedByUser() ) {
		k = (*target_it).returnState();
		(*evolver_it).Pik0 = rate_matrix->Pij.at((*evolver_it).returnCategory()).at(from(i)+to(k));
		if ((*evolver_it).Pik0 == 0) (*evolver_it).Pik0=1e-10;	/// Prevent division by zero. ///
		Pik_inverse = 1.0/(*evolver_it).Pik0;
		for (j = 0; j < numStates; j++) {
			if (i != j) {
				/// With the introduction of the rate_matrix passing into this function, there no longer needs to be a choice for which matrix to use.
				Qij = rate_matrix->Qij.at(from(i)+to(j));
				Qij *= rate_matrix->catRate.at((*evolver_it).returnCategory());
				(*evolver_it).Pjk0.at(j) = rate_matrix->Pij.at((*evolver_it).returnCategory()).at(from(j)+to(k));
				if (profile) 
					cerr << "state: " << stateCharacters.at((*evolver_it).returnState()) 
						 << "->" << stateCharacters.at(j) << "  s_a: " << sum_away << " += " 
						 << Qij << " * " 
						 << (*evolver_it).Pjk0.at(j) 
						 << " * " << Pik_inverse << endl;	//XOUT
				sum_away += Qij * ( (*evolver_it).Pjk0.at(j) * Pik_inverse );
				(*forward_site_sum_away) += Qij;
				(*forward_away_it) = Qij;
			} else (*forward_away_it) = 0;
			//
			// Add this to the site to the rate away.
			(*seq_away_it) = sum_away;
			++seq_away_it;
			++forward_away_it;
			//
			//////////
		}
	} else {
		//////////
		/// i is independent from j, but not from k. Pulling these statements out of the
		/// for loop over j saves 3/4 of the computation time.
		//////////
		double numerator = 0, denominator = 0, denominator_inverse;
		for (k = 0; k < numStates; k++) {
			Pik = rate_matrix->Pij.at((*evolver_it).returnCategory()).at(from(i)+to(k));
			Lk  = (*target_it).L_i.Li_xi_.at(k);
			denominator += Pik * Lk;
		}
		(*evolver_it).Pik0 = denominator;
		denominator_inverse = 1.0 / (*evolver_it).Pik0;
		//
		// Calculate the rates away given that the end state is the conditional likelihoods.
		for (j = 0; j < numStates; j++) {
			if (i != j) {
				if (order_3_markov) Qij = (*evolver_it).site_rate_away.at(j);
				else Qij = rate_matrix->Qij.at(from(i)+to(j));
				Qij *= rate_matrix->catRate.at((*evolver_it).returnCategory());
				/////
				// Calculate appropriate term:
				/////
				numerator = 0;
				for (k = 0; k < numStates; k++) {
					Lk  = (*target_it).L_i.Li_xi_.at(k);
					//Pjk = branch->rates->Pij.at((*evolver_it).returnCategory()).at(from(j)+to(k));						
					Pjk = rate_matrix->Pij.at((*evolver_it).returnCategory()).at(from(j)+to(k));						
					numerator += Qij * Pjk * Lk;
					(*evolver_it).Pjk0.at(j) += Pjk * Lk;
				}
				//
				// Inverse of denominator calculated before loops, done only 1 time.
				sum_away += Qij * (*evolver_it).Pjk0.at(j) * denominator_inverse;
				//
				/////
				(*forward_site_sum_away) += Qij;
				(*forward_away_it) = Qij;
			} else {
				//////////
				/// Not going away from site, so 0.
				//////////
				(*forward_away_it) = 0;
			}
			//
			// Add this to the site to the rate away.
			(*seq_away_it) = sum_away;
			++seq_away_it;
			++forward_away_it;
			//
			//////////
		}
	}

	(*cumulative_site_sum_away) += sum_away;
}

double 
TNode::calculateEndpointRateAwayFromSequence(
											 TTree *tree,
											 TNode *k_0,		// Target sequence. //
											 double T,
											 double at_dt,
											 int event_site
											)
{
	double site_sum_away = 0, cumulative_site_sum_away = 0, forward_site_sum_away = 0;
	setRateAway(TIME_RELATIVE_STEPS);	// Resets site_sum_away, allocates e_QijDt.
	unsigned int end_site;

	// THE -I OPTION
	//if (acc == independent_sites) rasmus_independent_proposals = true;		rasmus  = 0
	//else if (acc == QdPc) { Qd = Pc = true; nij = false; }					QdPc	= 1
	//else if (acc == QdP) { Qd = true; Pc = false; nij = false; }				QdP		= 2
	//else if (acc == QPc) { Qd = false; Pc = true; nij = false; }				QPc		= 3
	//else if (acc == QP) { Qd = false; Pc = false; nij = false; }				QP		= 4
	//else if (acc == Nij) { Qd = true; Pc = false; nij = true; }				QdN		= 5

	//////////
	/// Steps to follow:
	/// * Qij is set (for all delta_t, I think...)
	/// * if root_i diff than leaf_i
	///   - Set Pij for T-x*dt.
	//////////
	if (Pc || Qd) {	// Anything that needs Q to help build the structure needs to perform this step. //
		//////////
		/// Dependent sites. Recalculate Q matrices for positions at the beginning of the branch or
		/// affected by a change for subsequent calculation of the site-specific transition matrices.
		//////////
		set_site_window(tree->dep.front()->context.return_order(), &event_site, &end_site);
		site_specific_Qmat(tree, event_site, end_site);
	} else {
		//////////
		/// Independent sites: Copy the global model in for each site.
		//////////
		for (vector<Site>::iterator it = k_0->seq_evo.begin(); it != k_0->seq_evo.end(); ++it) 
			(*it).e_QijDt = k_0->branch->rates->Qij;
	}

	vector<Site>::iterator target_it = k_0->seq_evo.begin();
	RateMatrix temp_rates(*branch->rates);	// Copy constructor. This decl. keeps constant all data that is global, changing site-specific data in loop.
	int seq_pos = 0;
	if (nij) {
		vector<double> nij_row_sum (numStates, 0);
		int i = 0;
		vector<double>::iterator jt = nij_row_sum.begin();
		for (vector<double>::iterator it = branch->nij.begin(); it != branch->nij.end(); ++it, ++i) {
			if (i == numStates) { ++jt; i = 0; }
			(*jt) += (*it);
		}

		vector<double>::iterator pt = temp_rates.Pij.at(0).begin();	// Using nij, cannot use categories... yet. //
		jt = nij_row_sum.begin();
		i = 0;
		for (vector<double>::iterator it = branch->nij.begin(); it != branch->nij.end(); ++it, ++pt, ++i) {
			if (i == numStates) { ++jt; i = 0; }
			(*pt) = (*it) / (*jt);
		}

		pt = temp_rates.Pij.at(0).begin();
		i = 0;
		cerr << endl << "Transition Probabilities using nij." << endl;
		for (; pt != temp_rates.Pij.at(0).end(); ++pt, ++i) {
			if (i % numStates == 0) cerr << endl;
			cerr << (*pt) << " ";
		}
		//exit(0);
	} else if (!Pc) {	// Independent sites, need to scale the branch length for dependent sites.
		// setPij uses the Root and Cijk matrices, already set up for branch->rates
		branch->rates->setPij(seq_evo.front(), branch->S*(T-at_dt), nodeEnv->rateHetero);	// setPij dependes on the root Matrix
		// Transfer independent rates into the temp_rates holding the dependent Qij's.
		vector<vector<double> >::iterator ppt = temp_rates.Pij.begin(); 
		for (vector<vector<double> >::iterator PPt = branch->rates->Pij.begin(); PPt != branch->rates->Pij.end(); ++PPt, ++ppt) {
			vector<double>::iterator pt = (*ppt).begin();
			for(vector<double>::iterator Pt = (*PPt).begin(); Pt != (*PPt).end(); ++Pt, ++pt) {
				(*pt) = (*Pt);
			}
		}
	}

	for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it, ++target_it, ++seq_pos) {
		if (Qd && Pc) {
			/// Calculating and passing the site-specific transition probabilities.
			temp_rates.Qij = (*it).e_QijDt;							// Entire Qij matrix. For exponentiating.
			temp_rates.SetupMatrix(false);							// Set Cijk matrix, for exponentiating.
			temp_rates.setPij((*it), (T-at_dt), nodeEnv->rateHetero); // Exponentiation.
			calculateEndpointRateAwayFromSite(it, target_it, &cumulative_site_sum_away, &forward_site_sum_away, &temp_rates);
		} else if (!Qd && Pc) {
			temp_rates.Qij = (*it).e_QijDt;							// Entire Qij matrix. For exponentiating.
			temp_rates.SetupMatrix(false);							// Set Cijk matrix, for exponentiating.
			temp_rates.setPij((*it), (T-at_dt), nodeEnv->rateHetero); // Exponentiation.
			temp_rates.Qij = branch->rates->Qij;					// Reset dependent rates to independent rates.
			calculateEndpointRateAwayFromSite(it, target_it, &cumulative_site_sum_away, &forward_site_sum_away, &temp_rates);
		} else if (Qd && !Pc) {	// Qd set, !Pc set above to branch->rates so calculation of matrix does not happen every step.
			// Fast option
			calculateEndpointRateAwayFromSite(it, target_it, &cumulative_site_sum_away, &forward_site_sum_away, &temp_rates);
		} else {
			/// Passing global probabilities. This route is a great deal faster.
			calculateEndpointRateAwayFromSite(it, target_it, &cumulative_site_sum_away, &forward_site_sum_away, branch->rates);
		}
	}
	
	// Set values needed for calculating Forward and EPC log probabilities.
	evolvingSequence->Qidot_k__T__ = cumulative_site_sum_away;
	evolvingSequence->Qidot = forward_site_sum_away;
	cerr << " Qi.|k(t): " << evolvingSequence->Qidot_k__T__ << "         Qi: " << evolvingSequence->Qidot << endl; //XOUT

	return cumulative_site_sum_away;
}

void 
TNode::site_specific_Qmat(
						  TTree *tree,
						  int start, 		// Start position for calculating Q matrix.
						  unsigned int end	// End position
						 )
{
	int seq_pos = 0;
	int i,j;
	vector<unsigned int> power_1 (numStates, 0);
	for (i = 0; i < numStates; i++) power_1[i]=i*numStates; 
	vector<Site>::iterator it;
	vector<double>::iterator jt;
	short actual_state;
	string event = "XX";
	double tau_ij;

	seq_pos = start;
 	double row_total = 0;
	for (it = seq_evo.begin()+start; it != seq_evo.begin()+end; ++it, ++seq_pos) {
		actual_state = (*it).returnState();  // keep track of real state.
		event.at(0) = stateCharacters.at(actual_state);
		for (i = 0; i < numStates; ++i) {
			event.at(1) = stateCharacters.at(i);
			(*it).setState(i);
			tree->dep.front()->context.reset_sequence_indices(this, seq_pos, event);
			for (j = 0; j < numStates; ++j) {
				if (j != i) {
					tau_ij
					= TauIJ3(
							 tree, 
							 (*it).return_lookup_table_environment_index(), 
							 (*it).return_lookup_table_sequence_index(),
							 i,
							 j
							);
					if (tau_ij != 1) {
						(*it).e_QijDt.at(power_1[i]+j)
						= branch->rates->pi.at(j)
						  *
						  ( log(tau_ij)/(1 - 1.0/tau_ij) );
					} else (*it).e_QijDt.at(power_1[i]+j) = 0.25; 
					row_total += (*it).e_QijDt.at(power_1[i]+j);
				}
			}
			(*it).e_QijDt.at(power_1[i]+i) = -row_total;
			row_total = 0;
			event.at(0) = event.at(1);
		}

		(*it).setState(actual_state);	// Reset position to actual state.
		event.at(1) = stateCharacters.at(actual_state);
		tree->dep.front()->context.reset_sequence_indices(this, seq_pos, event);

		//////////
		/// Qij now holds the appropriate values for the site rate away for the position. Set it here.
		//////////
		j = 0;
		for (jt = (*it).site_rate_away.begin(); jt < (*it).site_rate_away.end(); ++jt, ++j) {
			(*jt) = (*it).e_QijDt.at(power_1[actual_state]+j);
			if ( (*jt) < 0 ) (*jt) = 0;
		}
		// Make CDF
		for (jt = (*it).site_rate_away.begin()+1; jt < (*it).site_rate_away.end(); ++jt, ++j) (*jt) += (*(jt-1));

//		cerr << "Position " << seq_pos << " Qmat: ";
//		for (jt = (*it).e_QijDt.begin(); jt != (*it).e_QijDt.end(); ++jt)
//			cerr << (*jt) << " ";
//		cerr << endl;
	}
}

int TNode::locateSite(
					  double *rate_away_locator,
					  int simulation_type
					)
{
	int site_no = 0;

	for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it, ++site_no) {
		*rate_away_locator -= (*it).site_rate_away.back();
		if (*rate_away_locator < 0) {
			*rate_away_locator += (*it).site_rate_away.back();
			return site_no;
		}
	}

	//////////
	/// The return statement is in the loop above. If we get here, we did not find a valid site
	/// in the sequence, i.e., site_no < 0 or site_no > sequence_length.
	//////////
	cerr << "Could not find the rate_away locator. This should be impossible. TNode::locateSite." << endl;
	abort();
}

bool TNode::EndpointCheck( void )
{
	for (vector<double>::iterator it = seq_evo.at(0).L_i.Li_xi_.begin(); it != seq_evo.at(0).L_i.Li_xi_.end(); ++it) {
		if ((*it) == 1) return true;
	}
	return false;
}

////////////////////
////// Returns true if the node is a descendant of the listed ancestral node.
////////////////////
bool 
TNode::isDescendant(
						 TTree *tree,
						 TNode *ancestor
						)
{
	//////////
	///
	//////////
	// If this is true, then the branch is ancestral (even if it is the root).
	if (branch0 == ancestor) return true;
	//
	// If the branch is the same as the root, then we are done looking, and node is not a descendant
	else if (branch0 == tree->root) return false;
	//
	// Otherwise, look at the next node up in the tree.
	else return branch0->isDescendant(tree, ancestor);

}

size_t 
TNode::findForwardNumberOfConstrained(
											 size_t fromSite, 
											 size_t numSites
											) 
{
	int num_unconstrained = 0;
	varSite *curr_varSite[2], *new_varSite[2];
	int num_chomped[2];
	
	num_chomped[0] = num_chomped[1] = 1;
	curr_varSite[0] 
	= new_varSite[0] 
	= seq_evo.at(fromSite).motif.active_properties.indel->del->my_sequence_template_varSite;
	curr_varSite[1] 
	= new_varSite[1] 
	= seq_evo.at(fromSite).motif.active_properties.indel->del->my_motif_varSite;

																						 // short circuit... use later //
	for (vector<Site>::iterator site_it = seq_evo.begin()+fromSite; site_it != seq_evo.end() && num_unconstrained < numSites; ++site_it, num_unconstrained++) {
		new_varSite[0] = (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite;
		new_varSite[1] = (*site_it).motif.active_properties.indel->del->my_motif_varSite;
		if (curr_varSite[0] != new_varSite[0]) {
			num_chomped[0] = 1;
			curr_varSite[0] = new_varSite[0];
		} 
		if (curr_varSite[1] != new_varSite[1]) {
			num_chomped[1] = 1;
			curr_varSite[1] = new_varSite[1];
		}
																	// This is for one_site_varSite. //
		if (curr_varSite[0]->member_set.size() - num_chomped[0] < curr_varSite[0]->min || curr_varSite[0]->min == curr_varSite[0]->max) {
			break;
		}
		if (curr_varSite[1]->member_set.size() - num_chomped[1] < curr_varSite[0]->min || curr_varSite[1]->min == curr_varSite[1]->max) {
			break;
		}

		num_chomped[0]++;
		num_chomped[1]++;
	}

	return num_unconstrained;
}

size_t 
TNode::findBackwardNumberOfConstrained(
											  size_t fromSite, 
											  size_t numSites
											 ) 
{
	int num_unconstrained = 0;
	varSite *curr_varSite[2], *new_varSite[2];
	int num_chomped[2];

	curr_varSite[0] 
	= new_varSite[0] 
	= seq_evo.at((seq_evo.size()-1)-fromSite).motif.active_properties.indel->del->my_sequence_template_varSite;
	curr_varSite[1] 
	= new_varSite[1] 
	= seq_evo.at((seq_evo.size()-1)-fromSite).motif.active_properties.indel->del->my_motif_varSite;
	num_chomped[0] = num_chomped[1] = 1;
																						 // short circuit... use later //
	for (vector<Site>::reverse_iterator site_rit = seq_evo.rbegin()+fromSite; site_rit != seq_evo.rend() && num_unconstrained < numSites; ++site_rit, num_unconstrained++) {
		new_varSite[0] = (*site_rit).motif.active_properties.indel->del->my_sequence_template_varSite;
		new_varSite[1] = (*site_rit).motif.active_properties.indel->del->my_motif_varSite;
		if (curr_varSite[0] != new_varSite[0]) {
			num_chomped[0] = 1;
			curr_varSite[0] = new_varSite[0];
		} 
		if (curr_varSite[1] != new_varSite[1]) {
			num_chomped[1] = 1;
			curr_varSite[1] = new_varSite[1];
		}
																	// This is for one_site_varSite. //
		if (curr_varSite[0]->member_set.size() - num_chomped[0] < curr_varSite[0]->min || curr_varSite[0]->min == curr_varSite[0]->max) {
			break;
		}
		if (curr_varSite[1]->member_set.size() - num_chomped[1] < curr_varSite[0]->min || curr_varSite[1]->min == curr_varSite[1]->max) {
			break;
		}
		num_chomped[0]++;
		num_chomped[1]++;
	}

	return num_unconstrained;
}

void 
TNode::addGeneral_varSites()
{
	variable_region_list.clear();
	one_site_varSite = new varSite(1,1,this);
	unconstrained_varSite = new varSite(0,numeric_limits<int>::max(),this);
}

void 
TNode::report_varSites()
{
	size_t total_varSites = 0;

	cout << "Address\tminl\tmaxl\tcurrent_num_members " << endl;
	for (list<varSite*>::iterator it = variable_region_list.begin(); it != variable_region_list.end(); ++it) {
		cout << (*it) << "\t";
		cout << (*it)->min << "\t" << (*it)->max << "\t" << (*it)->member_set.size() << endl;
		total_varSites += (*it)->member_set.size();
	}

	cout << "The total number of varSites is: " << total_varSites << endl;
}

void 
TNode::report_on_sites()
{
	for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); ++site_it) {
		cerr << "L_ins_ (st=" << (*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite
		     << " m=" << (*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite << ")   ";
		cerr << "R_ins_ (st=" << (*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite
		     << " m=" << (*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite << ") " << endl;
	}
}

varSite *
TNode::Generic_varSite()
{
	return one_site_varSite;
}

void 
TNode::inheritCategories(
					   		  TNode *anc
					  		 )
{
	vector<Site>::iterator anc_it = anc->seq_evo.begin(), des_it = seq_evo.begin();
	for ( ; anc_it != anc->seq_evo.end(); ++anc_it, ++des_it)
		(*des_it).setCategory( (*anc_it).returnCategory() );
}

void 
TNode::Print_Substitution_Properties()
{
	vector<Site>::iterator jt = branch1->seq_evo.begin();
	vector<double> cat_avg (4,0);
	vector<int> num_cat (4,0);
	for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it, ++jt) {
		printf(
			   "%c %c %3d %8d\n",
			   stateCharacters.at((*it).returnState()), 
			   stateCharacters.at((*jt).returnState()),
			   (*it).returnCategory(),
			   ( 
			     (
			      (*it).numSubstAway
			     ) 
			     ? (*it).numSubstAway 
			     : (*jt).numSubstAway 
			   )
		);
		num_cat.at((*it).returnCategory())++;
		cat_avg.at((*it).returnCategory()) += 
			( 
			 (
			  (*it).numSubstAway
			 ) 
			 ? (*it).numSubstAway 
			 : (*jt).numSubstAway 
			);
	}

	for (int i = 0; i < 4; i++) {
		cout << i << " " << cat_avg.at(i) / num_cat.at(i) << endl;
	}
}

//////////
/// Copies varSites from ancestor to descendant.
//////////
void 
TNode::Inherit_varSites() 
{
	size_t i = 0;

	anc->one_site_varSite->descendant_equiv = one_site_varSite;
	anc->unconstrained_varSite->descendant_equiv = unconstrained_varSite;
	for (list<varSite*>::iterator it = anc->variable_region_list.begin(); it != anc->variable_region_list.end(); ++it, i++) {
		if (i > 1) {
			variable_region_list.push_back(new varSite(*it,this));
		}
	}

	list<varSite*>::iterator anc_it = anc->variable_region_list.begin();
	for (list<varSite*>::iterator it = variable_region_list.begin(); it != variable_region_list.end(); ++it, ++anc_it) {
		if (anc_it == anc->variable_region_list.end()) {
			cerr << "something's wrong" << endl;
			exit(EXIT_FAILURE);
		}
	}
}

//////////
/// Copy the motif sites. Does this happen by default when copying the sequence from the ancestor?
//////////
void 
TNode::InheritMotifSites()
{	
	//////////
	/// Create new positions for each. Setting active sites will copy new things into each pos.
	//////////
	// SITE
	vector<Site>::iterator des_site_it = seq_evo.begin();
	for (vector<Site>::iterator anc_site_it = anc->seq_evo.begin(); 
								anc_site_it != anc->seq_evo.end(); 
								++anc_site_it, ++des_site_it) 
		(*des_site_it).motif.copy(&((*anc_site_it).motif));
	//
	// postProcess
	//
	this->Site_postProcess();
	//
	// SITE
	//
	seq_evo.front().motif.site_props.front()->indel.del->my_motif_varSite 
	= seq_evo.front().motif.site_props.front()->indel.L_ins_->my_motif_varSite_left
	= seq_evo.front().motif.site_props.front()->indel.L_ins_->my_motif_varSite_right
	= NULL;
	seq_evo.front().motif.site_props.back()->indel.del->my_sequence_template_varSite 
	= seq_evo.front().motif.site_props.back()->indel.L_ins_->my_sequence_template_varSite_left
	= seq_evo.front().motif.site_props.back()->indel.L_ins_->my_sequence_template_varSite_right
	= NULL;
	seq_evo.front().motif.site_props.front()->indel.del->my_sequence_template_varSite 
	= seq_evo.front().motif.site_props.front()->indel.L_ins_->my_sequence_template_varSite_left
	= seq_evo.front().motif.site_props.front()->indel.L_ins_->my_sequence_template_varSite_right
	= one_site_varSite;
	seq_evo.front().motif.site_props.back()->indel.del->my_motif_varSite 
	= seq_evo.front().motif.site_props.back()->indel.L_ins_->my_motif_varSite_left
	= seq_evo.front().motif.site_props.back()->indel.L_ins_->my_motif_varSite_right
	= one_site_varSite;
	// New ActiveProps:
	evolvingSequence->setActiveProps();
	//cerr << "Ancestor: " << endl;
	//anc->Print_Active_Properties();
	//cerr << "Descendant: " << endl;
	//Print_Active_Properties();
	//exit(0);
}

void 
TNode::setInvariableArrayPos()
{
	// Nearly all set, but still have to represent the invariable array constraints, applies as motif.
	if (nodeEnv->invariableSites) { 
		vector<Site>::iterator sit = seq_evo.begin();
		for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); ++site_it) {
			bool last_site = ( (site_it == seq_evo.end()-1) ? true : false );
			if ( !(*site_it).motif.active_properties.fromMotif ) {
				// Need to make current position invariable, and change motif pointers of the
				// previous and next position (if they exist).
				switch ( (*site_it).returnInvariableState() ) {
				case NO_CONSTRAINT:
					//////////
					/// Do nothing
					//////////
					break;
				case INVARIABLE:
					//////////
					/// Set to (1,1), make subst equal to only current aa.
					//////////
					(*site_it).motif.active_properties.indel->del->my_motif_varSite
					= (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right
					= (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left
					= one_site_varSite;
					
					(*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite 
					= (*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite 
					= unconstrained_varSite;

					(*site_it).motif.active_properties.subst->setInvariable((*site_it).returnState());
					break;
				case NO_INDEL:
					//////////
					/// set to (1,1), check next positions of invar to see if R_ins_->on_site needs set to (1,1)
					//////////
					(*site_it).motif.active_properties.indel->del->my_motif_varSite
					= (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right
					= (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left
					= one_site_varSite;

					if ( !last_site ) {
						//////////
						/// If either 2 or 3, this subsequece will be size 1 only throughout entire run.
						//////////
						if ( (*(site_it+1)).returnInvariableState() == NO_INDEL || 
							 (*(site_it+1)).returnInvariableState() == INVAR_AND_NOINDEL ) {
							(*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite 
							= one_site_varSite;
						}
					} // Else we don't set anything.
					break;
				case INVAR_AND_NOINDEL:
					//////////
					/// Do both INVAR & INDEL.
					//////////
					(*site_it).motif.active_properties.indel->del->my_motif_varSite
					= (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right
					= (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left
					= one_site_varSite;
	
					if ( !last_site ) {
						//////////
						/// As above.
						//////////
						if ( (*(site_it+1)).returnInvariableState() == NO_INDEL || 
							 (*(site_it+1)).returnInvariableState() == INVAR_AND_NOINDEL ) {
							(*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite 
							= one_site_varSite;
						}
					} // Else we don't set anything.	
					(*site_it).motif.active_properties.subst->setInvariable((*site_it).returnState());
					break;
				default:
					cerr << "What kind of invariable array position is " << (*site_it).returnInvariableState() << endl;
					exit(EXIT_FAILURE);
					break;
				}
			}
		}
	}
}

bool 
TNode::isVarSite(
					  varSite *chk_varSite
					 )
{
	if (chk_varSite == NULL) return true;

	for (list<varSite*>::iterator it = variable_region_list.begin(); it != variable_region_list.end(); ++it)
		if (chk_varSite == (*it) ) return true;

	return false;
}

void 
TNode::report(
				   TTree *tree
				  ) 
{
	cerr << "Node " << mytipNo << ":" << endl << "  ";
	for (vector<bool>::iterator it = bipartition.begin(); it != bipartition.end(); ++it)
		cerr << (*it);
	cerr << "  DistanceFromRoot: " << DistanceFromRoot << " trsDFR: " << trDistanceFromRoot << endl;
}

void 
TNode::FULL_REPORT()
{
	cout << "anc node: " << anc << endl;
	cout << "des node: " << this << endl;

	list<varSite*>::iterator anc_vS = anc->variable_region_list.begin();
	list<varSite*>::iterator des_vS = variable_region_list.begin();
	for (; des_vS != variable_region_list.end(); des_vS++, anc_vS++) {
		cout << setw(10) << (*anc_vS) << setw(4) << (*anc_vS)->min << setw(8) << (*anc_vS)->max << setw(10) << (*anc_vS)->descendant_equiv << setw(5) << (*anc_vS)->member_set.size();
		cout << setw(14) << "des " << (*des_vS) << setw(4) << (*des_vS)->min << setw(8) << (*des_vS)->max << setw(5) << (*des_vS)->member_set.size() << endl;
	}

	vector<Site>::iterator anc_it = anc->seq_evo.begin();
	vector<Site>::iterator des_it = seq_evo.begin();	
	for ( ; anc_it != anc->seq_evo.end(); ++anc_it, ++des_it) {
		cout << "subst: ";
		cout << setw(10) << &(*anc_it).motif.active_properties 
			 << setw(21) 
			 << (*anc_it).motif.active_properties.subst->report_bitset();
		cout << setw(10) << &(*des_it).motif.active_properties 
			 << setw(21) 
			 << (*des_it).motif.active_properties.subst->report_bitset();
		cout << endl;
	}
}

void 
Branch::report()
{
	cout << "Branch Lengths:" << endl;
	cout << "  length0: " << length0 << " trs: " << branch0_time_relative_length << endl;
	cout << "  length1: " << length1 << " trs: " << branch1_time_relative_length << " max_path: " << branch1_max_path << endl;
	cout << "  length2: " << length2 << " trs: " << branch2_time_relative_length << " max_path: " << branch2_max_path << endl;
	cout << "  perturbation: " << perturbation << endl;
	rates->fullReport();
}

void 
Branch::update_nij(
				   short from_state, 
				   short to_state,
				   short target_state
				  )
{
	nij.at(from_state*numStates+target_state)--;
	nij.at(to_state*numStates+target_state)++;
}

void
Branch::print_nij()
{
	int num_print = 0;
	cerr << "nij: ";
	for (vector<double>::iterator it = nij.begin(); it != nij.end(); ++it, ++num_print) {
		if (num_print % numStates == 0) cerr << endl;
		cerr << (*it) << " ";
	}
}

void 
TNode::clearMotifSites() 
{
	string seq = "";
	for (vector<Site>::iterator site_it = anc->seq_evo.begin(); site_it != anc->seq_evo.end(); ++site_it)
		seq += (*site_it).returnState();

	//////////
	/// Rebuild the evolutionaryAttributes array as a non-constrained array with the stored seq
	/// as the initial state.
	//////////
	constructUnconstrainedSequence(seq);
}

string 
TNode::output_sequence(
							  string::iterator start, 
							  string::iterator end
							 )
{
	string seqout = "";
	for (string::iterator it = start; it != end; ++it) {
		seqout += stateCharacters[*it];
	}
	return seqout;
}

string 
TNode::printSequence(bool print_nucleotide_sequence)
{
	string return_sequence;
	return_sequence.clear();
	for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it) {
		if ((*it).isSet())
			if ( print_nucleotide_sequence )
				return_sequence += stateCharacters.at((*it).returnState());
			else return_sequence += (*it).returnState();
		else
			return_sequence += "X";
	}

	return return_sequence;
}

void 
TNode::printGammaCategories()
{
	short cat = 0;
	for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it) {
		cat = (*it).returnCategory();
		if (cat < 0 || cat > branch->rates->num_categories) cerr << "X";
		else cerr << cat;
	}
	cerr << endl;
}

//////////
/// REMOVE_OBJECTS
//////////
void 
TNode::Remove_Objects(
						   bool root_node
						  )
{
	size_t i = 0, j = 0;
	delete seq_evo.front().motif.active_properties.indel->L_ins_;
	for (vector<Site>::iterator site_it = seq_evo.begin(); 
							    site_it != seq_evo.end(); 
							    ++site_it,i++) 
	{
		j = 0;
		for (list<siteProperties*>::iterator it2 = (*site_it).motif.site_props.begin(); 
											 it2 != (*site_it).motif.site_props.end(); 
											 ++it2, ++j) 
		{
			///////////
			/// General sites become activeProps. If I run across a general site, I don't want to
			/// delete the activeProps (b/c much non-hilarity ensues...)
			///////////
			if ( it2 == (*site_it).motif.site_props.begin() ) {
				delete (*it2)->indel.L_ins_;
			} else {
				if ( !root_node )
					delete (*it2)->indel.L_ins_;
			}
			delete (*it2)->indel.R_ins_;
			delete (*it2)->indel.del;
			delete (*it2);
		}
		(*site_it).motif.site_props.clear();
		delete (*site_it).motif.active_properties.subst;
		delete (*site_it).motif.active_properties.indel->R_ins_;
		delete (*site_it).motif.active_properties.indel->del;
		delete (*site_it).motif.active_properties.indel;
	}
	Remove_varSites();
}

void 
TNode::Remove_varSites()
{
	// variable sites list.
	for (list<varSite*>::iterator it = variable_region_list.begin(); it != variable_region_list.end(); ++it) {
		delete (*it);
	}
	variable_region_list.clear();	

}

void 
TNode::Remove_iz_Objects ( )
{
	for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); ++site_it) 
		(*site_it).motif.site_props.clear();

	Remove_varSites();
}

void 
TNode::Site_postProcess()
{
	// Set all pointers for the motifSite array and deeper structures.
	for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); ++it)
		(*it).motif.setNode(this);
}

void 
TNode::Print_Active_Properties()
{
	int i = 0;
	cerr << endl << "BIPARTITION: " << endl; report(); cerr << endl;
	for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); ++site_it, i++) {
		cerr << i << " " << (short)(*site_it).returnState() 
			 << " fromMotif:  ";
		if (&(*site_it).motif.active_properties.fromMotif != NULL) {
			cerr << (&(*site_it).motif.active_properties.fromMotif);
		} else { cerr << "NULL"; }
		cerr << " fromTemplate:  ";
 	 	if (&(*site_it).motif.active_properties.fromTemplate != NULL) {
			cerr << &(*site_it).motif.active_properties.fromTemplate;
		} else cerr << "NULL";
		cerr << " site:  " << (&(*site_it)) 
			 << " motif: " << (&(*site_it).motif) 
			 << " subst: " << (&(*site_it).motif.active_properties.subst)
			 << " indel: " << (&(*site_it).motif.active_properties.indel) << endl << "  ";

		cerr << "subst_bitstring: " << (*site_it).motif.active_properties.subst->report_bitset()
			 << endl << "  ";
		if ( (*site_it).motif.active_properties.indel->L_ins_ != NULL ) {
			cerr << "L_ins_ " 
				 << (*site_it).motif.active_properties.indel->L_ins_
				 << " ";
			if ( (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left != NULL)
				cerr << "stL("
					 << (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left->min
					 << ","
					 << (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left->max
					 << ") ";
			else {
				cerr << "-----------------------------------------NO mstv_l " << endl;
				abort();
			}
			if ( (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right != NULL)
				cerr << "stR("
					 << (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right->min
					 << ","
					 << (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right->max
					 << ") ";
			else {
				cerr << "-----------------------------------------NO mstv_r " << endl;
				abort();
			}
			if ( (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left != NULL)
				cerr << "mL("
					 << (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left->min
					 << ","
					 << (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left->max
					 << ") ";
			else {
				cerr << "-----------------------------------------NO mmv_l " << endl;
				abort();
			}
			if ( (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right != NULL)
				cerr << "mR("
					 << (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right->min
					 << ","
					 << (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right->max
					 << ") ";
			else {
				cerr << "-----------------------------------------NO mmv_r " << endl;
				abort();
			}
			if ( (*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite != NULL && 
				 isVarSite((*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite) )
				cerr << "osst("
					 << (*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite->min 
					 << "," 
					 << (*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite->max 
					 << ") ";
			if ( (*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite != NULL && 
				 isVarSite((*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite) ) 
				cerr << "osm(" 
					 << (*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite->min 
					 << "," 
					 << (*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite->max 
					 << ") ";
		} else cerr << "NONE   ";

		cerr << endl << "  ";

		if ( (*site_it).motif.active_properties.indel->del != NULL) {
			cerr << "del    " 
				 << (*site_it).motif.active_properties.indel->del
				 << " ";
			if ( (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite != NULL) {
			 	cerr << "st("
			 		 << (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite->min 
			 		 << ","
			 		 << (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite->max
			 		 << ")["
					 << (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite->member_set.size()
					 << "] ";
				//////////
				/// ERROR CHECKING: if varSite has more than the maximum occupancy or less than
				/// the least, die here.
				//////////
				if (
					(
					 (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite->min
					 > (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite->member_set.size()
					 ||
					 (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite->max
					 < (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite->member_set.size()
					)
					&&
					( 
					 (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite->min != 1
					 &&
					 (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite->max != 1
					)
				   )
				{
					cerr << "Illegal boundary constraints!!!" << endl;
					exit(EXIT_FAILURE);
				}

			} else {
				cerr << "-----------------------------------------NO mstv in del " << endl;
				abort();
			}
			if ( (*site_it).motif.active_properties.indel->del->my_motif_varSite != NULL)
			 	cerr << "m("
			 		 << (*site_it).motif.active_properties.indel->del->my_motif_varSite->min 
			 		 << ","
			 		 << (*site_it).motif.active_properties.indel->del->my_motif_varSite->max
			 		 << ") ";
			else {
				cerr << "-----------------------------------------NO mmv in del " << endl;
				abort();
			}
		} else 
			cerr << "NONE"
				 << "    ";

		cerr << endl << "  ";


		if ( (*site_it).motif.active_properties.indel->R_ins_ != NULL ) {
			cerr << "R_ins_ " 
				 << (*site_it).motif.active_properties.indel->R_ins_
				 << " ";
			if ( (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left != NULL)
				cerr << "stL("
				     << (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left->min
				     << ","
				     << (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left->max
					 << ") ";
			if ( (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right != NULL)
				cerr << "stR("
					 << (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right->min
					 << ","
					 << (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right->max
					 << ") ";
			if ( (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left != NULL)
				cerr << "mL("
					 << (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left->min
					 << ","
					 << (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left->max
					 << ") ";
			if ( (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right != NULL)
				cerr << "mR("
					 << (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right->min
					 << ","
					 << (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right->max
					 << ") ";
			if ( (*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite != NULL && 
				 isVarSite((*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite) )
				cerr << "osst("
					 << (*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite->min 
					 << "," 
					 << (*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite->max 
					 << ") ";
			if ( (*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite != NULL && 
				 isVarSite((*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite) )
				cerr << "osm(" 
					 << (*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite->min 
					 << "," 
					 << (*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite->max 
					 << ") ";
		} else cerr << "NONE   ";

		cerr << endl;

		cerr << "    Inactives: " << endl;
		for (list<siteProperties*>::iterator jt = (*site_it).motif.site_props.begin(); jt != (*site_it).motif.site_props.end(); ++jt) {
			cerr << "      " << (*jt) 
				 << " indel: " << (&(*jt)->indel) 
				 << " del: " << (*jt)->indel.del << " "	
				 << ( ((*jt)->indel.del->my_sequence_template_varSite != NULL) ? "s" : "-" )
				 << ( ((*jt)->indel.del->my_motif_varSite != NULL) ? "m" : "-" )
				 << " L_ins_: " << (*jt)->indel.L_ins_
				 << " R_ins_: " << (*jt)->indel.R_ins_ << endl;
		}
		cerr << endl;
	}
}

//////////
/// This is a patch function. At some point, the Sites in the evolving sequence maintain the
/// original ancestral varSite. Although this does not crash the program, it modifies the 
/// template constraints, causing unwanted output. This function resets all of the 
//////////
void 
TNode::setSitePointers(
							string from_function
						   )
{
	varSite *prev_mstv;
	cerr << "TNode::setSitePointers from " << from_function << endl;
	for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); ++site_it) {
		//////////
		/// For each insertion & deletion pointer, check to make sure it is in this TNode's
		/// varSites.
		//////////
		//
		// -----------------
		// 	Deletion:: varSite *my_sequence_template_varSite;
		prev_mstv = (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite;
		(*site_it).motif.active_properties.indel->del->my_sequence_template_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->del->my_sequence_template_varSite,
			"(*site_it).motif.active_properties.indel->del->my_sequence_template_varSite"
		);
		if(prev_mstv != (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite" << endl;
		//
		//  Deletion:: varSite *my_motif_varSite;
		prev_mstv = (*site_it).motif.active_properties.indel->del->my_motif_varSite;
		(*site_it).motif.active_properties.indel->del->my_motif_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->del->my_motif_varSite,
			"(*site_it).motif.active_properties.indel->del->my_motif_varSite"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->del->my_motif_varSite)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->del->my_motif_varSite" << endl;
		//
		// -----------------
		// Insertion::R_ins_ varSite *my_sequence_template_varSite_left;
		prev_mstv = (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left;
		(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left,
			"(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left" << endl;
		//
		// Insertion::R_ins_ varSite*my_sequence_template_varSite_right;
		prev_mstv = (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right;
		(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right,
			"(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right" << endl;
		//
		// Insertion::R_ins_ varSite *my_motif_varSite_left; 
		prev_mstv = (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left;
		(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left,
			"(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left" << endl;
		//
		// Insertion::R_ins_ varSite *my_motif_varSite_right;
		prev_mstv = (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right;
		(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right,
			"(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right" << endl;
		//
		// Insertion::R_ins_ varSite *on_site_motif_varSite; 
		prev_mstv = (*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite;
		(*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite
		= checkSitePointer(
			(*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite,
			"(*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite" << endl;
		//
		// Insertion::R_ins_ varSite *on_site_sequence_template_varSite;	
		prev_mstv = (*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite;
		(*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite,
			"(*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite" << endl;
		//
		// -----------------
		// Insertion::L_ins_ varSite *my_sequence_template_varSite_left;
		prev_mstv = (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left;
		(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left,
			"(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left" << endl;
		//
		// Insertion::L_ins_ varSite*my_sequence_template_varSite_right;
		prev_mstv = (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right;
		(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right,
			"(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right" << endl;
		//
		// Insertion::L_ins_ varSite *my_motif_varSite_left; 
		prev_mstv = (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left;
		(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left,
			"(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left" << endl;
		//
		// Insertion::L_ins_ varSite *my_motif_varSite_right;
		prev_mstv = (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right;
		(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right,
			"(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right" << endl;
		//
		// Insertion::L_ins_ varSite *on_site_motif_varSite; 
		prev_mstv = (*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite;
		(*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite,
			"(*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite)
			cerr << "Adjusted (*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite" << endl;
		//
		// Insertion::L_ins_ varSite *on_site_sequence_template_varSite;
		prev_mstv = (*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite;
		(*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite,
			"(*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite"
		);
		if (prev_mstv != (*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite)
		 cerr << "Adjusted (*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite" << endl;
		//
		//////////
	}
}

varSite*
TNode::checkSitePointer(
								 varSite *site_ptr, 
								 string message
								)
{
	if (! isVarSite(site_ptr) )
		if (anc->isVarSite(site_ptr)) 
			//////////
			/// We have discovered that the current varSite is contained in the ancestral node.
			/// Thus, we reset this varSite to the equivalent varSite in the descendant node.
			/////////
			return site_ptr->descendant_equiv;
//			site_ptr = site_ptr->descendant_equiv;
		else {
			cerr << message << endl;
			abort();
		}

	return site_ptr;
}

TTree::TTree()
{
	InitTree();
}

void 
TTree::InitTree()
{	
	treeEnv.clear();
	root=NULL;
	nodeList.clear(); numNodes=0;
	names.clear();
	numTips=0;
	totalLength=0.0;
	rooted=0;
	lengths=-1;
	global_alignment = new globalArray();
	epc_root_is_set = false;
}

void 
TTree::constructRootMotifSites(
									list<inMotif*>& motif_specs, 
									string& root_sequence, 
									string& invariable_in
								   ) 
{
	size_t found;
	vector<siteProperties*> site_properties;
	string::iterator sit, sit2;
	vector<short> dot_additions2template (root_sequence.size(),0);

	//////////
	/// "Immortal" link, in a way. If this is a -, causes a lot of pain and agony.
	//////////
	if (root_sequence.at(0) == '-') root_sequence.at(0) = 'A';
	//////////
	/// If any '.'s are placed on HMM-like match states, they will get here. This'll blast 'em.
	//////////
	for (list<inClade*>::iterator it = treeEnv.begin(); it != treeEnv.end(); ++it) {
		for (list<inMotif*>::iterator it2 = (*it)->my_motifs.begin(); it2 != (*it)->my_motifs.end(); ++it2) {
			found = (*it2)->sitemap.find_first_of(".");
			while (found != string::npos) {
				for (list<inClade*>::iterator it4 = treeEnv.begin(); it4 != treeEnv.end(); it4++) {
					for (list<inMotif*>::iterator it3 = (*it4)->my_motifs.begin(); it3 != (*it4)->my_motifs.end(); ++it3) {
						vector<short>::iterator rm = dot_additions2template.begin()+found;
						if ( (*it3)->sitemap.at(found) == '*' || (*it3)->sitemap.at(found) == '.') {
							(*it3)->sitemap.erase(found,1);
							if ( (*it3)->isTemplate() ) {
								dot_additions2template.erase(rm);
							}
						} else {
							if ( (*it3)->isTemplate() ) {
								(*it3)->sitemap.erase(found,1);
								dot_additions2template.erase(rm);
								dot_additions2template.at(found) = 1;
							} else {
								cerr << "Cannot have a '.' in the motif specifications overlapping with another motif (template overlap is acceptable)" << endl;
								exit(EXIT_FAILURE);
							}
						}
					}
				}
				root_sequence.erase(found,1);
				invariable_in.erase(found,1);
				found = (*it2)->sitemap.find_first_of(".");
			}
		}
	}

	//////////
	/// Build the sites with which we will play 
	//////////
	root->evolvingSequence = new Sequence(root, root_sequence.size());
	root->evolvingSequence->init(root, root_sequence, treeEnv.front());

	//////////
	/// Step 1: Relate all motifs to the root sequence.
	//////////
	size_t dot_additions_sum = 0;
	for (vector<short>::iterator it = dot_additions2template.begin(); it != dot_additions2template.end(); ++it) 
		dot_additions_sum += (*it);

	size_t motif_start_position, post_motif_size, post_motif, shouldbezero;
	for (list<inClade*>::iterator it = treeEnv.begin(); it != treeEnv.end(); ++it) {
		for (list<inMotif*>::iterator it2 = (*it)->my_motifs.begin(); it2 != (*it)->my_motifs.end(); ++it2) {
			int relationship = testRelation( (*it)->bipartition, root->bipartition);
			site_properties = (*it2)->enumerateRegEx(root);
			motif_start_position = (*it2)->sitemap.find_first_not_of("*");
			if (motif_start_position == string::npos) {
				cerr << "TTree::constructRootMotifSites: Cannot find the beginning of motif '" << (*it2)->marker << "'. This error may be generated when motif markers do not match between the specification file and the root sequence input file." << endl;
				exit(EXIT_FAILURE);
			}
			post_motif = (*it2)->sitemap.find_first_of("*", motif_start_position);
			if (post_motif == string::npos) post_motif = root_sequence.size();			
			shouldbezero = root_sequence.size() - (site_properties.size() + motif_start_position + (root_sequence.size()-post_motif));
			if ( (*it2)->isTemplate() ) shouldbezero += dot_additions_sum;
			if (shouldbezero != 0) {
				cerr << "Motif clade \'" << (*it)->clade_name << "\', marker \'" << (*it2)->marker << "': Regular expression size (" << site_properties.size() << ") does not agree with motif size (" << post_motif-motif_start_position << ")." << endl;
				(*it2)->report();
				exit(EXIT_FAILURE);
			}

			if ( (*it2)->isTemplate() ) {
				vector<siteProperties*>::iterator sPit = site_properties.begin();
				for (vector<short>::iterator rm = dot_additions2template.begin(); rm != dot_additions2template.end(); rm++) {
					if ( (*rm) ) {
						delete (*sPit)->indel.L_ins_;
						delete (*sPit)->indel.R_ins_;
						delete (*sPit)->indel.del;
						delete (*sPit);
						site_properties.erase(sPit);
					}
					++sPit;
				}
			}

			//////////
			/// Set the vector to the motif site we are updating.
			//////////
			vector<Site>::iterator m2Sit = root->seq_evo.begin() + motif_start_position;
			vector<short>::iterator rm = dot_additions2template.begin();

			// Get all of the siteProps to match the motifSites. Set active props later.
			for (vector<siteProperties*>::iterator it3 = site_properties.begin(); it3 != site_properties.end(); ++it3) {
				if (relationship == EXACT) {
					if ( (*it2)->isTemplate() && (*rm) ) {
						///////////
						/// This is the case that the motif has '.', but template does not.
						///////////
						(*m2Sit).push_siteProps(*it3);
					} else {
						(*m2Sit).push_siteProps(*it3);
					}
				} else if (relationship == NODE_IS_ANC) {
					///////////
					/// Push back property to the site prop. Only set dormant if not already an active site.
					///////////
					(*m2Sit).push_siteProps(*it3);
				} else {
					cerr << "Root node can be nothing but EXACT or NODE_IS_ANC" << endl;
					exit(EXIT_FAILURE);
				}
				rm++;
				++m2Sit;
			}
		}
	}

	root->Site_postProcess();
	root->seq_evo.front().motif.site_props.push_back(new siteProperties());
	root->seq_evo.front().motif.site_props.back()->indel.createObjects();
	root->seq_evo.front().motif.site_props.front()->indel.del->my_motif_varSite 
	= root->seq_evo.front().motif.site_props.front()->indel.L_ins_->my_motif_varSite_left
	= root->seq_evo.front().motif.site_props.front()->indel.L_ins_->my_motif_varSite_right
	= NULL;
	root->seq_evo.front().motif.site_props.back()->indel.del->my_sequence_template_varSite 
	= root->seq_evo.front().motif.site_props.back()->indel.L_ins_->my_sequence_template_varSite_left
	= root->seq_evo.front().motif.site_props.back()->indel.L_ins_->my_sequence_template_varSite_right
	= NULL;
	root->seq_evo.front().motif.site_props.front()->indel.del->my_sequence_template_varSite 
	= root->seq_evo.front().motif.site_props.front()->indel.L_ins_->my_sequence_template_varSite_left
	= root->seq_evo.front().motif.site_props.front()->indel.L_ins_->my_sequence_template_varSite_right
	= root->one_site_varSite;
	root->seq_evo.front().motif.site_props.back()->indel.del->my_motif_varSite 
	= root->seq_evo.front().motif.site_props.back()->indel.L_ins_->my_motif_varSite_left
	= root->seq_evo.front().motif.site_props.back()->indel.L_ins_->my_motif_varSite_right
	= root->one_site_varSite;
	root->seq_evo.front().motif.site_props.push_front(new siteProperties());
	root->seq_evo.front().motif.site_props.front()->indel.createObjects();
	root->evolvingSequence->setActiveProps();

	//////////
	/// Remove sites that are '-', populate varSite list of pointers to motifSites that belong to
	/// it. This will effectively remove the necessity of having the relateMotifSites function.
	//////////
	for (list<inMotif*>::iterator it2 = motif_specs.begin(); it2 != motif_specs.end(); ++it2) {
		found = root_sequence.rfind("-");
		while (found != string::npos) {
			(*it2)->sitemap.erase(found,1);
			found = root_sequence.rfind("-",found-1);
		}
	}

	sit = root_sequence.begin();
	sit2 = invariable_in.begin();
	int site_number = 0, num = 0;
	for (vector<Site>::iterator site_it = root->seq_evo.begin(); 
								site_it != root->seq_evo.end(); 
								++site_it) 
	{
		if ( (*site_it).returnState() == '-') {
			if (site_it != root->seq_evo.end()-1 ) {
				(*site_it).motif.Site_deleteMerge( (&(*(site_it-1)).motif), (&(*(site_it+1)).motif), false);
			} else {
				(*site_it).motif.Site_deleteMerge( (&(*(site_it-1)).motif), NULL, false);
			}
			root->seq_evo.erase(site_it--);
			root_sequence.erase(sit--);
			invariable_in.erase(sit2--);
		} else site_number++;
		++sit, ++sit2;
	}

	//cerr << "TTree::constructRootMotifSites: root_sequence:" << endl << root_sequence << endl;
	//cerr << endl << endl << "TTree::constructRootMotifSites: ROOT" << endl << endl;
	//root->Print_Active_Properties();
	//cerr << endl << endl << "TTree::constructRootMotifSites: ENDROOT" << endl << endl;
	//exit(0);
}

//////////
/// Similar, non-global function like constructRootMotifSites, used primarily to create a sequence
/// that is pseudogene'd.
//////////
void
TNode::constructUnconstrainedSequence(
										   string& sequence
										  ) 
{
	//////////
	/// Build the sites with which we will play with
	//////////
	evolvingSequence = new Sequence(this, sequence.size());
	evolvingSequence->init(this, sequence, nodeEnv);
    //
	Site_postProcess();
    //
	seq_evo.front().motif.site_props.push_back(new siteProperties());
	seq_evo.front().motif.site_props.back()->indel.createObjects();
	seq_evo.front().motif.site_props.front()->indel.del->my_motif_varSite 
	= seq_evo.front().motif.site_props.front()->indel.L_ins_->my_motif_varSite_left
	= seq_evo.front().motif.site_props.front()->indel.L_ins_->my_motif_varSite_right
	= NULL;
	seq_evo.front().motif.site_props.back()->indel.del->my_sequence_template_varSite 
	= seq_evo.front().motif.site_props.back()->indel.L_ins_->my_sequence_template_varSite_left
	= seq_evo.front().motif.site_props.back()->indel.L_ins_->my_sequence_template_varSite_right
	= NULL;
	seq_evo.front().motif.site_props.front()->indel.del->my_sequence_template_varSite 
	= seq_evo.front().motif.site_props.front()->indel.L_ins_->my_sequence_template_varSite_left
	= seq_evo.front().motif.site_props.front()->indel.L_ins_->my_sequence_template_varSite_right
	= one_site_varSite;
	seq_evo.front().motif.site_props.back()->indel.del->my_motif_varSite 
	= seq_evo.front().motif.site_props.back()->indel.L_ins_->my_motif_varSite_left
	= seq_evo.front().motif.site_props.back()->indel.L_ins_->my_motif_varSite_right
	= one_site_varSite;
	seq_evo.front().motif.site_props.push_front(new siteProperties());
	seq_evo.front().motif.site_props.front()->indel.createObjects();
	// Set the properties.
	evolvingSequence->setActiveProps();
	//
	//////////

//	cerr << endl << endl << "des_cUS" << endl << endl;
//	Print_Active_Sites();
//	cerr << endl << endl << "ENDdes_cUS" << endl << endl;
//	exit(0);
}

list<inMotif*> 
TTree::DrawMotifs(
								 string& root_sequence, 
								 double proportion_motif
								) 
{
	vector<inMotif*> prosite_library;
	vector<string> prosite_motifs;
	list<string> motif_data;
	vector<bool> occupied (root_sequence.size(), false);
	char beginning_marker = 'Z';
	size_t motif_size, motif_max;
	vector<int> acceptable_positions;
	bool exclude_last_position = false;

	// Extract the prosite motifs, place them into prosite_library. //
	Populate_Prosite_Motifs(prosite_motifs);
	for (vector<string>::iterator it = prosite_motifs.begin(); it != prosite_motifs.end(); ++it) {	
		motif_data = split((*it), ";");
		motif_data.pop_back(); 		// Remove a null element at end.
		list<string>::iterator jt = motif_data.begin();
		prosite_library.push_back(new inMotif((*jt), (*(++jt)), tips.size()));
	}

	// Start to cover the sequence with motifs.
	double proportion_sequence_covered = 0;
	int rand_motif;
	list<siteRegEx*> curr_regex;
	list<inMotif*> motifs_to_place;
	size_t total_motifs = 0, sequence_size = root_sequence.size();
	size_t attempts;
	while (attempts != 20) {
		attempts = 0;
		do {
			// Choose motif from library
			rand_motif = (double)rndu() * (prosite_library.size()-1);
			curr_regex = prosite_library.at(rand_motif)->parseRegEx();

			if (curr_regex.front()->last_motif_site_optional) {
				double test_value = (double)rndu();
				// If last site is optional, pop the back element if it is found. //
				if (((double)rndu() < test_value) ? 1 : 0 ) {
					curr_regex.pop_back();
					exclude_last_position = true;	// Need to remove library position, but only do if motif accepted
				} 
			}
			
			// Determine motif size. This includes randomly drawing the size of x(a,b) sites. //
			motif_size = motif_max = 0;
			for (list<siteRegEx*>::iterator jt = curr_regex.begin(); jt != curr_regex.end(); ++jt) {
				motif_max += (*jt)->allowable_characters.size();			
				motif_size += (*jt)->sites_occupied;
			}

			size_t site_num = 1;
			acceptable_positions.clear();
			if (motif_size < occupied.size()) {	// Motifs larger than sequence = infinite loop. //
				bool acceptable = true;
				if (curr_regex.front()->N_term_motif) {
					vector<bool>::iterator kt = occupied.begin()+1;
					for (size_t i = 0; i < motif_size; i++) {
						if ( (*(kt+i)) ) 
							acceptable = false; 
					}
					if (acceptable) {
						acceptable_positions.push_back(1);
					}					
				} else if (curr_regex.front()->C_term_motif) {
					for (vector<bool>::iterator kt = occupied.end()-motif_size; kt != occupied.end(); ++kt) {
						if (*kt) acceptable = false;
					}
					if (acceptable) {
						acceptable_positions.push_back(root_sequence.size()-motif_size);
						if (exclude_last_position) prosite_library.at(rand_motif)->removeLastPosition();
					}
				} else {
					for (vector<bool>::iterator kt = occupied.begin()+1; kt != occupied.end()-motif_size; ++kt, site_num++) {
						acceptable = true;
						for (size_t i = 0; i < motif_size; i++) {
							if ( (*(kt+i)) ) 
								acceptable = false; 
						}
						if (acceptable) acceptable_positions.push_back(site_num);
					}
				}
				attempts++;
			} else attempts = 20;
		} while ( proportion_sequence_covered + ((double)motif_size / sequence_size) > proportion_motif && attempts < 20);

		// Place the motifs.
		if (attempts != 20 && acceptable_positions.size() > 0) { 	// A suitable motif was found, so place it.
			Place_Motif(curr_regex, beginning_marker--, motifs_to_place, prosite_library.at(rand_motif), root_sequence, acceptable_positions.at((double)rndu()*(acceptable_positions.size()-1)));
			proportion_sequence_covered += ((double)motif_size / sequence_size);
		}

		// Since the root sequence size may change after motifs with x(a,b)'s, reconstruct the
		// occupied vector, setting all previous motifs to true.
		occupied.assign(root_sequence.size(), false);
		size_t first, last;
		for (list<inMotif*>::iterator it = motifs_to_place.begin(); it != motifs_to_place.end(); ++it) {
			first = (*it)->sitemap.find_first_not_of("*");
			last  = (*it)->sitemap.find_first_of("*", first+1);
			if (last == string::npos) last = (*it)->sitemap.size();
			for (size_t i = first; i < last; i++) occupied.at(i) = true;
		}
	}
	
	for (vector<inMotif*>::reverse_iterator rit = prosite_library.rbegin(); rit != prosite_library.rend(); ++rit) 
		delete *rit;

	return motifs_to_place;
}

//////////
/// Place motif into random sequence.
//////////
// - Creates the inMotif specs for the motif.
// - Places dashes in the root sequence for x(a,b) regexes, and pads those sites with '*' in previously placed motifs
// - Changes the root sequence to match the motif placed.
//////////
void 
TTree::Place_Motif(
					    list<siteRegEx*>& motif, 
					    char mark, 
					    list<inMotif*>& placed_motifs, 
					    inMotif *thisMotif, 
					    string& root_sequence, 
					    size_t start_site
					   ) 
{
	string sitemap (root_sequence.size(), '*');

	//////////
	/// Main loop. Cycle through the siteRegEx'es and make changes.
	//////////
	size_t sequence_site = start_site;
	for (list<siteRegEx*>::iterator it = motif.begin(); it != motif.end(); ++it) {
		// Going over all of the occupied sites.
		size_t num_added = 0;
		for (list<string>::iterator jt = (*it)->allowable_characters.begin();
									jt != (*it)->allowable_characters.end() && num_added < (*it)->sites_occupied; 
									++jt, sequence_site++, num_added++) 
		{
			root_sequence.at(sequence_site) = (*jt).at((double)rndu()*((*jt).size()-1));
			sitemap.at(sequence_site) = mark;
		}

		//////////
		/// Now pad root sequence with dashes for unfilled x(a,b) sites.
		//////////
		for (size_t j = num_added; j < (*it)->allowable_characters.size(); j++, sequence_site++) {
			root_sequence.insert(sequence_site, "-");
			sitemap.insert(sequence_site, 1, mark);
			for (list<inMotif*>::iterator jt = placed_motifs.begin(); jt != placed_motifs.end(); ++jt) {
				(*jt)->sitemap.insert(sequence_site, "*");
			}
		}

	}

	//////////
	/// Create new inMotif, copying the input inMotif.
	//////////
	placed_motifs.push_back(new inMotif(thisMotif, sitemap, mark));
}

TNode*
TTree::AddNode() 
{
	TNode *node = new TNode(this);
	nodeList.push_back(node);
	
	return node;
}

void 
TNode::resetSequence( TNode *node )
{
	vector<Site>::iterator my, it;
	my = seq_evo.begin();
	it = node->seq_evo.begin();
	
	for (; my != seq_evo.end(); ++my, ++it) (*my).setState((*it).returnState());
}

TNode*
TTree::copyNode( TNode *copy_node )
{
	TNode *node = new TNode(this);
	delete node->branch;		// When using this function, always making scratch seq that will use original seq branch.

	string sequence = copy_node->printSequence(false);
	// Sets up activeProperties, etc.
	node->nodeEnv = copy_node->nodeEnv;
	node->bipartition = copy_node->bipartition;
	node->branch = copy_node->branch;
	node->branch->rates->Qij = copy_node->branch->rates->Qij;
	node->anc = copy_node->anc;
	node->constructUnconstrainedSequence(sequence);
	node->inheritCategories(copy_node);

	return node;
}

inClade*
TTree::AddClade(
						 string& label
						) 
{
	inClade *clade = new inClade(label, this);
	treeEnv.push_back(clade);

	return clade;
}

TNode*
TTree::ReadTip(
					  string& tree_str, 
					  int *pos, 
					  char ch, 
					  TTree *tree, 
					  bool first_partition
					 )
{
	string name;
	size_t found;
	
	TNode *node = AddNode();
	
	found = tree_str.find_first_of(":,)",*pos);
	name = tree_str.substr(*pos, found - *pos);

	if(found - *pos > MAX_NAME_LEN) {
		cerr << "Name '" << name << "' longer than maximum allowed length." << endl;
		exit(EXIT_FAILURE);
	}

	*pos=found;	

	if (first_partition) {
		node->tipNo=tree->numTips;
		node->node_name = name;
		names.push_back(name);
	} else {
		size_t i = 0;
		while (i < names.size() && name.compare(names.at(i)) != 0)
			i++;
		if (i == names.size()) {
			sprintf(treeErrorMsg, "Taxon names in trees for different partitions do not match.");
			return NULL;
		}
		node->tipNo=i;
	}

	tree->tips.push_back(node);
	tree->numTips++;

	found = tree_str.find_first_of(":,)",*pos);
	*pos = found-1;

	if (found == tree_str.npos) {
		sprintf(treeErrorMsg, "Unexpected end of file");
		return NULL;
	}

	return node;
}

TNode*
TTree::ReadNode(
					   string& tree_str, 
					   int *pos, 
					   TTree *tree, 
					   int detectPolytomies, 
					   bool first_partition
					  )
{
	TNode *node, *node2;
	char ch;
	int found;

	try {
		node = AddNode();
	} catch (exception& e) {
		cerr << "Error creating TNode node (TTree::ReadNode): " << e.what() << endl;
		exit(EXIT_FAILURE);
	}

	if ((node2=ReadBranch(tree_str, pos, tree, first_partition))==NULL)
		return NULL;

	node->branch1=node2;
	node2->branch0=node;
	node->branch->length1=node2->branch->length0;

	ReadUntil(tree_str, pos, ',', "Comma");

	if (treeError)
		return NULL;
	if ((node2=ReadBranch(tree_str, pos, tree, first_partition))==NULL)
		return NULL;
	node->branch2=node2;
	node2->branch0=node;
	node->branch->length2=node2->branch->length0;
	
	found = tree_str.find_first_of(":,);", ++(*pos));
	*pos=found;
	ch = tree_str.at(found);
	
	if (detectPolytomies && ch==',') {
		fprintf(stderr, "This tree contains nodes which aren't bifurcations. Resolve the node\n");
		fprintf(stderr, "with zero branch lengths to obtain correct results. This can be done\n");
		fprintf(stderr, "with a program called TreeEdit: http://evolve.zoo.ox.ac.uk/software/TreeEdit\n");
		exit(EXIT_FAILURE);
	}

	if ((size_t)*pos == tree_str.npos) {
		sprintf(treeErrorMsg, "Unexpected end of file");
		return NULL;
	}
	(*pos)--;
	
	return node;
}

TNode*
TTree::ReadBranch(
						 string& tree_str, 
						 int *pos, 
						 TTree *tree, 
						 bool first_partition
						)
{

	char ch;
	double len, param=0.0;
	TNode *node;
	int found;
	inClade *clade;

	ch=ReadToNextChar(tree_str, pos);

	if (ch=='(') {	// is a node
		node=ReadNode(tree_str, pos, tree, 1, first_partition);
		ReadUntil(tree_str, pos, ')', "Closing bracket");
		if (treeError)
			return NULL;
	} else {		// is a tip
		node=ReadTip(tree_str, pos, ch, tree, first_partition);
	}
	
	ch=ReadToNextChar(tree_str, pos);

	if (isalpha(ch)) {	// Specifying clade:
		string clade_label = "";
		while(isalpha(ch) || isdigit(ch) || ch == '_') {
			clade_label += ch;
			ch=ReadToNextChar(tree_str,pos);
		}
		clade = AddClade(clade_label);
		if(clade_label.size() <= MAX_NAME_LEN) {
			node->clade_label.assign(clade_label);
		} else {
			cerr << "Clade name " << clade_label << " has too many characters. Please shorten label." << endl;
			exit(EXIT_FAILURE);
		}
	}
	
	if (ch==':') {
		if (tree->lengths==0) {
			sprintf(treeErrorMsg, "Some branches don't have branch lengths");
			return NULL;
		} else 
			tree->lengths=1;

		if (sscanf(&tree_str.at(++(*pos)), "%lf", &len)!=1) {
			sprintf(treeErrorMsg, "Unable to read branch length");
			return NULL;
		}

		found = tree_str.find_first_not_of(".0192837465",*pos);
		*pos = found-1;

		ch=ReadToNextChar(tree_str,pos);

		if (ch=='[') {
			if (sscanf(&tree_str.at(++(*pos)), "%lf", &param)!=1) {
				sprintf(treeErrorMsg, "Unable to read branch parameter");
				return NULL;
			}
			ReadUntil(tree_str, pos, ']', "Close square bracket");
		} else { (*pos)--; }
	} else {
		if (tree->lengths==1) {
			sprintf(treeErrorMsg, "Some branches don't have branch lengths");
			return NULL;
		} else 
			tree->lengths=0;
	
		len=0.0;
		(*pos)--;
	}

	node->branch->length0=len;
	node->branch->param=param;
	
	return node;
}	

void 
TTree::ReadTree(
					 string& tree_str, 
					 int *pos, 
					 TTree *tree, 
					 bool first_partition
					)
{
	TNode *P;
	size_t found;

	treeError=0;
	tree->numNodes=0;
	tree->numTips=0;
	tree->rooted=1;
	tree->lengths=-1;

    if (tree_str.at(*pos) !='(' || (tree->root=ReadNode(tree_str, pos, tree, 0, first_partition))==NULL) {
		fprintf(stderr, "Error reading tree: %s.\n", treeErrorMsg);
		exit(EXIT_FAILURE);
    }

    found = tree_str.find(",);",++(*pos));

    if (tree_str.at(*pos) == ',') {		
		//////////
		/// iSG currently is not equipped to handle unrooted trees. This may change in the future.
		//////////
		cerr << PROGRAM_NAME << VERSION_NUMBER << " currently does not support unrooted trees." << endl;
		cerr << "Note that you may make the unrooted tree into a rooted tree by adding a zero length subtree for any non-bifurcating subtree in Newick format. EX: (Taxon1:a, Taxon2:b, Taxon3:c); --> (Taxon1:a, (Taxon2:b, Taxon3:c):0.0);" << endl << endl;
		for (size_t i = 0; i <= *pos; i++) {
			cerr << tree_str.at(i);
		}
		cerr << "<--- Trifurcation occurs here." << endl;
		exit(EXIT_FAILURE);
    }

    tree->totalLength=0.0;

    if (tree->rooted) {
		P=tree->root;
		while (P!=NULL) {
	    	tree->totalLength+=P->branch->length0;
	    	P=P->branch1;
		}
    }
    
    for (list<inClade*>::iterator sit = treeEnv.begin(); sit != treeEnv.end(); ++sit) {
    	for (int i = 0; i < names.size(); i++) {
    		if (((names.at(i)).compare((*sit)->clade_name) == 0)) {
				cerr << "Error parsing tree: Clade name '" << names.at(i) << "' is the name of a taxon." << endl;
				exit(EXIT_FAILURE);
			}
    	}
    }
}

void 
TTree::report_clades() 
{
	cerr << "Clades: " << endl;
	for (list<inClade*>::iterator it = treeEnv.begin(); it != treeEnv.end(); ++it) {
		cerr << (*it)->clade_name << " reporting for duty.";
		if (!(*it)->environment_name.empty())
			cerr << " But you can call me " << (*it)->environment_name << ".";
		cerr << endl;
	}
}

char 
ReadToNextChar(
					string& tree_str, 
					int *pos
				   )
{
	int found;	

	found = tree_str.find_first_not_of(" \t\n",++(*pos));
	*pos = found;

	return tree_str.at(*pos);
}

void 
ReadUntil(
			   string& tree_str, 
			   int *pos, 
			   char stopChar, 
			   const char *what
			  )
{
	size_t found;
	string delims = "(,:);";
	delims += stopChar;

	found = tree_str.find_first_of(delims,++(*pos));
	*pos = found;

	if (found == tree_str.npos || tree_str.at(*pos)!=stopChar) {
		sprintf(treeErrorMsg, "%s missing", what);
		treeError=1;
	}
}

void 
TTree::report_branches()
{
	for (list<TNode*>::iterator it = nodeList.begin(); it != nodeList.end(); ++it) {
		if ( (*it)->tipNo != -1 ) cerr << "----------TIP---------" << endl;
		cerr << "length0:                      " << " " << (*it)->branch->length0 << endl;
		cerr << "length1:                      " << " " << (*it)->branch->length1 << endl;
		cerr << "length2:                      " << " " << (*it)->branch->length2 << endl;
		cerr << "param:                        " << " " << (*it)->branch->param << endl;

		cerr << "branch1_time_relative_length: " << " " << (*it)->branch->branch1_time_relative_length << endl;
		cerr << "branch2_time_relative_length: " << " " << (*it)->branch->branch2_time_relative_length << endl;
		cerr << "branch0_time_relative_length: " << " " << (*it)->branch->branch0_time_relative_length << endl;

		cerr << "branch1_max_path:             " << " " << (*it)->branch->branch1_max_path << endl;
		cerr << "branch2_max_path:             " << " " << (*it)->branch->branch2_max_path << endl;
		
		cerr << "perturbation:                 " << " " << (*it)->branch->perturbation << endl;

		cerr << endl;
	}
}

void 
TTree::report_sequences( void )
{
	for (list<TNode*>::iterator it = nodeList.begin(); it != nodeList.end(); ++it) {
		cerr << "tipNo: " << (*it)->tipNo 
			 << " mytipNo: " << (*it)->mytipNo 
			 << " anc mytipNo: " << (*it)->anc->mytipNo << endl;
		for (vector<Site>::iterator jt = (*it)->seq_evo.begin(); jt != (*it)->seq_evo.end(); ++jt) {
			if ((*jt).returnState() == -1) cerr << "X";
			else if ((*jt).returnState() < stateCharacters.size()) cerr << stateCharacters.at((*jt).returnState());
			else cerr << "TTree::report_sequences -> sequence state \"" << (*jt).returnState() << "\" set." << endl;
		}
		cerr << endl;
	}
}

void 
globalArray::Print()
{
	size_t site_number = 0;
//	for (vector<insertSite>::iterator it = insert_sites.begin(); it != insert_sites.end(); ++it, site_number++) {
//		for (list<siteModifier>::iterator jt = (*it).modifiers.begin(); jt != (*it).modifiers.end(); ++jt) 
//			cout << "***** " << (*jt).action << ", " << (*jt).indelNo << ", " << (*jt).fromAnc << endl;
//		cout << site_number << ": " << (*it).action << ", ";
//		cout << (*it).indelNo << ", " << (*it).fromAnc << endl;
//	}

	for (vector<insertSite>::iterator it = insert_sites.begin(); it != insert_sites.end(); ++it) {
		for (list<siteModifier>::iterator jt = (*it).modifiers.begin(); jt != (*it).modifiers.end(); ++jt) {
			if ( (*jt).fromAnc == -1) cerr << "-";
			else cerr << (*jt).action;
			cerr << " ";
			if ( (*jt).fromAnc == -1) cerr << "-";
			else cerr << (*jt).fromAnc;
			cerr << " ";
			if ( (*jt).indelNo == -1) cerr << "-";
			else cerr << (*jt).indelNo;
			cerr << endl;
		}
		if ( (*it).fromAnc == -1) cerr << "-";
		else cerr << (*it).action;
		cerr << " ";
		if ( (*it).fromAnc == -1) cerr << "-";
		else cerr << (*it).fromAnc;
		cerr << " ";
		if ( (*it).indelNo == -1) cerr << "-";
		else cerr << (*it).indelNo;
		cerr << endl;
	}
	cout << endl;
//	for (vector<insertSite>::iterator it = insert_sites.begin(); it != insert_sites.end(); ++it) {
//		for (list<siteModifier>::iterator jt = (*it).modifiers.begin(); jt != (*it).modifiers.end(); ++jt) {
//			if ( (*jt).fromAnc == -1) cout << "-";
//			else cout << (*jt).fromAnc;
//		}
//		if ( (*it).fromAnc == -1) cout << "-";
//		else cout << (*it).fromAnc;
//	}
//	cout << endl;
//	for (vector<insertSite>::iterator it = insert_sites.begin(); it != insert_sites.end(); ++it) {
//		for (list<siteModifier>::iterator jt = (*it).modifiers.begin(); jt != (*it).modifiers.end(); ++jt) {
//			if ( (*jt).indelNo == -1) cout << "-";
//			else cout << (*jt).indelNo;
//		}
//		if ( (*it).indelNo == -1) cout << "-";
//		else cout << (*it).indelNo;
//	}
//	cout << endl;
	
}

vector<insertSite>::iterator 
globalArray::locateEvent(
						 TNode *des, 
						 size_t position
						)
{
	int at_pos = 0;
	bool isSite;
	vector<insertSite>::iterator return_value;

//cerr << "Locating event at position: " << position << " from " << des->printBipartition() << endl;
//	cerr << "globalArray::locateEvent" << endl; Print();

	for (vector<insertSite>::iterator it = insert_sites.begin(); it != insert_sites.end(); ++it) {
		isSite = true;
//		cerr << "site " << at_pos << "    ";
		if ( isAnc(des, (*it).fromAnc) ) {
			for (list<siteModifier>::iterator jt = (*it).modifiers.begin(); jt != (*it).modifiers.end(); ++jt) {
				if ( (*jt).action == 'd' ) {	// Site has potentially been removed.
					if ( isAnc(des, (*jt).fromAnc) ) {
//						cerr << "Is not a site.";
						isSite = false;
					}
				}
			}
		} else isSite = false;
//		if (isSite) cerr << "is a site.";
//		cerr << endl;
		if (at_pos == position) return it;
		if (isSite) at_pos++;
	}

//cerr << "Event is at the end of the sequence." << endl;
	
	return insert_sites.end();
}

//////////
/// For MCMC, the arrays get huge keeping track of substitution events, even if the events are no
/// longer needed (e.g., from a previous MCMC cycle that was rejected, etc). This finds the minimum
/// eventID of the current path, and deletes all other events in the globalArray that are of a smaller
/// eventID.
//////////
void 
globalArray::cleanArray( list<eventTrack*>& events )
{
	// Save substitutions that exist in the current path.
	int num_saved = 0;
	for (list<eventTrack*>::iterator ht = events.begin(); ht != events.end(); ++ht) {
		for (vector<insertSite>::iterator it = insert_sites.begin(); it != insert_sites.end(); ++it) {
			for (list<siteModifier>::iterator jt = (*it).modifiers.begin(); jt != (*it).modifiers.end(); ++jt) {
				if ( (*jt).indelNo == (*ht)->ID || (*jt).action != 's') {
					(*jt).save = true;
				}
			}
		}
	}

	// Remove all substitutions that have not been saved.
	int site = 0, mod = 0;
	for (vector<insertSite>::iterator it = insert_sites.begin(); it != insert_sites.end(); ++it, ++site) {
		mod = 0;
		for (list<siteModifier>::iterator jt = (*it).modifiers.begin(); /*mod < mod_size*/ jt != (*it).modifiers.end(); ++jt, ++mod ) {
			if ( !(*jt).save ) {
				(*it).modifiers.erase(jt);
				// Major issues if I try to make compensatory changes to jt for erasure. 
				// Setting to begin should not make a big difference in run-time.
				jt = (*it).modifiers.begin();
			} 
		}
		
		// Reset all save flags back to false before advancing to next site in globalArray.
		for (list<siteModifier>::iterator jt = (*it).modifiers.begin(); jt != (*it).modifiers.end(); ++jt) 
			(*jt).save = false;
	}
}

bool 
globalArray::isDeleted(
							TNode *des, 
							insertSite insert_site
						   )
{
	for (list<siteModifier>::iterator jt = insert_site.modifiers.begin(); jt != insert_site.modifiers.end(); ++jt) {
		if ( (*jt).action == 'd' ) {	// Site has potentially been removed.
			if ( isAnc(des, (*jt).fromAnc) ) return false;
		}
	}
	return true;
}

