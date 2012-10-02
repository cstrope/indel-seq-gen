#include "motif.h"

extern vector<short> acgt;
extern vector<int> acagatcgctgt;

vector<int> chosenACGT (20, 0);
vector<int> choseACGT (20, 0);
vector<int> virtualACGT (20,0);
vector<int> value_check (528,0);

////////////////////
////// SEQUENCE
////////////////////
//////////
/// Initialize the created Sequence objects.
//////////
void Sequence::init(
					TNode *node,
				    string initial_state, 
				    inClade *env, 
				    motifSite *site_type, 
				    varSite *template_varSite, 
				    varSite *motif_varSite
				   )
{
	string::size_type i = 0;
	my_node = node;
	for (vector<Site>::iterator it = evolutionaryAttributes.begin(); 
								it != evolutionaryAttributes.end(); 
								++it) 
	{
		//////////
		/// Initialize items that cannot be done in the initialization list.
		//////////
		// Need to create motifSites for each of the <Site>s
		// Cannot do in initialization list -> All motifSites site_props are the same, causes
		// seg faults for each deletion of a new site.
		(*it).InitializeMotif(site_type, template_varSite, motif_varSite);
		//
		// Set the state
		(*it).setState(initial_state.at(i++));
		//
		// Set the initial rates
		switch(env->rateHetero) {
		case GammaRates:
			(*it).setGamma(rndgamma(env->gammaShape) / env->gammaShape);
			break;
		case DiscreteGammaRates:
			(*it).setCategory((int)(rndu()*env->numCats));
			(*it).setGamma(0);		// Flag saying that categories are being used. //
			break;
		case CodonRates:
			break;
		case NoRates:
			break;
		default:
			cerr << "Invalid rate category: " << env->rateHetero << endl;
			exit(EXIT_FAILURE);
		}
		//
		//////////
	}
}

void 
Sequence::ConvertRootSequence() 
{
	int i = 0;
	char this_state;
	size_t putative_state;
	string stateCharacters;
	if(isNucModel) stateCharacters = "ACGT";
	else stateCharacters = "ARNDCQEGHILKMFPSTWYV";

	//////////
	/// Put the root sequence into the evolving sequence structure.
	//////////
	for (vector<Site>::iterator it = evolutionaryAttributes.begin(); it != evolutionaryAttributes.end(); ++it, i++) {
		if ( (*it).returnState() == 'U' && isNucModel) (*it).setState('T');
		this_state = (*it).returnState();
		putative_state = stateCharacters.find(this_state);
		if (putative_state > stateCharacters.size()) {
			cerr << "Invalid character '" << this_state << "' at position " << i << " (" << putative_state << "/" << numStates << ")." << endl;
			exit(EXIT_FAILURE);
		}
		(*it).setState(putative_state);
	}
}

int 
Sequence::compare_sequence (
								Sequence *seq
							   )
{
	int num_diff = 0;

	vector<Site>::iterator jt = seq->evolutionaryAttributes.begin();
	for (vector<Site>::iterator it = evolutionaryAttributes.begin(); it != evolutionaryAttributes.end(); ++it, ++jt)
		if ( (*it).returnState() != (*jt).returnState() ) 
			num_diff++;

	return num_diff;
}

void 
Sequence::print_sequence( )
{
	for (vector<Site>::iterator it = evolutionaryAttributes.begin(); it != evolutionaryAttributes.end(); ++it) {
		if ((*it).isSet())
			cerr << stateCharacters[(*it).returnState()];
		else
			cerr << "X";
	}
	cerr << endl;
}

void Sequence::printSequenceRateAway()
{
	for (vector<Site>::iterator it = evolutionaryAttributes.begin(); it != evolutionaryAttributes.end(); ++it) {
		if ((*it).isSet()) {
			cerr << stateCharacters[(*it).returnState()] << ":   ";
			(*it).printSiteRateAway();
		} else {
			cerr << "X";
		}
		cerr << endl;
	}

}

//////////
/// Each time we calculate the rate away from a sequence, the rates for all possible changes is agnostic as to
/// whether the change will result in a stop codon being formed. This function checks such changes, and whenever
/// one possible stop codon can be formed, zeroes out the offending substitution.
/// - This function assumes that we are working with coding-data (block size 3), and that there is no sequence that
///   is not an order of 3 large. 
//////////
double
Sequence::zeroStopCodons (
						  int event_site,
						  int end_site
						 )
{
	int codon[3];
	// If we end up zeroing some possible transition, this needs to be reflected in the overall sequence probability.
	// For this, we return the amount that the sequence needs to be adjusted by.
	double probability_adjustment = 0;
	// First, set the event site to the nearest codon boundary.
	event_site -= 3 + (event_site%3);
	if (event_site < 0) event_site = 0;

	// and the end site similarly.
	vector<Site>::iterator end;
	end_site += 3 + (3-(end_site%3));		// ? End site at the end of a codon should be a 2... but if less than end site?
	if (end_site > evolutionaryAttributes.size()) {
		end = evolutionaryAttributes.end(); cerr << "***" << endl;
	} else {
		end = evolutionaryAttributes.begin()+end_site; cerr << "!!!" << endl;
	}

	cerr << "now event_site=" << event_site << " and end_site=" << end_site << endl;

	int i = event_site;
	for (vector<Site>::iterator it = evolutionaryAttributes.begin()+event_site; it != end; it+=3, i+=3) {
cerr << "site " << i << endl;
		codon[0] = (*it).returnState();
		codon[1] = (*(it+1)).returnState();
		codon[2] = (*(it+2)).returnState();
		
		// All stop codons start with T, this checks for that.
		if (codon[0] == 3) {		// T
			// This checks to see if any of the other positions match a stop codon. If so, then this
			// position is 1 change away from becoming a stop codon. So, we zero out the offending
			// change... make it impossible to happen.
			if (codon[1] == 0) {
				probability_adjustment += (*(it+1)).site_rate_away.at(0);
				(*(it+1)).site_rate_away.at(0) = 0;
			} else if (codon[1] == 2) {
				probability_adjustment += (*(it+1)).site_rate_away.at(2);
				(*(it+1)).site_rate_away.at(2) = 0;
			} else if (codon[2] == 0) {
				probability_adjustment += (*(it+2)).site_rate_away.at(0);
				(*(it+2)).site_rate_away.at(0) = 0;
			} else if (codon[2] == 2) {
				probability_adjustment += (*(it+2)).site_rate_away.at(2);
				(*(it+2)).site_rate_away.at(2) = 0;
			} // else this codon is fine.
			// Now check to see if second and third codon positions are stoppers.
		} else if (codon[1] == 0) {
			if (codon[2] == 2 || codon[2] == 0) {
				probability_adjustment += (*it).site_rate_away.at(3);
				(*it).site_rate_away.at(3) = 0;
			}
		} else if (codon[1] == 2) {
			if (codon[2] == 0) {
				probability_adjustment += (*it).site_rate_away.at(3);
				(*it).site_rate_away.at(3) = 0;
			}
		}
	}
	cerr << "EXITING ZEROSTOPCODONS" << endl;
	return probability_adjustment;
}

void 
Sequence::setActiveProps(
						 bool insertion
						)
{
	size_t current_site = 1;
	varSite *right_motif_vS, *left_motif_vS, *on_motif_vS;

	//////////
	/// my_node same as anc, root build. Need to build the entire sequence. Otherwise, it is an
	/// insert sequence.
	//////////
	if (!insertion) {
		current_site = 1;
		//////////
		/// Set the first site
		//////////
		evolutionaryAttributes.front().motif.active_properties.setFirstMotifPos(my_node);
		assert(my_node->anc != NULL);
		vector<Site>::iterator anc_it = my_node->anc->seq_evo.begin()+1;
		for (vector<Site>::iterator it = evolutionaryAttributes.begin()+1; 
									it != evolutionaryAttributes.end(); 
									++it, current_site++, ++anc_it) {
			//////////
			/// Set the active_properties for this site:
			//////////
			// Make the containers:
			(*it).motif.active_properties.indel = new Indel();
			(*it).motif.active_properties.indel->createObjects();
			(*it).motif.active_properties.subst = new Substitution();
			//
			// Set R_ins_ of the previous position to the L_ins_ of the current position (first set by setFirstMotifPos).
			// * delete the first R_ins_ (unnecessary for all other active_properties setting)
			// * Keep track of varSite pointers of deleted position: they will need to be reset
			//   after new position is made (only for post-motif inheritance?) May also need to
			//   do for between-sequence-template positions.
			delete (*(it-1)).motif.active_properties.indel->R_ins_;
			(*(it-1)).motif.active_properties.indel->R_ins_
			= (*it).motif.active_properties.indel->L_ins_;
			//
			// Set R_ins_, subst, del
			(*it).motif.active_properties.set_properties(
														 (*it).motif.site_props, 
														 (*anc_it).motif.site_props,
														 my_node,
														 &(*(it-1)).motif.active_properties,
														 current_site
														);
			//
			///////////
		}
		//////////
		/// Set the last R_ins_ site.
		//////////
		evolutionaryAttributes.back().motif.active_properties.setLastMotifPos();

		//////////
		/// Finishing: Add all deletion sites to the correct variable_region_list member site.
		//////////
		for (list<varSite*>::iterator it2 = my_node->variable_region_list.begin(); 
									  it2 != my_node->variable_region_list.end(); 
									  ++it2) 
		{
			for (vector<Site>::iterator site_it = evolutionaryAttributes.begin(); 
									  	site_it != evolutionaryAttributes.end(); 
									  	++site_it) 
			{
				// Each site should belong to 2 varSites: one for template, one for motif.
				if ( (*site_it).motif.active_properties.indel->del->my_sequence_template_varSite == (*it2)) 
					(*it2)->add2Member((*site_it).motif.active_properties.indel->del);
				if ( (*site_it).motif.active_properties.indel->del->my_motif_varSite == (*it2))
					(*it2)->add2Member((*site_it).motif.active_properties.indel->del);
			}
		}
	} else {
		//////////
		/// Insertion site. Need to set all sites to the motif sites of the (prev?) and next sites.
		//////////
		//////////
		/// Even though some of these will be reset later on, set all varSite pointers.
		//////////
		for (vector<Site>::iterator it = evolutionaryAttributes.begin(); 
									it != evolutionaryAttributes.end(); 
									++it) 
		{
			(*it).motif.active_properties.subst = new Substitution();
			(*it).motif.active_properties.indel = new Indel();
			(*it).motif.active_properties.indel->createObjects();

			if (it != evolutionaryAttributes.begin()) {
				delete (*it).motif.active_properties.indel->L_ins_;
				(*it).motif.active_properties.indel->L_ins_
				= (*(it-1)).motif.active_properties.indel->R_ins_;
			}

			varSite *mstv = NULL, *mmv = NULL;
			for (list<siteProperties*>::iterator jt  = (*it).motif.site_props.begin();
												 jt != (*it).motif.site_props.end();
												 ++jt) 
			{
				if ((*jt)->indel.del->my_sequence_template_varSite) {
					mstv = (*jt)->indel.del->my_sequence_template_varSite;
				}
				if ((*jt)->indel.del->my_motif_varSite) {
					mmv  = (*jt)->indel.del->my_motif_varSite;
				}
			}

			(*it).motif.active_properties.indel->del->my_motif_varSite 
			= (*it).motif.active_properties.indel->R_ins_->my_motif_varSite_left 
			= (*it).motif.active_properties.indel->R_ins_->my_motif_varSite_right 
			= mmv;
			(*it).motif.active_properties.indel->del->my_sequence_template_varSite 
			= (*it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left 
			= (*it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right 
			= mstv;

			//////////
			/// Take care of on_site pointers.
			//////////
			if (mmv->min == 0) 
				(*it).motif.active_properties.indel->del->my_motif_varSite = mmv;
			if (mstv->min == 0) 
				(*it).motif.active_properties.indel->del->my_sequence_template_varSite = mstv;

			for (list<varSite*>::iterator it2 =  my_node->variable_region_list.begin(); 
										  it2 != my_node->variable_region_list.end(); 
									  	  ++it2) 
			{
				if ( (*it).motif.active_properties.indel->del->my_sequence_template_varSite == (*it2)) {
					cerr << "Sequence::setActiveProps: Adding new site to site with membership: ";
				 	cerr << "st("
			 		 << (*it).motif.active_properties.indel->del->my_sequence_template_varSite->min 
			 		 << ","
			 		 << (*it).motif.active_properties.indel->del->my_sequence_template_varSite->max
			 		 << ")["
					 << (*it).motif.active_properties.indel->del->my_sequence_template_varSite->member_set.size()
					 << "] ";
					if ((*it).motif.active_properties.indel->del->my_sequence_template_varSite->member_set.size()+1 > (*it).motif.active_properties.indel->del->my_sequence_template_varSite->max) {
						(*it2)->add2Member((*it).motif.active_properties.indel->del);
					} else {
						(*it2)->add2Member((*it).motif.active_properties.indel->del);
					}
					cerr << "NOW: ";
				 	cerr << "st("
			 		 << (*it).motif.active_properties.indel->del->my_sequence_template_varSite->min 
			 		 << ","
			 		 << (*it).motif.active_properties.indel->del->my_sequence_template_varSite->max
			 		 << ")["
					 << (*it).motif.active_properties.indel->del->my_sequence_template_varSite->member_set.size()
					 << "] ";
				}
				if ( (*it).motif.active_properties.indel->del->my_motif_varSite == (*it2))
					(*it2)->add2Member((*it).motif.active_properties.indel->del);
			}
		}
	}
}

double
Sequence::forward_rate_away_from_sequence (
										   Branch *branch,
										   int event_site,
										   int end_site
										  )
{
	double forward_sum_away = 0;

	for (vector<Site>::iterator it = evolutionaryAttributes.begin()+event_site; 
		 it != evolutionaryAttributes.begin()+end_site && it != evolutionaryAttributes.end(); 
		 ++it
		) 
	{
		forward_sum_away += (*it).forward_rate_away_from_site(branch);
	}

	Qidot_k__T__ = -1;	// Forward simulation doesn't do this.

	return forward_sum_away;
}


////////////////////
////// LIKELIHOOD
////////////////////
void Likelihood::calculateStateLikelihood (
									  	   TNode *node,
									  	   int category,		/// For pruning to infer rate category of site
									  	   int position,
									  	   double branch_length_scalar
							   			  )
{
	size_t xi, xj;
	vector<double> Li_branch1 (numStates,0);
	vector<double> Li_branch2 (numStates,0);
	vector<double>::iterator it, jt, my;
	TNode *which_bifurcation;
	
	//////////
	/// BRANCH 1
	//////////
	which_bifurcation = node->branch1;

//	node->branch->rates->setPij(
//								node->seq_evo.at(position), 
//								node->branch->length1 * branch_length_scalar, 
//								node->nodeEnv->rateHetero
//							   );
	//
	// Cycle over all values of current sequence
	xi = 0;
	for (it = Li_xi_.begin(); it != Li_xi_.end(); ++it, xi++) {
		xj = 0;
		// Then over all values of next sequence
		for (jt = which_bifurcation->seq_evo.at(position).L_i.Li_xi_.begin(); jt != which_bifurcation->seq_evo.at(position).L_i.Li_xi_.end(); ++jt, xj++) {
			Li_branch1.at(xi)
//			+= node->branch->rates->Pij.at(category).at(from(xi)+to(xj)) 
			+= which_bifurcation->branch->rates->Pij.at(category).at(from(xi)+to(xj)) 
			   * (*jt);
		}

		/// Set the 1st term of the likelihood eq.
		(*it) = Li_branch1.at(xi);
	}

	//////////
	/// BRANCH 2
	//////////
	which_bifurcation = node->branch2;
//	node->branch->rates->setPij(
//								node->seq_evo.at(position), 
//								node->branch->length2 * branch_length_scalar, 
//								node->nodeEnv->rateHetero
//							   );
	//
	// Cycle over all values of current sequence
	xi = 0;
	for (it = Li_xi_.begin(); it != Li_xi_.end(); ++it, xi++) {
		xj = 0;
		// Then over all values of next sequence
		for (jt = which_bifurcation->seq_evo.at(position).L_i.Li_xi_.begin(); jt != which_bifurcation->seq_evo.at(position).L_i.Li_xi_.end(); ++jt, xj++) {
			Li_branch2.at(xi)
//			+= node->branch->rates->Pij.at(category).at(from(xi)+to(xj)) 
			+= which_bifurcation->branch->rates->Pij.at(category).at(from(xi)+to(xj)) 
			   * (*jt);
		}

		/// Multiply the first term (set in `BRANCH 1' loop) by the likelihood eq. 2nd term.
		(*it) *= Li_branch2.at(xi);
	}
}

vector<double> Likelihood::calculateCatLikelihood (
									  			   TNode *node,
									  			   int category,		/// For pruning to infer rate category of site
									  			   int position,
									  			   vector<double> b1_Li,
									  			   vector<double> b2_Li
							   				      )
{
	size_t xi, xj;
	vector<double> Li_branch1 (numStates,0);
	vector<double> Li_branch2 (numStates,0);
	vector<double>::iterator jt, my;
	vector<double> my_Li (numStates,0);

	//////////
	// Calculate the likelihood of the observed data given the category and the model:
	// Calculate P( X_1,X_2,X_3 | C=j,\theta )
	//             = ·_i P( X_1,X_2,X_3,X_R=i | C=j,\theta )
	// 			   = ·_i P( X_1,X_2,X_3 | X_R=i,C=j,\theta ) x P( X_R=i | C=j,\theta )
	// 			   = ·_i P( X_1,X_2,X_3 | X_R=i,C=j,\theta ) x ¹_i
	//
	// BRANCH 1
	//
//cerr << "Node: " << node->printBipartition() << " BL: " << node->branch->length1 << " " << node->branch->length2 << endl;
//cerr << "Before (branch1): " << endl;
//cerr << "node->branch->rates->setPij(," << node->branch->length1*1 << ",)" << endl;
//node->branch1->branch->rates->printPij();
//	node->branch->rates->setPij(
//	node->branch1->branch->rates->setPij(
//								node->seq_evo.at(position), 
//								node->branch->length1 * 1, 
//								node->nodeEnv->rateHetero
//							   );
//cerr << "After (branch1): " << endl;
//node->branch1->branch->rates->printPij();
	xi = 0;
	for (my = my_Li.begin(); my != my_Li.end(); xi++, my++) {
		xj = 0;
		// Then over all values of next sequence
		for (jt = b1_Li.begin(); jt != b1_Li.end(); ++jt, xj++) {
			Li_branch1.at(xi)
//			+= node->branch->rates->Pij.at(category).at(from(xi)+to(xj)) 
			+= node->branch1->branch->rates->Pij.at(category).at(from(xi)+to(xj)) 
			   * (*jt);
		}

		/// Set the 1st term of the likelihood eq.
		(*my) = Li_branch1.at(xi);
	}
	//
	// BRANCH 2
	//
//cerr << "Before (branch2): " << endl;
//node->branch2->branch->rates->printPij();
//	node->branch->rates->setPij(
//	node->branch2->branch->rates->setPij(
//								node->seq_evo.at(position), 
//								node->branch->length2 * 1, 
//								node->nodeEnv->rateHetero
//							   );
//cerr << "After (branch2): " << endl;
//node->branch->rates->printPij();
//node->branch2->branch->rates->printPij();
//cerr << endl << endl;
	xi = 0;
	for (my = my_Li.begin(); my != my_Li.end(); xi++, my++) {
		xj = 0;
		// Then over all values of next sequence
		for (jt = b2_Li.begin(); jt != b2_Li.end(); ++jt, xj++) {
			Li_branch2.at(xi)
//			+= node->branch->rates->Pij.at(category).at(from(xi)+to(xj)) 
			+= node->branch2->branch->rates->Pij.at(category).at(from(xi)+to(xj)) 
			   * (*jt);
		}

		/// Multiply the first term (set in `BRANCH 1' loop) by the likelihood eq. 2nd term.
		(*my) *= Li_branch2.at(xi);
	}

	//cerr << "          Category: " << category << "    " << endl;
	//cerr << "b1:  ";
	//for(jt = b1_Li.begin(); jt != b1_Li.end(); ++jt)
	//	cerr << "  " << (*jt);
	//cerr << endl;
	//cerr << "b2:  ";
	//for(jt = b2_Li.begin(); jt != b2_Li.end(); ++jt)
	//	cerr << "  " << (*jt);
	//cerr << endl << "my:  ";
	//for (my = my_Li.begin(); my != my_Li.end(); my++) {
	//	cerr << (*my) << " "; 
	//}
	//cerr << endl;

	return my_Li;
}


////////////////////
////// SITE
////////////////////
void Site::InitializeMotif(
						   motifSite *site_type, 
						   varSite *template_varSite, 
						   varSite *motif_varSite
						  ) 
{
	if (template_varSite == NULL && motif_varSite == NULL) {
		motif = motifSite();
		motif.site_props.clear();
		motif.site_props.push_back(new siteProperties());
		motif.site_props.front()->indel.createObjects();
	} else {
		motif = motifSite();
		motif.site_props.push_back(new siteProperties());
		motif.site_props.back()->indel.createObjects();
		if (template_varSite != NULL)
			motif.site_props.push_back(new siteProperties(site_type, template_varSite, NULL));
		if (motif_varSite != NULL) 
			motif.site_props.push_back(new siteProperties(site_type, NULL            , motif_varSite));

		motif.active_properties.fromMotif = site_type->active_properties.fromMotif;
		motif.active_properties.fromTemplate = site_type->active_properties.fromTemplate;
	}
}

void Site::setSiteRateAway (
							double site_width,
							double gamma,
							RateMatrix *rates
						   )
{
	int base = 0;
	double sum_rate_away = 0;

	site_rate_away.assign(numStates,0);
	if (e_QijDt.empty()) e_QijDt.assign(numStates_squared,0);

	//////////
	/// Set uniformization bin:
	//////////
	//  For each i­j, Pij = (Qij*gamma)/bin_max
	//  For i=j       Pii = 1-\sum_{i­j}Pij
	// where bin_max = max{Qii}*max{gamma_j}*BUFFER, where gamma_j is the maximum gamma value over
	// the entire sequence.
	sum_rate_away = base = 0;
	for (vector<double>::iterator it = site_rate_away.begin(); it != site_rate_away.end(); ++it, base++)
		if (state != base) {
			(*it) = (rates->Qij[state*numStates+base] * gamma) / site_width;
			sum_rate_away += (*it);
		}
	site_rate_away.at(state) = 1 - sum_rate_away;
	//
	//////////

	//////////
	/// Have Qij values for each individual transition. Make bin additive.
	//////////
	for (vector<double>::iterator it = site_rate_away.begin()+1; it != site_rate_away.end(); ++it)
		(*it) += (*(it-1));
}

double Site::setSiteRateAway(
							 vector<double>& Qij,
							 RateMatrix *rates
						    ) 
{
	short base = 0;
	double Pij, prev_Pij = 0;

	//////////
	/// Initialize the site_rate_away.
	/// * Assign drops all elements previously assigned to site_rate_away, and replaces it with
	///   the specified elements.
	//////////
	site_rate_away.assign(numStates,0);
	if (e_QijDt.empty()) e_QijDt.assign(numStates_squared,0);

	//////////
	/// Calculate the rate away for a site.
	//////////
	for (vector<double>::iterator it = site_rate_away.begin(); it != site_rate_away.end(); ++it, base++) {
		if (state == base) Pij = 0;
		else Pij = Qij[state*numStates+base];
		(*it) = prev_Pij + Pij * rates->catRate.at(returnCategory());
		prev_Pij = (*it);
	}
	//
	//////////

	//cerr << "Category rate: " << rates->catRate.at(returnCategory()) << endl;
	//for (vector<double>::iterator it = site_rate_away.begin(); it != site_rate_away.end(); ++it)
	//	cerr << "  " << (*it);
	//cerr << "    " << site_rate_away.back() << endl;

	return site_rate_away.back();
}

bool Site::doSubstitution (
						   double value, 
						   int step_type, 
						   TNode *node, 
						   int *codon_position, 
						   string& event
						  )
{
	int base = 0, newState = -1;
	bool set_state = false;
	char from = stateCharacters.at(state), to;

	if (motif.active_properties.subst->siteInvariable()) return false; // Substituion failed.
	for (vector<double>::iterator it = site_rate_away.begin(); it != site_rate_away.end(); ++it, base++) {
		if (value < (*it)) {
			newState = base;
			break;
		}
	}

	if (newState == -1) {																					//XOUT
		cerr << "Site::doSubstitution: Failed to find acceptable state! " << endl;							//XOUT
		cerr << "  key: " << value << endl << "  rate away from site: " << site_rate_away.back() << endl;	//XOUT
		exit(EXIT_FAILURE);																					//XOUT
	}																										//XOUT
	
	//////////
	/// Make sure that site can accept state.
	//////////
	if ( !(motif.active_properties.subst->substitution_bitstring.test(state)) )
		set_state = true;
	else if (motif.active_properties.subst->substitution_bitstring.test(newState) )
		set_state = true;
	else return false;
	
	if (set_state) {
		if (step_type == UNIFORMIZATION) {
			if (newState == state) { 
				virtual_events++; 
				acagatcgctgt.at(state*4+state)++;
				virtualACGT.at(state)++;
				return false;
			} else {
				//////////
				/// If codon rates, then set the position. Decision as to whether the position is
				/// a stop codon or not is done in substitute.
				//////////
				acagatcgctgt.at(state*4+newState)++;
				chosenACGT.at(state)++;
				if (node->nodeEnv->rateHetero == CodonRates) *codon_position = newState;
				else state = newState;
				choseACGT.at(state)++;
				to = stateCharacters.at(newState);
				event.push_back(from); event.push_back(to);
				actual_events++;
				setSiteRateAway(
								node->rate_away_site_width,
								(
								 (gammaRates)
								 ? gammaRates
								 : node->branch->rates->catRate.at(returnCategory())
								),
								node->branch->rates
							   );
				acgt.at(newState)++;
				numSubstAway++;
			}
		} else {
			if (node->nodeEnv->rateHetero == CodonRates) *codon_position = newState;
			else {
				state = newState;
			}
			to = stateCharacters.at(newState);
			event.push_back(from); 
			event.push_back(to);
			numSubstAway++;
		}
	} else return false;
	return true;
}

bool siteDependencies::isDependentInteraction()
{
	return 0;
}

void Site::printSiteRateAway()
{
	size_t base = 0;
	for (vector<double>::iterator it = site_rate_away.begin(); it != site_rate_away.end(); ++it, base++) {
		cerr << stateCharacters.at(base) << "->";
		if (base == state) cerr << "v";
		cerr << (*it) << " ";
	}
}

bool Site::isDefinedByUser ( void )
{
	size_t newState = 0;
	//////////
	/// If this is defined by user at input, then there will be a 1 in one location, with all else
	/// zeroes.
	//////////
	// One issue that needs to be taken care of is: If Category 0 matrix is I, then the site will
	// appear to be defined by the user, and will act accordingly. This is an issue in calculate-
	// RateAwayFromEndpointSequence(), since the state will be set as a -1. The purpose of newState
	// is to do exactly that.
	//////////
	bool a_one = false;
	for (vector<double>::iterator it = L_i.Li_xi_.begin(); it != L_i.Li_xi_.end(); ++it, newState++) {
		if ( (*it) == 1 ) {
			a_one = true;
			//setState(newState);
		} else if ( (*it) != 0 ) return false;
	}
	if ( a_one ) return true;
	else {
		cerr << "Site::isDefinedByUser()" << endl
			 << "  sum (L_i.Li_xi_) is either greater than 1 or something really strange is going on."
			 << endl;
		for (vector<double>::iterator it = L_i.Li_xi_.begin(); it != L_i.Li_xi_.end(); ++it, newState++) {
			cerr << (*it) << "  ";
		}
		cerr << endl;
		exit(EXIT_FAILURE);
	}
}

double
Site::forward_rate_away_from_site(
								  Branch *branch
						   		 )
{
	short i;
	bool profile=false;
	double Qij;
	double forward_sum_away = 0;

	i = returnState();
	forward_rate_away.assign(numStates, 0);
	for (int j = 0; j < numStates; j++) {
		Qij = 0;
		if (order_3_markov) Qij = site_rate_away.at(j); 
		else Qij = branch->rates->Qij.at(from(i)+to(j));
		Qij *= branch->rates->catRate.at(returnCategory());
		if ( j != i ) {
			forward_rate_away.at(j) = Qij;
			forward_sum_away += Qij;
		} else forward_rate_away.at(j) = 0;
	}

	if (order_3_markov || Human_Data_simulation) {
		for (vector<double>::iterator jt = site_rate_away.begin()+1; jt != site_rate_away.end(); ++jt) 
			(*jt) += (*(jt-1));
	} else {
		vector<double>::iterator jt, it = forward_rate_away.begin();
		site_rate_away.assign(numStates, 0);
		it = forward_rate_away.begin();
		jt = site_rate_away.begin();
		(*jt) = (*it);
		++jt; ++it;
		for (; jt != site_rate_away.end(); ++jt, ++it) (*jt) = (*(jt-1))+(*it);
	}
	
	return forward_sum_away;
}


////////////////////
////// INMOTIF
////////////////////
bool inMotif::isTemplate()
{
	if (isdigit(marker)) return true;
	else return false;
}

void inMotif::report() 
{
	cerr << "  name =        " << name << endl;
	cerr << "  marker =      " << marker << endl;
	cerr << "  bipartition = "; 
	for(vector<bool>::iterator it = bipartition.begin(); it != bipartition.end(); ++it)
		cerr << (*it);
	cerr << endl;
	cerr << "  regex =       " << regex << endl;
	cerr << "  sitemap =     " << sitemap << endl;
}

vector<siteProperties*> inMotif::enumerateRegEx(
												TNode *node
											   )
{
	size_t found, found2, found3;
	string site;
	list<string> return_sites;
	list<string> regex_sites;
	list<string> tmp;
	string my_stateCharacters;
	string reverse_sites;
	string x_length;
	int motif_min, motif_max;
	vector<siteProperties*> new_site_props;
	siteProperties *newSite;
	varSite *newVar;
	bool isTemplate;

	regex_sites = split(regex, "-");
	my_stateCharacters.clear();
	for (int i = 0; i < numStates; i++) my_stateCharacters += stateCharacters.at(i);
	return_sites.clear();

	for (list<string>::iterator it = regex_sites.begin(); it != regex_sites.end(); ++it) {
		site.clear();
		switch ((*it).at(0)) {
		case '[':
			found = (*it).find("]", 1);
			if (found == string::npos) {
				cerr << "Regex position starts with '[', but has no ending ']'." << endl;
				exit(EXIT_FAILURE);
			}

			// Check characters in brackets. //
			site = (*it).substr(1,found-1);
			found = site.find_first_not_of(my_stateCharacters);
			if (found != string::npos) {
				cerr << "Unrecognized character in motif regular expression: " << (*it) << endl;
				exit(EXIT_FAILURE);
			}

			// Check to see if it is a multi-site, e.g., [LIV](3) //
			found2 = (*it).find("(", found);
			if (found2 != string::npos) {
				found3 = (*it).find(")", found2+1);
				x_length = (*it).substr(found2+1,found3-found2-1);
				int xlen = atoi(x_length.c_str());

				newVar = new varSite(motif_min, motif_max, node);

				for (int i = 0; i < xlen; i++) {
					newSite = new siteProperties(new_site_props, newVar, this, sequence_template);
					new_site_props.push_back(newSite);
					return_sites.push_back(site);
				}
			} else {
				newSite = new siteProperties(new_site_props, node->Generic_varSite(), this, sequence_template);
				new_site_props.push_back(newSite);
				return_sites.push_back(site);
			}
			break;
		case '{':
			found = (*it).find("}", 1);
			if (found == string::npos) {
				cerr << "Regex position starts with '{', but has no ending '}'." << endl;
				exit(EXIT_FAILURE);
			}
			reverse_sites = (*it).substr(1, found-1);
			found = reverse_sites.find_first_not_of(my_stateCharacters);
			if (found != string::npos) {
				cerr << "Unrecognized character in motif regular expression: " << (*it) << endl;
				exit(EXIT_FAILURE);
			}
			
			for (string::iterator it2 = my_stateCharacters.begin(); it2 != my_stateCharacters.end(); ++it2) {
				found = reverse_sites.find((*it2));
				if (found == string::npos) site += (*it2);
			}

			found2 = (*it).find("(", found);
			if (found2 != string::npos) {
				found3 = (*it).find(")", found2+1);
				x_length = (*it).substr(found2+1,found3-found2-1);
				int xlen = atoi(x_length.c_str());

				newVar = new varSite(motif_min, motif_max, node);

				for (int i = 0; i < xlen; i++) {
					newSite = new siteProperties(new_site_props, newVar, this, sequence_template);
					new_site_props.push_back(newSite);
					return_sites.push_back(site);
				}
			} else {
				newSite = new siteProperties(new_site_props, node->Generic_varSite(), this, sequence_template);
				new_site_props.push_back(newSite);
				return_sites.push_back(site);
			}
			break;
		case 'x':
			if ((*it).size() == 1) {
				newSite = new siteProperties(new_site_props, node->Generic_varSite(), this, sequence_template);
				new_site_props.push_back(newSite);
				site = my_stateCharacters;	// just 'x'
				return_sites.push_back(my_stateCharacters);
			} else {
				tmp = split((*it), ",");
				if (tmp.size() > 1) {
					tmp.front().erase(0,2);
					tmp.back().erase(tmp.back().size()-1, 1);
					motif_min = atoi(tmp.front().c_str());
					motif_max = atoi(tmp.back().c_str());

					newVar = new varSite(motif_min, motif_max, node);

					// Fill out motif_min number of sites.
					for (size_t i = 0; i < motif_max; i++) {
						newSite = new siteProperties(new_site_props, newVar, this, sequence_template);
						new_site_props.push_back(newSite);
						return_sites.push_back(my_stateCharacters);
					}
				} else {
					size_t found2;
					(*it).erase(0,2);
					found2 = (*it).find(")");
					x_length = (*it).substr(0,found2);
					int xlen = atoi(x_length.c_str());

					newVar = new varSite(xlen, xlen, node);

					for (int i = 0; i < xlen; i++) {
						newSite = new siteProperties(new_site_props, newVar, this, sequence_template);
						new_site_props.push_back(newSite);
						return_sites.push_back(my_stateCharacters);
					}
				}
			}
			break;
		default:
			found = my_stateCharacters.find((*it).at(0));
			if (found != string::npos) {
				newSite = new siteProperties(new_site_props, node->Generic_varSite(), this, sequence_template);
				new_site_props.push_back(newSite);
				site.push_back((*it).at(0));
				return_sites.push_back(site);
			} else {
				cerr << "Unrecognized format for regular expression: " << (*it) << endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
	}

	// Set bit string.
	vector<siteProperties*>::iterator it = new_site_props.begin();
	for (list<string>::iterator it3 = return_sites.begin(); 
								it3 != return_sites.end(); 
								++it3, ++it) 
	{
		bitset<20> regex_site;
		regex_site.set();
		for (string::iterator it4 = (*it3).begin(); it4 != (*it3).end(); it4++) {
			found = my_stateCharacters.find((*it4));
			if (found != string::npos) regex_site.flip(found);
		}
		regex_site.flip();
		(*it)->subst.setSiteBits(regex_site);
	}

	return new_site_props;
}

void inMotif::removeLastPosition() 
{
	// If last position is optional (e.g., [G>] in motif), and it is excluded from the motif,
	// this routine removes it from the inMotif prosite_library so that the regex size will
	// agree with the regex placed.
	list<string> regex_sites;
	string fixed_regex = "";
	regex_sites = split(regex, "-");
	regex_sites.pop_back();
	for (list<string>::iterator it = regex_sites.begin(); it != regex_sites.end(); ++it) {
		fixed_regex += (*it);
		fixed_regex += "-";
	}
	fixed_regex.erase(fixed_regex.find_last_of("-"),1);
	regex = fixed_regex;
}

list<siteRegEx*> inMotif::parseRegEx()
{
	size_t found, found2, found3;
	string site;
	list<string> return_sites;
	list<string> regex_sites;
	list<string> tmp;
	string my_stateCharacters;
	string reverse_sites;
	string x_length;
	int motif_min, motif_max;
	bool isTemplate;
	list<siteRegEx*> return_list;
	siteRegEx *newSite;
	bool N_term = false, C_term = false, last_site_optional = false;

	found = regex.find("<", 0);
	if (found != string::npos) { 
		N_term = true; 
		regex.erase(found, 1); 
	}
	found = regex.find(">", 0);
	if (found != string::npos) {
		C_term = true;
		if (found != regex.size()-1) // ">" comes inside of [], e.g. [G>], meaning G is optional. //
			last_site_optional = true;
		regex.erase(found, 1);		
	}
	
	regex_sites = split(regex, "-");
	my_stateCharacters.clear();
	for (int i = 0; i < numStates; i++) my_stateCharacters += stateCharacters.at(i);
	return_sites.clear();
	return_list.clear();

	for (list<string>::iterator it = regex_sites.begin(); it != regex_sites.end(); ++it) {
		site.clear();
		switch ((*it).at(0)) {
		case '[':
			newSite = new siteRegEx(*it, 1, N_term, C_term, last_site_optional);
			found = (*it).find("]", 1);
			if (found == string::npos) {
				cerr << "Regex position starts with '[', but has no ending ']'." << endl;
				exit(EXIT_FAILURE);
			}

			// Check characters in brackets. //
			site = (*it).substr(1,found-1);
			found = site.find_first_not_of(my_stateCharacters);
			if (found != string::npos) {
				cerr << "Unrecognized character in motif regular expression: " << (*it) << endl;
				exit(EXIT_FAILURE);
			}

			// Check to see if it is a multi-site, e.g., [LIV](3) //
			found2 = (*it).find("(", found);
			if (found2 != string::npos) {
				found3 = (*it).find(")", found2+1);
				x_length = (*it).substr(found2+1,found3-found2-1);
				int xlen = atoi(x_length.c_str());

				for (int i = 0; i < xlen; i++) newSite->allowable_characters.push_back(site);
				newSite->sites_occupied = xlen;
			} else {
				newSite->allowable_characters.push_back(site);
			}
			return_list.push_back(newSite);
			break;
		case '{':
			newSite = new siteRegEx(*it, 1, N_term, C_term, last_site_optional);
			found = (*it).find("}", 1);
			if (found == string::npos) {
				cerr << "Regex position starts with '{', but has no ending '}'." << endl;
				exit(EXIT_FAILURE);
			}
			reverse_sites = (*it).substr(1, found-1);
			found = reverse_sites.find_first_not_of(my_stateCharacters);
			if (found != string::npos) {
				cerr << "Unrecognized character in motif regular expression: " << (*it) << endl;
				exit(EXIT_FAILURE);
			}
			
			for (string::iterator it2 = my_stateCharacters.begin(); it2 != my_stateCharacters.end(); ++it2) {
				found = reverse_sites.find((*it2));
				if (found == string::npos) site += (*it2);
			}

			found2 = (*it).find("(", found);
			if (found2 != string::npos) {
				found3 = (*it).find(")", found2+1);
				x_length = (*it).substr(found2+1,found3-found2-1);
				int xlen = atoi(x_length.c_str());

				for (int i = 0; i < xlen; i++) newSite->allowable_characters.push_back(site);
				newSite->sites_occupied = xlen;
			} else {
				newSite->allowable_characters.push_back(site);
			}
			return_list.push_back(newSite);
			break;
		case 'x':
			newSite = new siteRegEx(*it, 1, N_term, C_term, last_site_optional);
			if ((*it).size() == 1) {
				site = my_stateCharacters;	// just 'x'
				newSite->allowable_characters.push_back(site);
			} else {
				tmp = split((*it), ",");
				if (tmp.size() > 1) {
					tmp.front().erase(0,2);
					tmp.back().erase(tmp.back().size()-1, 1);
					motif_min = atoi(tmp.front().c_str());
					motif_max = atoi(tmp.back().c_str());

					newSite->sites_occupied = (double)rndu() * (motif_max-motif_min) + motif_min;
					// Fill out motif_max number of sites.
					for (size_t i = 0; i < motif_max; i++) newSite->allowable_characters.push_back(my_stateCharacters);
				} else {
					size_t found2;
					(*it).erase(0,2);
					found2 = (*it).find(")");
					x_length = (*it).substr(0,found2);
					int xlen = atoi(x_length.c_str());

					for (int i = 0; i < xlen; i++) newSite->allowable_characters.push_back(my_stateCharacters);
					newSite->sites_occupied = xlen;
				}
			}
			return_list.push_back(newSite);
			break;
		default:
			newSite = new siteRegEx(*it, 1, N_term, C_term, last_site_optional);
			found = my_stateCharacters.find((*it).at(0));
			if (found != string::npos) {
				site.push_back((*it).at(0));
				newSite->allowable_characters.push_back(site);
			} else {
				cerr << "Unrecognized format for regular expression: " << (*it) << endl;
				exit(EXIT_FAILURE);
			}
			return_list.push_back(newSite);
			break;
		
		}
	}

	// Set bit string.
	for (list<string>::iterator it3 = return_sites.begin(); it3 != return_sites.end(); ++it3) {
		bitset<20> regex_site;
		regex_site.set();
		for (string::iterator it4 = (*it3).begin(); it4 != (*it3).end(); it4++) {
			found = my_stateCharacters.find((*it4));
			if (found != string::npos) regex_site.flip(found);
		}
		regex_site.flip();
	}

	return return_list;
}

bitset<20> inMotif::getRegExValues(
								   size_t which_site
								  )
{
	list<string> regex_sites;
	string my_stateCharacters;
	size_t at_current_site = 0;
	bitset<20> return_values;
	string site, reverse_sites;
	size_t found;
	list<string> tmp;
	string x_length;
	int motif_min, motif_max;
	bool done = false;

	regex_sites = split(regex, "-");
	my_stateCharacters.clear();
	for (int i = 0; i < numStates; i++) my_stateCharacters += stateCharacters.at(i);
	return_values.reset();

	for (list<string>::iterator it = regex_sites.begin(); it != regex_sites.end() && !done; ++it) {
		site.clear();
		switch ((*it).at(0)) {
		case '[':
			at_current_site++;
			if (at_current_site == which_site) {

				found = (*it).find("]", 1);
				if (found == string::npos) {
					cerr << "Regex position starts with '[', but has no ending ']'." << endl;
					exit(EXIT_FAILURE);
				}
				site = (*it).substr(1,found-1);
				found = site.find_first_not_of(my_stateCharacters);
				if (found != string::npos) {
					cerr << "Unrecognized character in motif regular expression: " << (*it) << endl;
					exit(EXIT_FAILURE);
				} else {
					// Set return_values
					for (string::iterator sit = site.begin(); sit != site.end(); ++sit) {
						found = my_stateCharacters.find_first_of(*sit);
						return_values.set(found);
					}
				}
				done = true;
			}
			break;
		case '{':
			at_current_site++;
			if (at_current_site == which_site) {
				found = (*it).find("}", 1);
				if (found == string::npos) {
					cerr << "Regex position starts with '{', but has no ending '}'." << endl;
					exit(EXIT_FAILURE);
				}
				reverse_sites = (*it).substr(1, found-1);
				found = reverse_sites.find_first_not_of(my_stateCharacters);
				if (found != string::npos) {
					cerr << "Unrecognized character in motif regular expression: " << (*it) << endl;
					exit(EXIT_FAILURE);
				}
				for (string::iterator it2 = my_stateCharacters.begin(); it2 != my_stateCharacters.end(); ++it2) {
					found = reverse_sites.find((*it2));
					if (found == string::npos) site += (*it2);
				}
				for (string::iterator sit = site.begin(); sit != site.end(); ++sit) {
					found = my_stateCharacters.find_first_of(*it);
					return_values.set(found);
				}
				done = true;
			}
			break;
		case 'x':
			if ((*it).size() == 1) {
				at_current_site++;
				if (at_current_site == which_site) {
					return_values.set();
					done = true;
				}
			} else {
				tmp = split((*it), ",");
				if (tmp.size() > 1) {
					tmp.front().erase(0,2);
					tmp.back().erase(tmp.back().size()-1, 1);
					motif_min = atoi(tmp.front().c_str());
					motif_max = atoi(tmp.back().c_str());

					if (motif_min == motif_max) {
						if (at_current_site + motif_min >= which_site) {
							return_values.set();
							done = true;
						} else at_current_site += motif_min;
					} else if (at_current_site + motif_max >= which_site) {
						done = true;
					} else at_current_site += motif_max;
				} else {
					size_t found2;
					(*it).erase(0,2);
					found2 = (*it).find(")");
					x_length = (*it).substr(0,found2);
					int xlen = atoi(x_length.c_str());

					if (at_current_site + xlen >= which_site) {
						return_values.set();
						done = true;
					} else at_current_site += xlen;
				}
			}
			break;
		default:
			at_current_site++;
			if (at_current_site == which_site) {
				found = my_stateCharacters.find((*it).at(0));
				if (found != string::npos) return_values.set(found);
				else {
					cerr << "Unrecognized format for regular expression: " << (*it) << endl;
					exit(EXIT_FAILURE);
				}
				done = true;
			}
			break;
		}
	}
	
	return return_values;
}


////////////////////
////// MOTIFSITE
////////////////////
motifSite::motifSite()
{
	site_props.clear();
}

void motifSite::copy(
					 motifSite* site2copy
					) 
{
	site_props.clear();
	for (list<siteProperties*>::iterator it = site2copy->site_props.begin(); it != site2copy->site_props.end(); ++it) {
		site_props.push_back(new siteProperties(*it,this));
		if (site_props.back()->indel.del == NULL)
			site_props.back()->indel.createObjects();
	}
}

void motifSite::setNode(
					    TNode *node
					   )
{
	my_node = node;
}

void motifSite::Site_deleteMerge(
								 motifSite *prev, 
								 motifSite *next, 
								 bool deletion
								)
{
	///////////
	/// Briefly, variables are as such (with respect to this):
	///  prev      this      next
	///      \    /    \    /
	///       \  /      \  /
	///        1 2       3 4
	///       L_ins_    R_ins_
	///////////

	if (next != NULL && prev != NULL) {
		next->active_properties.indel->L_ins_ 
		= prev->active_properties.indel->R_ins_;
		next->active_properties.indel->L_ins_->my_sequence_template_varSite_right 
		= next->active_properties.indel->del->my_sequence_template_varSite;
		next->active_properties.indel->L_ins_->my_motif_varSite_right 
		= next->active_properties.indel->del->my_motif_varSite;
	} else {
		if (next == NULL) {
			active_properties.indel->L_ins_->my_sequence_template_varSite_left
			= active_properties.indel->L_ins_->my_sequence_template_varSite_right
			= prev->active_properties.indel->del->my_sequence_template_varSite;
			active_properties.indel->L_ins_->my_motif_varSite_left
			= active_properties.indel->L_ins_->my_motif_varSite_right
			= prev->active_properties.indel->del->my_motif_varSite;
		} else if (prev == NULL) cerr << "front end being deleted." << endl;
	}

	if (active_properties.indel->del->my_sequence_template_varSite != my_node->one_site_varSite
		&& active_properties.indel->del->my_sequence_template_varSite != my_node->unconstrained_varSite)
		active_properties.indel->del->my_sequence_template_varSite->removeFromMember(active_properties.indel->del);
	if (active_properties.indel->del->my_motif_varSite != my_node->one_site_varSite
		&& active_properties.indel->del->my_motif_varSite != my_node->unconstrained_varSite)
		active_properties.indel->del->my_motif_varSite->removeFromMember(active_properties.indel->del);

	deleteSiteObjects(deletion);
}

void motifSite::deleteSiteObjects(
								  bool deletion
								 )
{
	for (list<siteProperties*>::iterator it2 = site_props.begin(); it2 != site_props.end(); ++it2) {
		if (deletion)
			delete (*it2)->indel.R_ins_;
		delete (*it2)->indel.del;
		delete (*it2);
	}

	delete active_properties.subst;
	delete active_properties.indel->R_ins_;
	delete active_properties.indel->del;
	delete active_properties.indel;

	site_props.clear();
}


////////////////////
////// ACTIVEPROPERTIES
////////////////////
void activeProperties::setFirstMotifPos(
									    TNode *node
									   )
{
	//////////
	/// Beginning should be a (1,1) varSite.
	//////////
	indel = new Indel();
	indel->createObjects();
	subst = new Substitution();

	//////////
	/// Need to manually set the first position insertion sites, since the connecting routine
	/// checks only those sites that 
	//////////
	indel->del->my_sequence_template_varSite 
	= indel->del->my_motif_varSite 
	= indel->L_ins_->my_sequence_template_varSite_left
	= indel->L_ins_->my_motif_varSite_left
	= indel->L_ins_->my_sequence_template_varSite_right
	= indel->L_ins_->my_motif_varSite_right
	= node->one_site_varSite;
}

void activeProperties::setLastMotifPos()
{
	indel->R_ins_->my_sequence_template_varSite_left
	= indel->R_ins_->my_sequence_template_varSite_right
	= indel->del->my_sequence_template_varSite;

	indel->R_ins_->my_motif_varSite_right
	= indel->R_ins_->my_motif_varSite_left
	= indel->del->my_motif_varSite;
	
	indel->R_ins_->on_site_sequence_template_varSite
	= ( (indel->del->my_sequence_template_varSite->min == 0)
		? indel->del->my_sequence_template_varSite
		: NULL
	  );

	indel->R_ins_->on_site_motif_varSite
	= ( (indel->del->my_motif_varSite->min == 0)
		? indel->del->my_motif_varSite
		: NULL
	  );
}

void activeProperties::set_properties(
									  list<siteProperties*>& site_props, 
									  list<siteProperties*>& anc_site_props,
									  TNode *node, 
									  activeProperties *prev_site, 
									  size_t site_number
									 )
{
	int relationship;
	bool constrainTemplate=false, constrainMotif=false;
	short template_type = NO_RELATION, motif_type = NO_RELATION;
	siteProperties *template_props = NULL; 
	siteProperties *anc_template_props = NULL;
	siteProperties *motif_props = NULL;
	siteProperties *anc_motif_props = NULL;
	siteProperties *anc_props = NULL;

	//////////
	/// Find the template.
	//////////
	list<siteProperties*>::iterator anc_it = anc_site_props.begin();
	for (list<siteProperties*>::iterator it = site_props.begin(); it != site_props.end(); ++it, ++anc_it) {
		if ( (*it)->fromTemplate != NULL) {	
			relationship = testRelation ((*it)->fromTemplate->bipartition, node->bipartition);
			if (relationship == EXACT) {
				if (constrainTemplate) {
					cerr << "Can only set one template per position per subtree." << endl;
					exit(EXIT_FAILURE);
				}
				constrainTemplate = true;
				template_type = EXACT;
				template_props = (*it);
			} else if (relationship == NODE_IS_ANC) {  // Is a template, but relationship is not exact.
				if (template_type != EXACT && template_type != NODE_IS_DES) {
					constrainTemplate = true;
					template_type = NODE_IS_ANC;
					template_props = (*it);
				}
			} else if (relationship == NODE_IS_DES) {
				if (template_type != EXACT) {
					constrainTemplate = true;
					template_type = NODE_IS_DES;
					template_props = (*it);
					anc_template_props = (*anc_it);
				}
			} else {		// NO_RELATIONSHIP
				template_type = NO_RELATION;
			}
		}
	}

	if (constrainTemplate) {
		setProps(
		  	     TEMPLATE, 
				 template_props,
				 prev_site,
				 template_type,
				 ( (node == node->anc)	/// Root build.
				   ? NULL
				   : &(anc_template_props->indel)
				 )
				);
		fromTemplate = template_props->fromTemplate;
	} else {
		indel->del->my_sequence_template_varSite
		= indel->L_ins_->my_sequence_template_varSite_right
		= node->unconstrained_varSite;
		indel->L_ins_->my_sequence_template_varSite_left
		= prev_site->indel->del->my_sequence_template_varSite;
		indel->L_ins_->on_site_sequence_template_varSite
		= indel->del->my_sequence_template_varSite->set_on_site_ptr(prev_site->indel->del->my_sequence_template_varSite);
		fromTemplate = NULL;
	}

	//////////
	/// Find the motif
	//////////
	list<siteProperties*>::iterator motif_anc_it = anc_site_props.begin();
	for (list<siteProperties*>::iterator it = site_props.begin(); it != site_props.end(); ++it, ++motif_anc_it) {
		if ( (*it)->fromMotif != NULL) {		// is a motif.
			relationship = testRelation ((*it)->fromMotif->bipartition, node->bipartition);
			if (relationship == EXACT) {
				if (constrainMotif) {
					cerr << "Can only set one motif per position per subtree." << endl;
					exit(EXIT_FAILURE);
				}	
				constrainMotif = true;
				motif_type = EXACT;
				motif_props = (*it);
			} else if (relationship == NODE_IS_ANC) { 	// is a motif, relationship is not exact...
				if (motif_type != EXACT && motif_type != NODE_IS_DES) {
					constrainMotif = true;
					motif_type = NODE_IS_ANC;
					motif_props = (*it);
				}
cerr << "NODE_IS_ANC" << endl;
			} else if (relationship == NODE_IS_DES) {
				if (motif_type != EXACT) { 
					constrainMotif = true;
					motif_type = NODE_IS_DES;
					motif_props = (*it);
					anc_motif_props = (*motif_anc_it);
				}
			} else {		// NO_RELATIONSHIP
				motif_type = NO_RELATION;
cerr << "NO_RELATION" << endl;
			}
		}
	}

	if (constrainMotif) {
		setProps(
				 MOTIF,
				 motif_props,
				 prev_site,
				 motif_type,
				 ( (node == node->anc)
				   ? NULL
				   : &(anc_motif_props->indel)
				 )
				);
		fromMotif = motif_props->fromMotif;
	} else {
		indel->del->my_motif_varSite
		= indel->L_ins_->my_motif_varSite_right
		= node->unconstrained_varSite;
		indel->L_ins_->my_motif_varSite_left
		= prev_site->indel->del->my_motif_varSite;
		indel->L_ins_->on_site_motif_varSite
		= indel->del->my_motif_varSite->set_on_site_ptr(prev_site->indel->del->my_motif_varSite);
		fromMotif = NULL;
	}
}

//////////
/// SETPROPS
//////////
// This function sets the properties of the insertion and deletion sites pointers.
// * my_motif_varSitexxx
// * my_sequence_template_varSitexxx
// * on_site_motif_varSitexxx
// * on_site_sequence_template_varSitexxxx
//////////
void activeProperties::setProps(
								short type, 
								siteProperties *props, 
								activeProperties *prev_site, 
								short relation, 
								Indel *anc_indel
							   )
{
	if (type == TEMPLATE) {
		switch (relation) {
		case EXACT:
			indel->set(							///////////
					   &props->indel, 			/// * Indel properties to set.
					   prev_site->indel,		/// * Previous site's indel properties
					   type,					/// * Type of props (template, here)
					   false, 					/// * Is node descendant?
					   ( (anc_indel == NULL) 	/// * Is there an ancestor?
						 ? NULL 				///   - NO: node is root node
						 : anc_indel 			///   - YES: May need to inherit properties.
					   )						///
					  );						//////////
			break;
		case NODE_IS_DES:
			assert( anc_indel );	// Makes sure that anc_indel is not null. If so, aborts.
			indel->set(
					   anc_indel,
					   prev_site->indel,
					   type,
					   true,
					   anc_indel
					  );
			break;
		case NODE_IS_ANC:
			indel->set(							///////////
					   &props->indel, 			/// * Indel properties to set.
					   prev_site->indel,		/// * Previous site's indel properties
					   type,					/// * Type of props (template, here)
					   false, 					/// * Is node descendant?
					   ( (anc_indel == NULL) 	/// * Is there an ancestor?
						 ? NULL 				///   - NO: node is root node
						 : anc_indel 			///   - YES: May need to inherit properties.
					   )						///
					  );						//////////
			cerr << "activeProperties::setProps TEMPLATE NODE_IS_ANC" << endl; exit(0);
			break;
		case NO_RELATION:
			//////////
			/// DO NOTHING
			//////////
			break;
		default:
			cerr << "Illegal relationship: " << relation << endl;
			exit(EXIT_FAILURE);
			break;
		}	
	} else {	//// Motif.
		switch (relation) {
		case EXACT:
			subst->copy(&props->subst);
			indel->set(						///////////
					   &props->indel, 			/// * Indel properties to set.
					   prev_site->indel,		/// * Previous site's indel properties
					   type,					/// * Type of props (motif, here)
					   false, 					/// * Is node descendant?
					   ( (anc_indel == NULL) 	/// * Is there an ancestor?
						 ? NULL 				///   - NO: node is root node
						 : anc_indel 			///   - YES: May need to inherit properties.
					   )						///
					  );						//////////
			break;
		case NODE_IS_DES:
			assert( anc_indel );	// Makes sure that anc_indel is not null. If so, aborts.
			subst->copy(&props->subst);
			indel->set(
					   anc_indel,
					   //&props->indel,
					   prev_site->indel,
					   type,
					   true,
					   anc_indel
					  ); /// Part of problem 2.
			break;
		case NODE_IS_ANC:
			//////////
			/// Virtually same as EXACT, but substitution constraints not yet in effect.
			//////////
			indel->set(						///////////
					   &props->indel, 			/// * Indel properties to set.
					   prev_site->indel,		/// * Previous site's indel properties
					   type,					/// * Type of props (motif, here)
					   false, 					/// * Is node descendant?
					   ( (anc_indel == NULL) 	/// * Is there an ancestor?
						 ? NULL 				///   - NO: node is root node
						 : anc_indel 			///   - YES: May need to inherit properties.
					   )						///
					  );						//////////
			break;
		case NO_RELATION:
			//////////
			/// DO NOTHING
			//////////
			break;
		default:
			cerr << "Illegal relationship: " << relation << endl;
			exit(EXIT_FAILURE);
			break;
		}	
	}
	
}

void activeProperties::report() 
{
	cout << this << ": " << endl;
	
	cout << "  L_ins_ " << indel->L_ins_ << " mstv[LEFT]: ";
	if (indel->L_ins_->my_sequence_template_varSite_left == NULL) cout << "NULL";
	else cout << indel->L_ins_->my_sequence_template_varSite_left->min 
			  << " " << indel->L_ins_->my_sequence_template_varSite_left->max << "  ";
	cout << " mstv[RIGHT]: ";
	if (indel->L_ins_->my_sequence_template_varSite_right == NULL) cout << "NULL";
	else cout << indel->L_ins_->my_sequence_template_varSite_right->min << " " 
			  << indel->L_ins_->my_sequence_template_varSite_right->max << "  ";
	cout << " mmv[LEFT]: ";
	if (indel->L_ins_->my_motif_varSite_left == NULL) cout << "NULL";
	else cout << indel->L_ins_->my_motif_varSite_left->min << " " 
			  << indel->L_ins_->my_motif_varSite_left->max << "  ";
	cout << " mmv[RIGHT]: ";
	if (indel->L_ins_->my_motif_varSite_right == NULL) cout << "NULL";
	else cout << indel->L_ins_->my_motif_varSite_right->min << " " 
			  << indel->L_ins_->my_motif_varSite_right->max << "  ";
	cout << endl; 

	cout << "  R_ins_ " << indel->R_ins_ << " mstv[LEFT]: ";
	if (indel->R_ins_->my_sequence_template_varSite_left == NULL) cout << "NULL";
	else cout << indel->R_ins_->my_sequence_template_varSite_left->min << " " 
			  << indel->R_ins_->my_sequence_template_varSite_left->max << "  ";
	cout << " mstv[RIGHT]: ";
	if (indel->R_ins_->my_sequence_template_varSite_right == NULL) cout << "NULL";
	else cout << indel->R_ins_->my_sequence_template_varSite_right->min << " " 
			  << indel->R_ins_->my_sequence_template_varSite_right->max << "  ";
	cout << " mmv[LEFT]: ";
	if (indel->R_ins_->my_motif_varSite_left == NULL) cout << "NULL";
	else cout << indel->R_ins_->my_motif_varSite_left->min << " " 
			  << indel->R_ins_->my_motif_varSite_left->max << "  ";
	cout << " mmv[RIGHT]: ";
	if (indel->R_ins_->my_motif_varSite_right == NULL) cout << "NULL";
	else cout << indel->R_ins_->my_motif_varSite_right->min << " " 
			  << indel->R_ins_->my_motif_varSite_right->max << "  ";

	cout << endl << "  del " << indel->del << " ";
	if (indel->del->my_sequence_template_varSite == NULL) cout << "NULL";
	else 
		cout << indel->del->my_sequence_template_varSite->min << " " 
			 << indel->del->my_sequence_template_varSite->max << " mmv: ";
	if (indel->del->my_motif_varSite == NULL) cout << "NULL";
	else 
		cout << indel->del->my_motif_varSite->min << " " 
			 << indel->del->my_motif_varSite->max << " mmv: ";
	
	cout << endl << "  ";
	subst->report_bitset();

	cout << endl;
}


////////////////////
////// SITE PROPERTIES
////////////////////
siteProperties::siteProperties(
							   vector<siteProperties*> prev_site, 
							   varSite *memberOfSite, 
							   inMotif *infromMotif, 
							   bool isTemplate
							  ) 
			   : subst (Substitution()),
			     indel(Indel(prev_site)),
				 fromMotif
				 	(
				 	  (
				 	    (isTemplate)
						? NULL
						: infromMotif
					  )
					),
				 fromTemplate
				 	(
				 	  (
				 	    (isTemplate)
						? infromMotif
						: NULL
					  )
					)
{
	// Set up indel sites.
	if (prev_site.empty()) {
		indel.L_ins_->setMembership(memberOfSite, isTemplate, RIGHT);
		indel.R_ins_->setMembership(memberOfSite, isTemplate, LEFT);
	} else {
		indel.R_ins_->setMembership(memberOfSite, isTemplate, LEFT);
		// Also need to set the previous site's membership to this spot.
		prev_site.back()->indel.L_ins_->setMembership(memberOfSite, isTemplate, RIGHT);
	}
	indel.del->setMembership(memberOfSite, isTemplate);
}

void siteProperties::copy(
						  siteProperties *stuff, 
						  short which, 
						  short type, 
						  motifSite *prev, 
						  siteProperties *anc, 
						  bool last_site
						 )
{
	bitset<20> general;
	general.set();
	Indel *toPass = NULL;

	if (anc != NULL) toPass = &anc->indel;

	// Substitutions are set ONLY in motif regions.
	if (which == TEMPLATE) {
		if (type != NO_RELATION) {
			fromTemplate = stuff->fromTemplate;
		}
		if (type == EXACT) {
			indel.copy(&stuff->indel, prev, true /*isTemplate*/, last_site, false/*isDescendant*/, toPass);	// True is the isTemplate
		} else if (type == NODE_IS_DES) {
			if (anc != NULL) {
				indel.copy(&anc->indel, prev, true, last_site, true, toPass);
			} else {
				cerr << "Should this be possible???" << endl;
				exit(EXIT_FAILURE);
			}
		} else if (type == NODE_IS_ANC) {
			indel.copy(&stuff->indel, prev, true, last_site, false, toPass);
		} else {
			// No relationship, do nothing.
		}
	} else if (which == MOTIF) {
		if (type != NO_RELATION) fromMotif = stuff->fromMotif;
		if (type == EXACT) {
			subst.copy(&stuff->subst);
			indel.copy(&stuff->indel, prev, false, last_site, false, toPass);
		} else if (type == NODE_IS_DES) {
			subst.copy(&anc->subst);
			indel.copy(&anc->indel,prev, false,last_site, true, toPass);
		} else if (type == NODE_IS_ANC) {
			subst.setSiteBits(general);
			indel.copy(&stuff->indel, prev, false, last_site, false, toPass);
		} else {
			// No relationship, do nothing
		}
	} else {
		// This is for Pseudogene'd item. If pseudogene'd then we do not want to copy anything
		indel.copy(&stuff->indel, prev, false, last_site, false, NULL);
		subst.setSiteBits(general);
	}
}

////////////////////
////// Indel
////////////////////
Indel::Indel(
			 vector<siteProperties*> prev_site
			) 
	:  L_ins_
		(
		  (
		    (prev_site.empty())
		    ? (new Insertion())
		    : prev_site.back()->indel.R_ins_
		  )
		),
  	   R_ins_ (new Insertion()),
  	   del (new Deletion())
{ }

// Should not auto-create insertion and deletions?
Indel::Indel(
			 varSite *in_template_varSite, 
			 varSite *in_motif_varSite
			)
	:  L_ins_ (new Insertion(in_template_varSite, in_motif_varSite)),
	   R_ins_ (new Insertion(in_template_varSite, in_motif_varSite)),
	   del    (new Deletion(in_template_varSite, in_motif_varSite))
{ }

// same?
Indel::Indel(
			 Indel *site2copy, 
			 motifSite *anc, 
			 bool isTemplate
			)
	:  L_ins_ (new Insertion(site2copy->L_ins_, anc, isTemplate, NULL)),
	   R_ins_ (new Insertion(site2copy->R_ins_, anc, isTemplate, NULL)),
	   del    (new Deletion(site2copy->del, anc, isTemplate))
{ }

Indel::Indel() 
	:  L_ins_ (NULL),
	   R_ins_ (NULL),
	   del (NULL)
{ }

void Indel::createObjects()
{
	L_ins_ = new Insertion();
	R_ins_ = new Insertion();
	del = new Deletion();
}

void Indel::copy(
				 Indel *site2copy, 
				 motifSite *prev_site, 
				 bool isTemplate, 
				 bool last_site, 
				 bool isDescendant, 
				 Indel *anc_indel
				)
{
	Insertion *anc_L_ins_ = NULL, *anc_R_ins_ = NULL;
	Deletion  *anc_del    = NULL;
	
	if (anc_indel != NULL) {
		if (anc_indel->L_ins_ != NULL) anc_L_ins_ = anc_indel->L_ins_;
		if (anc_indel->R_ins_ != NULL) anc_R_ins_ = anc_indel->R_ins_;
		if (anc_indel->del    != NULL) anc_del    = anc_indel->del;
	}

	L_ins_->copy(site2copy->L_ins_, isTemplate, isDescendant, anc_L_ins_);
	if (last_site) {
		R_ins_ = new Insertion();
		R_ins_->copy(site2copy->R_ins_,isTemplate, isDescendant, anc_R_ins_); // B/c otherwise NULL. maybe not smart to do.
		R_ins_->my_sequence_template_varSite_right 
		= R_ins_->my_sequence_template_varSite_left;
		R_ins_->my_motif_varSite_right 
		= R_ins_->my_motif_varSite_left;
	}
	del->copy(site2copy->del, isTemplate, isDescendant, anc_del);
}

void Indel::report()
{
	cout << "indel (des): " << setw(10) << this << endl;
	if (L_ins_ != NULL) {
		cout << "st L_ins_=" << setw(10) << L_ins_; 
		if (L_ins_->my_sequence_template_varSite_left != NULL) {
			cout << setw(10) << L_ins_->my_sequence_template_varSite_left;
			cout << "(" << setw(2) << L_ins_->my_sequence_template_varSite_left->min; 
			cout << setw(6) << L_ins_->my_sequence_template_varSite_left->max << ")"; 
		} else cerr << "NULL";
		if (L_ins_->my_sequence_template_varSite_right != NULL) {
			cout << setw(10) << L_ins_->my_sequence_template_varSite_right;
			cout << "(" << setw(2) << L_ins_->my_sequence_template_varSite_right->min; 
			cout << setw(6) << L_ins_->my_sequence_template_varSite_right->max << ")"; 
		} else cerr << "NULL";

		if (L_ins_->my_motif_varSite_left != NULL) {
			cout << "m L_ins_=" << setw(10) << L_ins_; 
			cout << setw(10) << L_ins_->my_motif_varSite_left;
			cout << "(" << setw(2) << L_ins_->my_motif_varSite_left->min; 
			cout << setw(6) << L_ins_->my_motif_varSite_left->max << ")"; 
		} else cerr << "NULL";
		if (L_ins_->my_motif_varSite_right != NULL) {
			cout << setw(10) << L_ins_->my_motif_varSite_right;
			cout << "(" << setw(2) << L_ins_->my_motif_varSite_right->min; 
			cout << setw(6) << L_ins_->my_motif_varSite_right->max << ")"; 
		} else cerr << "NULL";
	}

	if (R_ins_ != NULL) {
		cout << "st R_ins_=" << setw(10) << L_ins_; 
		if (R_ins_->my_sequence_template_varSite_left != NULL) {
			cout << setw(10) << R_ins_; 
			cout << setw(10) << R_ins_->my_sequence_template_varSite_left;
			cout << "(" << setw(2) << R_ins_->my_sequence_template_varSite_left->min; 
			cout << setw(6) << R_ins_->my_sequence_template_varSite_left->max << ")"; 
		} else cerr << "NULL";
		if (R_ins_->my_sequence_template_varSite_right != NULL) {
			cout << setw(10) << R_ins_->my_sequence_template_varSite_right;
			cout << "(" << setw(2) << R_ins_->my_sequence_template_varSite_right->min; 
			cout << setw(6) << R_ins_->my_sequence_template_varSite_right->max << ")"; 
			cout << endl;
		} else cerr << "NULL";

		if (R_ins_->my_motif_varSite_left != NULL) {
			cout << "m R_ins_=" << setw(10) << L_ins_; 
			cout << setw(10) << R_ins_; 
			cout << setw(10) << R_ins_->my_motif_varSite_left;
			cout << "(" << setw(2) << R_ins_->my_motif_varSite_left->min; 
			cout << setw(6) << R_ins_->my_motif_varSite_left->max << ")"; 
		} else cerr << "NULL";
		if (L_ins_->my_motif_varSite_right != NULL) {
			cout << setw(10) << R_ins_->my_motif_varSite_right;
			cout << "(" << setw(2) << R_ins_->my_motif_varSite_right->min; 
			cout << setw(6) << R_ins_->my_motif_varSite_right->max << ")"; 
			cout << endl;
		} else cerr << "NULL";
	}

	if (del != NULL) {
		cout << "st:m del=" << setw(10) << del; 
		if (del->my_sequence_template_varSite != NULL) {
			cout << setw(10) << del->my_sequence_template_varSite;
			cout << "(" << setw(2) << del->my_sequence_template_varSite->min;
			cout << setw(6) << del->my_sequence_template_varSite->max;
		} else cerr << "NULL";

		if (del->my_motif_varSite != NULL) {
			cout << setw(10) << del->my_motif_varSite;
			cout << "(" << setw(2) << del->my_motif_varSite->min;
			cout << setw(6) << del->my_motif_varSite->max;
		} else cerr << "NULL";
	}
}

void Indel::set(
				Indel *site2copy, 
				Indel *prev_site, 
				short type, 
				bool isDescendant, 
				Indel *anc_indel
			   )
{
	Deletion  *anc_del = NULL;
	if (anc_indel != NULL) 
		if (anc_indel->del != NULL) 
			anc_del = anc_indel->del;
	
	//////////
	/// Del must be before the insertion info, b/c L_ins_ will ALWAYS use the del mstv and mmv
	/// structures. More reliable, less bookwork.
	//////////
	del->copy(site2copy->del, type, isDescendant, anc_del);
	L_ins_->L_ins_copy(site2copy->del, prev_site->del, type);
}


////////////////////
////// Insertion
////////////////////
Insertion::Insertion(
					 Insertion *site2copy, 
					 motifSite *anc, 
					 bool isTemplate, 
					 Insertion *anc_copy
					)
{
	// Copying ancestral, pointing to current struct.
	if (isTemplate) {
		if (site2copy->my_sequence_template_varSite_left == NULL) return;
		my_sequence_template_varSite_left = site2copy->my_sequence_template_varSite_left;
		my_sequence_template_varSite_right = site2copy->my_sequence_template_varSite_right;
		on_site_sequence_template_varSite = site2copy->on_site_sequence_template_varSite;

		if (anc_copy != NULL) {
			my_sequence_template_varSite_left 
			= anc_copy->my_sequence_template_varSite_left->descendant_equiv;
			if (anc_copy->on_site_sequence_template_varSite != NULL)
				on_site_sequence_template_varSite 
				= anc_copy->on_site_sequence_template_varSite->descendant_equiv;			
			else on_site_sequence_template_varSite = NULL;
			if (site2copy->my_sequence_template_varSite_right == NULL) { 
				cerr << "null." << endl; 
				return; 
			}
			my_sequence_template_varSite_right 
			= anc_copy->my_sequence_template_varSite_right->descendant_equiv;

		}
	} else {
		if (site2copy->my_motif_varSite_left == NULL) return;
		my_motif_varSite_left 
		= site2copy->my_motif_varSite_left;
		my_motif_varSite_right 
		= site2copy->my_motif_varSite_right;
		on_site_motif_varSite = site2copy->on_site_motif_varSite;

		if (anc_copy != NULL) {
			my_motif_varSite_left 
			= anc_copy->my_motif_varSite_left;
			if (anc_copy->on_site_motif_varSite != NULL)
				on_site_motif_varSite 
				= anc_copy->on_site_motif_varSite->descendant_equiv;
			else on_site_motif_varSite = NULL;
			if (site2copy->my_motif_varSite_right == NULL) { 
				cerr << "null." << endl; 
				return; 
			}
			my_motif_varSite_right 
			= anc_copy->my_motif_varSite_right;
		}
	}
}

void Insertion::copy(
					 Insertion *site2copy, 
					 bool isTemplate, 
					 bool isDescendant, 
					 Insertion *anc_copy
					)
{
	if (isTemplate) {
		if (site2copy->my_sequence_template_varSite_left == NULL) return;
		if (isDescendant) {
			my_sequence_template_varSite_left 
			= anc_copy->my_sequence_template_varSite_left->descendant_equiv;
			if (anc_copy->on_site_sequence_template_varSite != NULL)
				on_site_sequence_template_varSite 
				= anc_copy->on_site_sequence_template_varSite->descendant_equiv;
			else on_site_sequence_template_varSite = NULL;
			if (anc_copy->my_sequence_template_varSite_right != NULL)
				my_sequence_template_varSite_right 
				= anc_copy->my_sequence_template_varSite_right->descendant_equiv;
		} else {
			my_sequence_template_varSite_left 
			= site2copy->my_sequence_template_varSite_left;
			my_sequence_template_varSite_right 
			= site2copy->my_sequence_template_varSite_right;
			on_site_sequence_template_varSite = site2copy->on_site_sequence_template_varSite;
		}
	} else {
		if (site2copy->my_motif_varSite_left == NULL) return;
		if (isDescendant) {
			my_motif_varSite_left 
			= anc_copy->my_motif_varSite_left->descendant_equiv;
			if (anc_copy->on_site_motif_varSite != NULL) {
				on_site_motif_varSite 
				= anc_copy->on_site_motif_varSite->descendant_equiv;
			} else on_site_motif_varSite = NULL;
			if (anc_copy->my_motif_varSite_right != NULL)
				my_motif_varSite_right 
				= anc_copy->my_motif_varSite_right->descendant_equiv;
		} else {
			my_motif_varSite_left 
			= site2copy->my_motif_varSite_left;
			my_motif_varSite_right 
			= site2copy->my_motif_varSite_right;		
			on_site_motif_varSite = site2copy->on_site_motif_varSite;
		}
	}
}

//////////
/// This version of copy is utilized when an insertion is made (i.e., during the evolution run, 
/// rather than the evolutionary setup).
//////////
void Insertion::copy(
					 Insertion *site2copy
					)
{
	my_sequence_template_varSite_left 
	= site2copy->my_sequence_template_varSite_left;
	my_sequence_template_varSite_right 
	= site2copy->my_sequence_template_varSite_right;
	on_site_sequence_template_varSite = site2copy->on_site_sequence_template_varSite;
	my_motif_varSite_left 
	= site2copy->my_motif_varSite_left;
	my_motif_varSite_right 
	= site2copy->my_motif_varSite_right;
	on_site_motif_varSite = site2copy->on_site_motif_varSite;
}

void Insertion::setMembership(
							  varSite *memberOfSite, 
							  bool isTemplate, 
							  size_t side
							 )
{
	if (isTemplate) {
		if (side == LEFT) my_sequence_template_varSite_left = memberOfSite;
		else my_sequence_template_varSite_right = memberOfSite;
		if (memberOfSite->min == 0) on_site_sequence_template_varSite = memberOfSite;
	} else {
		if (side == LEFT) my_motif_varSite_left = memberOfSite;
		else my_motif_varSite_right = memberOfSite;
		if (memberOfSite->min == 0) on_site_motif_varSite = memberOfSite;
	}
}

bool Insertion::insertionAllowed(
								 int size
								)
{
	bool profile=false;

	if (profile) {
	//	cerr << "Insertion::insertionAllowed: INSERTION_STATS: (insert size = " << size << ")" << endl;
	//	if (on_site_sequence_template_varSite != NULL)
	//		cerr << "Insertion::insertionAllowed:   osst: " << on_site_sequence_template_varSite->min << " " << on_site_sequence_template_varSite->max << " " << on_site_sequence_template_varSite->member_set.size() << endl;
	//	if (on_site_motif_varSite != NULL)
	//		cerr << "Insertion::insertionAllowed:   osm:  " << on_site_motif_varSite->min << " " << on_site_motif_varSite->max << " " << on_site_motif_varSite->member_set.size() << endl;
	//	cerr << "Insertion::insertionAllowed:   mstl: " << my_sequence_template_varSite_left->min << " " << my_sequence_template_varSite_left->max << " " << my_sequence_template_varSite_left->member_set.size() << endl;
	//	cerr << "Insertion::insertionAllowed:   mstr: " << my_sequence_template_varSite_right->min << " " << my_sequence_template_varSite_right->max << " " << my_sequence_template_varSite_right->member_set.size() << endl;
	//	cerr << "Insertion::insertionAllowed:   mml:  " << my_motif_varSite_left->min << " " << my_motif_varSite_left->max << " " << my_motif_varSite_left->member_set.size() << endl;
	//	cerr << "Insertion::insertionAllowed:   mmr:  " << my_motif_varSite_right->min << " " << my_motif_varSite_right->max << " " << my_motif_varSite_right->member_set.size() << endl;
	}

	//////////
	/// on_site_xxx set only if the minimum in the motif or sequence template is zero.
	//////////
	// If there is both a template and motif covering this site with a zero minimum
	if (
		 on_site_sequence_template_varSite != NULL && 
		 on_site_motif_varSite != NULL
	   ) 
	{
		if (
			 on_site_sequence_template_varSite->insertion(size) && 
			 on_site_motif_varSite->insertion(size)
		   ) 
		{
			if (profile) cerr << "1true" << endl;
			return true;
		}
	}
	//
	// Just a motif
	if (
		on_site_sequence_template_varSite == NULL && 
		on_site_motif_varSite != NULL
	   ) {
		// Motif will accept, but need to check if template left or right will accept an insertion.
		if (
			 on_site_motif_varSite->insertion(size) && 		// Will on_site motif accept?
		     (
		       my_sequence_template_varSite_left->insertion(size) || // WIll either st accept?
		       my_sequence_template_varSite_right->insertion(size)
		     )
		   ) 
		{
			if (profile) cerr << "2true" << endl;
			return true;
		}
	} 
	//
	// Just a template
	if (
		 on_site_sequence_template_varSite != NULL && 
		 on_site_motif_varSite == NULL
	   ) 
	{
		if (
			 on_site_sequence_template_varSite->insertion(size) && 
			 (
			   my_motif_varSite_left->insertion(size) || 
			   my_motif_varSite_right->insertion(size)
			 )
		   ) 
		{
			if (profile) cerr << "3true" << endl;
			return true;
		}
	}
	//
	// Neither motif nor template have a minimum size of zero. In this case, need to check
	// the M or T in order to see if the maximum-current_num is large enough to handle the insert.
	// Left:
	if ( my_sequence_template_varSite_left->insertion(size)
		 && my_motif_varSite_left->insertion(size) ) {
		 if (profile) cerr << "4true" << endl;
		 return true;
	}
	//
	// Right:
	if ( my_sequence_template_varSite_right->insertion(size) 
		 && my_motif_varSite_right->insertion(size) ) {
		 if (profile) cerr << "5true" << endl;
		 return true;
	}
	//
	// No insert can be fit in this site.
	if (profile) cerr << "false" << endl;
	return false;
	//
	//////////
}

void Insertion::L_ins_copy(
						   Deletion *site2copy, 
						   Deletion *prev_site, 
						   short type
						  )
{
	//////////
	/// Because we are setting del before L_ins_, no longer need to know if this is a descendant
	/// node! Del will contain the necessary descendant_equiv information, and is much more
	/// reliable than working with insertion sites, since there is 1 del/site.
	//////////
	if (type == TEMPLATE) {
		if (site2copy->my_sequence_template_varSite == NULL) { cerr << "mstv null??" << endl; return; }
		//////////
		/// Set the Insertion pointers to the Deletion of the actual site.
		//////////
		// * pointer to the left -> prev_site.
		my_sequence_template_varSite_left 
		= prev_site->my_sequence_template_varSite;
		//
		// * pointer to the right -> this site.
		my_sequence_template_varSite_right 
		= site2copy->my_sequence_template_varSite;
		//
		// * on site -> (if either site == 0) ? set : NULL
		//   where set is min(a,b) if both are to set.
		on_site_sequence_template_varSite 
		= site2copy->my_sequence_template_varSite->set_on_site_ptr(prev_site->my_sequence_template_varSite);
		//
		//////////
	} else {
		if (site2copy->my_motif_varSite == NULL) return;
		//////////
		/// Set the (L_ins_) Insertion pointers to the Deletion of the actual site.
		//////////
		// * pointer to the left -> prev_site.
		my_motif_varSite_left 
		= prev_site->my_motif_varSite;
		//
		// * pointer to the right -> this site.
		my_motif_varSite_right 
		= site2copy->my_motif_varSite;
		//
		// * on site -> (if either site == 0) ? set : NULL
		//   where set is min(a,b) if both are to set.
		on_site_motif_varSite 
		= site2copy->my_motif_varSite->set_on_site_ptr(prev_site->my_motif_varSite);
		//
		//////////
	}
}


////////////////////
////// DELETION
////////////////////
void Deletion::copy(
					Deletion *site2copy, 
					bool isTemplate, 
					bool isDescendant, 
					Deletion *anc_copy
				   )
{
	if (isTemplate) { 
		if (isDescendant) my_sequence_template_varSite = anc_copy->my_sequence_template_varSite->descendant_equiv;
		else my_sequence_template_varSite = site2copy->my_sequence_template_varSite;
	} else { 
		if (isDescendant) 
			my_motif_varSite 
			= anc_copy->my_motif_varSite->descendant_equiv;
		else 
			my_motif_varSite 
			= site2copy->my_motif_varSite;
	}
}

void Deletion::setMembership(
							 varSite *memberOfSite, 
							 bool isTemplate
							)
{
	if (isTemplate) my_sequence_template_varSite = memberOfSite;
	else my_motif_varSite = memberOfSite;
}

bool Deletion::deletionAllowed(
							   int size
							  )
{

	//cerr << "deletionAllowed size = " << size;
	//cerr << "  mstv: " << my_sequence_template_varSite->min << " " << my_sequence_template_varSite->max << " [" << my_sequence_template_varSite->member_set.size();
	//cerr << "]  mmv:  " << my_motif_varSite->min << " " << my_motif_varSite->max << " [" << my_motif_varSite->member_set.size() << "]       ";

	if ( my_sequence_template_varSite->deletion(size) && my_motif_varSite->deletion(size) ) {
		//cerr << "true" << endl;
		return true;
	}

	//cerr << "false" << endl;
	return false;
}

void Deletion::copy(
					Deletion *site2copy, 
					short type, 
					bool isDescendant, 
					Deletion *anc_copy
				   )
{
	if (type == TEMPLATE) { 
		if (isDescendant) {
			my_sequence_template_varSite 
			= anc_copy->my_sequence_template_varSite->descendant_equiv;
		} else 
			my_sequence_template_varSite 
			= site2copy->my_sequence_template_varSite;
	} else { 
		if (isDescendant) {
			my_motif_varSite 
			= anc_copy->my_motif_varSite->descendant_equiv;
		} else 
			my_motif_varSite 
			= site2copy->my_motif_varSite;
	}
}


////////////////////
////// SUBSTITUTION
////////////////////
bool Substitution::siteInvariable()
{
	if (substitution_bitstring.count() == 1) return true;
	else return false;
}

void Substitution::setInvariable(
								 char residue
								)
{
	substitution_bitstring.reset();
	substitution_bitstring.set(residue);
}

void Substitution::copy(
						Substitution *site2copy
					   )
{
	substitution_bitstring.set();
	substitution_bitstring &= site2copy->substitution_bitstring;
}

string Substitution::report_bitset()
{
	string report = "";
	for (int i = 0; i < numStates; i++)
		if (substitution_bitstring[i]) report += stateCharacters[i];

	return report;
}

void Substitution::setSiteBits(
							   bitset<20> bitstring
							  )
{
	substitution_bitstring.set();
	substitution_bitstring &= bitstring;
}


////////////////////
/////// VAR SITE
////////////////////
void varSite::add2Member (
						  Deletion *newMember
						 )
{
	//cerr << "Adding member: " << newMember << endl;
	member_set.insert(newMember);
}

void varSite::removeFromMember (
								Deletion *thisMember
							   )
{
	member_set.erase(thisMember);
}

bool varSite::insertion(
					    int size
					   )
{
	if (member_set.size() + size > max) return false; 
	if (max < member_set.size() + size) return false;
	return true;
}

bool varSite::deletion(
					   int size
					  )
{
	if (min == max) return false;
	if (min > member_set.size() - size) return false;

	return true;
}

bool varSite::isUnconstrained() 
{
	if (min == 0 && max == ISG_INT_MAX) return true;
	return false;
}

void varSite::print_members()
{
	cerr << "Members (set): " << member_set.size() << "   " << min << " " << max << endl;
	for (set<Deletion*>::iterator it = member_set.begin(); it != member_set.end(); ++it) {
		cerr << (*it) << endl;
	}


}

varSite *varSite::set_on_site_ptr(
								  varSite *prev_site
								 )
{
	//////////
	/// * on site -> (if either site == 0) ? set : NULL
	///   where set is min(a,b) if both are to set.
	//////////
	if (prev_site->min == 0) 
		if(min == 0)
			return (
			    	 (prev_site->max < max)
			  		 ? prev_site
					 : this
			  	   );
		else
			return prev_site;
	else if (min == 0)
		return this;
	else
		return NULL;

}
