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

#include "evolve.h"

#define FSL_FAIL -1

using namespace std;

static double 	NeighborRates[400];
int 			indel_fill_aas;
double 			freqRate[MAX_RATE_CATS];
extern int		num_insert, num_delete, len_insert, len_delete;
int 			num_subst = 0, rejected_subst = 0, num_ins = 0, num_del = 0;
extern vector<int> value_check;
double			global_sum_Pr_subst = 0;
double	Pr_subst_site_sum_I, Pr_subst_site_sum_D;
double	total_diff_Pr_subst_site_sum = 0;
bool round1 = false;

vector<int>		cat_chosen (4,0);
vector<short>	acgt (20, 0);
vector<int>		acagatcgctgt (400, 0);
extern int changed_site;
extern int prev_state;

// functions 

void Create_Global_Arrays(
						  TTree *tree, 
						  int inNumSites
						 ) 
{
	//////////
	/// Root sequence, coming from root node.
	//////////
	tree->global_alignment->insert_sites.clear();
	tree->global_alignment->insert_sites.assign(inNumSites, insertSite('i',-1,-1));
}

void SetCategories(
				   TNode *node, 
				   int inNumSites, 
				   seqGenOptions *options
				  )
{
	if (node->anc != NULL) {
		stringstream message;
		switch (node->nodeEnv->rateHetero) {
		case GammaRates:
			// GammaRates->GammaRates, NoRates->GammaRates
			if (node->anc->nodeEnv->rateHetero != GammaRates && node->anc->nodeEnv->rateHetero != NoRates) {
				message << "Rate heterogeneity must remain gamma rates in clade parameters: Alpha parameter can be changed." << endl;
				message << "  Cannot change " << node->anc->nodeEnv->report_rateHetero() << " to Gamma Rates.";
				options->SpoolWarnings(message.str());
				node->nodeEnv->rateHetero = node->anc->nodeEnv->rateHetero;
				if (node->anc->nodeEnv->rateHetero == DiscreteGammaRates) 
					node->nodeEnv->numCats = node->anc->nodeEnv->numCats;
				// Changed the rateHetero of the current node, reprocess.
				SetCategories(node, inNumSites, options);
				return;
			}
			break;
		case DiscreteGammaRates:
			// DiscreteGammaRates->DiscreteGammaRates, NoRates->DiscreteGammaRates
			if (node->anc->nodeEnv->rateHetero != DiscreteGammaRates && node->anc->nodeEnv->rateHetero != NoRates) {
				message << "Rate heterogeneity must remain discrete gamma rates in clade parameters.";
				options->SpoolWarnings(message.str());
				node->nodeEnv->rateHetero = node->anc->nodeEnv->rateHetero;
				if (node->anc->nodeEnv->rateHetero == GammaRates)
					node->nodeEnv->gammaShape = node->anc->nodeEnv->gammaShape;
				SetCategories(node, inNumSites, options);
				return;
			}
			break;
		case CodonRates:
			// CodonRates->CodonRates
			if (node->anc->nodeEnv->rateHetero != CodonRates) {
				message << "Cannot change from codon rates in lineages";
				options->SpoolWarnings(message.str());
				node->nodeEnv->rateHetero = node->anc->nodeEnv->rateHetero;
				SetCategories(node, inNumSites, options);
				return;
			}
			break;
		case NoRates:
			// NoRates->NoRates, GammaRates->NoRates, DiscreteGammaRates->NoRates
			if (node->anc->nodeEnv->rateHetero == CodonRates) {
				message << "Cannot change from Uniform Rates to Codon Rates in lineages.";
				options->SpoolWarnings(message.str());
				node->nodeEnv->rateHetero = CodonRates;
				node->nodeEnv->catRate[0] = node->anc->nodeEnv->catRate[0];
				node->nodeEnv->catRate[1] = node->anc->nodeEnv->catRate[1];
				node->nodeEnv->catRate[2] = node->anc->nodeEnv->catRate[2];
			} else if (node->anc->nodeEnv->rateHetero == NoRates) {
				return;
			} else {
				node->nodeEnv->rateHetero = node->anc->nodeEnv->rateHetero;
				SetCategories(node, inNumSites, options);
				return;
			}
			break;
		}
	}

	if (node->nodeEnv->rateHetero != node->anc->nodeEnv->rateHetero || node->mytipNo == -1) {
		// If not returned from above, we're cool to change the rates.
	    if (node->nodeEnv->rateHetero==CodonRates) {
			double sumRates=node->nodeEnv->catRate[0]+node->nodeEnv->catRate[1]+node->nodeEnv->catRate[2];
			if (sumRates!=3.0) {
		    	node->nodeEnv->catRate[0]*=3.0/sumRates;
		    	node->nodeEnv->catRate[1]*=3.0/sumRates;
		    	node->nodeEnv->catRate[2]*=3.0/sumRates;
			}
	    } else if (node->nodeEnv->rateHetero==GammaRates) {
			for (size_t i=0; i<inNumSites; i++) {
				node->seq_evo.at(i).setGamma(rndgamma(node->nodeEnv->gammaShape) / node->nodeEnv->gammaShape);
			}
	    } else if (node->nodeEnv->rateHetero==DiscreteGammaRates) {
			DiscreteGamma(freqRate, node->nodeEnv->catRate, node->nodeEnv->gammaShape, node->nodeEnv->gammaShape, node->nodeEnv->numCats, 0);
			for (size_t i=0; i<inNumSites; i++) node->seq_evo.at(i).setCategory((int)(rndu()*node->nodeEnv->numCats));
//			for (size_t i=0; i<inNumSites; i++) node->seq_evo.at(i).setCategories((int)(rndu()*node->nodeEnv->numCats));
	    }
	} else {
		//////////
		/// Need to make sure everything is inherited.
		//////////
		int i;
		switch (node->nodeEnv->rateHetero) {
		case DiscreteGammaRates:
			node->nodeEnv->gammaShape = node->anc->nodeEnv->gammaShape;
			node->nodeEnv->numCats = node->anc->nodeEnv->numCats;
			for (i=0; i < node->nodeEnv->numCats; i++) 
				node->nodeEnv->catRate[i] = node->anc->nodeEnv->catRate[i];
			break;
		case GammaRates:
			node->nodeEnv->gammaShape = node->anc->nodeEnv->gammaShape;
			break;
		case CodonRates: 
			node->nodeEnv->catRate[0] = node->anc->nodeEnv->catRate[0];
			node->nodeEnv->catRate[1] = node->anc->nodeEnv->catRate[1];
			node->nodeEnv->catRate[2] = node->anc->nodeEnv->catRate[2];
			break;
		default:
			break;
		}
	}
}

char SetState(
			  double *P, 
			  string& caller
			 )
{
    char j = 0;
    double r;
	bool done = false;
	int num_rounds = 0;
	double *original_P = P;
	
	do {
		P = original_P;
		num_rounds++;
	    r=rndu();
	    for (j=0; r>=(*P) && j<numStates; j++) P++;
		if (j < numStates) {
			done = true;
		}
		if (num_rounds > 250) {
			cerr << "Received value " << (double)(*P) << " from " << caller << ": could not find a suitable value to return in SetState." << endl;
			cerr << "Value r = " << r << " vs value *P = " << *original_P << " that was received." << endl;
			exit(EXIT_FAILURE);
		}
	} while (!done);

    return j;
}

short IsInvariable(
				   double proportionInvariable
				  )
{
    double r;

    r=rndu();
    if (r < proportionInvariable) return INVARIABLE;
    else return 0;
}

void RandomSequence(
					char *seq, 
					int inNumSites, 
					seqGenOptions *options,
					RateMatrix *rates
				   )
{
    char *P;
	string calling_routine = "RandomSequence";
	double addFreq[numStates];

	addFreq[0] = rates->pi.at(0);
	int i = 1;
	for (vector<double>::iterator it = rates->pi.begin()+1; it != rates->pi.end(); ++it, i++)
		addFreq[i] = addFreq[i-1]+(*it);

    P=seq;
    for (size_t i=0; i<inNumSites; i++) {
		*P=SetState(addFreq, calling_routine);
		P++;
    }

	if (options->default_rateHetero == CodonRates) {
		int codon[3];
		for (size_t i = 0; i < (int)(inNumSites/3); i++) {
			codon[0] = seq[i*3];
			codon[1] = seq[i*3+1];
			codon[2] = seq[i*3+2];
			while (Stop_Codon(codon)) {
				codon[0] = seq[i*3] = SetState(addFreq, calling_routine);
				codon[1] = seq[i*3+1] = SetState(addFreq, calling_routine);
				codon[2] = seq[i*3+2] = SetState(addFreq, calling_routine);
			}
		}
	}
}

string insertFillSequence(
						  int insertSize,
						  RateMatrix *rates
						 )
{
    vector<double> bayes_probabilities_for_next (numStates, 0);
    double r;
	string insert_sequence (insertSize+1, '\0');
	string calling_routine;
	double addFreq[numStates];

	addFreq[0] = rates->pi.at(0);
	int i = 1;
	for (vector<double>::iterator it = rates->pi.begin()+1; it != rates->pi.end(); ++it, i++)
		addFreq[i] = addFreq[i-1]+(*it);
	
	if (!isNucModel) {
		calling_routine.assign("insertFillSequence, isNucModel");
	    insert_sequence.at(0) = SetState(addFreq, calling_routine);
	    for(size_t i=1; i<insertSize+1; i++) {
			bayes_probabilities_for_next.at(0)=(rates->pi.at(0)*NeighborRates[(int)insert_sequence.at(i-1)*numStates])/(rates->pi.at(insert_sequence.at(i-1)));
	        for(size_t j=1;j<numStates;j++) {
		    	bayes_probabilities_for_next.at(j)=bayes_probabilities_for_next.at(j-1);
		    	bayes_probabilities_for_next.at(j)+=(rates->pi.at(j)*NeighborRates[(int)insert_sequence.at(i-1)*numStates+j])/(rates->pi.at(insert_sequence.at(i-1)));
			}

			r=rndu();
			r*=bayes_probabilities_for_next.at(numStates-1);
			size_t j=0;
			while(r>bayes_probabilities_for_next.at(j)) j++; 
			if (j >= numStates) {
				cerr << "Illegal character found in iFS: " << j << endl;
				exit(EXIT_FAILURE);
			}
			insert_sequence.at(i) = j;
	    }
	} else {
		for (size_t i = 0; i < insertSize+1; i++) {
			calling_routine.assign("insertFillSequence, isNucModel");
			insert_sequence.at(i) = SetState(addFreq, calling_routine);
			if (insert_sequence.at(i) >= numStates) {
				cerr << "illegal character found for insertSize in iFS: " << insert_sequence.at(i) << endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	return insert_sequence;
}

bool Stop_Codon(
				int *codon
			   ) 
{
        size_t T = stateCharacters.find("T"), 
        	   A = stateCharacters.find("A"), 
        	   G = stateCharacters.find("G");
                                
        if (codon[0] == T && codon[1] == A && codon[2] == G) return true;
        if (codon[0] == T && codon[1] == G && codon[2] == A) return true;
        if (codon[0] == T && codon[1] == A && codon[2] == A) return true;
                                
        return false;
}

int calcTrials(
			   TNode *des, 
			   int *E_in_, 
			   int *E_del_
			  )
{
	*E_in_ = *E_del_ = 0;
	if(des->nodeEnv->invariableSites) {
		for (vector<Site>::iterator it = des->seq_evo.begin(); it != des->seq_evo.end(); ++it) {
			//////////
			/// Calculations or ins/del probabilities.
			//////////
			// Insertion calculations:
			if ( (*it).motif.active_properties.indel->L_ins_->insertionAllowed() ) (*E_in_)++;
			//
			// Deletion  calculations:		
			if ( (*it).motif.active_properties.indel->del->deletionAllowed() ) (*E_del_)++;
			//
			// Last site right needs to be checked.
			if ( it+1 == des->seq_evo.end() )
				if ( (*it).motif.active_properties.indel->R_ins_->insertionAllowed() ) (*E_in_)++;
			//
			//////////
		}
	} else {
		*E_in_ = des->seq_evo.size()+1;
		*E_del_ = des->seq_evo.size();
	}

	// C&B are half n half.
	if (des->nodeEnv->P_ins_ == 0 && des->nodeEnv->P_del_ == 0) 
		*E_in_ = *E_del_ = (*E_in_+*E_del_) / 2;

	return *E_in_ + *E_del_;
}

int  Find_Action(
				 double len, 
				 double P_insert_, 
				 double P_delete_, 
				 int CnB_divisor, 
				 int action_test
				) 
{
    if( !P_insert_ && !P_delete_) {
		//////////
		/// 021210:
		/// len comes in as the length of the branch / CnB_divisor. Undo this by multiplying the
		/// len by CnB_divisor.
		//////////
        double P_indel_ = ( 0.0224 - 0.0219 * exp(-0.01168*((len*CnB_divisor)*100))) / (double)CnB_divisor;
		if((double)rndu() < P_indel_) {
            if((double)rndu() < 0.5) return INSERT;
            else return DELETE;
        } else {
        	return NO_ACTION;
		}
    } else {
		if (action_test == INSERT) {
			if ((double)rndu() < P_insert_*len) return INSERT;
		} else  {
			double chkval = (double)rndu();
			if (chkval < P_delete_*len) {
				return DELETE;
			}
		}
    } 
    return NO_ACTION;
} 

void Find_Motif_Positions(
						  TTree *tree, 
						  TNode *des, 
						  motifSite *this_site, 
						  int indel_size, 
						  bool back, 
						  varSite **in_template_varSite, 
						  varSite **in_motif_varSite
						 )
{
	bool st_L, m_L, st_R, m_R;

	st_L = this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left->insertion(indel_size);
	m_L  = this_site->active_properties.indel->R_ins_->my_motif_varSite_left->insertion(indel_size);
	st_R = this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right->insertion(indel_size);
	m_R  = this_site->active_properties.indel->R_ins_->my_motif_varSite_right->insertion(indel_size);

	bool os_st = false, os_m = false;
	if (this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite)
		os_st = this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite->insertion(indel_size);
	if (this_site->active_properties.indel->R_ins_->on_site_motif_varSite)
		os_m = this_site->active_properties.indel->R_ins_->on_site_motif_varSite->insertion(indel_size);

	// For the randomization of finding which cases I try. Each time I try a case, I set that bit
	// to a '1' to avoid infinite loop (if it is even truly a possibility).
	bitset<6> tried_case (string("000000"));

	double random = (double)rndu(), random2;
	if ( back ) {
		*in_template_varSite = this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left;
		*in_motif_varSite = this_site->active_properties.indel->R_ins_->my_motif_varSite_left;
	} else {
		bool found_site = false;
		int num_rounds = 0;

		// Randomize the position of the varSite where insertion occurs.
		while(!found_site) {
			found_site = true;
			random2 = (double)rndu();
			size_t rand = (size_t)(random2 * 6.0);

			if(num_rounds++ > 250) {
				des->Print_Active_Properties(); 
				cerr << "=====case left undealt with in Insert routine:" << endl;
				cerr << "st_L " << st_L << " st_R " << st_R << " m_L " << m_L << " m_R " << m_R << " os_st " << os_st << " os_m " << os_m << endl;
				exit(EXIT_FAILURE);
			}

			// If already tested, no need to re-invent the wheel (i.e., it still won't work).
			if (!tried_case.test(rand)) {
				tried_case.set(rand, true);
				switch (rand) {

				case 0:
//					cerr << "case 0" << endl;
					if (st_L) {
						*in_template_varSite 
						= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left;
						if (m_L && os_m)
							if (random < 0.5) 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->my_motif_varSite_left;
							else 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_motif_varSite;
						else if (os_m) 
							*in_motif_varSite 
							= this_site->active_properties.indel->R_ins_->on_site_motif_varSite;
						else if (m_L) 
							*in_motif_varSite 
							= this_site->active_properties.indel->R_ins_->my_motif_varSite_left;
						else found_site = false;
					} else found_site = false;

					break;

				case 1:
//					cerr << "case 1" << endl;
					if (st_R) {
						*in_template_varSite 		
						= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right;
						if (m_R && os_m)
							if (random < 0.5) 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->my_motif_varSite_right;
							else 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_motif_varSite;
						else if (os_m) 
							*in_motif_varSite 
							= this_site->active_properties.indel->R_ins_->on_site_motif_varSite;
						else if (m_R) 
							*in_motif_varSite 
							= this_site->active_properties.indel->R_ins_->my_motif_varSite_right;
						else found_site = false;
					} else found_site = false;

					break;

				case 2:
//					cerr << "case 2" << endl;
					if (m_R) {
						*in_motif_varSite 		
						= this_site->active_properties.indel->R_ins_->my_motif_varSite_right;
						if (st_R && os_st) {
							if (random < 0.5) {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right;
							} else {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite;
							}
						} else if (st_R) {	
							*in_template_varSite 
							= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right;
						} else if (os_st) {
							*in_template_varSite 
							= this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite;
						} else {
							found_site = false;
						}
					} else {
						found_site = false;
					}

					break;

				case 3:
//					cerr << "case 3" << endl;
					if (m_L) {
						*in_motif_varSite 
						= this_site->active_properties.indel->R_ins_->my_motif_varSite_left;
						if (st_L && os_st)
							if (random < 0.5) 
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left;
							else 
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite;
						else if (st_L) 
							*in_template_varSite 
							= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left;
						else if (os_st) 
							*in_template_varSite 
							= this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite;
						else found_site = false;
					} else found_site = false;

					break;

				case 4:
//					cerr << "case 4" << endl;
					if (os_st) {
						*in_template_varSite 
						= this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite;
						if (m_L && m_R && os_m) 
							if (random < 0.333) 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->my_motif_varSite_left;
							else if (random < 0.666) 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->my_motif_varSite_right;
							else 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_motif_varSite;
						else if (m_L && m_R) 
							if (random < 0.5) 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->my_motif_varSite_left;
							else 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->my_motif_varSite_right;
						else if (m_L && os_m) 
							if (random < 0.5) 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->my_motif_varSite_left;
							else 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_motif_varSite;
						else if (m_R && os_m) 
							if (random < 0.5) 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->my_motif_varSite_right;
							else 
								*in_motif_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_motif_varSite;
						else if (m_L) 
							*in_motif_varSite 
							= this_site->active_properties.indel->R_ins_->my_motif_varSite_left;
						else if (m_R) 
							*in_motif_varSite 
							= this_site->active_properties.indel->R_ins_->my_motif_varSite_right;
						else if (os_m) 
							*in_motif_varSite 
							= this_site->active_properties.indel->R_ins_->on_site_motif_varSite;
						else found_site = false;
					} else found_site = false;

					break;

				case 5:
//					cerr << "case 5" << endl;
					if (os_m) {
						*in_motif_varSite 
						= this_site->active_properties.indel->R_ins_->on_site_motif_varSite;
						if (st_L && st_R && os_st) {
							if (random < 0.333) {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left;
							} else if (random < 0.666) {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right;
							} else {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite;
							}
						} else if (st_L && st_R) {	
							if (random < 0.5) {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left;
							} else {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right;
							}
						} else if (st_L && os_st) {
							if (random < 0.5) {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left;
							} else {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite;
							}
						} else if (st_R && os_st) {
							if (random < 0.5) {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right;
							} else {
								*in_template_varSite 
								= this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite;
							}
						} else if (st_L) {
							*in_template_varSite 
							= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left;
						} else if (st_R) {		
							*in_template_varSite 
							= this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right;
						} else if (os_st) {
							*in_template_varSite 
							= this_site->active_properties.indel->R_ins_->on_site_sequence_template_varSite;
						} else {
							found_site = false;
						}
					} else found_site = false;

					break;

				default:
					// should be impossible.
					des->evolvingSequence->print_sequence();
					cerr << this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left->insertion(indel_size) 
						 << " (" 
						 << this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left->min 
						 << "," 
						 << this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left->max 
						 << ") " 
						 << this_site->active_properties.indel->R_ins_->my_motif_varSite_left->insertion(indel_size) 
						 << "  (" 
						 << this_site->active_properties.indel->R_ins_->my_motif_varSite_left->min 
						 << "," 
						 << this_site->active_properties.indel->R_ins_->my_motif_varSite_left->max 
						 << ")  " 
						 << this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right->insertion(indel_size) 
						 << " (" 
						 << this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right->min 
						 << "," 
						 << this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_right->max 
						 << ") " 
						 << this_site->active_properties.indel->R_ins_->my_motif_varSite_right->insertion(indel_size) 
						 << " (" 
						 << this_site->active_properties.indel->R_ins_->my_motif_varSite_right->min 
						 << "," 
						 << this_site->active_properties.indel->R_ins_->my_motif_varSite_right->max 
						 << ")" 
						 << endl 
						 << "members.size() = " 
						 << this_site->active_properties.indel->R_ins_->my_sequence_template_varSite_left->member_set.size() 
						 << " " 
						 << this_site->active_properties.indel->R_ins_->my_motif_varSite_left->member_set.size() 
						 << " Wrong insertion site????" 
						 << endl;

					exit(EXIT_FAILURE);
					break;
				}
			} else found_site = false;
		}
	}
}

int Insert(
		   TTree *tree, 
		   TNode *des, 
		   int indel_size
		  ) 
{
    int start_aa, final_aa;
    int indel_aa_in_globals = 0;
    int position_in_globals = 0;
	bool profile=true;

	//////////
	/// When converting from sequence to evolvingSequence, need to change these functions to the
	/// proper evolvingSequence data. For now, ignoring, try to make evolvingSequence equivalent.
	//////////
	final_aa=des->seq_evo.size();
	cerr << "Evolve::Insert: Just before F_S_L..." << endl;
	start_aa=Find_Start_Location(des,&indel_size,des->seq_evo.size(),INSERT);

	if (start_aa == FSL_FAIL) return 0;

	// For insertion, working on half sites, which in turn means we have to go to the next site
	// to represent an insertion correctly. 
	start_aa++;

	tree->global_alignment->insert_sites.insert(
												tree->global_alignment->locateEvent(des, start_aa),
												indel_size,
												insertSite('i', des->mytipNo, eventNo)
											   );

	//tree->global_alignment->Print();

	if (profile) cerr << endl << endl << "INSERT size " << indel_size << " at " << start_aa << endl;

	//////////
	/// Insert into existing sequence.
	//////////
	// Obtain environment parameters:
	varSite *in_template_varSite, *in_motif_varSite;
	vector<Site>::iterator site_it = des->seq_evo.begin();
	bool Site_back = ( ( des->seq_evo.size() == start_aa ) ? true : false );
	Find_Motif_Positions(
						 tree, 
						 des, 
						 (&(*(site_it+start_aa-1)).motif), 
						 indel_size, 
						 Site_back, 
						 &in_template_varSite, 
						 &in_motif_varSite
						);

	if (profile) {
		cerr << "MOTIF (" << in_motif_varSite 
		     << ") TEMPLATE(" << in_template_varSite << ")" 
		     << " NODE(" << des << ")"
		     << endl;
		cerr << "Ancestral and Descendant varSites: " << endl;

		cerr << "Sequence::setActiveProps: Adding new site to site with membership: ";
		cerr << "st("
 			 << in_template_varSite->min 
			 << ","
			 << in_template_varSite->max
			 << ")["
			 << in_template_varSite->member_set.size()
			 << "] ";

	    //list<varSite*>::iterator des_it =  des->variable_region_list.begin(); 
		//cerr << "PRE-INSERT:" << endl;
		//cerr << "anc length: " << des->anc->seq_evo.size() << "  des length: " << des->seq_evo.size() << endl;
		//for (list<varSite*>::iterator anc_it =  des->anc->variable_region_list.begin();
		//							  anc_it != des->anc->variable_region_list.end(); 
		//						  	  ++anc_it,
		//						  	  ++des_it) 
		//{
		//	cerr << (*anc_it) << "->" << (*anc_it)->descendant_equiv << " " << (*anc_it)->min << "," << (*anc_it)->max << " [" << (*anc_it)->member_set.size() << "]" << "\t\t"
		//		 << (*des_it) << " " << (*des_it)->min << "," << (*des_it)->max << " [" << (*des_it)->member_set.size() << "]" <<endl;
		//}
	}

	if (!des->isVarSite(in_template_varSite)) {
		if (des->isVarSite(in_template_varSite->descendant_equiv)) {
			in_template_varSite = in_template_varSite->descendant_equiv;
		} else {
			cerr << "template->descendant_equiv is NOT a varSite" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (!des->isVarSite(in_motif_varSite)) { 
		if (des->isVarSite(in_motif_varSite->descendant_equiv)) {
			in_motif_varSite = in_motif_varSite->descendant_equiv;
		} else {
			cerr << "motif->descendant_equiv is NOT a varSite" << endl;
			exit(EXIT_FAILURE);
		}
	}

	if (in_template_varSite == NULL || in_motif_varSite == NULL) {
		cerr << "Unsuccessful in locating correct insertion position." << endl;
		exit(EXIT_FAILURE);
	}
	//
	// Create insert Sequence
	//  * Passing the clade environment, so no need to worry about rates later in the function.
	//  * Passing the motif environment, so no need to worry about integrating sites (except for 
	//    L_ins_ in first position and R_ins_ in last position, which need to be connected to
	//    R_ins_ and L_ins_ of the site before and site after, respectively).
	if (profile) cerr << "INSERT! size " << indel_size << "(" << des->seq_evo.size() << ")" << endl;
	auto_ptr<Sequence> insert_sites ( new Sequence(des, indel_size)	);
	insert_sites->init(
					   des,
			  		   insertFillSequence(indel_size, des->branch->rates), 
			  		   des->nodeEnv,
					   (
					     (Site_back)
					     ? &(*(site_it+start_aa-1)).motif
					     : &(*(site_it+start_aa)).motif	
					   ),	
			  		   in_template_varSite, 
			  		   in_motif_varSite
			 		  );
	//
	// Set the properties.
	insert_sites->setActiveProps(true);
	//
	// Connect inner activeProps
	//  * Remove L_ins_ from front, R_ins_ from back
	//  * Make subsequent sites L_ins_ and R_ins_ point to same activeProp.
	//  * Point the deletion objects to the correct mstv and mmv locations.
	//
	// Integrate insertion sequence into current sequence
	//  * Get insert location
	//  * Connect front L_ins_ to previous location
	//  * Connect back R_ins_ to next location.
	vector<Site>::iterator insert_location;
	if (start_aa >= des->seq_evo.size()) insert_location = des->seq_evo.end();
	else insert_location = des->seq_evo.begin()+start_aa;
	//
	// Connect front
	//	
	// Now a vector, so cannot simply use sequence member functions
	//   * The sequence position (initialized to inserted sequence)
	//   * The invariable position, initialized to 0 (which is what we want, so above invar can be removed.)
	//
	// Fig:
    //       <--INSERT HERE
    //    1     2              i1    i2 
	//    S     S              S     S
	//   / \   / \            / \   / \  
	//  /   \ /   \          /   \ /   \
	// I  D  I  D  I        I  D  I  D  I
	//
	// Step IS:
	//  (IS1) Delete i1->L_ins_
	//  (IS2) i1->L_ins_ =  1->R_ins_
	//  (IS3)  2->L_ins_ = i2->R_ins_
	//
	// (IS1)
	delete insert_sites->evolutionaryAttributes.front().motif.active_properties.indel->L_ins_;
	//
	// (IS2) / (IS3)
	if (Site_back) {
		//////////
		/// Set the relevant L_ins_ and R_ins_ pointers to the appropriate spots.
		//////////
		//
		// (IS2)
		insert_sites->evolutionaryAttributes.front().motif.active_properties.indel->L_ins_
		= des->seq_evo.back().motif.active_properties.indel->R_ins_;
		//
		// (IS3) R_ins_ should be fine.
		// Add sites to the end of the sequence.
		for (vector<Site>::iterator site_it =  insert_sites->evolutionaryAttributes.begin();
									site_it != insert_sites->evolutionaryAttributes.end();
									++site_it) 
		{
			des->evolvingSequence->evolutionaryAttributes.push_back(*site_it);
		}
		//
		//////////
	} else {
		//////////
		/// Connect the front and back.
		//////////
		// (IS2)
		insert_sites->evolutionaryAttributes.front().motif.active_properties.indel->L_ins_
		= (*(insert_location-1)).motif.active_properties.indel->R_ins_;
		//
		// (IS3)
		(*(insert_location)).motif.active_properties.indel->L_ins_
		= insert_sites->evolutionaryAttributes.back().motif.active_properties.indel->R_ins_;
		//
		// Insert sites into the correct position.
		des->evolvingSequence->evolutionaryAttributes.insert(
															 insert_location, 
															 insert_sites->evolutionaryAttributes.begin(),
															 insert_sites->evolutionaryAttributes.end()
															);
		//
		//////////
	}
	//
	//////////
	des->Site_postProcess();

	//////////
	/// This will point the deletion and insertion pointers appropriately. The insertion routine
	/// has a bad habit of adding new members to the deletion structure of the position following the
	/// the insertion, instead of the position that it will be a representative of.
	//////////
	// SEE FIGURE ABOVE: Insert routine places insertion into deletion element of position 2. Using
	// this function will set it to the deletion element of position 1. Very important for template
	// and motif constraints. Without this function, we can exceed the specified min and max constraints
	// of template and motif elements.
	des->setSitePointers("Insert");
	//
	//////////

	if (profile) {
		//cerr << "INSERT size " << indel_size << " at " << start_aa << endl;
		//des->Print_Active_Properties();
		//cerr << "ENDINSERT" << endl;
		//exit(0);
	}
	return 1;
}

int Delete(
		   TTree *tree,
		   TNode *des, 
		   int *indel_size
		  ) 
{
    int start_aa, final_aa;
    int indel_aa_in_globals = 0; 
    int position_in_globals; 
    int num_indels_in_globals; 
    int temp_in_globals; 
    int num_rounds = 0; 


    final_aa=des->seq_evo.size();
    start_aa=Find_Start_Location(des,indel_size,des->seq_evo.size(),DELETE);

	if (start_aa == FSL_FAIL) return 0;

	cerr << "Evolve::Delete: sequence before = " << des->printSequence() << endl;
	cerr << "Evolve::Delete: *********************DELETE size " << *indel_size << " at " << start_aa << endl;

	//////////
	/// For multi-residue deletions, note that we will always track to start_aa when we delete 
	/// the residue at start_aa, it will no longer be part of the current sequence, and thus will 
	/// not be counted during our routine. Therefore, we start at the furthest residue, and delete
	/// backwards.
	//////////
	for (int x = *indel_size-1; x >= 0; x--) 
		(*(tree->global_alignment->locateEvent(
											   des,
											   start_aa+x
											  ))).modifiers.push_back(siteModifier('d',des->mytipNo, eventNo));

	vector<Site>::iterator site_it = des->seq_evo.begin()+start_aa;
	for (size_t num_dels = 0; num_dels != *indel_size; num_dels++) {
		if (site_it != des->seq_evo.end()-1 )
			(*site_it).motif.Site_deleteMerge( (&(*(site_it-1)).motif), (&(*(site_it+1)).motif), true );
		else
			(*site_it).motif.Site_deleteMerge( (&(*(site_it-1)).motif), NULL, true);
		des->seq_evo.erase(site_it);
	}

	cerr << "Evolve::Delete: sequence after event " << eventNo << " = " << des->printSequence() << endl;

	//cerr << "Evolve::Delete: END DELETE" << endl;

	return 1;
} 

int Find_Start_Location(
					    TNode *des, 
					    int *indel_size, 
					    int final_aa, 
					    int action
					   ) 
{
    int position;
    int i,j; 
    int ok; 
    int numAcceptablePositions = 0;
    vector<bool> acceptable_positions (final_aa+2, false);
	vector<bool>::iterator bit = acceptable_positions.begin();
	size_t first_constrained = 0, last_constrained = 0;

	des->setSitePointers("Find_Start_Location");
    if(action == INSERT) {
		bit = acceptable_positions.begin();
		numAcceptablePositions = 0;
		position = 0;
		for (vector<Site>::iterator site_it = des->seq_evo.begin(); site_it != des->seq_evo.end(); ++site_it, ++bit, position++) {
			if ( (*site_it).motif.active_properties.indel->R_ins_->insertionAllowed(*indel_size) ) {
				(*bit) = true;
				numAcceptablePositions++;
			}
		}
 	} else { 
		// Two arrays, if a deletion of indel_size is acceptable for each structure, then this pos
		// is an acceptable pos. If not, no dice.
 		vector<bool> st_accept (final_aa+2, false);
 		vector<bool> m_accept (final_aa+2, false);
		varSite *var_site, *var_site2, *current_varSite;
		size_t members, min, current_indel_size, num_dels_in_current_varSite, max;
		size_t number_of_unconstrained;

		// For Gillespie: Deletions can take portions of the sequence before and after the 
		// sequence being simulated. This finds the first and last positions in the sequence
		// that are constrained (i.e., cannot be deleted), and then calculates the number of
		// sub-sequence positions that can be taken by a deletion that overlaps the borders of
		// the sequence being simulated (Since first position is always invariant, start loop
		// after this position. It will make little difference to the overall routine.
		// THIS IS ONLY OK FOR A GUIDE TREE FILE THAT IS ASSUMED TO BE GENES, NOT DOMAINS.
		// Sets the maximum number of positions that can be deleted for the template in N-term.
		if (des->seq_evo.size() > 1) {	// Make sure there is sequence left... Otherwise, KABOOM! //
			// N-TERMINUS: (start at 1 because first position ALWAYS constrained)
			first_constrained = des->findForwardNumberOfConstrained(1, *indel_size);	// Function to find first constrained site
			// C-TERMINUS: (start at 0 [last in rbegin] position because constraints not always imposed)
			last_constrained = des->findBackwardNumberOfConstrained(0,*indel_size);
		}

		bool success;
		vector<bool>::iterator pos;
		pos = st_accept.begin();
		for (vector<Site>::iterator site_it = des->seq_evo.begin(); site_it != des->seq_evo.end(); ++site_it, pos++) {
			// Template array fillage.
			current_varSite = NULL;
			current_indel_size = 0;
			success = true;
			for (vector<Site>::iterator site_it2 = site_it; site_it2 != des->seq_evo.end() && current_indel_size != *indel_size; ++site_it2, current_indel_size++) {
				// If we are in a new varSite, need to reflect that in all data members.
				if ( (*site_it2).motif.active_properties.indel->del->my_sequence_template_varSite != current_varSite) {	
					current_varSite = (*site_it2).motif.active_properties.indel->del->my_sequence_template_varSite;
					members = current_varSite->member_set.size();
					min = current_varSite->min;
					max = current_varSite->max;
					num_dels_in_current_varSite = 0;
				}


				// Too many dels in this site, less than minimum acceptable size.
				num_dels_in_current_varSite++;
				if (members - num_dels_in_current_varSite < min || max == min ) success = false;
			}
			// Near the end of the sequence.
			if (current_indel_size < *indel_size) success = false;
			(*pos) = success;
		}

		pos = m_accept.begin();
		for (vector<Site>::iterator site_it = des->seq_evo.begin(); site_it != des->seq_evo.end(); ++site_it, pos++) {
			// Motif array fillage.
			current_varSite = NULL;
			current_indel_size = 0;
			success = true;
			for (vector<Site>::iterator site_it2 = site_it; site_it2 != des->seq_evo.end() && current_indel_size != *indel_size; ++site_it2, current_indel_size++) {
				// If we are in a new varSite, need to reflect that in all data members.
				if ( (*site_it2).motif.active_properties.indel->del->my_motif_varSite != current_varSite) {	
					current_varSite = (*site_it2).motif.active_properties.indel->del->my_motif_varSite;
					members = current_varSite->member_set.size();
					min = current_varSite->min;
					max = current_varSite->max;
					num_dels_in_current_varSite = 0;
				}

				// Too many dels in this site, less than minimum acceptable size.
				num_dels_in_current_varSite++;
				if (members - num_dels_in_current_varSite < min || max == min ) success = false;
			}
			// Near the end of the sequence.
			if (current_indel_size < *indel_size) success = false;
			(*pos) = success;
		}

		vector<bool>::iterator it2 = m_accept.begin();
		bit = acceptable_positions.begin();
		for (vector<bool>::iterator it = st_accept.begin(); it != st_accept.end(); ++it, ++it2, ++bit) {
			(*bit) = (*it) && (*it2);
			if (*bit) numAcceptablePositions++;
		}
    } 

    if(numAcceptablePositions < 0) {
		fprintf(stderr,"evolve.c, Find_Start_Locations:\nNegative number of acceptable positions.\n");
		fprintf(stderr,"numAcceptablePositions = %d\n",numAcceptablePositions);
		exit(EXIT_FAILURE);
    }

    ok=0; 
    i=0; 
    while(!ok && i<final_aa) { 
        if(acceptable_positions.at(i)) { ok=1; } 
        i++; 
    } 

    if(i >= final_aa) {
		return FSL_FAIL;
    } else {
    	if (action == INSERT) i=(int)(rndu()*(double)numAcceptablePositions+1);
		else i = (int)(rndu()*(double)(numAcceptablePositions+1+first_constrained+last_constrained));
	    position=0; j=0;
	    while(!acceptable_positions.at(j) && j < acceptable_positions.size()) j++;
	    while(position < i && j < final_aa) {
	        if(acceptable_positions.at(j)) position++;
	        if(position < i) j++;
	    } 
		if (action == DELETE) {
			if (i > position) {	// Deletion occurred in either before main seq or after main seq. //
				if (position + first_constrained >= i) { // Occurring at the beginning of seq.
					int x = i - position;
					j = 1;
					*indel_size = x;
				} else if (position + (first_constrained + last_constrained) < i) {
					cerr << "Error in the routine for finding indels at N- and C-terminus." << endl;
					exit(EXIT_FAILURE);
				} else {
					int x = i - (position + first_constrained);
					j = des->seq_evo.size()-x;
					*indel_size = x;
				}
			}
		}
	}

	cerr << "Sequence template accept: ";
	for (vector<bool>::iterator it = acceptable_positions.begin(); it != acceptable_positions.end(); ++it) {
		cerr << (*it);
	}
	cerr << endl;
	cerr << "Evolve::Find_Start_Location: chose site " << j << endl;

    return j; 
}
                     
int Find_Indel_Size(
					vector<double>& dist, 
					int max_indel
				   ) 
{
    int size=0;   
    int i=1;  
    double sum=0;
    double rand; 
     
    rand=rndu(); 
    if(rand==0) rand=rand+0.00001; 
     
    while(i <= max_indel && size == 0) { 
        sum+=dist.at(i);
        if(sum > rand) size=i; 
        i++; 
    } 
    return size;
}        
     
int  isAnc(
		   TNode *thisNode, 
		   int anc
		  ) 
{
//cerr << "Evolve::isAnc: mytipNo = " << thisNode->mytipNo << endl;
    if(thisNode -> mytipNo == anc) {
        return 1; 
    } else if(thisNode -> mytipNo == -1) {
        return 0; 
    } else { 
        return isAnc(thisNode->anc, anc);
    } 
} 

//////////
/// This changes all environment parameters for Nodes (no motif changes).
//////////
void initializeCladeParameters(
							   TTree *tree, 
							   TNode *node
							  ) 
{
	inClade *env;
	if (node->anc == NULL) env = tree->treeEnv.front();
	else env = node->anc->nodeEnv;

	if (!(node->nodeEnv->set_in_clade[__indelFlag__])) {
		node->nodeEnv->indelFlag = env->indelFlag;
		if (env->indelFlag) {
			if (!(node->nodeEnv->set_in_clade[__P_ins___]))
				node->nodeEnv->P_ins_ = env->P_ins_;

			if (!(node->nodeEnv->set_in_clade[__P_del___]))
				node->nodeEnv->P_del_ = env->P_del_;

			if (!(node->nodeEnv->set_in_clade[__maxIndel__]))
				node->nodeEnv->maxIndel = env->maxIndel;

			if (!(node->nodeEnv->set_in_clade[__insert_lengthDistribution__]))
				node->nodeEnv->insert_lengthDistribution = env->insert_lengthDistribution;

			if (!(node->nodeEnv->set_in_clade[__delete_lengthDistribution__])) 
					node->nodeEnv->delete_lengthDistribution = env->delete_lengthDistribution;
		}
	}

	//////////
	/// NOTE: if values2Export2Freq is set, then the proper values will already be in
	/// node->nodeEnv. Below, we are checking only if the model is set and not the
	/// frequencies. If the model is set, then we set the model frequencies.
	//////////
	if (!(node->nodeEnv->set_in_clade[__values2Export2Freq__])) {
		// Model frequencies should be used if the model is set.
		if (node->nodeEnv->set_in_clade[__model__] && !isNucModel) {
			node->nodeEnv->values2Export2Freq.assign(numStates,0);
//			SetModelFreqs(node->nodeEnv->model, node->nodeEnv);
			node->branch->rates->SetAAFreqs(node->nodeEnv->model, node->nodeEnv);
		}
		else node->nodeEnv->values2Export2Freq = env->values2Export2Freq;
	}

	if (!(node->nodeEnv->set_in_clade[__model__])) 
		node->nodeEnv->model = env->model;

	if (!(node->nodeEnv->set_in_clade[__invariableSites__]))
		node->nodeEnv->invariableSites = env->invariableSites;

	if (!(node->nodeEnv->set_in_clade[__proportion_invariable__]))
		node->nodeEnv->proportion_invariable = env->proportion_invariable;

	if (!(node->nodeEnv->set_in_clade[__rateHetero__])) {
		node->nodeEnv->rateHetero=env->rateHetero;
		node->nodeEnv->gammaShape = env->gammaShape;
		node->nodeEnv->numCats = env->numCats;
		for (int i = 0; i < MAX_RATE_CATS; i++) 
			node->nodeEnv->catRate[i]=env->catRate[i];
	}
	
	if (!(node->nodeEnv->set_in_clade[__constraintChange__])) 
		node->nodeEnv->constraintChange = env->constraintChange;
}

//////////
/// Primarily deals with invariable sites and site rates. 
//////////
// Other "jobs":
//  * PSEUDOGENE's lineages
//  * Sets site categories for clade changes.
//////////
void initializeBranchRun(
						 inTree *iTree, 
						 TNode *node, 
						 int inNumSites, 
						 seqGenOptions *options
						) 
{
//	SetModel(node->nodeEnv->model, node->nodeEnv, options);
	node->branch->rates->SetModel(node->nodeEnv);

	if (node->nodeEnv->constraintChange == PSEUDOGENE) {
		///////////
		/// Remove all constraints from lineage. It's free-flowing.
		//////////
		// Think about removing this whole "invariableSites" flag. Is it really necessary???
		//////////
		if (node->nodeEnv->invariableSites) {
			node->nodeEnv->invariableSites = false;
			for (vector<Site>::iterator it = node->seq_evo.begin(); it != node->seq_evo.end(); ++it)
				(*it).setInvariableState(NO_CONSTRAINT);
			node->nodeEnv->proportion_invariable = 0;
		}
		switch (node->nodeEnv->rateHetero) {
		case GammaRates:
			node->nodeEnv->gammaShape = 1;
			break;
		case DiscreteGammaRates:
			node->nodeEnv->gammaShape = 1;
			break;
		default:
			break;
		}
		node->nodeEnv->rateHetero = NoRates;
    } else if (node->nodeEnv->invariableSites) {
        for(int i=0; i<inNumSites; i++) {
			if ( (iTree->randomInvariableAssignment) ) {
				// No invariable sites, therefore set the proportion_invariable.
				node->seq_evo.at(i).setInvariableState(IsInvariable(node->nodeEnv->proportion_invariable));
			} else if (node->anc->nodeEnv->invariableSites && node->anc->nodeEnv->proportion_invariable != node->nodeEnv->proportion_invariable) {
				//////////
				/// Preserves invariable regions.
				//////////
				// More invariable sites needed. Add some.
				//////////
				double ratio, r;
				if (node->anc->nodeEnv->proportion_invariable < node->nodeEnv->proportion_invariable) {
					ratio = node->nodeEnv->proportion_invariable - node->anc->nodeEnv->proportion_invariable;
					//////////
					/// Coin tossing to find additional sites that will be held invariable (1 in invar array)
					//////////
					if (node->anc->seq_evo.at(i).returnInvariableState() == NO_CONSTRAINT) {
						r = rndu();
						if (r < ratio) node->seq_evo.at(i).setInvariableState(INVARIABLE);
						else node->seq_evo.at(i).setInvariableState(NO_CONSTRAINT);
					} else {
						if (node->anc->seq_evo.at(i).returnInvariableState() == NO_INDEL 
							|| node->anc->seq_evo.at(i).returnInvariableState() == INVAR_AND_NOINDEL)
							node->seq_evo.at(i).setInvariableState(node->anc->seq_evo.at(i).returnInvariableState());
						else node->seq_evo.at(i).setInvariableState(NO_CONSTRAINT);
					}
				//////////
				/// Fewer invariable sites in des, remove some of invariable sites from anc.
				//////////
				} else if (node->anc->nodeEnv->proportion_invariable > node->nodeEnv->proportion_invariable) {
					ratio = node->nodeEnv->proportion_invariable / node->anc->nodeEnv->proportion_invariable;
					//////////
					/// Coin toss to find states to remove invariability from.
					//////////
					if (node->anc->seq_evo.at(i).returnInvariableState() == INVARIABLE) {
						r = rndu();
						if (r < ratio) node->seq_evo.at(i).setInvariableState(NO_CONSTRAINT);
						else node->seq_evo.at(i).setInvariableState(INVARIABLE);
					} else node->seq_evo.at(i).copyInvariableState(node->anc->seq_evo.at(i));
				} else node->seq_evo.at(i).copyInvariableState(node->anc->seq_evo.at(i));
			} else {
				//////////
				/// Do I need this after the copy constructor???
				//////////
				node->seq_evo.at(i).setInvariableState(node->anc->seq_evo.at(i).returnInvariableState());
			}
        }
    } 

	if (node->nodeEnv->constraintChange != PSEUDOGENE)
		SetCategories(node, inNumSites, options);
}

int testRelation(
				 vector<bool>& clade, 
				 vector<bool>& node
				)
{
	bool all_bits_same = true;
	bool some_bits_same = false;

	if(clade.size() != node.size()) {
		cerr << "Error: Bipartitions in testRelation are not the same size: " << clade.size() << " vs. " << node.size() << endl;
		abort();
	}
	vector<bool>::iterator it2 = node.begin();
	for (vector<bool>::iterator it = clade.begin(); it != clade.end(); ++it, ++it2) {
		if ( (*it) != (*it2) ) all_bits_same = false;
		if ( (*it) && (*it2) ) some_bits_same = true;
	}
	if (all_bits_same) return EXACT;

	if (some_bits_same) {
		vector<bool>::iterator it2 = node.begin();
		for (vector<bool>::iterator it = clade.begin(); it != clade.end(); ++it, ++it2) {
			if ( (*it) && !(*it2) ) return NODE_IS_DES;
			if ( !(*it) && (*it2) ) return NODE_IS_ANC;
		}
	} else return NO_RELATION;
}

////////// InDel deprecated, removed (version 1980, from codon_models/src svn repository) /////////

void initializeMotifPositions(
							  TNode *des
							 ) 
{
	//////////
	/// Inherit the varSites from the ancestor. Does not deal with motifs.
	//////////
	des->Inherit_varSites();

	if ( !(des->nodeEnv->constraintChange == PSEUDOGENE) ) {
		des->evolvingSequence = new Sequence(des->anc->evolvingSequence, des);
		des->InheritMotifSites();
		des->setInvariableArrayPos();
	} else des->clearMotifSites();
	//des->setSitePointers("InitializeMotifPositions");

	//cerr << endl << endl << "Evolve::initializeMotifPositions: DES" << endl << endl;
	//des->Print_Active_Properties();
	//cerr << endl << endl << "Evolve::initializeMotifPositions ENDDES" << endl << endl;
	//des->evolvingSequence->print_sequence();
	//exit(0);
}

bool substitute (
				 TTree *tree, 
				 TNode *des, 
				 int *event_site,
				 int simulation_step_type, 
				 string& event, 
				 bool track,
				 double rateAwayFromSequence
				)
{
	double key;
	int mid;

	if (simulation_step_type == UNIFORMIZATION) {	// UNIFORMIZATION
		//////////
		/// Uniformization approach.
		//////////
		//
		key = (double)rndu() * (des->seq_evo.size() * des->rate_away_site_width);
		mid = key / des->rate_away_site_width;
		key = (key - (mid * des->rate_away_site_width)) / des->rate_away_site_width;
		//////////
	} else if (simulation_step_type == INDEPENDENT_ENDPOINT_CONDITIONED) {
		// Independent sites, special case of Uniformization, where bin_size is set to 1. We also
		// specifically know which site is under consideration, since the state for independence
		// is the site, rather than the sequence, which is different from uniformization and dependence.
		// rateAwayFromSequence points to the beginning of the next position under consideration. Subtract
		// 1 to make it point to the beginning of the current position.
		mid = *event_site;
		key = (double)rndu();
		prev_state = des->seq_evo.at(mid).returnState();
	} else {
		//////////
		/// For TRS, simulation_step_type is the site under consideration.
		//////////
		key = (double)rndu() * rateAwayFromSequence;
//		cerr << "rateAwayFromSequence: " << rateAwayFromSequence << " key: " << key << "  ";
		mid = des->locateSite(&key, simulation_step_type);
		prev_state = des->seq_evo.at(mid).returnState();
//		cerr << "mid: " << mid << " remainder: " << key << endl;
//		des->evolvingSequence->print_sequence();
	}
	if (des->nodeEnv->rateHetero == CodonRates) {
		int codon[3];
		int *return_value;
		vector<Site>::iterator posit = des->seq_evo.begin() + mid;
		if (mid-(3-tree->my_iTree->codon_offset) < 0 || mid + (3-tree->my_iTree->codon_offset) >= des->seq_evo.size()) {
			; // Punt, since I do not keep track of other partitions
		} else if ( (mid+tree->my_iTree->codon_offset) % 3 == 0 ) {
			return_value = &codon[0];
			des->seq_evo.at(mid).doSubstitution(key, simulation_step_type, des, return_value, event);
			codon[1] = (*(posit+1)).returnState();
			codon[2] = (*(posit+2)).returnState();
		} else if ( (mid+tree->my_iTree->codon_offset) % 3 == 1 ) {
			codon[0] = (*(posit-1)).returnState();
			return_value = &codon[1];
			des->seq_evo.at(mid).doSubstitution(key, simulation_step_type, des, return_value, event);
			codon[2] = (*(posit+1)).returnState();
		} else if ( (mid+tree->my_iTree->codon_offset) % 3 == 2 ) {
			codon[0] = (*(posit-2)).returnState();
			codon[1] = (*(posit-1)).returnState();
			return_value = &codon[2];
			des->seq_evo.at(mid).doSubstitution(key, simulation_step_type, des, return_value, event);
		} else {
			cerr << "Huh?? how can there be more than 3 codon positions?" << endl;
			exit(EXIT_FAILURE);
		}
		if ( !Stop_Codon(codon) && track) {
			event.clear();
			event.push_back( stateCharacters.at( (*posit).returnState() ) );
			(*posit).setState(*return_value);
			event.push_back( stateCharacters.at( (*posit).returnState() ) );
			(*(tree->my_iTree->my_tree->global_alignment->locateEvent(des,mid))).modifiers.push_back(siteModifier('s',des->mytipNo, eventNo));
			return true;
		} else return false;
	}
	
	// Commented out tracking because it stops changes from being made during the evolution to equilibrium.
	if (des->seq_evo.at(mid).doSubstitution(key, simulation_step_type, des, NULL, event)) {
		cat_chosen.at(des->seq_evo.at(mid).returnCategory())++;
		(*(tree->global_alignment->locateEvent(des,mid))).modifiers.push_back(siteModifier('s',des->mytipNo, eventNo));
		*event_site = mid;
		changed_site = mid;
		return true;
	} else return false;
}
