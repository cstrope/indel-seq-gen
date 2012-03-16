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
#define PRINT_ID 0

using namespace std;

int				node_num = 0;
int				total_events = 0;
int				change_events = 0;
static double 	NeighborRates[400];
int 			indel_fill_aas;
double 			freqRate[MAX_RATE_CATS];
static double 	*matrix[MAX_RATE_CATS];
static double 	*cvector;
int 			indelNo; 
extern int		num_insert, num_delete, len_insert, len_delete;
int 			num_subst = 0;
bool 			insert_debug = false;
bool 			delete_debug = false;
extern int 		num_inserted_positions;
extern int		num_deleted_positions;
string			simulation_status = "";

// functions 

void CreateMatrix() 
{
    for (int i = 0; i < MAX_RATE_CATS; i++) {
		try {
			matrix[i] = new double [numStates * numStates];
		} catch (exception& e) {
			cerr << "Error allocating matrix: " << e.what() << endl;
			exit(EXIT_FAILURE);
		}
    }
	if (!isNucModel) Set_Neighbor_Rates(NeighborRates,indel_fill_aas);
	try {
		cvector = new double [numStates];
	} catch (exception& e) {
		cerr << "Error allocating cvector: " << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

void DestroyMatrix() 
{
    for (int i = 0; i < MAX_RATE_CATS; i++) {
		try {
			delete [] matrix[i];
		} catch (exception& e) {
			cerr << "Error de-allocating matrix: " << e.what() << endl;
			exit(EXIT_FAILURE);
		}
    }
	try {
		delete [] cvector;
	} catch (exception& e) {
		cerr << "Error de-allocating cvector: " << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

void Create_Global_Arrays(TTree *tree, int inNumSites) 
{
	//////////
	/// Root sequence, coming from root node.
	//////////
	tree->global_arrays.clear();
	tree->global_arrays.assign(inNumSites+1, globalAlignment('i',-1,-1));
}

void SetCategories(TNode *node, int inNumSites, seqGenOptions *options)
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
		for (size_t i=0; i<inNumSites; i++) node->seq_evo.at(i).setCategories((int)(rndu()*node->nodeEnv->numCats));
    }

}

char SetState(double *P, string& caller)
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
//		if (!done) cerr << "orig_P = " << *original_P << "    r = " << r << endl; 
	} while (!done);

    return j;
}

short IsInvariable(double proportionInvariable)
{
    double r;

    r=rndu();
    if (r < proportionInvariable) return INVARIABLE;
    else return 0;
}

void RandomSequence(char *seq, int inNumSites, seqGenOptions *options)
{
    char *P;
	string calling_routine = "RandomSequence";

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

string insertFillSequence(int insertSize)
{
    vector<double> bayes_probabilities_for_next (numStates, 0);
    double r;
	string insert_sequence (insertSize+1, '\0');
	string calling_routine;
	
	if (!isNucModel) {
		calling_routine.assign("insertFillSequence, isNucModel");
	    insert_sequence.at(0) = SetState(addFreq, calling_routine);
	    for(size_t i=1; i<insertSize+1; i++) {
			bayes_probabilities_for_next.at(0)=(freq[0]*NeighborRates[(int)insert_sequence.at(i-1)*numStates])/(freq[insert_sequence.at(i-1)]);
	        for(size_t j=1;j<numStates;j++) {
		    	bayes_probabilities_for_next.at(j)=bayes_probabilities_for_next.at(j-1);
		    	bayes_probabilities_for_next.at(j)+=(freq[j]*NeighborRates[(int)insert_sequence.at(i-1)*numStates+j])/(freq[insert_sequence.at(i-1)]);
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

void MutateSequence(inTree *iTree, TNode *des, double inlen, double motif_correction, list<eventTrack*> *events)
{
    int j, cat, pos;
	vector<Site>::size_type i;
    char result, prev_state;
    double len, in_value;
	bool do_sub = false;
	string calling_routine;

    switch (des->nodeEnv->rateHetero) {
	case GammaRates:
		for (i = 0; i < des->seq_evo.size(); i++) {
			pos = i;
			vector<Site>::iterator posit = des->seq_evo.begin()+pos;
			len = inlen + motif_correction;
			in_value = des->seq_evo.at(pos).returnGamma() * len;
			SetVector(cvector, (*posit).returnState(), in_value);
		    if (des->nodeEnv->invariableSites) {
				if ( !des->seq_evo.at(pos).motif.active_properties.subst->siteInvariable() ) {
					calling_routine.assign("MutateSequence, GammaRates, invariableSites");
					result=SetState(cvector, calling_routine);
					// If not part of the motif, then change regardless of result value.
					if ( !(des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test((*posit).returnState())) ) {
						(*posit).setState(result);
					// Otherwise, if it is part of the motif, and it is another accepting value, accept.
					} else if (des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test(result) ) {
						(*posit).setState(result);
					}
					//else cerr << pos << ": Rejected " << stateCharacters[result] << " at position accepting " << des->seq_evo.at(pos).motif.active_properties.subst->report_bitset() << endl;
		    	}

		    } else {
				calling_routine.assign("MutateSequence, GammaRates");
				prev_state = (*posit).returnState();
				result=SetState(cvector, calling_routine);
				if (result != prev_state) change_events++; total_events++;
				if ( !(des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test((*posit).returnState())) ) {
					(*posit).setState(result);
				} else if (des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test(result) ) {
					(*posit).setState(result);
				}
			}
	    }
	    break;

	case DiscreteGammaRates:
	    for (j=0; j < des->nodeEnv->numCats; j++)
	    	SetMatrix(matrix[j], des->nodeEnv->catRate[j]*(inlen+motif_correction));
	    if (des->nodeEnv->invariableSites) {
			for (i = 0; i < des->seq_evo.size(); i++) {
				pos = i;
				vector<Site>::iterator posit = des->seq_evo.begin()+pos;
				if ( !des->seq_evo.at(pos).motif.active_properties.subst->siteInvariable() ) {
					calling_routine.assign("MutateSequence, DiscreteGammaRates, invariableSites");
					result=SetState(
									matrix[des->seq_evo.at(pos).returnCategories()]+((*posit).returnState()*numStates), 
									calling_routine
								   );
					if ( !(des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test((*posit).returnState())) ) {
						(*posit).setState(result);
					} else if (des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test(result) ) {
						(*posit).setState(result);
					}
				}
			}
	    } else {
			for (i = 0; i < des->seq_evo.size(); i++) {
				pos = i;
				vector<Site>::iterator posit = des->seq_evo.begin()+pos;
				calling_routine.assign("MutateSequence, DiscreteGammaRates");
				result=SetState(
								matrix[des->seq_evo.at(pos).returnCategories()]+((*posit).returnState()*numStates),
								calling_routine
							   );
				if ( !(des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test((*posit).returnState())) ) {
					(*posit).setState(result);
				} else if (des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test(result) ) {
					(*posit).setState(result);
				}
			}
	    }
	    break;

	case CodonRates:
	    for (j=0; j < 3; j++)
	    	SetMatrix(matrix[j], des->nodeEnv->catRate[j]*inlen);
		for (i = 0; i < des->seq_evo.size(); i++) {
			pos = i;
			vector<Site>::iterator posit = des->seq_evo.begin()+pos;
			cat=(pos+iTree->codon_offset)%3;
			calling_routine.assign("MutateSequence, CodonRates");
			result=SetState(
							matrix[cat]+((*posit).returnState() * numStates), 
							calling_routine
						   );
			if ( !(des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test((*posit).returnState())) ) 
				do_sub=true;
			else if (des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test(result) )
				do_sub = true;


			if (do_sub) {
				// Check if stop codon is result of substitution.
				int codon[3];
				if (pos - (3-iTree->codon_offset) < 0 || pos + (3-iTree->codon_offset) >= des->seq_evo.size()) {
					; // Punt, since I do not keep track of other partitions (what about subsequences??)
				} else if ( (pos+iTree->codon_offset) % 3 == 0 ) {
					codon[0] = result;
					codon[1] = (*(posit+1)).returnState();
					codon[2] = (*(posit+2)).returnState();
				} else if ( (pos+iTree->codon_offset) % 3 == 1 ) {
					codon[0] = (*(posit-1)).returnState();
					codon[1] = result;
					codon[2] = (*(posit+1)).returnState();
				} else if ( (pos+iTree->codon_offset) % 3 == 2 ) {
					codon[0] = (*(posit-2)).returnState();
					codon[1] = (*(posit-1)).returnState();
					codon[2] = result;
				} else {
					cerr << "Huh?? how can there be more than 3 codon positions?" << endl;
					exit(EXIT_FAILURE);
				}
				if ( !Stop_Codon(codon) ) {
					(*posit).setState(result);
				}
			}
	    }
	    break;

	case NoRates:
	    SetMatrix(matrix[0], inlen+motif_correction);
	    if (des->nodeEnv->invariableSites) {
			for (i = 0; i < des->seq_evo.size(); i++) {
				pos = i;
				vector<Site>::iterator posit = des->seq_evo.begin()+pos;
				if ( !des->seq_evo.at(pos).motif.active_properties.subst->siteInvariable() ) {
					calling_routine.assign("MutateSequence, NoRates, invariableSites");
					result=SetState(
									matrix[0]+((*posit).returnState() * numStates), 
									calling_routine
								   );
					if ( !(des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test((*posit).returnState())) ) {
						(*posit).setState(result);
					} else if (des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test(result) ) {
						(*posit).setState(result);
					}
				}
			}
	    } else {
			for (i = 0; i < des->seq_evo.size(); i++) {
				pos = i;
				vector<Site>::iterator posit = des->seq_evo.begin()+pos;
				prev_state = (*posit).returnState();
				calling_routine.assign("MutateSequence, NoRates");
		    	result=SetState(
		    					matrix[0]+(des->seq_evo.at(pos).returnState()*numStates), 
		    					calling_routine
		    				   );
				if (result != prev_state) change_events++; total_events++;
				if ( !(des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test((*posit).returnState())) ) {
					(*posit).setState(result);
				} else if (des->seq_evo.at(pos).motif.active_properties.subst->substitution_bitstring.test(result) ) {
					(*posit).setState(result);
				}

			}
	    }
	    break;
    }
}

bool Stop_Codon(int *codon) 
{
        size_t T = stateCharacters.find("T"), 
        	   A = stateCharacters.find("A"), 
        	   G = stateCharacters.find("G");
                                
        if (codon[0] == T && codon[1] == A && codon[2] == G) return true;
        if (codon[0] == T && codon[1] == G && codon[2] == A) return true;
        if (codon[0] == T && codon[1] == A && codon[2] == A) return true;
                                
        return false;
}

int calcTrials(TNode *des, int *E_in_, int *E_del_)
{
	*E_in_ = *E_del_ = 0;
	if(des->nodeEnv->invariableSites) {
		for (vector<Site>::iterator it = des->seq_evo.begin(); it != des->seq_evo.end(); it++) {
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

int  Find_Action(double len, double P_insert_, double P_delete_, int CnB_divisor, int action_test) 
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

void Find_Motif_Positions(TTree *tree, TNode *des, motifSite *this_site, int indel_size, bool back, varSite **in_template_varSite, varSite **in_motif_varSite)
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

int Insert(TTree *tree, TNode *des, int indel_size) 
{
    int start_aa, final_aa;
    int indel_aa_in_globals = 0;
    int position_in_globals = 0;

	//////////
	/// When converting from sequence to evolvingSequence, need to change these functions to the
	/// proper evolvingSequence data. For now, ignoring, try to make evolvingSequence equivalent.
	//////////
    final_aa=des->seq_evo.size();
    start_aa=Find_Start_Location(des,&indel_size,des->seq_evo.size(),INSERT);

	if (start_aa == FSL_FAIL) return 0;
	// For insertion, working on half sites, which in turn means we have to go to the next site
	// to represent an insertion correctly. 
	start_aa++;

	//cerr << endl << endl << "INSERT size " << indel_size << " at " << start_aa << endl;

    if(start_aa > final_aa) {
		fprintf(stderr,"start_aa=%d\t+\tindel_size=%d\t>\tfinal_aa=%d\n",start_aa,indel_size,final_aa);
		exit(EXIT_FAILURE);
    }

    while(indel_aa_in_globals <= start_aa) { 
        if(position_in_globals > tree->global_arrays.size()) { 
            fprintf(stderr,"INSERT: Over g_a_end.  indel_aa_in_globals = %d,\tstart_aa = %d\n", indel_aa_in_globals,start_aa);
			cerr << "tree->global_arrays.size() " << tree->global_arrays.size() << " position_in_globals " << position_in_globals << endl;
            exit(EXIT_FAILURE); 
        } 

        if(isSpot(tree,&position_in_globals,des)) { indel_aa_in_globals++; } 
        position_in_globals += 1; 
    } 

    while(indel_aa_in_globals <= start_aa+indel_size) indel_aa_in_globals++;

	tree->global_arrays.insert(
							   tree->global_arrays.begin()+position_in_globals, 
							   indel_size, 
							   globalAlignment('i', des->mytipNo, indelNo)
							  );

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
//	cerr << "MOTIF (" << in_motif_varSite 
//	     << ") TEMPLATE(" << in_template_varSite << ")" 
//	     << " NODE(" << des << ")"
//	     << endl;
//	cerr << "Ancestral and Descendant varSites: " << endl;
//   list<varSite*>::iterator des_it =  des->variable_region_list.begin(); 
//	cerr << "PRE-INSERT:" << endl;
//	cerr << "anc length: " << des->anc->seq_evo.size() << "  des length: " << des->seq_evo.size() << endl;
//	for (list<varSite*>::iterator anc_it =  des->anc->variable_region_list.begin();
//								  anc_it != des->anc->variable_region_list.end(); 
//							  	  anc_it++,
//							  	  des_it++) 
//	{
//		cerr << (*anc_it) << "->" << (*anc_it)->descendant_equiv << " " << (*anc_it)->min << "," << (*anc_it)->max << " [" << (*anc_it)->member_set.size() << "]" << "\t\t"
//			 << (*des_it) << " " << (*des_it)->min << "," << (*des_it)->max << " [" << (*des_it)->member_set.size() << "]" <<endl;
//	}
	if (!des->isVarSite(in_template_varSite)) {
		cerr << "template not varsite." << endl; 
		if (des->isVarSite(in_template_varSite->descendant_equiv)) {
			cerr << "template->descendant_equiv is a varSite." << endl;
			in_template_varSite = in_template_varSite->descendant_equiv;
		} else {
			cerr << "template->descendant_equiv is NOT a varSite" << endl;
			exit(1);
		}
	}
	if (!des->isVarSite(in_motif_varSite)) { 
		cerr << "motif not varsite." << endl;  
		if (des->isVarSite(in_motif_varSite->descendant_equiv)) {
			cerr << "motif->descendant_equiv is a varSite." << endl;
			in_motif_varSite = in_motif_varSite->descendant_equiv;
		} else {
			cerr << "motif->descendant_equiv is NOT a varSite" << endl;
			exit(1);
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
	//cerr << "INSERT! size " << indel_size << "(" << des->seq_evo.size() << ")" << endl;
	auto_ptr<Sequence> insert_sites ( new Sequence(des, indel_size)	);
	insert_sites->init(
					   des,
			  		   insertFillSequence(indel_size), 
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
	num_inserted_positions += indel_size;
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
									site_it++) 
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
	des->setSitePointers();

//	cerr << "INSERT size " << indel_size << " at " << start_aa << endl;
//	des->Print_Active_Properties();
//	cerr << "ENDINSERT" << endl;
//	exit(0);
	return 1;
}

int Delete(TTree *tree,TNode *des, int *indel_size) 
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

 	// Find corresponding position in globals...
    position_in_globals = 0;
    num_indels_in_globals = 0;
 
    while(indel_aa_in_globals <= start_aa) { 
        if(position_in_globals > tree->global_arrays.size()) {
            printf("DELETE: Over g_a_end.  indel_aa_in_globals = %d,\tstart_aa = %d\tfinal_aa = %d\tdelete_size = %d\n", 
		            indel_aa_in_globals,start_aa,final_aa,*indel_size);
		    cerr << "position_in_globals = " << position_in_globals 
		    	 << " tree->global_arrays.size() = " << tree->global_arrays.size() << endl;
            exit(EXIT_FAILURE);
        }

        if(isSpot(tree,&position_in_globals,des)) { indel_aa_in_globals++; }
        position_in_globals += 1;
    }

    num_indels_in_globals = position_in_globals;
    while(indel_aa_in_globals <= start_aa+*indel_size) {
        if(isSpot(tree,&num_indels_in_globals,des)) { indel_aa_in_globals++; num_rounds++; }
        num_indels_in_globals += 1; 
    } 

    temp_in_globals = position_in_globals;
    for(size_t i = position_in_globals; i < num_indels_in_globals+num_rounds; i++) { 
        temp_in_globals = i;
        if(isSpot(tree,&temp_in_globals,des)) { 
			tree->global_arrays.insert(
									   tree->global_arrays.begin()+i,
									   1,
									   globalAlignment('d', des->mytipNo, indelNo)
									  );
            temp_in_globals += 1;
        }    
        i = temp_in_globals;
    }

    if(delete_debug) Print_Globals(tree); 

	//cerr << "*********************DELETE size " << *indel_size << " at " << start_aa << endl;

	num_deleted_positions += *indel_size;

	vector<Site>::iterator site_it = des->seq_evo.begin()+start_aa;
	for (size_t num_dels = 0; num_dels != *indel_size; num_dels++) {
		if (site_it != des->seq_evo.end()-1 )
			(*site_it).motif.Site_deleteMerge( (&(*(site_it-1)).motif), (&(*(site_it+1)).motif), true );
		else
			(*site_it).motif.Site_deleteMerge( (&(*(site_it-1)).motif), NULL, true);
		des->seq_evo.erase(site_it);
	}

	return 1;
} 

int Find_Start_Location(TNode *des, /*TTree *tree,*/ int *indel_size, int final_aa, int action) 
{
    int position;
    int i,j; 
    int ok; 
    int numAcceptablePositions = 0;
    vector<bool> acceptable_positions (final_aa+2, false);
	vector<bool>::iterator bit = acceptable_positions.begin();
	size_t first_constrained = 0, last_constrained = 0;

	des->setSitePointers();
    if(action == INSERT) {
		bit = acceptable_positions.begin();
		numAcceptablePositions = 0;
		for (vector<Site>::iterator site_it = des->seq_evo.begin(); site_it != des->seq_evo.end(); site_it++, bit++) {
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
		// after this position. It will make litte difference to the overall routine.
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
		for (vector<Site>::iterator site_it = des->seq_evo.begin(); site_it != des->seq_evo.end(); site_it++, pos++) {
			// Template array fillage.
			current_varSite = NULL;
			current_indel_size = 0;
			success = true;
			for (vector<Site>::iterator site_it2 = site_it; site_it2 != des->seq_evo.end() && current_indel_size != *indel_size; site_it2++, current_indel_size++) {
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
		for (vector<Site>::iterator site_it = des->seq_evo.begin(); site_it != des->seq_evo.end(); site_it++, pos++) {
			// Motif array fillage.
			current_varSite = NULL;
			current_indel_size = 0;
			success = true;
			for (vector<Site>::iterator site_it2 = site_it; site_it2 != des->seq_evo.end() && current_indel_size != *indel_size; site_it2++, current_indel_size++) {
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
		for (vector<bool>::iterator it = st_accept.begin(); it != st_accept.end(); it++, it2++, bit++) {
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
		num_template_violations++;
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
					//cerr << "Error in the routine for finding indels at N- and C-terminus." << endl;
					exit(EXIT_FAILURE);
				} else {
					int x = i - (position + first_constrained);
					j = des->seq_evo.size()-x;
					*indel_size = x;
				}
			}
		}
	}

    return j; 
}
                     
int Find_Indel_Size(vector<double>& dist, int max_indel) 
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
     
int  isAnc(TNode *thisNode, int anc) 
{
    if(thisNode -> mytipNo == anc) {
        return 1; 
    } else if(thisNode -> mytipNo == -1) {
        return 0; 
    } else { 
        return isAnc(thisNode->anc, anc);
    } 
} 
             
int  isSpot(TTree *tree, int *pos, TNode *curr_taxa) 
{
    int new_pos;
	int save_pos = *pos;

    while(tree->global_arrays.at(*pos).action != 'i' && *pos < tree->global_arrays.size())
    	*pos+=1;

    if(!isAnc(curr_taxa,tree->global_arrays.at(*pos).fromAnc)) { return 0; }
    else { 
        new_pos = *pos-1;
        while(new_pos > 0 && tree->global_arrays.at(new_pos).action != 'i') { new_pos--; }
         
        new_pos++;
             
        if(new_pos == *pos) { return 1; } 
        else {
            while(tree->global_arrays.at(new_pos).action == 'd') {
                if(isAnc(curr_taxa,tree->global_arrays.at(new_pos).fromAnc)) { return 0; } 
                new_pos++;
            } 
            return 1; 
        } 
    } 

} 

//////////
/// This changes all environment parameters for Nodes (no motif changes).
//////////
void initializeCladeParameters(TTree *tree, TNode *node) 
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
			SetModelFreqs(node->nodeEnv->model, node->nodeEnv);
		}
		else {
			node->nodeEnv->values2Export2Freq = env->values2Export2Freq;
		}
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
void initializeBranchRun(inTree *iTree, TNode *node, double *en_invar_scale, int inNumSites, seqGenOptions *options) 
{
	SetModel(node->nodeEnv->model, node->nodeEnv);

	if (node->nodeEnv->constraintChange == PSEUDOGENE) {
		///////////
		/// Remove all constraints from lineage. It's free-flowing.
		//////////
		// Think about removing this whole "invariableSites" flag. Is it really necessary???
		//////////
		if (node->nodeEnv->invariableSites) {
			node->nodeEnv->invariableSites = false;
			for (vector<Site>::iterator it = node->seq_evo.begin(); it != node->seq_evo.end(); it++)
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
			if ( !(iTree->randomInvariableAssignment) ) {
				//////////
				/// Do I need this after the copy constructor???
				//////////
				node->seq_evo.at(i).setInvariableState(node->anc->seq_evo.at(i).returnInvariableState());
			} else if (node->anc->nodeEnv->invariableSites) {
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
				// No invariable sites, therefore set the proportion_invariable.
				node->seq_evo.at(i).setInvariableState(IsInvariable(node->nodeEnv->proportion_invariable));
			}
        }

		//////////
		/// Invalidates above.
		//////////
		if ( !(iTree->randomInvariableAssignment) ) {
			int num_invar = 0;
			for (int i = 0; i < inNumSites; i++) {
				if (node->seq_evo.at(i).returnInvariableState() == INVARIABLE 
					|| node->seq_evo.at(i).returnInvariableState() == INVAR_AND_NOINDEL)	
					num_invar++;
			}
			node->nodeEnv->proportion_invariable = (double)num_invar / (double)inNumSites;
		}
    } 

	if (node->nodeEnv->constraintChange != PSEUDOGENE)
		SetCategories(node, inNumSites, options);
}

//////////
/// Maybe get rid of this function, instead "retrying" whenever a substitution fails, rather
/// than scaling up the substitution rate for each site.
//////////
void initializeTreeScalingParameters(inTree *iTree, double *scale, double *invar_scale) 
{
	// Calc scale
	*scale = iTree->partitionRate;

	// Calc invar_scale.
	if (iTree->my_tree->treeEnv.front()->invariableSites) {
		if (iTree->randomInvariableAssignment) {
			if (iTree->rootSeqType == RANDOM) {
				for (int j=0; j<iTree->partitionLength; j++) {
   					iTree->my_tree->root->seq_evo.at(j).setInvariableState(IsInvariable(iTree->proportion_invariable)); 
				}
   				for (int j=iTree->partitionLength; j<iTree->partitionLength; j++) 
   					iTree->my_tree->root->seq_evo.at(j).setInvariableState(0);
			} else {
				cerr << "Didn't think I could get to this spot:: InitializeTreeParameters, iTree->randomInvariableAssignment, iTree->rootSeqType != RANDOM." << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			if (!(iTree->rootSeqType == RANDOM)) {
				for (int j=0; j<iTree->partitionLength; j++) 
					if (iTree->my_tree->root->seq_evo.at(j).returnInvariableState() == INVARIABLE 
						|| iTree->my_tree->root->seq_evo.at(j).returnInvariableState() == INVAR_AND_NOINDEL)
						(iTree->proportion_invariable)++;
				iTree->proportion_invariable /= iTree->partitionLength;
			} else {
				cerr << "Didn't think I could get to this spot:: InitializeTreeParameters, !(iTree->randomInvariableAssignment), iTree->rootSeqType == RANDOM." << endl;
				exit(EXIT_FAILURE);
			}
		}
		*invar_scale = (1.0 - iTree->my_tree->treeEnv.front()->proportion_invariable);
	} else *invar_scale = 1;
}

void calculateMotifScale(inTree *iTree, TNode *des, seqGenOptions *options, double subst_len, double *motif_scale)
{
	size_t current_position;

	current_position = 0;

	if (des->seq_evo.size() > 1)
		for (vector<Site>::iterator site_it = des->seq_evo.begin()+1; site_it != des->seq_evo.end(); site_it++, current_position++)
			*motif_scale += CalculateMotifScale(des, current_position, subst_len, (*site_it).motif.active_properties.subst->substitution_bitstring);

	*motif_scale /= des->seq_evo.size();
}

double CalculateMotifScale(TNode *des, int position, double subst_len, bitset<20>& motif_characters)
{
	double scale = 0;
	double position_score, in_val;
	int cat;
	int num_pos = 0;

	switch (des->nodeEnv->rateHetero) {
	case GammaRates:
		for (int i = 0; i < numStates; i++) {
			position_score = 0;
			if (motif_characters.test(i)) {
				num_pos++;
				in_val = des->seq_evo.at(position).returnGamma() * subst_len;
				SetVector(cvector, i, in_val);
				for (int j = 0; j < numStates; j++)
					if (!motif_characters.test(j))
						scale += ( (j == (numStates-1)) ? (cvector[j+1]-cvector[j]) : (1-cvector[j]) );
			}
		}
	    break;

	case DiscreteGammaRates:
		cat = des->seq_evo.at(position).returnCategories();
		SetMatrix(matrix[cat], des->nodeEnv->catRate[cat]*subst_len);
		for (int i = 0; i < numStates; i++) {
			position_score = 0;
			if (motif_characters.test(i)) {
				num_pos++;
				for (int j = 0; j < numStates; j++)
					if (!motif_characters.test(j))
						scale += ( (j == (numStates-1)) ? (1-matrix[cat][i*numStates+j]) : (matrix[cat][i*numStates+j+1]-matrix[cat][i*numStates+j]) );
			}
		}
	    break;

	case CodonRates:
		cat = position % 3;
		SetMatrix(matrix[cat], des->nodeEnv->catRate[cat]*subst_len);
		for (int i = 0; i < numStates; i++) {
			position_score = 0;
			if (motif_characters.test(i)) {
				num_pos++;
				for (int j = 0; j < numStates; j++)
					if (!motif_characters.test(j))
						scale += ( (j == (numStates-1)) ? (1-matrix[cat][i*numStates+j]) : (matrix[cat][i*numStates+j+1]-matrix[cat][i*numStates+j]) );
			}
		}
	    break;

	case NoRates:
	    SetMatrix(matrix[0], subst_len);
		for (int i = 0; i < numStates; i++) {
			position_score = 0;
			if (motif_characters.test(i)) {
				num_pos++;
				for (int j = 0; j < numStates; j++)
					if (!motif_characters.test(j))
						scale += ( (j == (numStates-1)) ? (1-matrix[0][i*numStates+j]) : (matrix[0][i*numStates+j+1]-matrix[0][i*numStates+j]) );
			}
		}
	    break;
	} 

	return scale / (double)num_pos;
}

int testRelation(vector<bool>& clade, vector<bool>& node)
{
	bool all_bits_same = true;
	bool some_bits_same = false;

	if(clade.size() != node.size()) {
		cerr << "Error: Bipartitions in testRelation are not the same size: " << clade.size() << " vs. " << node.size() << endl;
		abort();
	}
	vector<bool>::iterator it2 = node.begin();
	for (vector<bool>::iterator it = clade.begin(); it != clade.end(); it++, it2++) {
		if ( (*it) != (*it2) ) all_bits_same = false;
		if ( (*it) && (*it2) ) some_bits_same = true;
	}
	if (all_bits_same) return EXACT;

	if (some_bits_same) {
		vector<bool>::iterator it2 = node.begin();
		for (vector<bool>::iterator it = clade.begin(); it != clade.end(); it++, it2++) {
			if ( (*it) && !(*it2) ) return NODE_IS_DES;
			if ( !(*it) && (*it2) ) return NODE_IS_ANC;
		}
	} else return NO_RELATION;
}

void InDel(TTree *tree, TNode  *des, double len, list<eventTrack*> *events, int CnB_divisor, seqGenOptions *options) 
{
    int action = INSERT;
    int indel_size;
    int i;
	int E_in_ = 0;
	int E_del_ = 0;
	int total_trials;

	total_trials = calcTrials(des,&E_in_,&E_del_);		// Calculate rounds by examining indel objects.

	bool success;
	for(i=0; i<total_trials; i++) {
		if (des->nodeEnv->rateHetero == CodonRates) {
			if ((double)rndu() < (double)E_in_ / ((double)E_in_+E_del_))
				action = Find_Action(len,des->nodeEnv->P_ins_ / 3.,des->nodeEnv->P_del_ / 3., CnB_divisor, INSERT); 
			else
				action = Find_Action(len,des->nodeEnv->P_ins_ / 3.,des->nodeEnv->P_del_ / 3., CnB_divisor, DELETE); 
		} else {
			if ((double)rndu() < (double)E_in_ / ((double)E_in_+E_del_))
				action = Find_Action(len,des->nodeEnv->P_ins_,des->nodeEnv->P_del_, CnB_divisor, INSERT); 
			else
				action = Find_Action(len,des->nodeEnv->P_ins_,des->nodeEnv->P_del_, CnB_divisor, DELETE); 
		}
		
		if(action == INSERT || action == DELETE) {
			success = false;
	    	if(action == INSERT) {
				des->setSitePointers();
	        	indel_size=Find_Indel_Size(des->nodeEnv->insert_lengthDistribution,des->nodeEnv->maxIndel);	
				if ( (success=Insert(tree,des,indel_size)) ) {
					des->setSitePointers();
					num_insert++;
					len_insert += indel_size;
				}
	    	} else {
				des->setSitePointers();
	        	indel_size=Find_Indel_Size(des->nodeEnv->delete_lengthDistribution,des->nodeEnv->maxIndel);
				if ( (des->seq_evo.size()-indel_size) > 0)
		    		if ( (success=Delete(tree,des,&indel_size)) ) {
						des->setSitePointers();
		    			num_delete++;
						len_delete += indel_size;
					}
	    	}
			if (success && traceEvents) {
				eventTrack *new_event;// = new eventTrack(INSERT | DELETE);
				new_event = new eventTrack(
										   indelNo,
										   action,
										   des->atEpochTime,
										   des->bipartition,
										   indel_size
										  );
				(*events).push_back(new_event);
				indelNo++;
			}
		}
    }
	if (PRINT_ID) {
		cout << des->atEpochTime << " "; 
		cout << num_insert << " " << len_insert << " " << num_delete << " " << len_delete << endl;
	}
}

void initializeMotifPositions(TNode *des) 
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
	des->setSitePointers();

	//cerr << endl << endl << "DES" << endl << endl;
	//des->Print_Active_Properties();
	//cerr << endl << endl << "ENDDES" << endl << endl;
	//des->evolvingSequence->print_sequence();
	//exit(0);
}

void EvolveNode(inTree *iTree, TNode *anc, TNode *des, double scale, double invar_scale, int inNumSites, list<eventTrack*> *events, seqGenOptions *options)
{ 
    double subst_len, indel_len;   
    double en_scale = scale, en_invar_scale = invar_scale;
	double motif_scale = 0;

	ofstream step_out;
	string trs_outfile;

	node_num++;
	if (!options->quiet && node_num % 100 == 0) cerr << "Node " << node_num << endl;
    des->anc=anc;
	des->atEpochTime = anc->atEpochTime;

	initializeMotifPositions(des);
	initializeCladeParameters(iTree->my_tree, des);
	initializeBranchRun(iTree, des, &en_invar_scale, anc->seq_evo.size(), options);

	if (options->simulation_step_type == GILLESPIE) {
		des->numEventsToSimulate = 1;
		subst_len = des->branch->anc_length() * (en_scale / en_invar_scale);
		indel_len = des->branch->anc_length() * en_scale;
	} else {
	    subst_len=des->BL_step_size * (en_scale / en_invar_scale);
		indel_len=des->BL_step_size * en_scale;
	}

	calculateMotifScale(iTree, des, options, subst_len, &motif_scale);
	motif_scale /= des->numEventsToSimulate;

	if (options->simulation_step_type == GILLESPIE) {
		gillespieIndel(iTree, des, indel_len, events, options);
		MutateSequence(iTree, des, subst_len, motif_scale, events);
	} else {
		for (int i = 0; i < des->numEventsToSimulate; i++) {
			des->atEpochTime += iTree->my_tree->epoch_step_size;
	    	MutateSequence(iTree, des, subst_len, motif_scale, events);
	    	if (des->nodeEnv->indelFlag)
		   		InDel(iTree->my_tree,des,indel_len,events, des->numEventsToSimulate, options);
		}
		if (output_trs_example) step_out.close();
	}

//	cerr << "node num " << node_num << endl;
//	des->Print_Active_Properties();
//	cerr << "end node num " << node_num << endl;

    if (des->tipNo==-1) { 
    	des->setSitePointers();
        EvolveNode(iTree, des, des->branch1, scale, invar_scale, inNumSites, events, options);
    	EvolveNode(iTree, des, des->branch2, scale, invar_scale, inNumSites, events, options);
    } else des->calcMotifAccuracy();	// calculator for pub (both iSGv1 & iSGv2)

	des->Remove_Objects();
} 

void EvolveSequences(inTree *iTree, list<eventTrack*> *events, seqGenOptions *options)
{
	double scale, invar_scale;
	num_subst = 0;
	initializeTreeScalingParameters(iTree, &scale, &invar_scale);

	EvolveNode(iTree, iTree->my_tree->root, iTree->my_tree->root->branch1, scale, invar_scale, iTree->partitionLength, events, options);
    EvolveNode(iTree, iTree->my_tree->root, iTree->my_tree->root->branch2, scale, invar_scale, iTree->partitionLength, events, options);
    if (!iTree->my_tree->rooted) {
		cerr << "Unrooted tree??" << endl; 
		exit(EXIT_FAILURE);
        EvolveNode(iTree, iTree->my_tree->root, iTree->my_tree->root->branch0, scale, invar_scale, iTree->partitionLength, events, options);
	}

	//////////
	/// Remove objects from root sequence.
	//////////
	// The following line of code blows iSG up when a motif sequence is followed by a RANDOM,
	// sequence with the -1 option unset.
	if ( iTree->rootSeqType != RANDOM || !iTree->randomInvariableAssignment )
		iTree->my_tree->root->Remove_Objects(true);
	//
	// Regardless of the motif situation, the general varSites still need to be removed. This is
	// done automatically if the if-stmt is true.
	else iTree->my_tree->root->Remove_varSites();
	//
	//////////
} 

void Print_Globals(TTree *tree) 
{ 
    stringstream getlen;
    int fromAnc_len, indelNo_len;
   
    getlen.str("");
	//////////
	/// Print global array positions.
	//////////
    getlen.str("");
	// action
    for(int i=1; i<tree->global_arrays.size(); i++) {
        if(tree->global_arrays.at(i).fromAnc == -1) { fprintf(stderr,"-"); }
        else {
			if(tree->global_arrays.at(i).fromAnc != -1) {
        		getlen << tree->global_arrays.at(i).fromAnc;
        		fromAnc_len = (getlen.str()).size();
        		getlen.str("");
        	} else fromAnc_len = 1;
        	if(tree->global_arrays.at(i).indelNo != -1) {
        		getlen << tree->global_arrays.at(i).indelNo;
        		indelNo_len = (getlen.str()).size();
				getlen.str("");
			} else indelNo_len = 1;

			int numspace = (fromAnc_len > indelNo_len) ? fromAnc_len : indelNo_len;
			if(numspace > 1) 
				for(int j = 1; j < numspace; j++) fprintf(stderr, " ");

            fprintf(stderr,"%c",tree->global_arrays.at(i).action);
        }    
    } 

	// fromAnc
    fprintf(stderr,"\n"); 
    for(int i=1; i<tree->global_arrays.size(); i++) {
        if(tree->global_arrays.at(i).fromAnc == -1) { fprintf(stderr,"-"); }
        else { 
        	if(tree->global_arrays.at(i).indelNo != -1) {
        		getlen << tree->global_arrays.at(i).indelNo;
        		indelNo_len = (getlen.str()).size();
        		getlen.str("");
			} else indelNo_len = 1;

			if(indelNo_len > 1) 
				for(int j = 1; j < indelNo_len; j++) fprintf(stderr, " ");

        	fprintf(stderr,"%d",tree->global_arrays.at(i).fromAnc); 
        }
    }
    fprintf(stderr,"\n");

	// indelNo
	for(int i=1; i<tree->global_arrays.size(); i++) {
        if(tree->global_arrays.at(i).fromAnc == -1) { fprintf(stderr,"-"); }
        else { 
        	if(tree->global_arrays.at(i).indelNo != -1) {
        		getlen << tree->global_arrays.at(i).fromAnc;
        		fromAnc_len = (getlen.str()).size();
        		getlen.str("");
			} else fromAnc_len = 1;

			if(fromAnc_len > 1) 
				for(int j = 1; j < fromAnc_len; j++) fprintf(stderr, " ");

        	fprintf(stderr,"%d",tree->global_arrays.at(i).indelNo); 
        }
    }
	fprintf(stderr,"\n");

} 

double calcGillespieLambda(TNode *des, double *lambda_ins, int *ins_L, double *lambda_del, int *del_L)
{
	double lambda = 0;
	double del_freq_by_size=0, del_freq=0, u_del=0;
	int i = 0;
	for (vector<double>::iterator it = des->nodeEnv->delete_lengthDistribution.begin(); it != des->nodeEnv->delete_lengthDistribution.end(); it++, i++) {
		del_freq += (*it);
		del_freq_by_size += (*it)*i;
	}

	if (del_freq == 0) u_del = 0;
	else u_del = del_freq_by_size / del_freq;

	*lambda_ins = des->nodeEnv->P_ins_;
	*lambda_del = des->nodeEnv->P_del_;

	calcTrials(des, ins_L, del_L);
	lambda = *lambda_ins*(*ins_L) + *lambda_del*(*del_L) + *lambda_del*u_del - *lambda_del + *lambda_ins;

	return lambda;
}
               
void gillespieIndel(inTree *iTree, TNode *des, double indel_len, list<eventTrack*> *events, seqGenOptions *options)
{
	// Gillespie style indels. 
	int ins_L, del_L;
	double lambda_ins, lambda_del, lambda_T, u_del, exp_mean;
	bool success;
	int indel_size = 0, action = INSERT;

   	if (des->nodeEnv->indelFlag) {
		//////////
		/// The Chang and Benner model does not reconcile with the Gillespie model. Run
		/// the C&B model as is done in iSGv1.
		//////////
		if (des->nodeEnv->P_ins_ == 0 && des->nodeEnv->P_del_ == 0) {
			InDel(iTree->my_tree, des, indel_len, events, 1, options);
		} else {
			//////////
			// Indel routine	
			//  * Should codonRates affect the probability of an indel occurring?
			//////////
			lambda_T = calcGillespieLambda(des, &lambda_ins, &ins_L, &lambda_del, &del_L);
			exp_mean = 1.0/lambda_T;
			des->setSitePointers();
			for (double dt = rand_exp(exp_mean); dt <= indel_len; dt += rand_exp(exp_mean)) {
				success = false;
				if ( (double)rndu() < ((lambda_ins*(double)ins_L + lambda_ins) / lambda_T) ) {	// Insertion.
		        	indel_size=Find_Indel_Size(des->nodeEnv->insert_lengthDistribution,des->nodeEnv->maxIndel);	
					action = INSERT;
					des->setSitePointers();
					if ( (success=Insert(iTree->my_tree,des,indel_size)) ) {
						num_insert++;
						len_insert += indel_size;
					}
				} else { // Deletion.
		        	indel_size=Find_Indel_Size(des->nodeEnv->delete_lengthDistribution,des->nodeEnv->maxIndel);
					action = DELETE;
					if ( (des->seq_evo.size()-indel_size) > 0) {
						des->setSitePointers();
			    		if ( (success=Delete(iTree->my_tree,des,&indel_size)) ) {
			    			num_delete++;
							len_delete += indel_size;
						}
					}
				}
				if (success && traceEvents) {
					eventTrack *new_event;
					//////////
					/// Calculate the relative time that the event occurs at.
					//////////
					// * anc->trDistanceFromRoot: The scaled distance from the root
					// * (dt/indel_len): Percentage of the branch that has been simulated
					// * des->branch0_time_relative_length: Time rel length of the branch being simulated.
					// * iTree->global_max_path_length: Scalar for the sum of above terms to set time of occurrence between 0 and 1.
					//////////
					des->atEpochTime = (des->anc->trDistanceFromRoot+(dt/indel_len)*des->branch->branch0_time_relative_length) / iTree->global_max_path_length;
					new_event = new eventTrack(
											   indelNo,
											   action,
											   des->atEpochTime,
											   des->bipartition,
											   indel_size
											  );
					(*events).push_back(new_event);
					indelNo++;
				}
	
				lambda_T = calcGillespieLambda(des, &lambda_ins, &ins_L, &lambda_del, &del_L);
				exp_mean = 1.0/lambda_T;
			}
		}
	}
}
