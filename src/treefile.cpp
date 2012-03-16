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

#include "treefile.h"
#include "prosite_motif.h"

char treeErrorMsg[256];
int treeError;
extern int run_no;

TNode *avail=NULL;
long usedAvail=0;
long usedMalloc=0;
extern int num_ms_siteProps;

/* functions */

TNode::TNode(TTree *tree) 
	  : evolvingSequence(NULL),
	    branch0(NULL),
	    branch1(NULL),
	    branch2(NULL),
	    anc(NULL),
	    trDistanceFromRoot(0.0),
	    DistanceFromRoot(0.0),
	    tipNo(-1),
		numEventsToSimulate(0),
		BL_step_size(0),
	    atEpochTime(0.0),
	    nodeEnv(NULL),
		branch (new Branch()),
	    mytipNo(-1)
{
	clade_label.clear();
	tree->numNodes++;
	bipartition.clear();
	addGeneral_varSites();
} 

size_t TNode::findForwardNumberOfConstrained(size_t fromSite, size_t numSites) 
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
	for (vector<Site>::iterator site_it = seq_evo.begin()+fromSite; site_it != seq_evo.end() && num_unconstrained < numSites; site_it++, num_unconstrained++) {
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

size_t TNode::findBackwardNumberOfConstrained(size_t fromSite, size_t numSites) 
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
	for (vector<Site>::reverse_iterator site_rit = seq_evo.rbegin()+fromSite; site_rit != seq_evo.rend() && num_unconstrained < numSites; site_rit++, num_unconstrained++) {
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

void TNode::addGeneral_varSites()
{
	variable_region_list.clear();
	one_site_varSite = new varSite(1,1,this);
	unconstrained_varSite = new varSite(0,numeric_limits<int>::max(),this);
}

void TNode::report_varSites()
{
	size_t total_varSites = 0;

	cout << "Address\tminl\tmaxl\tcurrent_num_members " << endl;
	for (list<varSite*>::iterator it = variable_region_list.begin(); it != variable_region_list.end(); it++) {
		cout << (*it) << "\t";
		cout << (*it)->min << "\t" << (*it)->max << "\t" << (*it)->member_set.size() << endl;
		total_varSites += (*it)->member_set.size();
	}

	cout << "The total number of varSites is: " << total_varSites << endl;
}

void TNode::report_on_sites()
{
	for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); site_it++) {
		cerr << "L_ins_ (st=" << (*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite
		     << " m=" << (*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite << ")   ";
		cerr << "R_ins_ (st=" << (*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite
		     << " m=" << (*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite << ") " << endl;
	}
}

varSite *TNode::Generic_varSite()
{
	return one_site_varSite;
}

void TNode::report(TTree *tree) 
{
	for (vector<bool>::iterator it = bipartition.begin(); it != bipartition.end(); it++)
		cerr << (*it);

	if (tipNo != -1 && tree != NULL)
		cerr << " \"" << tree->names.at(tipNo) << "\"";
}

//////////
/// Copies varSites from ancestor to descendant.
//////////
void TNode::Inherit_varSites() 
{
	size_t i = 0;

	anc->one_site_varSite->descendant_equiv = one_site_varSite;
	anc->unconstrained_varSite->descendant_equiv = unconstrained_varSite;
	for (list<varSite*>::iterator it = anc->variable_region_list.begin(); it != anc->variable_region_list.end(); it++, i++) {
		if (i > 1) {
			variable_region_list.push_back(new varSite(*it,this));
		}
	}

	list<varSite*>::iterator anc_it = anc->variable_region_list.begin();
	for (list<varSite*>::iterator it = variable_region_list.begin(); it != variable_region_list.end(); it++, anc_it++) {
		if (anc_it == anc->variable_region_list.end()) {
			cerr << "something's wrong" << endl;
			exit(EXIT_FAILURE);
		}
	}
}

//////////
/// Copy the motif sites. Does this happen by default when copying the sequence from the ancestor?
//////////
void TNode::InheritMotifSites()
{	
	//////////
	/// Create new positions for each. Setting active sites will copy new things into each pos.
	//////////
	// SITE
	vector<Site>::iterator des_site_it = seq_evo.begin();
	for (vector<Site>::iterator anc_site_it = anc->seq_evo.begin(); 
								anc_site_it != anc->seq_evo.end(); 
								anc_site_it++, des_site_it++) 
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
//	cerr << "Ancestor: " << endl;
//	anc->Print_Active_Properties();
//	cerr << "Descendant: " << endl;
//	Print_Active_Properties();
//	exit(0);
}

void TNode::setInvariableArrayPos()
{
	// Nearly all set, but still have to represent the invariable array constraints, applies as motif.
	if (nodeEnv->invariableSites) { 
		vector<Site>::iterator sit = seq_evo.begin();
		for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); site_it++) {
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

bool TNode::isVarSite(varSite *chk_varSite)
{
	if (chk_varSite == NULL) return true;

	for (list<varSite*>::iterator it = variable_region_list.begin(); it != variable_region_list.end(); it++)
		if (chk_varSite == (*it) ) return true;

	return false;
}

void TNode::FULL_REPORT()
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
	for ( ; anc_it != anc->seq_evo.end(); anc_it++, des_it++) {
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

void TNode::clearMotifSites() 
{
	string seq = "";
	for (vector<Site>::iterator site_it = anc->seq_evo.begin(); site_it != anc->seq_evo.end(); site_it++)
		seq += (*site_it).returnState();

	//////////
	/// Rebuild the evolutionaryAttributes array as a non-constrained array with the stored seq
	/// as the initial state.
	//////////
	constructUnconstrainedSequence(seq);
}

string TNode::output_sequence(string::iterator start, string::iterator end)
{
	string seqout = "";
	for (string::iterator it = start; it != end; it++) {
		seqout += stateCharacters[*it];
	}
	return seqout;
}

//////////
/// REMOVE_OBJECTS
//////////
void TNode::Remove_Objects(bool root_node)
{
	size_t i = 0;

	delete seq_evo.front().motif.active_properties.indel->L_ins_;
	
	for (vector<Site>::iterator site_it = seq_evo.begin(); 
							    site_it != seq_evo.end(); 
							    site_it++,i++) 
	{
		for (list<siteProperties*>::iterator it2 = (*site_it).motif.site_props.begin(); 
											 it2 != (*site_it).motif.site_props.end(); 
											 it2++) 
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

void TNode::Remove_varSites()
{
	// variable sites list.
	for (list<varSite*>::iterator it = variable_region_list.begin(); it != variable_region_list.end(); it++) {
		delete (*it);
	}
	variable_region_list.clear();	

}

void TNode::Site_postProcess()
{
	// Set all pointers for the motifSite array and deeper structures.
	for (vector<Site>::iterator it = seq_evo.begin(); it != seq_evo.end(); it++)
		(*it).motif.setNode(this);
}

//////////
/// Old. Used for iSGv2 pub.
//////////
void TNode::calcMotifAccuracy()
{
//	for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); site_it++) {
//		if ( (*site_it).motif.activeProps->fromMotif != NULL ) {
//			total_motif_positions++;
//			if ( (*site_it).motif.activeProps->subst.substitution_bitstring.test((*site_it).returnState()) ) {
//				correct_motif_positions++;
//			}		
//		}
//	}
}

TTree::TTree()
{
	InitTree();
}

void TTree::InitTree()
{	
	treeEnv.clear();
	root=NULL;
	nodeList.clear(); numNodes=0;
	names.clear();
	numTips=0;
	totalLength=0.0;
	rooted=0;
	lengths=-1;
	global_arrays.clear();
}

void TNode::Print_Active_Properties()
{
	int i = 0;
	cerr << endl << "BIPARTITION: " << endl; report(); cerr << endl;
	for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); site_it++, i++) {
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
		for (list<siteProperties*>::iterator jt = (*site_it).motif.site_props.begin(); jt != (*site_it).motif.site_props.end(); jt++) {
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
void TNode::setSitePointers()
{
	for (vector<Site>::iterator site_it = seq_evo.begin(); site_it != seq_evo.end(); site_it++) {
		//////////
		/// For each insertion & deletion pointer, check to make sure it is in this TNode's
		/// varSites.
		//////////
		//
		// -----------------
		// 	Deletion:: varSite *my_sequence_template_varSite;
		(*site_it).motif.active_properties.indel->del->my_sequence_template_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->del->my_sequence_template_varSite,
			"(*site_it).motif.active_properties.indel->del->my_sequence_template_varSite"
		);
		//
		//  Deletion:: varSite *my_motif_varSite;
		(*site_it).motif.active_properties.indel->del->my_motif_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->del->my_motif_varSite,
			"(*site_it).motif.active_properties.indel->del->my_motif_varSite"
		);
		//
		// -----------------
		// Insertion::R_ins_ varSite *my_sequence_template_varSite_left;
		(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left,
			"(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_left"
		);
		//
		// Insertion::R_ins_ varSite*my_sequence_template_varSite_right;
		(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right,
			"(*site_it).motif.active_properties.indel->R_ins_->my_sequence_template_varSite_right"
		);
		//
		// Insertion::R_ins_ varSite *my_motif_varSite_left; 
		(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left,
			"(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_left"
		);
		//
		// Insertion::R_ins_ varSite *my_motif_varSite_right;
		(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right,
			"(*site_it).motif.active_properties.indel->R_ins_->my_motif_varSite_right"
		);
		//
		// Insertion::R_ins_ varSite *on_site_motif_varSite; 
		(*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite
		= checkSitePointer(
			(*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite,
			"(*site_it).motif.active_properties.indel->R_ins_->on_site_motif_varSite"
		);
		//
		// Insertion::R_ins_ varSite *on_site_sequence_template_varSite;	
		(*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite,
			"(*site_it).motif.active_properties.indel->R_ins_->on_site_sequence_template_varSite"
		);
		//
		// -----------------
		// Insertion::L_ins_ varSite *my_sequence_template_varSite_left;
		(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left,
			"(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_left"
		);
		//
		// Insertion::L_ins_ varSite*my_sequence_template_varSite_right;
		(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right,
			"(*site_it).motif.active_properties.indel->L_ins_->my_sequence_template_varSite_right"
		);
		//
		// Insertion::L_ins_ varSite *my_motif_varSite_left; 
		(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left,
			"(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_left"
		);
		//
		// Insertion::L_ins_ varSite *my_motif_varSite_right;
		(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right,
			"(*site_it).motif.active_properties.indel->L_ins_->my_motif_varSite_right"
		);
		//
		// Insertion::L_ins_ varSite *on_site_motif_varSite; 
		(*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite,
			"(*site_it).motif.active_properties.indel->L_ins_->on_site_motif_varSite"
		);
		//
		// Insertion::L_ins_ varSite *on_site_sequence_template_varSite;
		(*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite
		= checkSitePointer( 
			(*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite,
			"(*site_it).motif.active_properties.indel->L_ins_->on_site_sequence_template_varSite"
		);
	}
}

varSite *TNode::checkSitePointer(varSite *site_ptr, string message)
{
	if (! isVarSite(site_ptr) )
		if (anc->isVarSite(site_ptr)) 
			//////////
			/// We have discovered that the current varSite is contained in the ancestral node.
			/// Thus, we reset this varSite to the equivalent varSite in the descendant node.
			/////////
			return site_ptr->descendant_equiv;
//			site_ptr = site_ptr->descendant_equiv;
		else crash(message);

	return site_ptr;
}

void crash(string message)
{
	cerr << "Invalid pointer to a varSite from descendant or ancestor: " << message << endl;
	exit(EXIT_FAILURE);
}

void TTree::constructRootMotifSites(list<inMotif*>& motif_specs, string &root_sequence, string &invariable_in) 
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
	for (list<inClade*>::iterator it = treeEnv.begin(); it != treeEnv.end(); it++) {
		for (list<inMotif*>::iterator it2 = (*it)->my_motifs.begin(); it2 != (*it)->my_motifs.end(); it2++) {
			found = (*it2)->sitemap.find_first_of(".");
			while (found != string::npos) {
				for (list<inClade*>::iterator it4 = treeEnv.begin(); it4 != treeEnv.end(); it4++) {
					for (list<inMotif*>::iterator it3 = (*it4)->my_motifs.begin(); it3 != (*it4)->my_motifs.end(); it3++) {
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
	for (vector<short>::iterator it = dot_additions2template.begin(); it != dot_additions2template.end(); it++) 
		dot_additions_sum += (*it);

	size_t motif_start_position, post_motif_size, post_motif, shouldbezero;
	for (list<inClade*>::iterator it = treeEnv.begin(); it != treeEnv.end(); it++) {
		for (list<inMotif*>::iterator it2 = (*it)->my_motifs.begin(); it2 != (*it)->my_motifs.end(); it2++) {
			int relationship = testRelation( (*it)->bipartition, root->bipartition);
			site_properties = (*it2)->enumerateRegEx(root);
			motif_start_position = (*it2)->sitemap.find_first_not_of("*");
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
					sPit++;
				}
			}

			//////////
			/// Set the vector to the motif site we are updating.
			//////////
			vector<Site>::iterator m2Sit = root->seq_evo.begin() + motif_start_position;
			vector<short>::iterator rm = dot_additions2template.begin();

			// Get all of the siteProps to match the motifSites. Set active props later.
			for (vector<siteProperties*>::iterator it3 = site_properties.begin(); it3 != site_properties.end(); it3++) {
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
				m2Sit++;
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
	for (list<inMotif*>::iterator it2 = motif_specs.begin(); it2 != motif_specs.end(); it2++) {
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
								site_it++) 
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
		sit++, sit2++;
	}

//	cerr << "root_sequence:" << endl << root_sequence << endl;
//	cerr << endl << endl << "ROOT" << endl << endl;
//	root->Print_Active_Properties();
//	cerr << endl << endl << "ENDROOT" << endl << endl;
//	exit(0);
}

//////////
/// Similar, non-global function like constructRootMotifSites, used primarily to create a sequence
/// that is pseudogene'd.
//////////
void TNode::constructUnconstrainedSequence(string &sequence) 
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

list<inMotif*> TTree::DrawMotifs(string &root_sequence, double proportion_motif) 
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
	for (vector<string>::iterator it = prosite_motifs.begin(); it != prosite_motifs.end(); it++) {	
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
			for (list<siteRegEx*>::iterator jt = curr_regex.begin(); jt != curr_regex.end(); jt++) {
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
					for (vector<bool>::iterator kt = occupied.end()-motif_size; kt != occupied.end(); kt++) {
						if (*kt) acceptable = false;
					}
					if (acceptable) {
						acceptable_positions.push_back(root_sequence.size()-motif_size);
						if (exclude_last_position) prosite_library.at(rand_motif)->removeLastPosition();
					}
				} else {
					for (vector<bool>::iterator kt = occupied.begin()+1; kt != occupied.end()-motif_size; kt++, site_num++) {
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
		for (list<inMotif*>::iterator it = motifs_to_place.begin(); it != motifs_to_place.end(); it++) {
			first = (*it)->sitemap.find_first_not_of("*");
			last  = (*it)->sitemap.find_first_of("*", first+1);
			if (last == string::npos) last = (*it)->sitemap.size();
			for (size_t i = first; i < last; i++) occupied.at(i) = true;
		}
	}
	
	for (vector<inMotif*>::reverse_iterator rit = prosite_library.rbegin(); rit != prosite_library.rend(); rit++) 
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
void TTree::Place_Motif(list<siteRegEx*>& motif, char mark, list<inMotif*>& placed_motifs, inMotif *thisMotif, string &root_sequence, size_t start_site) 
{
	string sitemap (root_sequence.size(), '*');

	//////////
	/// Main loop. Cycle through the siteRegEx'es and make changes.
	//////////
	size_t sequence_site = start_site;
	for (list<siteRegEx*>::iterator it = motif.begin(); it != motif.end(); it++) {
		// Going over all of the occupied sites.
		size_t num_added = 0;
		for (list<string>::iterator jt = (*it)->allowable_characters.begin();
									jt != (*it)->allowable_characters.end() && num_added < (*it)->sites_occupied; 
									jt++, sequence_site++, num_added++) 
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
			for (list<inMotif*>::iterator jt = placed_motifs.begin(); jt != placed_motifs.end(); jt++) {
				(*jt)->sitemap.insert(sequence_site, "*");
			}
		}

	}

	//////////
	/// Create new inMotif, copying the input inMotif.
	//////////
	placed_motifs.push_back(new inMotif(thisMotif, sitemap, mark));
}

TNode *TTree::AddNode() 
{
	TNode *node = new TNode(this);
	nodeList.push_back(node);
	
	return node;
}

inClade *TTree::AddClade(string& label) 
{
	inClade *clade = new inClade(label, this);
	treeEnv.push_back(clade);

	return clade;
}

TNode *TTree::ReadTip(string tree_str, int *pos, char ch, TTree *tree, bool first_partition)
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

TNode *TTree::ReadNode(string tree_str, int *pos, TTree *tree, int detectPolytomies, bool first_partition)
{
	TNode *node, *node2;
	char ch;
	int found;

	try {
		node = AddNode();
	} catch (exception& e) {
		cerr << "Error creating TNode node (ReadTree): " << e.what() << endl;
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

TNode *TTree::ReadBranch(string tree_str, int *pos, TTree *tree, bool first_partition)
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

void TTree::ReadTree(string tree_str, int *pos, TTree *tree, bool first_partition)
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
    
    for (list<inClade*>::iterator sit = treeEnv.begin(); sit != treeEnv.end(); sit++) {
    	for (int i = 0; i < names.size(); i++) {
    		if (((names.at(i)).compare((*sit)->clade_name) == 0)) {
				cerr << "Error parsing tree: Clade name '" << names.at(i) << "' is the name of a taxon." << endl;
				exit(EXIT_FAILURE);
			}
    	}
    }
}

void TTree::report_clades() 
{
	cerr << "Clades: " << endl;
	for (list<inClade*>::iterator it = treeEnv.begin(); it != treeEnv.end(); it++) {
		cerr << (*it)->clade_name << " reporting for duty.";
		if (!(*it)->environment_name.empty())
			cerr << " But you can call me " << (*it)->environment_name << ".";
		cerr << endl;
	}
}

char ReadToNextChar(string tree_str, int *pos)
{
	int found;	

	found = tree_str.find_first_not_of(" \t\n",++(*pos));
	*pos = found;

	return tree_str.at(*pos);
}

void ReadUntil(string tree_str, int *pos, char stopChar, char *what)
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

void TTree::report_branches()
{
	for (list<TNode*>::iterator it = nodeList.begin(); it != nodeList.end(); it++) {
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