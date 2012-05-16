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

#include "trace.h"
#include "forward_simulation.h"

eventTrack::eventTrack(
					   int events_to_track
					  ) 
{
	init(events_to_track);
}

void eventTrack::init(
					  int events_to_track
					 ) 
{
	ID = -1;
	eventTime = -1;
	bipartitionInduced.clear();
	MSA_positions.clear();
	eventType = NO_EVENT;	// insert, delete, subst.
	size.clear();
}

bool isEvent ( string& str )
{
	if (str.empty()) return false;

	if (str.at(0) == '[' && str.at(1) != '-') {	// It is an event.
		str.erase(0,1);	// Erase the '['
		str.erase(str.size()-1);	// Erase the ']'
		return true;
	}
	
	return false;
}

vector<eventTrack*> sortEventsByTime( list<eventTrack*> *events ) 
{
	vector<eventTrack*> sorted_events;
	double curr_event_time = -1;
	double next_lowest_event_time = 100000;
	eventTrack *selected_event;
	bool success = false;
	double T;

	for (list<eventTrack*>::iterator it = (*events).begin(); it != (*events).end(); ++it) {
		next_lowest_event_time = 100000;
		success = false;
		for (list<eventTrack*>::iterator jt = (*events).begin(); jt != (*events).end(); ++jt) {
			if ( (*jt)->eventTime > curr_event_time && (*jt)->eventTime < next_lowest_event_time ) {
				next_lowest_event_time = (*jt)->eventTime;
				selected_event = (*jt);
				T = (*jt)->eventTime;
				success = true;
			}
		}
		curr_event_time = selected_event->eventTime;
		if (success) sorted_events.push_back(selected_event);
	}

	(*events).clear();
	for (vector<eventTrack*>::iterator it = sorted_events.begin(); it != sorted_events.end(); ++it) {
		(*events).push_back(*it);
	}


	return sorted_events;
}

void eventTrack::assign_Q(
						  TNode *node,
						  Site& event_site,
						  int action,
						  int number_of_differences_between_sequences
						 )
{
	bool profile=false;

	assert (action == SUBSTITUTION);

	//////////
	/// "Forward" rates away:
	/// * These are not summed during the rate away calculation, so we need to sum them for idot.
	/// * The current state of the site represents the change that happened, thus we just use the
	///   current state to retrieve the ij rates. 
	/// "End-point Conditioned" rates away:
	/// * These are summed, so we simply retrieve the back value for idot_k__T__. Need to subtract
	///   the previous value to get the if_k__T__.
	//////////

	if (profile) {
		cerr << "assign_Q::Forward rate away: ";
		for (vector<double>::iterator it = event_site.forward_rate_away.begin(); it != event_site.forward_rate_away.end(); ++it)
			cerr << (*it) << "  ";
		cerr << endl;
		cerr << "Site " << event_site.returnState() << " rate away chosen: " << event_site.forward_rate_away.at(event_site.returnState()) << endl;
	}

	Q.assign(
			 node->evolvingSequence->Qidot,	// idot
			 event_site.forward_rate_away.at(event_site.returnState()),	// ij
			 node->evolvingSequence->Qidot_k__T__,	// idot_k__T__
			 event_site.return_epc_ij(),	// ij_k__T__
			 ( (!event_site.Pjk0.empty()) ? event_site.Pjk0.at(event_site.returnState()) : 0 ),	// Pjk
			 ( (!event_site.Pjk0.empty()) ? event_site.Pik0 : 0 ),	// Pik
			 number_of_differences_between_sequences	// diff between i_z and k_0
			);
}

int eventTrack::Compute_MSA_Positions(
									  TTree *tree, 
									  int start_msa_pos
									 ) 
{
	int curr_pos = 0; 
	for (vector<insertSite>::iterator it = tree->global_alignment->insert_sites.begin(); it != tree->global_alignment->insert_sites.end(); ++it, curr_pos++) {
		for (list<siteModifier>::iterator jt = (*it).modifiers.begin(); jt != (*it).modifiers.end(); ++jt) {
			if (ID == (*jt).indelNo) {
				if (!MSA_positions.empty() && MSA_positions.back() != ID) continue;
				MSA_positions.push_back(curr_pos+start_msa_pos);
			}
		}
		if (ID == (*it).indelNo) MSA_positions.push_back(curr_pos+start_msa_pos);
	}

	return curr_pos;
}

string eventTrack::Print_Short_Event()
{
	stringstream short_event;
	
	if (eventType == INSERT) short_event << "I";
	else if (eventType == DELETE) short_event << "D";
	else if (eventType == SUBSTITUTION) short_event << "S";
	else if (eventType == FOSSIL) short_event << "F";
	else if (eventType == BRANCH_BEGIN) short_event << "B";		// Beginning of branch for MCMC. //
	else if (eventType == BRANCH_END) short_event << "E";		// END of branch for MCMC. //
	else short_event << "X";
	short_event << ",";
	if (eventTime!=-1) 
		short_event << eventTime;
	short_event << ",";
	for (vector<bool>::iterator it = bipartitionInduced.begin(); it != bipartitionInduced.end(); ++it)
		short_event << (*it);
	short_event << "," << size << ",";
	if (eventType != FOSSIL && eventType != BRANCH_BEGIN && eventType != NO_EVENT && eventType != BRANCH_END) {
		for (vector<int>::iterator it = MSA_positions.begin(); it != MSA_positions.end(); ++it) {
			if (it != MSA_positions.begin()) short_event << ":";
			short_event << (*it);
		}
	} else short_event << "X";
	
	return short_event.str();
}

string eventTrack::Print_Event() 
{
	stringstream event;

	event << "[" << ID << ",";
	event << Print_Short_Event();

	//////////
	/// Print out the elements of the Q matrix:
	///   1. Q_{i.}
	///   2. Q_{ij}
	///   3. Q_{i.|k}(t)
	///   4. Q_{ij|k}(t)
	//////////
	event << "," << Q.idot
		  << "," << Q.ij
		  << "," << Q.idot_k__T__
		  << "," << Q.ij_k__T__
		  << "," << Q.Pjk
		  << "," << Q.Pik
		  << "," << Q.i2k;
	event << "]" << endl;

	return event.str();
}

string eventTrack::write_path_event()
{
	stringstream event;
	
	event << eventTime << ",";
	if (eventType != FOSSIL && eventType != BRANCH_BEGIN && eventType != NO_EVENT && eventType != BRANCH_END) {
		for (vector<int>::iterator it = MSA_positions.begin(); it != MSA_positions.end(); ++it) {
			if (it != MSA_positions.begin()) event << ":";
			event << (*it);
		}
	} else event << "X";

	event << "," << size 
		  << "," << Q.idot 
		  << "," << Q.idot_k__T__ << endl;

	return event.str();
}

void 
eventTrack::ratesAway::assign (
									long double idot_in, 
									long double ij_in, 
									long double idot_k__T___in, 
									long double ij_k__T___in, 
									long double Pjk_in, 
									long double Pik_in, 
									int i2k_in
								   )
{
	idot = idot_in;
	ij = ij_in;
	idot_k__T__ = idot_k__T___in;
	ij_k__T__ = ij_k__T___in;
	Pjk = Pjk_in;
	Pik = Pik_in;
	i2k = i2k_in;
}

void 
eventTrack::ratesAway::print()
{
	cerr << "Q values: " << endl;
	cerr << "  Qi.: " << idot << endl;
	cerr << "  Qij: " << ij << endl;
	cerr << "  Qi.|k(T-t): " << idot_k__T__ << endl;
	cerr << "  Qij|k(T-t): " << ij_k__T__ << endl;
	cerr << "  Pjk: " << Pjk << endl;
	cerr << "  Pik: " << Pik << endl;
	cerr << "  ~#iz!=k0: " << i2k << endl;

}

eventTrack *branch_terminal_event(
								  int which_node, 
								  int type, 
								  double distance_from_root, 
								  vector<bool>& bipartition, 
								  double branch_start, 
								  double branch_end
								 )
{
	eventTrack *return_event;
	string str = to_string(-which_node);
	return_event = new eventTrack(
							      -1,			/// Not really an event, so make it flag-like. ///
							      type,			/// S,I,D,F,B,E
							      ( (type == BRANCH_END) ? branch_end : branch_start ),	/// Sum of paths from ancestors ///
							      bipartition,	/// Who are the ancestors of this node. ///
								  str,
							      branch_end-branch_start
							     ); 
	return return_event;
}
