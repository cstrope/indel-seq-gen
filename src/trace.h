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

#ifndef _TRACE_H_
#define _TRACE_H_

using namespace std;

#include "tree.h"
#include "motif.h"
#include <iostream>
#include <sstream>

//////////
/// Definitions
//////////
extern int eventNo;
extern bool forward_simulation;

class eventTrack : private Counter<eventTrack>
{
private:
	void init(int events_to_track);
public:
	using Counter<eventTrack>::howMany;
	double event_occurrence_branch_length;	// The length of the branch upon which the event occurred. For EPC calculations.

	int ID;						// Event counter
	vector<int> MSA_positions;	// Positions in the MSA that event appears
	vector<bool> bipartitionInduced;	// Tree edge of event.
	double eventTime;	// Relative time of event occurrence.
	int eventType;	// insert, delete, substitution, fossil
	string size;	// Event-specific: I/D-indel_length, S-from/to state, F-global_fossil_deposition_time

	class ratesAway : private Counter<ratesAway>
	{
	public:
		using Counter<ratesAway>::howMany;
		long double idot, ij, idot_k__T__, ij_k__T__, Pjk, Pik;
		int	i2k;		// Keeps track of the number of differences between i_z and k_0;
		ratesAway() :
			idot (0),
			ij (0),
			idot_k__T__ (0),
			ij_k__T__ (0),
			i2k (1),
			Pjk (0),
			Pik (0)
		{ }

		void assign (long double idot_in, long double ij_in, long double idot_k__T___in, long double ij_k__T___in, long double Pjk_in, long double Pik_in, int i2k_in);
		void print();
	} Q;

	eventTrack(int eventNo, int action, double atEpochTime, vector<bool>& bipartition, string& event, double BL)
			  : ID(eventNo), bipartitionInduced(bipartition), eventTime(atEpochTime),
			    eventType(action), size(event), event_occurrence_branch_length(BL)
	{ }
	// Below: Reading in substitution history. BL will be filled in during emulation.
	eventTrack(int eventNo, char action, double atEpochTime, vector<bool> bipartition, string event, 
			   int pos, double Qidot, double Qij, double Qidot_k__T__, double Qij_k__T__, double i2k, double Pjk, double Pik)
			  : ID(eventNo), bipartitionInduced(bipartition), eventTime(atEpochTime), size(event), event_occurrence_branch_length(-1)
	{ 
		eventType = ( (action == 'S') ? SUBSTITUTION : ( (action == 'I') ? INSERT : ( (action == 'D') ? DELETE : ( (action == 'F') ? FOSSIL : 42) ) ) );
		if (action == 'S') {
			action = SUBSTITUTION;
		} else if (action == 'D') {
			action = INSERT;
		} else if (action == 'I') {
			action = DELETE;
		} else if (action == 'F') {
			action = FOSSIL;
		} else action = 42;		/// Because it's always the answer.
		MSA_positions.push_back(pos);
		Q.assign (Qidot, Qij, Qidot_k__T__, Qij_k__T__, Pjk, Pik, i2k);
	}
	eventTrack(int events_to_track);

	string Print_Event();
	string write_path_event();
	string Print_Short_Event();
	void assign_Q(TNode *node, Site& event_site, int action, int number_of_differences_between_sequences);
	int Compute_MSA_Positions(TTree *tree, int start_msa_pos);
};

//////////
/// Functions
//////////
bool isEvent(string& str);
vector<eventTrack*> sortEventsByTime( list<eventTrack*> *events );
eventTrack *branch_terminal_event(int which_node, int type, double distance_from_root, vector<bool>& bipartition, double branch_start, double branch_end);

#endif
