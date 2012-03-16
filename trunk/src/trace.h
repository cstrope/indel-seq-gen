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

//////////
/// Definitions
//////////
#define NO_EVENT -1
#define INSERTION 0
#define DELETION 1
#define SUBSTITUTION 2

#include "tree.h"
#include <iostream>
#include <sstream>

class eventTrack : private Counter<eventTrack>
{
private:
	void init(int events_to_track);
public:
	using Counter<eventTrack>::howMany;
	int ID;
	vector<int> MSA_positions;
	vector<bool> bipartitionInduced;
	double eventTime;
	int eventType;	// insert, delete, subst.
	int size;
	int trackThisEvent;

	eventTrack(int indelNo, int action, double atEpochTime, vector<bool>& bipartition, int indel_size)
			  : ID(indelNo), bipartitionInduced(bipartition), eventTime(atEpochTime),
			    eventType(action), size(indel_size), trackThisEvent(INSERT | DELETE)
	{ }
	eventTrack(int events_to_track);
	string Print_Event();
	int Compute_MSA_Positions(TTree *tree, int start_msa_pos);
};
#endif
