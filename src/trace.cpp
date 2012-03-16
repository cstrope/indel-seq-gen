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

eventTrack::eventTrack(int events_to_track) 
{
	init(events_to_track);
}

void eventTrack::init(int events_to_track) 
{
	ID = -1;
	eventTime = -1;
	bipartitionInduced.clear();
	MSA_positions.clear();
	eventType = NO_EVENT;	// insert, delete, subst.
	trackThisEvent = events_to_track;
	size = 0;
}

int eventTrack::Compute_MSA_Positions(TTree *tree, int start_msa_pos) 
{
	int curr_pos = 0; 

	for (vector<globalAlignment>::iterator it = tree->global_arrays.begin(); it != tree->global_arrays.end(); it++) {
		if (ID == (*it).indelNo) MSA_positions.push_back(curr_pos+start_msa_pos);
		if ((*it).action == 'i') curr_pos++;
	}

	return curr_pos-1;
}

string eventTrack::Print_Event() 
{
	stringstream event;
	bool stop = false;

	event << "[" << ID << "," << ((eventType==0)?"I":"D") << ",";
	if (eventTime!=-1) 
		event << eventTime;
	event << ",";
	for (vector<bool>::iterator it = bipartitionInduced.begin(); it != bipartitionInduced.end(); it++)
		event << (*it);
	event << "," << size << ",";
	for (vector<int>::iterator it = MSA_positions.begin(); it != MSA_positions.end(); it++) {
		if (it != MSA_positions.begin()) event << ":";
		if (size != MSA_positions.size()) { cerr << (*it) << ": " << size << " vs " << MSA_positions.size() << endl; stop = true; }
		event << (*it);
	}
	
	if (stop) { cerr << size << " type: " << eventType << endl; exit(EXIT_FAILURE); }
	event << "]" << endl;

	return event.str();
}
