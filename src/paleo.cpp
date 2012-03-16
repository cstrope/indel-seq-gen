#include "paleo.h"

using namespace std;

int num_fossils = 0;

void paleobiology::doPaleontology(TTree *tree, list<eventTrack*> *events)
{
	//////////
	/// Want to do a minimum of FOSSIL_DEPOSITION_GRANULARITY checks. I think this does it...
	//////////
	paleo_step_size = 1.0 / FOSSIL_DEPOSITION_GRANULARITY;
	DepositFossils(tree->root->branch1, events);
    DepositFossils(tree->root->branch2, events);
}

void paleobiology::DepositFossils(TNode *des, list<eventTrack*> *events)
{
	double t_top, t_bottom;
	double percentage_of_branch, epoch_time;
	double RN;
	string global_time;

	t_top = des->anc->trDistanceFromRoot;
	t_bottom = des->trDistanceFromRoot;

	size_t i = 0;
	for (double dt = paleo_step_size; dt <= t_bottom-t_top; dt += paleo_step_size, i++) {
		RN = (double)rndu();
		if ( RN < paleo_step_size*(root_node_age / fossil_deposition_rate) ) {
			//////////
			/// Calculate the relative time that the event occurs at.
			//////////
			// * anc->trDistanceFromRoot: The scaled distance from the root
			// * (dt/indel_len): Percentage of the branch that has been simulated
			// * des->branch0_time_relative_length: Time rel length of the branch being simulated.
			// * iTree->global_max_path_length: Scalar for the sum of above terms to set time of occurrence between 0 and 1.
			//////////
			num_fossils++;
			epoch_time = (t_top + dt) / max_path_length;
			global_time = to_string(epoch_time * root_node_age);
			eventTrack *new_event;
			new_event = new eventTrack(
								       eventNo++,
									   FOSSIL,
									   epoch_time,
									   des->bipartition,
									   global_time,
									   t_bottom
									  );
  			(*events).push_back(new_event);
		} // else no fossil deposited in this step.
	}

    if (des->tipNo==-1) { 
        DepositFossils(des->branch1, events);
    	DepositFossils(des->branch2, events);
    }
}

template <class T> string to_string (const T& t)
{
	stringstream ss;
	ss << t;
	return ss.str();
}
