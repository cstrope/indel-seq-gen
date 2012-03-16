#ifndef _PALEO_H_
#define _PALEO_H_

using namespace std;

#include <vector>
#include <string>
#include <list>
#include "trace.h"
#include "seqGenOptions.h"
#include "evolve.h"

#define FOSSIL_DEPOSITION_GRANULARITY 1024

class seqGenOptions;
class inTree;
class eventTrack;
class TNode;
class TTree;

template <class T> string to_string (const T& t);

class paleobiology : private Counter<paleobiology>
{
	public:
		using Counter<paleobiology>::howMany;
		double paleo_step_size;
		double max_path_length;
	
		paleobiology() 
			: root_node_age (1.0), 
			  fossil_deposition_rate (0.0),
			  paleo_step_size (0.0)
		{ }
		paleobiology(double root_age, double fossil_deposition_rate)
			: root_node_age (root_age),
			  fossil_deposition_rate (fossil_deposition_rate)
		{ }
		void doPaleontology(TTree *tree, list<eventTrack*> *events);
		void DepositFossils(TNode *des, list<eventTrack*> *events);
	
	private:
		double root_node_age;
		double fossil_deposition_rate;
};



#endif
