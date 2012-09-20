#include "dependency.h"

#define INVERSE true

////////////////////
////// EPC simulation: Reads the Markov dependency from a file created by the forward simulation.
////////////////////

//////////
/// CONSTRUCTORS
//////////
// Forward 3rd-order MM simulation, Randomly assigned frequencies and transition probabilities taken to a power
Dependency::Dependency(
					   double dep_order,
					   double dependence_superscript,
					   string& outfile_name_root		// After generating 3MM, Write out file for EPC.
					  )
	: context(dep_order)
{
	cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) IN" << endl;
	context.allocate_lookup_context_vector();
	cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) LOOKUP CONTEXT VECTOR" << endl;
	context.generateDependencies(dependence_superscript);
	cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) GENERATE DEPENDENCIES" << endl;
	context.set_lookup_table();
	cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) SET LOOKUP TABLE " << endl;
	context.setOffset();
	cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) SET OFFSET " << endl;
	context.outputDependencies(outfile_name_root);
	cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) OUT" << endl;
}

// EPC Third-order Markov Model.
Dependency::Dependency(
					   double dep_order,
					   string& file
					  )
	: context(dep_order)
{
	cerr << "Point-> Dependency::Dependency(dep_order, &file) IN" << endl;
	context.allocate_lookup_context_vector();
	cerr << "Point-> Dependency::Dependency(dep_order, &file) allocate_lookup_context_vector" << endl;
	context.readDependencies(file);
	cerr << "Point-> Dependency::Dependency(dep_order, &file) readDependencies" << endl;
	context.set_lookup_table();
	cerr << "Point-> Dependency::Dependency(dep_order, &file) set_lookup_table" << endl;
	context.setOffset();
	cerr << "Point-> Dependency::Dependency(dep_order, &file) OUT" << endl;
}

// Neutral model setting for 3MM.
Dependency::Dependency(
					   double dep_order,
					   inClade *environment
					  )
	: context(dep_order)
{
	context.allocate_lookup_context_vector();
	context.set_neutral_lookup_vector(environment);
}

Dependency::Dependency(
					   double dep_order,	// Still useful...
					   int block_size,		// Dummy variable for distinguishing from prev decl., but the idea would work.
					   string& file			// File holding dependency counts
					  )
	: context(dep_order)
{
	cerr << "Point-> Dependency::Dependency(dep_order, block_size, file) 1" << endl;
	context.allocate_lookup_context_vector();
	cerr << "Point-> Dependency::Dependency(dep_order, block_size, file) 2" << endl;
	context.readDependencies(file);
	cerr << "Point-> Dependency::Dependency(dep_order, block_size, file) 3" << endl;
	//exit(0);
}

void contextDependence::set_neutral_lookup_vector(
												  inClade *environment
												 )
{
	LookUp *lookup;
	double RN, RN_inv;
	
	/// inClade, though passed, means nothing at the moment. Hard-coding JC. ///

	// Beginning of sequence. Set all to 1, hardcoding JC makes all cancel out.
	tuplet_pi.assign(index_position_multiplier.at(order), 1);

	// Set everything else to 1 also.
	for (vector<vector<LookUp*> >::iterator it = lookup_table.begin(); it != lookup_table.end(); ++it) {
		cerr << (*it).size() << endl;
		for (vector<LookUp*>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt) {
			lookup = new LookUp();
			lookup->value = 1;		// Doesn't matter what I set this to. It will cancel out anyways. //
			lookup->inverse = 1;
			(*jt) = lookup;
		}
	}
}

void 
contextDependence::set_lookup_table()
{
	LookUp *lookup;
	double total;
	int i;
	string sequence;
	double RN, RN_inv;

	//////////
	/// Use seed values to construct the lookup table values for patterns not involving the first <order> residues.
	//////////
	string substr_key;
	for (int j = 2; j <= order; j++) {
		for (i = 0; i < index_position_multiplier.at(order+j); i++) {
			sequence = lookup_table_sequence(i, order+j);
			//cerr << " sequence: " << sequence << endl;
			RN = 1;
			for (int k = 1; k <= j; k++) {
				substr_key = sequence.substr(sequence.size() - order - k , order + 1);
				substr_key += "*";
				//cerr << " -> substr_key: " << substr_key << endl;
				RN *= lookup_table2[substr_key]->value;
			}
			lookup = new LookUp();
			lookup->value = RN;
			RN_inv = 1.0 / RN;
			lookup->inverse = RN_inv;
			sequence += "*";	// Marker is still at the end. // 
			lookup_table2.insert( pair<string, LookUp*>(sequence, lookup) );

			//lookup_context.at(lookup_context_index(int2iSG_seq(i, order+j), true)) = lookup;


			lookup_table.at(order+j).at(i) = lookup;

			//cerr << "lookup_table[" << sequence << "] = v(" << lookup_table2[sequence]->value << ") i(" << lookup_table2[sequence]->inverse << ")";
			//cerr << "lookup_context[" << lookup_context_index(int2iSG_seq(i, order+j), true); 
			//cerr << "] = v(" << lookup_context.at(lookup_context_index(int2iSG_seq(i, order+j), true))->value << ") i("; 
			//cerr << lookup_context.at(lookup_context_index(int2iSG_seq(i, order+j), true))->inverse << ")" << endl;

			//cerr << "   " << order+j << ", " << i << "/" << lookup_table.at(order+j).size() << "       ";
			//cerr << lookup_table.at(order+j).at(i)->value << ", " << lookup_table.at(order+j).at(i)->inverse << endl;
		}
	}

	//////////
	/// Fill in the middle of sequence values
	/// Substitute (order+1) for all instances of j, do not add "*" at end of sequence.
	//////////
	for (i = 0; i < index_position_multiplier.at(order+order+1); i++) {
		sequence = lookup_table_sequence(i, order+order+1);
		//cerr << " sequence: " << sequence << endl;
		RN = 1;
		for (int k = 1; k <= order+1; k++) {
			substr_key = sequence.substr(sequence.size() - order - k , order + 1);
			substr_key += "*";
			//cerr << " -> substr_key: " << substr_key << endl;
			RN *= lookup_table2[substr_key]->value;
		}
		lookup = new LookUp();
		lookup->value = RN;
		RN_inv = 1.0 / RN;
		lookup->inverse = RN_inv;
		lookup_table2.insert( pair<string, LookUp*>(sequence, lookup) );

		//lookup_context.at(lookup_context_index(int2iSG_seq(i, order*2+1))) = lookup;

		lookup_table.at(order).at(i) = lookup;
//		cerr << "lookup_table[" << sequence << "] = v(" << lookup_table2[sequence]->value << ") i(" << lookup_table2[sequence]->inverse << ")" << endl;
//		cerr << "lookup_index[" << lookup_context_index(int2iSG_seq(i, order*2+1)) << "] = v(" << lookup_context.at(lookup_context_index(int2iSG_seq(i, order*2+1)))->value << ") i(" << lookup_context.at(lookup_context_index(int2iSG_seq(i, order*2+1)))->inverse << ")" << endl;
//		cerr << "   " << lookup_table.at(order).at(i)->value << "  " << lookup_table.at(order).at(i)->inverse << endl;
	}

	//////////
	/// Calculate the changes for the very first sequence position. This will be used for future sequences.
	/// * This position will use the transition frequencies for the last site. 
	//////////
	string tuplet;
	string sequence_5prime;
	for (i = 0; i < index_position_multiplier.at(order+1); i++) {
		sequence = lookup_table_sequence(i, order+1);
		tuplet = sequence.substr(0, order);
		RN = tuplet_pi.at(lookup_ktuplet_index(tuplet));
		sequence_5prime = "*"+sequence;
		//cerr << " sequence: " << sequence_5prime << endl;
		RN *= lookup_table2[sequence+"*"]->value;
		lookup = new LookUp();
		lookup->value = RN;
		RN_inv = 1.0 / RN;
		lookup->inverse = RN_inv;
		lookup_table2.insert( pair<string, LookUp*>(sequence_5prime, lookup) );

		//lookup_context.at(lookup_context_index(int2iSG_seq(i, order+1))) = lookup;

		lookup_table.at(0).at(i) = lookup;
		//cerr << "lookup_table[" << sequence_5prime << "] = v(" << lookup_table2[sequence_5prime]->value << ") i(" << lookup_table2[sequence_5prime]->inverse << ")  ";
		//cerr << "lookup_table[" << lookup_context_index(int2iSG_seq(i, order+1)) << "] = v(" << lookup_context.at(lookup_context_index(int2iSG_seq(i, order+1)))->value << ") i(" << lookup_context.at(lookup_context_index(int2iSG_seq(i, order+1)))->inverse << ")" << endl;
		//cerr << "     " << lookup_table.at(0).at(i)->value << ", " << lookup_table.at(0).at(i)->inverse << endl;
	}

	//////////
	/// Finish off the lookup table by adding changes to the 2nd, 3rd, ..., order positions.
	//////////
	for (int j = 2; j <= order; j++) {
		for (i = 0; i < index_position_multiplier.at(order+j); i++) {
			sequence = lookup_table_sequence(i, order+j);
			tuplet = sequence.substr(0, order);
			RN = tuplet_pi.at(lookup_ktuplet_index(tuplet));
			sequence_5prime = "*" + sequence;
			//cerr << " sequence: " << sequence_5prime << endl;
			
			/// Pre-calculated value for ktuplet for current sequence size - 1
			substr_key = "*";
			substr_key += sequence.substr(0, sequence.size()-1);
			//cerr << "lookup_table[" << substr_key << "] = "; cerr << lookup_table[substr_key]->value << endl;
			RN = lookup_table2[substr_key]->value;
			/// Pre-calculated transition probabilities for size order+1 sequence (for last position).
			substr_key = sequence.substr(sequence.size()-order-1, order+1);
			substr_key += "*";
			//cerr << " X lookup_table[" << substr_key << "] = " << lookup_table[substr_key]->value << endl;
			RN *= lookup_table2[substr_key]->value;

			lookup = new LookUp();
			lookup->value = RN;
			RN_inv = 1.0 / RN;
			lookup->inverse = RN_inv;
			lookup_table2.insert( pair<string, LookUp*>(sequence_5prime, lookup) );
			//lookup_context.at(lookup_context_index(int2iSG_seq(i, order+j))) = lookup;
			
			lookup_table.at(j-1).at(i) = lookup;
			//cerr << "lookup_table[" << sequence_5prime << "] = v(" << lookup_table2[sequence_5prime]->value << ") i(" << lookup_table2[sequence_5prime]->inverse << ")" << endl;
			//cerr << "lookup_table[" << lookup_context_index(int2iSG_seq(i, order+j)) << "] = v(" << lookup_context.at(lookup_context_index(int2iSG_seq(i, order+j)))->value << ") i(" << lookup_context.at(lookup_context_index(int2iSG_seq(i, order+j)))->inverse << ")" << endl;
			//cerr << "    " << j-1 << ", " << i << "       " << lookup_table.at(j-1).at(i)->value << ", " << lookup_table.at(j-1).at(i)->inverse << endl;
		}
	}

	//////////
	/// Following code is FOR TESTING ONLY: Get the range of values for XXXYXXX in order-3 markov model,
	/// i.e., P(Y|XXX) x P(X|YXX), etc.
	//////////
	//int gg = 0;
	//double minV = 1; double maxV = 0;
	//vector<LookUp*>::iterator lt_ptr = lookup_table.at(order).begin();
	//for (; lt_ptr != lookup_table.at(order).end(); ++lt_ptr, ++gg) {
	//	cout << gg << " " << (*lt_ptr)->value << endl;
	//	if ((*lt_ptr)->value > maxV) maxV = (*lt_ptr)->value;
	//	if ((*lt_ptr)->value < minV) minV = (*lt_ptr)->value;
	//}

	//cout << minV << " " << maxV << endl;

	//exit(0);
}

//////////
/// Functions
//////////
void contextDependence::readDependencies(string& file)
{
	ifstream is;
	char line[2096];
	int i = 0;
	string sequence;
	double val;
	LookUp *lookup;
	list<string> split_line;
	list<string>::iterator it;

	cerr << "Point-> contextDependence::readDependencies(file) reading " << file << endl;

	tuplet_pi.assign(index_position_multiplier.at(order), 0);
	vector<double>::iterator tup_it = tuplet_pi.begin();
	is.open(file.c_str());
	int num_lines = 0;
	while(is.good()) {
		// Read in line of values. Element 1 is seed ktuplet (K) probability. Following are Pr(N|K).
		is.getline(line, 2096);

		// Split into each element, as separated by spaces.
		split_line = split (line, " ");

		// last line might have no values, will cause problems.
		if (split_line.size() < 2) break;

		// Iterator for the list of elements, start at beginning.
		it = split_line.begin();

		// Set element 1 (ktuplet probability).
		(*tup_it) = atof((*it).c_str());

		// Set to the first transition probability.
		++tup_it, ++it;
		for (; it != split_line.end(); ++it, ++i) {
//			cerr << "i: " << i ;
			if (order_3_markov) sequence = lookup_table_sequence(i, order+1);
			else sequence = lookup_table_sequence(i, (order+1) * 3);	// Codons...
//			cerr << "  sequence: " << sequence;
			sequence += "*";
			val = atof((*it).c_str());
			lookup = new LookUp(val);
			if (order_3_markov) lookup_table2.insert( pair<string, LookUp*>(sequence, lookup) );
			//lookup_context.at(lookup_context_index(int2iSG_seq(i, order+1), true)) = lookup;
//			cerr << "     Element " << i << "/" << lookup_table.at(order+1).size() << endl;
			lookup_table.at(order+1).at(i) = lookup;
		}
	}
	is.close();

//	i = 0;
//	for (vector<LookUp*>::iterator dt = lookup_table.at(order+1).begin(); dt != lookup_table.at(order+1).end(); ++dt, ++i) {
//		cerr << "Lookup table element " << i << ": " << (*dt)->value << endl;
//	}
	//report_tuplet_pi();
}

void contextDependence::allocate_lookup_context_vector()
{
	vector<int>::iterator it;
	int lookup_context_size = 0;
	LookUp *lookup;
	LookUp *dummy;	// For allocating vector<LookUp*> below.
	lookup_table.assign( order*2+1, vector<LookUp*> (1, dummy) );

	cerr << "ORDER? " << order << endl;

	if (order_3_markov) {
		//////////
		/// General purpose vector that pre-multiplies the numStates. Used for setting offsets, for one. ///
		index_position_multiplier.assign(order*2+2, 1);
		for (it = index_position_multiplier.begin()+1; it != index_position_multiplier.end(); ++it)
			(*it) = (*(it-1))*numStates;
		/// Allocate lookup array
		for (it = index_position_multiplier.begin()+order+1; it != index_position_multiplier.end(); ++it) {
			lookup_context_size += 2*(*it);
		}
		lookup_context_size -= index_position_multiplier.at(order+4);
		//lookup_context.assign(lookup_context_size, lookup);
		//cerr << "lookup_context_size = " << lookup_context_size << endl;
		/// Setup vector that indexes the star of each order setting.
		///		For example, order 3 markov, this vector will hold:
		///    *4   *5    *6      7      4*      5*      6*     XX
		///     0 | 256 | 1280 | 5376 | 21760 | 22016 | 23040 | 27136
		/// Bookend values not really necessary, but are useful as placeholders (define beginning and ending of array).
		context_specific_index_offset.assign(order*2+1,0);
		vector<int>::iterator jt = context_specific_index_offset.begin()+1;
		for (it = index_position_multiplier.begin()+order+1; it != index_position_multiplier.end(); ++it, ++jt)
			(*jt) = (*(jt-1))+(*it);
		for (it = index_position_multiplier.begin()+order+1; it != index_position_multiplier.end(); ++it, ++jt)
			(*jt) = (*(jt-1))+(*it);

		// 2D vector lookup table.
		vector<vector<LookUp*> >::iterator xt = lookup_table.begin();
		for (it = index_position_multiplier.begin()+order+1; it != index_position_multiplier.end(); ++it, ++xt)
			(*xt).assign((*it), dummy);
		for (it = index_position_multiplier.begin()+order+1; it != index_position_multiplier.end()-1; ++it, ++xt)
			(*xt).assign((*it), dummy);
		int i = 0;

		int o=0;
		for (jt = context_specific_index_offset.begin(); jt != context_specific_index_offset.end(); ++jt, ++o) {
			cerr << "Order " << o << " offset: " << (*jt) << endl;
		} //exit(0);
	} else if (Human_Data_simulation) {
		// This is hard coded to make human example work. To generalize, need to include that expected order of the
		// data, and proceed similarly as above. This is for codons (triplets), where only 1 site will vary.
		
		cerr << "Point-> allocate_lookup_context_vector()" << endl;

		//////////
		/// General purpose vector that pre-multiplies the numStates. Used for setting offsets, for one. ///
		index_position_multiplier.assign(order*6+2, 1);
		for (it = index_position_multiplier.begin()+1; it != index_position_multiplier.end(); ++it)
			(*it) = (*(it-1))*numStates;
		/// Allocate lookup array
//		for (it = index_position_multiplier.begin()+order+1; it != index_position_multiplier.end(); ++it) {
//			lookup_context_size += 2*(*it);
//		}
//		lookup_context_size -= index_position_multiplier.at(order+4);
//		lookup_context.assign(lookup_context_size, lookup);

		
		// Beginning of sequence change probability:
		// i_h^1 i_h^2 i_h^3   i_{4}^1 i_{5}^2 i_{6}^3
		// We will assume that this is always ATG
		lookup_table.at(0).assign(pow(numStates, 6), dummy);
		lookup_table.at(1).assign(pow(numStates, 9), dummy);
		lookup_table.at(2).assign(pow(numStates, 6), dummy);
		
		cerr << "Sizes of each lookup table element:" << endl;
		cerr << "  0: " << lookup_table.at(0).size() << endl;
		cerr << "  1: " << lookup_table.at(1).size() << endl;
		cerr << "  2: " << lookup_table.at(2).size() << endl;
	} else {
		cerr << "Queer... should not reasonably get here in contextDependence::allocate_lookup_context_vector()." << endl;
		exit(EXIT_FAILURE);
	}

}

contextDependence::~contextDependence()
{
	vector<vector<LookUp*> >::iterator xt = lookup_table.begin();
	vector<LookUp*>::iterator yt;
	for (xt = lookup_table.begin(); xt != lookup_table.end(); ++xt) {
		for (yt = (*xt).begin(); yt != (*xt).end(); ++yt) {
			delete (*yt);
		}
	}
}

void 
contextDependence::reset_sequence_indices(
										  TNode *node, 
										  int event_site,
										  string event
										 )
{
	int start = event_site;
	unsigned int end;
	node->set_site_window(order, &start, &end);
	int event_offset = stateCharacters.find_first_of(event.at(1)) - stateCharacters.find_first_of(event.at(0));
	vector<int>::iterator bt;

//	cerr << "Changed site " << event_site << "->" << event << endl;
//	cerr << node->printSequence() << endl;
//	cerr << "Modifying sites " << start << "->" << end << endl;

//	cerr << "Previous indices: " << endl;
//	for (vector<Site>::iterator it = node->seq_evo.begin(); it != node->seq_evo.end(); ++it)
//		cerr << (*it).return_lookup_table_sequence_index() << " ";
//	cerr << endl;
	
	// MIDDLE: If in the middle of the sequence, the index position will change by the nature of the event 
	// (e.g., T->A is a 3 to a 0 in the stateCharacters array), multiplied by the index_position_multiplier,
	// i.e., if we are changing the index at the beginning of the string in order-3, the index will be
	// old_index+4^6*(-3)
	// END:    If at the end of the sequence, still start at ipm beginning, but the nature of the change
	// to the sequence indices is different in that after entering the last order+1 sites, the value removed
	// from the previous indices will remain constant. The Nth site will always have the multiplier 4^0, 
	// the (N-1)st site will always have the multiplier 4^1, and so on.
	// BEGINNING: If in the beginning of the sequence, the ipm index will have to start a little ways into the
	// ipm array, since a change to first site will have a specific multiplier > 1 (actually, it will be 
	// 4^order), and so on...
	bt = index_position_multiplier.begin();
	if (node->seq_evo.at(event_site).return_lookup_table_environment_index() < order) {
		bt += order-node->seq_evo.at(event_site).return_lookup_table_environment_index();
	}

	int mult;
	int at_site = start;
	for (vector<Site>::iterator it = node->seq_evo.begin()+start; at_site < end && it != node->seq_evo.end(); ++it, ++bt, ++at_site) {
		if (at_site >= node->seq_evo.size()-(order+1)) {
			for (; at_site < end && it != node->seq_evo.end(); ++it, ++at_site)
				(*it).set_lookup_table_sequence_index((*it).return_lookup_table_sequence_index()+(*bt)*event_offset);
			--it;
			continue;
		} else (*it).set_lookup_table_sequence_index((*it).return_lookup_table_sequence_index()+(*bt)*event_offset);
	}

//	cerr << "New indices: " << endl;
//	for (vector<Site>::iterator it = node->seq_evo.begin(); it != node->seq_evo.end(); ++it)
//		cerr << (*it).return_lookup_table_sequence_index() << " ";
//	cerr << endl;
}

void contextDependence::set_sequence_indices(TNode *node)
{
	vector<Site>::iterator it = node->seq_evo.begin();
	vector<short> seq;

	//cerr << node->printSequence();

	// Beginning <order> positions.
	for (int i = 0; i < order; ++i) seq.push_back((*(it+i)).returnState());
	seq.push_back((*(it+order)).returnState());
	(*it).set_lookup_table_environment_index(0);
	(*it).set_lookup_table_sequence_index(sequence_specific_index_offset(seq));
	//cerr << "  " << (*it).return_lookup_table_environment_index() << ", " << (*it).return_lookup_table_sequence_index() << endl;
	++it;
	int i = 1;
	for (i = 1; i < order; ++i, ++it) {
		(*it).set_lookup_table_environment_index(i);
		(*it).set_lookup_table_sequence_index( (*(it-1)).return_lookup_table_sequence_index()*numStates + (*(it+order)).returnState() );
		//cerr << "  " << (*it).return_lookup_table_environment_index() << ", " << (*it).return_lookup_table_sequence_index() << endl;
	}

	// Middle.
	for (; it != node->seq_evo.end()-order; ++it, ++i) {
		(*it).set_lookup_table_environment_index(order);
		(*it).set_lookup_table_sequence_index( ((*(it-1)).return_lookup_table_sequence_index()*numStates + (*(it+order)).returnState()) % index_position_multiplier.at(order*2+1));
		//cerr << "  " << (*it).return_lookup_table_environment_index() << ", " << (*it).return_lookup_table_sequence_index() << endl;
	}

	// End <order> positions.
	int j = order*2;
	int k = 1;
	for (; it != node->seq_evo.end(); ++it, ++i, ++k, --j) {
		(*it).set_lookup_table_environment_index(2*order+1-k);
		//cerr << (*(it-1)).return_lookup_table_sequence_index() << " - " << (*(it-order)).returnState()*index_position_multiplier.at(j) << endl;
		(*it).set_lookup_table_sequence_index((*(it-1)).return_lookup_table_sequence_index() - (*(it-order-1)).returnState()*index_position_multiplier.at(j));
		//cerr << "  " << (*it).return_lookup_table_environment_index() << ", " << (*it).return_lookup_table_sequence_index() << endl;
	}
}

void
contextDependence::setOffset()
{
	// Need to keep track of the offsets of neighboring patterns. 
	// ASSUMING WE ARE MID-SEQUENCE, given an order 3 markov dependence,
	// if the current nucleotide is an A, XXXAXXX with index indexA, then the neighboring sequences 
	// XXXCXXX, XXXGXXX, and XXXTXXX will be indexA+4^3, indexA+2*4^3, and indexA+3*4^3
	// For XXXCXXX, indices will be indexC-4^3, indexC+4^3, and indexC+2*4^3  for A,G,T
	// For XXXGXXX, indices will be indexG-2*4^3, indexG-4^3, and indexG+4^3  for A,C,T
	// For XXXTXXX, indices will be indexT-3*4^3, indexT-2*4^3, and indexT-2*4^3 for A,C,G
	// If we are not mid-sequence, then 4^3 needs to be changed to the appropriate value,
	// seq position 0: 4^3
	// seq position 1: 4^3
	// seq position 2: 4^3
	// seq position 3: 4^3 (order)
	// seq position N-2: 4^2
	// seq position N-1: 4^1
	// seq position N: 4^0
	int multiple_of_numStates = index_position_multiplier.at(order);
	// Makes array for the multiplier above. This array essentially removes the zero multiplier,
	// in this case, it will be -3 -2 -1  1  2  3
	vector<int> num_steps_from_current (order*2, 0);
	int start = -order;
	for (vector<int>::iterator it = num_steps_from_current.begin(); it != num_steps_from_current.end(); ++it, start++) {
		if (start == 0) start++;
		(*it) = start;
	}
	index_offset.assign(order*2+1, vector<int> (numStates*numStates, 0));
	vector<vector<int> >::iterator at = index_offset.begin();
	for (; at != index_offset.begin()+order+1; ++at) {
		for (int j = 0; j < numStates; ++j) {		// current == sequence i
			for (int k = 0; k < numStates; ++k) {		// proposed	== sequence j
				// e.g., current T, proposed A, should be -3...
				// current T -> iterator j = 3
				// proposed A -> iterator k = 0
				// Thus, to index num_steps_from_current, use k-j.
				// for j == k, who cares. TauIJ will not check reference since it calculates only sequences differing in one position.
				(*at).at(j*numStates+k) = multiple_of_numStates * (k-j);
			}
		}
	}
	int ipm = 0;
	for (at = index_offset.begin()+order+1; at != index_offset.end(); ++at) {
		multiple_of_numStates = index_position_multiplier.at(ipm++);
		for (int j = 0; j < numStates; ++j) {		// current == sequence i
			for (int k = 0; k < numStates; ++k) {		// proposed	== sequence j
				(*at).at(j*numStates+k) = multiple_of_numStates * (k-j);
			}
		}
	}

//	cerr << "Offset setup: " << endl;
//	int p = 0;
//	for (at = index_offset.begin(); at != index_offset.end(); ++at, ++p) {
//		cerr << "Level " << p << ": " << endl;
//		for (int j = 0; j < numStates; ++j) {		// current == sequence i
//			for (int k = 0; k < numStates; ++k) {		// proposed	== sequence j
//				cerr << (*at).at(j*numStates+k) << " ";
//			}
//			cerr << endl;
//		}
//	}
}

int contextDependence::getOffset(
								 int environment,
								 short i, 
								 short j
								)
{
	return index_offset.at(environment).at(i*numStates+j);
}

void contextDependence::generateDependencies(double dependence_strength_superscript)
{
	//////////
	/// Seed the random model with values. (Generator sequence positions)
	/// * First objects (order+1 length strings).
	/// * tuplet frequencies
	//////////
	string sequence;
	vector<double> tmp_val (numStates, 0);
	vector<string> sequences (numStates, "");
	double RN, RN_inv; 
	vector<double>::iterator tt;
	vector<string>::iterator st;
	int i, j;
	double total;
	LookUp *lookup;

	int generator_index = order+1;	// By indexing scheme. //
	for (i = 0; i < index_position_multiplier.at(order+1);) {
		tt = tmp_val.begin();
		st = sequences.begin();
		total = 0;
		for (j = 0; j < numStates; ++j, ++i, ++tt, ++st) {
			(*st) = lookup_table_sequence(i, order+1) + "*";
			//cerr << "sequence: " << sequence << "  " << (*st) << endl;
			(*tt) = pow(rndu(), dependence_strength_superscript);
			total += (*tt);
		}

		tt = tmp_val.begin();
		for (tt = tmp_val.begin(); tt != tmp_val.end(); ++tt) (*tt) /= total;

		st = sequences.begin();
		j = 0;
		for (tt = tmp_val.begin(); tt != tmp_val.end(); ++tt, ++st, ++j) {
			lookup = new LookUp();
			lookup->value = (*tt);
			lookup->inverse = 1.0/(*tt);
			lookup_table2.insert( pair<string, LookUp*>((*st), lookup) );

			// Place into lookup_context vector. //
			//lookup_context.at( lookup_context_index(int2iSG_seq(i+(j-numStates), generator_index), true) )= lookup;

			lookup_table.at(order+1).at(i+j-numStates) = lookup;
		}
	}

//	cerr << "Lookup table contains " << lookup_table2.size() << " elements:" << endl;
//	for (i = 0; i < index_position_multiplier.at(order+1); i++) {
//		sequence = lookup_table_sequence(i, order+1);
//		sequence += "*";
//		cerr << "lookup_table[" << sequence << "] = v(" << lookup_table2[sequence]->value << ") i(" << lookup_table2[sequence]->inverse << ")";
//		cerr << "vs lookup_context[" << context_specific_index_offset.at(generator_index) << "+" << sequence_specific_index_offset(int2iSG_seq(i, generator_index)) << " = ";		cerr << lookup_context_index(int2iSG_seq(i, generator_index), true) << "] = v(" 
//			 << lookup_context.at(lookup_context_index(int2iSG_seq(i, generator_index), true))->value << ") i("
//	 		 << lookup_context.at(lookup_context_index(int2iSG_seq(i, generator_index), true))->inverse << ")" << endl;
//
//		cerr << "   " << lookup_table.at(order*2).at(i)->value << ", " << lookup_table.at(order*2).at(i)->inverse << endl;
//	}

	//////////
	/// Fill in the values for the first <order> positions
	//////////
	// Set order-tuplets
	total = 0;
	double total2 = 0;
	tuplet_pi.assign(index_position_multiplier.at(order), 0);
	for (vector<double>::iterator it = tuplet_pi.begin(); it != tuplet_pi.end(); ++it) {
		RN = pow(rndu(), dependence_strength_superscript);
		total += RN;
		(*it) = RN;
	}
	i = 0;
	for (vector<double>::iterator it = tuplet_pi.begin(); it != tuplet_pi.end(); ++it, ++i) {
		(*it) /= total;
		total2+=(*it);
		sequence = lookup_table_sequence(i, order);
		//cerr << lookup_ktuplet_index(sequence) << " sequence " << sequence << "  tuplet_freq = " << (*it) << endl;
	}

	//cerr << "CHECK: total of tuplet freqs = " << total2 << endl;
}

int contextDependence::sequence_specific_index_offset(vector<short> sequence)
{
	int index = 0;
	vector<int>::iterator mt = index_position_multiplier.begin()+sequence.size() - 1;
	for (vector<short>::iterator it = sequence.begin(); it != sequence.end(); ++it, --mt) index += (*it)*(*mt);
	return index;
}

vector<short> 
contextDependence::int2iSG_seq(int idx, int sequence_length)
{
	vector<short> sequence (sequence_length, 0);
	vector<int>::iterator mt = index_position_multiplier.begin()+(sequence_length-1);
	short residue;
	for (int i = 0; i < sequence_length; i++) {
		residue = (short)idx / (*mt);
		sequence.at(i) = residue;
		idx -= residue * (*mt);
		--mt;
	}

	return sequence;
}

string contextDependence::lookup_table_sequence(int idx, int sequence_length)
{
	string sequence = "";
	vector<int>::iterator mt = index_position_multiplier.begin()+(sequence_length-1);
	short residue;
	for (int i = 0; i < sequence_length; i++) {
		residue = (short)idx / (*mt);
		sequence.push_back(stateCharacters.at(residue));
		idx -= residue * (*mt);
		--mt;
	}
	
	return sequence;
}

int contextDependence::lookup_ktuplet_index(string sequence)
{
	string::iterator it = sequence.begin();
	vector<int>::iterator mt = index_position_multiplier.begin()+order-1;
	int index = 0;
	for (; it != sequence.end(); ++it, --mt) {
		index += stateCharacters.find(*it) * (*mt);
	}
	return index;
}

void contextDependence::report_tuplet_pi()
{
	string sequence;
	double hi, lo; hi = lo = tuplet_pi.at(0);
	sequence = lookup_table_sequence(0, order);
	//cerr << sequence << ": " << tuplet_pi.at(0) << endl;
	for (int i = 1; i < index_position_multiplier.at(order); ++i) {
		sequence = lookup_table_sequence(i, order);
	//	cerr << sequence << ": " << tuplet_pi.at(i) << endl;
		if (tuplet_pi.at(i) > hi) hi = tuplet_pi.at(i);
		if (tuplet_pi.at(i) < lo) lo = tuplet_pi.at(i);
	}
	//cerr << "Tuplet pi range: (" << lo << ", " << hi << ")" << endl;
}

void contextDependence::outputDependencies(string& outfile_name_root)
{
	vector<double>::iterator it;
	int i = 0, j;
	string sequence;

	ofstream *dependency_out;
	string outfile_name = outfile_name_root;
	outfile_name += ".dep";
	dependency_out = new ofstream(outfile_name.c_str(), ios::trunc | ios::out);
	it = tuplet_pi.begin();
	for (it = tuplet_pi.begin(); it != tuplet_pi.end(); ++it) {
		*dependency_out << (*it);
		for (j = 0; j < numStates; ++j, ++i) {
			sequence = lookup_table_sequence(i, order+1);
			sequence += "*";
			//*dependency_out << " " << lookup_table2[sequence]->value;
			//*dependency_out << " " << lookup_context.at(lookup_context_index(int2iSG_seq(i, order+1), true))->value;
			*dependency_out << " " << lookup_table.at(order+1).at(i)->value;
		}
		*dependency_out << endl;
	}
	(*dependency_out).close();
}

int contextDependence::lookup_table_index(vector<short> sequence)
{
	int index = 0;
	vector<int>::reverse_iterator mt = index_position_multiplier.rbegin()+1;
	for (vector<short>::iterator it = sequence.begin(); it != sequence.end(); ++it, ++mt) index += (*it)*(*mt);
	return index;
}

