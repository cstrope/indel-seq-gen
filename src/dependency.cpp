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
	//cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) IN" << endl;
	context.allocate_lookup_context_vector();
	//cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) LOOKUP CONTEXT VECTOR" << endl;
	context.generateDependencies(dependence_superscript);
	//cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) GENERATE DEPENDENCIES" << endl;
	context.set_lookup_table();
	//cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) SET LOOKUP TABLE " << endl;
	context.setOffset();
	//cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) SET OFFSET " << endl;
	context.outputDependencies(outfile_name_root);
	//cerr << "Point-> Dependency::Dependency(dep_order, dependence_superscript, &file) OUT" << endl;
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
	//cerr << "Point-> Dependency::Dependency(dep_order, block_size, file) 1" << endl;
	context.allocate_lookup_context_vector();
	//cerr << "Point-> Dependency::Dependency(dep_order, block_size, file) 2" << endl;
	context.readDependencies(file);
	//cerr << "Point-> Dependency::Dependency(dep_order, block_size, file) 3" << endl;
	context.set_lookup_table();
	context.setOffset();
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
	if (order_3_markov) {
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
	} else {
		// Human first order markov model. //
		// Right now, lookup_table.at(2) is fully filled with 4096 sequences, as is the tuplet_pi.
		// Use these to do lookup_table.at(0) (tuplet_pi * 2codon prob.) and lookup_table.at(1) (middle.)
		// lookup_table.at(0);
		vector<double>::iterator tup_it = tuplet_pi.begin();
		for (int i = 0; i < 64; i++, ++tup_it) {
			for (int j = 0; j < 64; j++) {
				lookup = new LookUp((*tup_it) * lookup_table.at(2).at(i*64+j)->value); 
				// Assume P(ATG) = 1, therefore, there is 0 chance that it will change.
				// With these values, tau_ij = P(j)/P(i) * P_0(i)/P_0(j) = 1, so Rij = u\pi log(tau_ij)/1-1/tau_ij = 0; 
				lookup->value = numeric_limits<double>::min();
				lookup->inverse = numeric_limits<double>::max();
				lookup_table.at(0).at(i*64+j) = lookup;
			}
		}

		int i = 0;
		int lt_index_1, lt_index_2;
		for (vector<LookUp*>::iterator it = lookup_table.at(1).begin(); it != lookup_table.at(1).end(); ++it, ++i) {
			get_first_order_indices(i, &lt_index_1, &lt_index_2);
			lookup = new LookUp(lookup_table.at(2).at(lt_index_1)->value * lookup_table.at(2).at(lt_index_2)->value);
			(*it) = lookup;
		}
	}
}

void 
contextDependence::get_first_order_indices (int tuplet, int *index1, int *index2)
{
	vector<int> decompose_9tuplet (9, 0);

	// For ACAGGTTTG, this consists of ACAGGT and GGTTTG. This vector will be used to get the
	// nucleotide at each position of the string, then we can get the indices of the component
	// parts.
	vector<int>::iterator it = decompose_9tuplet.begin();
	vector<int>::reverse_iterator rit = index_position_multiplier.rbegin();
	for (; rit != index_position_multiplier.rend(); ++rit, ++it) {
		(*it) = tuplet / (*rit);		// Num integers that can go into i
		tuplet -= (*it)*(*rit);		// Subtract this amount, proceed.
	}

	(*index1) = 0;
	rit = index_position_multiplier.rbegin()+3;
	for (it = decompose_9tuplet.begin(); it != decompose_9tuplet.end()-3; ++it, ++rit) (*index1) += (*it)*(*rit);
	(*index2) = 0;
	rit = index_position_multiplier.rbegin()+3;
	for (it = decompose_9tuplet.begin()+3; it != decompose_9tuplet.end(); ++it, ++rit) (*index2) += (*it)*(*rit);
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

	if (order_3_markov) tuplet_pi.assign(index_position_multiplier.at(order), 0);
	else tuplet_pi.assign(numStates*numStates*numStates, 0);		// Codons, triplets.
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
			val = atof((*it).c_str());
			// If the read-in value is zero, it causes havoc with nan's in the tau_ij routine. Set to minimum double
			// value, instead. Need to make sure that this number can be squared!! This is because if we are going to 
			// multiply two "0"s, and these were both at numeric_limits<double>::min(), we would get underflow, and
			// the problem with zeroes would persist. So, set val to the square root of the minimum double value.
			if (val == 0) val = sqrt(numeric_limits<double>::min());
			lookup = new LookUp(val);
			if (order_3_markov) {
				sequence = lookup_table_sequence(i, order+1);
//				cerr << "  sequence: " << sequence;
				sequence += "*";
				lookup_table2.insert( pair<string, LookUp*>(sequence, lookup) );
			}
			//lookup_context.at(lookup_context_index(int2iSG_seq(i, order+1), true)) = lookup;
			//cerr << "     Element " << i << "/" << lookup_table.at(order+1).size();
			//cerr << "  VALUE: " << lookup->value << endl;
			lookup_table.at(order+1).at(i) = lookup;
			//if (lookup_table.at(order+1).at(0)->value < 0.25446 && lookup_table.at(order+1).at(0)->value > 0.254467) {
			//	cerr << lookup_table.at(order+1).at(0) << endl;
			//	exit(0);
			//}
		}
	}

	cerr << "Exited loop" << endl;
	is.close();

//	i = 0;
//	for (vector<LookUp*>::iterator dt = lookup_table.at(order+1).begin(); dt != lookup_table.at(order+1).end(); ++dt, ++i) {
//		cerr << "Lookup table element " << i << ": " << (*dt)->value << endl;
//	}
	//report_tuplet_pi();

	cerr << "Exiting function  " << lookup_table.at(order+1).size() << endl;
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
		index_position_multiplier.assign(order*6+3, 1);
		for (it = index_position_multiplier.begin()+1; it != index_position_multiplier.end(); ++it)
			(*it) = (*(it-1))*numStates;
		
		// Beginning of sequence change probability:
		// i_h^1 i_h^2 i_h^3   i_{4}^1 i_{5}^2 i_{6}^3
		// We will assume that this is always ATG
		lookup_table.at(0).assign(pow((double)numStates, 6.), dummy);
		lookup_table.at(1).assign(pow((double)numStates, 9.), dummy);
		lookup_table.at(2).assign(pow((double)numStates, 6.), dummy);
		
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
										  // Block size (for codons)
										  int event_site,
										  string event
										 )
{
	int start = event_site;
	unsigned int end;
	node->set_site_window(order, &start, &end);
	int event_offset = stateCharacters.find_first_of(event.at(1)) - stateCharacters.find_first_of(event.at(0));
	vector<int>::iterator bt;

	//cerr << "Changed site " << event_site << "->" << event << endl;
	//cerr << node->printSequence() << endl;
	//cerr << "Modifying sites " << start << "->" << end << endl;

	//cerr << "Previous indices: " << endl;
	//for (vector<Site>::iterator it = node->seq_evo.begin(); it != node->seq_evo.end(); ++it)
	//	cerr << (*it).return_lookup_table_sequence_index() << " ";
	//cerr << endl;

	if (Human_Data_simulation) {
		int block_size = 3;
		int codon1[3], codon2[3], codon3[3];	// Get the codons. Can translate them and multiply to get correct indices.
		int from_codon[3];
		vector<int> new_indices;
		vector<int>::iterator nt = new_indices.begin();
		// Again, assuming that a change cannot occur in the first 3 sites of the sequence, but still need to calculate
		// sequence indices for EPC (?).
		if (end - start == order*2*block_size) {
			if (start < order) {
				
			} else {
				// At the end of the sequence.
				for (int i = 0; i < block_size; ++i) codon1[i] = node->seq_evo.at(start+i).returnState();
				for (int i = 0; i < block_size; ++i) codon2[i] = node->seq_evo.at(start+block_size+i).returnState();
				from_codon[0] = codon2[0];
				from_codon[1] = codon2[1];
				from_codon[2] = codon2[2];
				int event_site_codon_position = event_site % block_size;
				from_codon[event_site_codon_position] = stateCharacters.find(event.at(0));
				//cerr << "codon 1: ";
				//for (int i = 0; i < block_size; ++i) cerr << stateCharacters.at(codon1[i]);
				//cerr << endl << "codon 2: ";
				//for (int i = 0; i < block_size; ++i) cerr << stateCharacters.at(from_codon[i]);
				//cerr << " --> ";
				//for (int i = 0; i < block_size; ++i) cerr << stateCharacters.at(codon2[i]);
				//cerr << endl;
				for (vector<Site>::iterator it = node->seq_evo.begin()+start; it != node->seq_evo.begin()+end; ++it) {
					new_indices.push_back((*it).return_lookup_table_sequence_index());
				}
				// Populate vector with older indices.
				vector<int> new_indices;
				for (vector<Site>::iterator it = node->seq_evo.begin()+start; it != node->seq_evo.begin()+end; ++it) {
					new_indices.push_back((*it).return_lookup_table_sequence_index());
				}
				int codon1_index = codon1[0]*16+codon1[1]*4+codon1[2];
				int codon2_index = codon2[0]*16+codon2[1]*4+codon2[2];
				int codon3_index = codon3[0]*16+codon3[1]*4+codon3[2];
				int from_codon_index = from_codon[0]*16+from_codon[1]*4+from_codon[2];
				nt = new_indices.begin();
				int index_change_base = codon2_index - from_codon_index;
				(*nt) += index_change_base; nt++;
				(*nt) += index_change_base; nt++;
				(*nt) += index_change_base; nt++;
				// No multiply by 64 here, since ne'er does more codon cometh into th' picture.
				// E.g. only 6 multipliers are used between the codons:
				// 1024 256 64 16 4 1
				// If one of the last three change, the same multiplier is used for both (codon1 and codon2) index modifications.
				(*nt) += index_change_base; nt++;
				(*nt) += index_change_base; nt++;
				(*nt) += index_change_base; nt++;
				//for (nt = new_indices.begin(); nt != new_indices.end(); ++nt) cerr << (*nt) << endl;

				nt = new_indices.begin();
				for (vector<Site>::iterator it = node->seq_evo.begin()+start; it != node->seq_evo.begin()+end; ++it, ++nt) 
					(*it).set_lookup_table_sequence_index(*nt);

				//for (vector<Site>::iterator it = node->seq_evo.begin(); it != node->seq_evo.end(); ++it) 
				//	cerr << (*it).return_lookup_table_sequence_index() << " ";
				//cerr << endl;
			}
		} else {
			// Middle of the sequence.
			for (int i = 0; i < block_size; ++i) codon1[i] = node->seq_evo.at(start+i).returnState();
			for (int i = 0; i < block_size; ++i) codon2[i] = node->seq_evo.at(start+block_size+i).returnState();
			for (int i = 0; i < block_size; ++i) codon3[i] = node->seq_evo.at(start+2*block_size+i).returnState();

			from_codon[0] = codon2[0];
			from_codon[1] = codon2[1];
			from_codon[2] = codon2[2];
			int event_site_codon_position = event_site % block_size;
			//cerr << "event_site_codon_position: " << event_site_codon_position << endl;
			from_codon[event_site_codon_position] = stateCharacters.find(event.at(0));

			//cerr << "codon 1: ";
			//for (int i = 0; i < block_size; ++i) cerr << stateCharacters.at(codon1[i]);
			//cerr << endl << "codon 2: ";
			//for (int i = 0; i < block_size; ++i) cerr << stateCharacters.at(from_codon[i]);
			//cerr << " --> ";
			//for (int i = 0; i < block_size; ++i) cerr << stateCharacters.at(codon2[i]);
			//cerr << endl << "codon 3: ";
			//for (int i = 0; i < block_size; ++i) cerr << stateCharacters.at(codon3[i]); 
			//cerr << endl;

			// Populate vector with older indices.
			for (vector<Site>::iterator it = node->seq_evo.begin()+start; it != node->seq_evo.begin()+end; ++it) {
				new_indices.push_back((*it).return_lookup_table_sequence_index());
			}

			int codon1_index = codon1[0]*16+codon1[1]*4+codon1[2];
			int codon2_index = codon2[0]*16+codon2[1]*4+codon2[2];
			int codon3_index = codon3[0]*16+codon3[1]*4+codon3[2];
			int from_codon_index = from_codon[0]*16+from_codon[1]*4+from_codon[2];
			nt = new_indices.begin();
			int index_change_base = codon2_index - from_codon_index;
			(*nt) += index_change_base; nt++;
			(*nt) += index_change_base; nt++;
			(*nt) += index_change_base; nt++;
			index_change_base *= 64;
			(*nt) += index_change_base; nt++;
			(*nt) += index_change_base; nt++;
			(*nt) += index_change_base; nt++;
			if (end != node->seq_evo.size()) index_change_base *= 64;
			(*nt) += index_change_base; nt++;
			(*nt) += index_change_base; nt++;
			(*nt) += index_change_base; nt++;

			//for (nt = new_indices.begin(); nt != new_indices.end(); ++nt) cerr << (*nt) << endl;

			nt = new_indices.begin();
			for (vector<Site>::iterator it = node->seq_evo.begin()+start; it != node->seq_evo.begin()+end; ++it, ++nt) 
				(*it).set_lookup_table_sequence_index(*nt);

			//for (vector<Site>::iterator it = node->seq_evo.begin(); it != node->seq_evo.end(); ++it) 
			//	cerr << (*it).return_lookup_table_sequence_index() << " ";
			//cerr << endl;
		}
	} else {
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
	}

	//cerr << "New indices: " << endl;
	//for (vector<Site>::iterator it = node->seq_evo.begin(); it != node->seq_evo.end(); ++it)
	//	cerr << (*it).return_lookup_table_sequence_index() << " ";
	//cerr << endl;
}

void contextDependence::set_sequence_indices(
											 TNode *node,
											 int block_size
											)
{
	vector<Site>::iterator it = node->seq_evo.begin();
	vector<short> seq;
	int i, site;

	if (order_3_markov) {
		// Beginning <order> positions.
		for (i = 0; i < order*block_size; ++i) seq.push_back((*(it+i)).returnState());
		seq.push_back((*(it+order)).returnState());
		(*it).set_lookup_table_environment_index(0);
		(*it).set_lookup_table_sequence_index(sequence_specific_index_offset(seq));
		++it;
		i = 1;
		for (i = block_size; i < order*block_size; ++i, ++it) {
			(*it).set_lookup_table_environment_index(i);
			// This is a shortcut that works well for order 1, but not so well for higher orders...
			// Can make a function to handle this.
			(*it).set_lookup_table_sequence_index( 
												  (*(it-1)).return_lookup_table_sequence_index()*numStates 
												  + (*(it+order)).returnState() 
												 );
		}

		// Middle.
		for (; it != node->seq_evo.end()-order*block_size; ++it, ++i) {
			(*it).set_lookup_table_environment_index(order);
			(*it).set_lookup_table_sequence_index( ((*(it-1)).return_lookup_table_sequence_index()*numStates + (*(it+order)).returnState()) % index_position_multiplier.at(order*2+1));
		}

		// End <order> positions.
		int j = order*2;
		int k = 1;
		for (; it != node->seq_evo.end(); ++it, ++i, ++k, --j) {
			(*it).set_lookup_table_environment_index(2*order+1-k);
			(*it).set_lookup_table_sequence_index((*(it-1)).return_lookup_table_sequence_index() - (*(it-order-1)).returnState()*index_position_multiplier.at(j));
		}
	} else {
		// Beginning order*block_size positions (beginning of the sequence only?)
		// The first position will be for when changes occur in the beginning of the sequence. These
		// types of changes will affect the first position as well as the downstream dependencies.
		// We represent the first position as order+block_size, and the downstream dependencies as
		// + block_size. The first block_size positions will be assigned an index into the lookup_table.
		for (i = 0; i < order*block_size+block_size; ++i) 
			seq.push_back((*(it+i)).returnState());

		//for (vector<short>::iterator q = seq.begin(); q != seq.end(); ++q)
		//	cerr << stateCharacters.at(*q);
		//cerr << endl;

		site = 0;
		for (it = node->seq_evo.begin(); it != node->seq_evo.begin()+block_size; ++it, ++site) {
			(*it).set_lookup_table_environment_index(0);
			(*it).set_lookup_table_sequence_index(sequence_specific_index_offset(seq));
			//cerr << "Site " << site << ":" << "  env->" << (*it).return_lookup_table_environment_index() << "  idx->" << (*it).return_lookup_table_sequence_index() << endl;
		}
		

		// Middle of sequence
		for (it = node->seq_evo.begin()+order*block_size; it != node->seq_evo.end()-block_size; it += block_size) {
			// First begin by moving the next <block_size> to the front of the sequence:
			// First pass through, the seq.size() will be 6 (for 1OMM, BS=3), so skip this as the first codon will still
			// be needed to find the correct index. Beyond that, need to move data to the front.
			i = 0;
			if (seq.size() >= 3*order*block_size) {
				// Move to the front of the array.
				// Why the 2 in the 2nd term? Because we are moving the first 2 codons to the front of the array,
				// then adding the third codon.
				for (vector<short>::iterator q = seq.begin(); q != seq.begin()+2*order*block_size; ++q, ++i)
					(*q) = (*(q+order*block_size));
				for (int x = 0; x < block_size; x++) seq.pop_back();
			}
			// Replace the back elements with the next 3 positions of the sequence
			for (i = block_size; i < 2*block_size; i++) 
				seq.push_back((*(it+i)).returnState());

			//for (vector<short>::iterator q = seq.begin(); q != seq.end(); ++q)
			//	cerr << stateCharacters.at(*q);
			//cerr << endl;

			// Multiply the two previous entries by 64, then add new block.
			for (i = 0; i < block_size; i++, site++) {
				(*(it+i)).set_lookup_table_environment_index(order);
				(*(it+i)).set_lookup_table_sequence_index( sequence_specific_index_offset(seq) );
				//cerr << "Site " << site << ":" << "  env->" << (*(it+i)).return_lookup_table_environment_index() << "  idx->" << (*(it+i)).return_lookup_table_sequence_index() << endl;
			}
		}

		// End of sequence. Need to move the 2 codons to the front, and remove the last 3 positions.
		for (vector<short>::iterator q = seq.begin(); q != seq.begin()+2*order*block_size; ++q, ++i)
			(*q) = (*(q+order*block_size));
		for (int x = 0; x < block_size; x++) seq.pop_back();
		for (; it != node->seq_evo.end(); it += block_size, ++site) {
			for (i = 0; i < block_size; ++i) {
				(*(it+i)).set_lookup_table_environment_index(2);
				(*(it+i)).set_lookup_table_sequence_index( sequence_specific_index_offset(seq) );
				//cerr << "Site " << site << ":" << "  env->" << (*(it+i)).return_lookup_table_environment_index() << "  idx->" << (*(it+i)).return_lookup_table_sequence_index() << endl;
			}
		}
	}
	
	// Report
	//i = 0;
	//cerr << node->printSequence() << endl;
	///for (it = node->seq_evo.begin(); it != node->seq_evo.end(); ++it, ++i) {
	//	cerr << "Site " << i << ":" << endl;
	//	cerr << "  env->" << (*it).return_lookup_table_environment_index() << endl;
	//	cerr << "  idx->" << (*it).return_lookup_table_sequence_index() << endl;
	//}
}

void
contextDependence::setOffset()
{
	vector<vector<int> >::iterator at;

	if (order_3_markov) {
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
		// Level 0: (Very beginning of sequence? AXXX)
		//  Offset to: A    C    G    T
		// 	Curr: A	   0    64   128  192 
		// 		  C	   -64  0    64   128 
		// 		  G	   -128 -64  0    64 
		// 		  T	   -192 -128 -64  0 
		// Level 1: 
		// 0 64 128 192 
		// -64 0 64 128 
		// -128 -64 0 64 
		// -192 -128 -64 0 
		// Level 2: 
		// 0 64 128 192 
		// -64 0 64 128 
		// -128 -64 0 64 
		// -192 -128 -64 0 
		// Level 3: 
		// 0 64 128 192 
		// -64 0 64 128 
		// -128 -64 0 64 
		// -192 -128 -64 0 
		// Level 4: 
		// 0 1 2 3 
		// -1 0 1 2 
		// -2 -1 0 1 
		// -3 -2 -1 0 
		// Level 5: 
		// 0 4 8 12 
		// -4 0 4 8 
		// -8 -4 0 4 
		// -12 -8 -4 0 
		// Level 6: 
		// 0 16 32 48 
		// -16 0 16 32 
		// -32 -16 0 16 
		// -48 -32 -16 0 
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
		at = index_offset.begin();
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
	} else if (Human_Data_simulation) {
		// Here is a little tougher case, since we are thinking with respect to codons.
		// Middle of sequence, codon position 1, sequence AAAAAAAAA for simplicity:
		// Starting at index_offset.at(1).at(0):
		//            AAAAAAAAA  AAACAAAAA AAAGAAAAA AAATAAAAA
		// AAAAAAAAA   0          1024      2048      3072
		// AAACAAAAA  -1024       0         1024      2048
		// AAAGAAAAA  -2048      -1024      0	      1024
		// AAATAAAAA  -3072      -2048     -1024	  0
		// codon position 2 varies:
		// Starting at index_offset.at(1).at(16):
		//            AAAAAAAAA  AAAACAAAA AAAAGAAAA AAAATAAAA
		// AAAAAAAAA   0          256       512       768
		// AAAACAAAA  -256        0         256       512
		// AAAAGAAAA  -512       -256       0	      256
		// AAAATAAAA  -768       -512      -256	      0
		// codon position 3 varies:
		// Starting at index_offset.at(1).at(32):
		//            AAAAAAAAA  AAAAACAAA AAAAAGAAA AAAAATAAA
		// AAAAAAAAA   0          64        128       192
		// AAAAACAAA  -64         0         64        128
		// AAAAAGAAA  -128       -64        0	      64
		// AAAAATAAA  -192       -128      -64	      0
		//
		// In summary, for environment 1 (which is middle of the sequence), the index offset takes
		// into account which codon position, and what the nature of the change is.
		// If codon position 1, offset is 4**5
		// If codon position 2, offset is 4**4
		// If codon position 3, offset is 4**3
		//
		// Thus, the index offset will have 3 environments, beginning of sequence, middle of sequence, and end of sequence
		// Within each environment, access will be obtained as:
		// (nature_of_change) * (3 - pow(numStates, sequence_site_number % 3)) * pow(numStates, 3)
		// e.g., codon position 2 changes from C->A:
		// -1 * (4^(3-2)) * 4^3 ---> -4^4 ---> -256
		// Checking above, the entry AAAACAAAA -> AAAAAAAAA, -256!

		// Make index_offset array:
		// There are 3 environments, Codons are order 1, so order*2+1 = 3. Perfect!
		// Each environment will have 3 codon positions with a 4*4 flattened matrix. So the environment
		// needs 3*4*4 entries.
		index_offset.assign(order*2+1, vector<int> (numStates*numStates*3));
		int block_size = 3;	// Codons, should be made a variable passed into this function at some time.
		for (int bs = 0; bs < block_size; bs++) {	// Number of items in a codon... Could be replaced by a variable "block_size"
			for (int seqi = 0; seqi < numStates; seqi++) {	// Nucleotide on sequence i, or "from" state
				for (int seqj = 0; seqj < numStates; seqj++) {	// Nucleotide on sequence j, or "to" state
					// I will be explicit in the environments here, but this could be replaced by a loop
					// over all environments with some coding magic.
					index_offset.at(1).at(
										    bs * numStates * numStates	// Each codon position has 16 flattened matrix positions, so this moves to the correct codon position.
										  + seqi*numStates
										  + seqj
										 ) 
										 =  (seqj-seqi)	// Nature of change
										  * pow((double)numStates,(double)(block_size-(bs+1)))	// codon position (0, 1, or 2) needs to be (1, 2, 3)... Conceptually, which is why bs+1
										  * pow((double)numStates,(double)(block_size))		// size of codons before current.;
										 ;
					// Interestingly, the begin environment is the SAME as middle. Makes sense... if you think about it.
					// To be explicit, mid, we change AAAXXXAAA the X's, so, the positions with 4^4, 4^5, and 4^6.
					// Likewise, for beginning, we change XXXAAA, again with positions with 4^4, 4^5, and 4^6.
					index_offset.at(0).at(
										    bs * numStates * numStates	// Each codon position has 16 flattened matrix positions, so this moves to the correct codon position.
										  + seqi*numStates
										  + seqj
										 ) 
										 =  (seqj-seqi)	// Nature of change
										  * pow((double)numStates,(double)(block_size-(bs+1)))	// codon position (0, 1, or 2) needs to be (1, 2, 3)... Conceptually, which is why bs+1
										  * pow((double)numStates,(double)block_size)		// size of codons before current.;
										 ;
					// The end of the sequence, however, we need to do differently.
					// We change AAAXXX, so codon position 0 is 4^2, codon pos 1 is 4^1, and codon pos 2 is 4^0.
					index_offset.at(2).at(
										    bs * numStates * numStates	// Each codon position has 16 flattened matrix positions, so this moves to the correct codon position.
										  + seqi*numStates
										  + seqj
										 ) 
										 =  (seqj-seqi)	// Nature of change
										  * pow((double)numStates,(double)(block_size-(bs+1)))	// codon position (0, 1, or 2) needs to be (1, 2, 3)... Conceptually, which is why bs+1
										 ;
				}
			}
		}
	}

	//cerr << "Offset setup: " << endl;
	//int p = 0;
	//for (at = index_offset.begin(); at != index_offset.end(); ++at, ++p) {
	//	cerr << "Level " << p << ": " << endl;
	//	for (vector<int>::iterator bt = (*at).begin(); bt != (*at).end(); ++bt) {
	//		cerr << (*bt) << " ";
	//	}
	//	cerr << endl;
//		for (int j = 0; j < numStates; ++j) {		// current == sequence i
//			for (int k = 0; k < numStates; ++k) {		// proposed	== sequence j
//				cerr << (*at).at(j*numStates+k) << " ";
//			}
//			cerr << endl;
//		}
	//}
	//exit(0);
}

int contextDependence::getOffset(
								 int environment,
								 int codon_position,
								 short i, 
								 short j
								)
{
	int offset;
	//cerr << "Point-> contextDependence::getOffset(environment=" << environment << " codon_position=" << codon_position << " i=" << i << " j=" << j << ")" << endl;

	return index_offset.at(environment).at(codon_position*numStates*numStates+i*numStates+j);
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

