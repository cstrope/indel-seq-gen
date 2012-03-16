// Dependency matrix inputs file.
// Each dependency should be described with a single matrix. Two items need to // be specified:
// (1) Matrix specification
//   a. Modifications to the global matrix
//     Nucleotide R matrix
//     \  ac ag at
//     ca \  cg ct
//     ga gc \  gt
//     ta tc tg \
//
//     Command: 
//     RMATMOD \ @ac @ag @at @ca \ @cg @ct @ga @gc \ @gt @ta @tc @tg \
//     
//     where '@' is the corresponding value in the Q matrix of the global Q.
//     For example, the CpG example, where CpGs change at 10X the normal rate
//     RMATMOD \ @ @ @ @ \ @ @*10 @ @ \ @ @ @ @ \
//     RMATMOD \ @ @ @ @ \ @ @ @*10 @ \ @ @ @ @ \
//
//   b. Input of a new matrix, specify parameters of new matrix, frequencies
//     Command:
//     NEWRMAT F84(kappa=2) 0.1 0.4 0.4 0.1
//
//     where F84 is the matrix, which requires kappa to be specified, and the
//     nucleotide frequencies are given to the matrix in the order
//       pi_A   pi_C   pi_G   pi_T
//
//     Note that only matrices that allow the frequencies to be specified are
//     accepted. Thus, to specify the JC69 matrix, use F81 0.25 0.25 0.25 0.25,
//     and likewise for other nucleotide (or AA, or even possibly codon, later)
//     matrices.
//
// (2) Rule for application of the dependency matrix
//   a. Markovian model
//     Command:
//     MARKOV(ORDER) X*|X* <PATTERN_FREQ>
//
//     where X is any residue and the X* term refers to the regular expression
//     of 0 or more X's. In this case, the X's refer to the pattern for which
//     the Markov rule takes effect. The symbol '|' thus refers to the position
//     for which we are calculating the rate. ORDER refes to the order of the
//     Markov process. Thus, this Markov model allows for dependency for 
//     preceding and proceeding patterns. Finally, PATTERN_FREQ is an optional
//     specification that supplies the stationary frequency of the X* pattern,
//     for cases in which the dependency is not applicable (see example b.).
//
//
// (*) EXAMPLES
// Following are examples for the usage of the items described above.
//
//   a. CpG model, where CpGs change at a rate of 10X faster than the normal
//      rate.
//     DEPENDENCY MARKOV(1) C| RMATMOD \@@@@\@@*10@@\@@@@\
//     DEPENDENCY MARKOV(1) |G RMATMOD \@@@@\@@@*10@\@@@@\
//
//     - DEPENDENCY marks the beginning of a new dependency rule. 
//
//     - MARKOV(1) C| specifies a 1st-order Markov process, where if the current
//     position is preceded by a C, apply the dependency matrix. In this case, 
//     if this is true, then @ct is multiplied by 10.
//
//     - MARKOV(1) |G specifies a 1st-order Markov process, where if the current
//     position is followed by a G, apply the dependency matrix. In this case,
//     if this is true, then @ga is multiplied by 10.
//
//     -- NOTE that if the C is not followed by a G, the rates for the position
//     remain the same as the global model.
//
//
//   b. 3rd-order Markov process applied globally on the sequence
//     DEPENDENCY MARKOV(3) AAA| NEWRMAT FREQ=0.01 F84(kappa=2) 0.1 0.3 0.5 0.1
//     DEPENDENCY MARKOV(3) AAC| NEWRMAT FREQ=0.02 F84(kappa=2) 0.2 0.2 0.2 0.4
//     DEPENDENCY MARKOV(3) AAG| NEWRMAT FREQ=0.01 F84(kappa=2) 0.1 0.1 0.1 0.7
//         .         .        .    .        .           .        .   .   .   .
//         .         .        .    .        .           .        .   .   .   .
//         .         .        .    .        .           .        .   .   .   .
//     DEPENDENCY MARKOV(3) TTT| NEWRMAT FREQ=0.05 F84(kappa=2) 0.3 0.3 0.3 0.1
//
//     - This specifies all 64 3mers of nucleotides. Each 3mer defines a
//     a dependency in this model, where the dependency is the transition
//     probabilities from an F84 matrix with kappa=2, and nucleotide freqs
//     pi_A pi_C pi_G pi_T. The first 3mer on the 5' side of the sequence is
//     assigned a probability as assigned by FREQ, otherwise the probability
//     of a site depends on the 3mer preceding it, which defines the new rate
//     matrix.
//

#ifndef _DEPENCENCY_H_
#define _DEPENDENCY_H_

#include "tree.h"
#include "model.h"

extern int changed_site;
extern bool fast_simulation;
extern bool optimize;

using namespace std;

class TNode;

class LookUp : private Counter<LookUp>
{
public:
	double 			value;
	double 			inverse;
	vector<double> 	tau_ij;
	vector<double> 	Qmat;
	using Counter<LookUp>::howMany;

	LookUp() : value(0), inverse(0) { }
	LookUp (double value) : value(value), inverse(1.0/value) { }
};

class contextDependence : private Counter<contextDependence>
{
public:
	vector<int>   index_position_multiplier;
	vector<vector<int> > index_offset;	// For TauIJ: calculate indices for neighboring residues.
	contextDependence(int dep_order) : order(dep_order) { }
	using Counter<contextDependence>::howMany;

	// Accessors
	double return_tuplet_pi(int index, bool inverse = false) { return ( (inverse) ? (1.0/tuplet_pi.at(index)) : tuplet_pi.at(index) ); }
	double return_lookup_table_value(string sequence) { return lookup_table2[sequence]->value; }
	double return_lookup_table_inverse(string sequence) { return lookup_table2[sequence]->inverse; }
	double return_order() { return order; }
	
	void 				set_sequence_indices(TNode *node);
	void 				reset_sequence_indices(TNode *node, int event_site, string subst_event);
	int					getOffset(int env, short i, short j);
	void				setOffset();
	void 				set_lookup_table();
	void				calculate_tau_ij();
	string				lookup_table_sequence(int idx, int sequence_length);
	int 				lookup_table_index(vector<short> sequence);
	int 				lookup_ktuplet_index(string sequence);
	void 				outputDependencies(string& outfile_name_root);
	void 				readDependencies(string& file);
	void				generateDependencies(double dependence_superscript);
	void				allocate_lookup_context_vector();
	vector<short>		int2iSG_seq(int index, int sequence_length);
	int 				sequence_specific_index_offset(vector<short> sequence);
	int					lookup_context_index(vector<short> sequence, bool sequence_end = false)
						{ return sequence_specific_index_offset(sequence)+context_specific_index_offset.at(sequence.size() - ( (!sequence_end) ? (order+1) : 0 ) ); }
	double 				markov_ratio(vector<short>& i, vector<short>& j, size_t sequence_position);
	double				lt_markov_ratio(double env_index, int i_seq_index, int j_seq_index);
	void				set_Qmat(TNode *node);
	double				TauIJ(TNode *node, int env_index, int i_seq_index, short residue_i, short residue_j);

	// Print functions
	void				report_tuplet_pi();

private:
	int order;		// The order of dependence (Markovian). //
	vector<vector<LookUp*> > lookup_table;
	map<string, LookUp*> 	lookup_table2;
	// Triplet frequencies
	vector<double>			tuplet_pi;
	vector<int>				context_specific_index_offset;	// see allocate_lookup_context_vector();
	vector<LookUp*>			lookup_context;
};

class Dependency : private Counter<Dependency>
{
public:
	vector<vector<double> > Qij_D;
	vector<vector<double> > pi;
	vector<vector<double> > pi_inv;
	vector<double>			triplet_pi;
	vector<double>			triplet_pi_inv;
	contextDependence		context;
	using Counter<Dependency>::howMany;

	Dependency(double dep_order, double dependence_superscript, string& output_files);
	Dependency(double dep_order, string& dependency_file);
	void buildDependenceStructures();
};

double inline contextDependence::lt_markov_ratio (
												  double env_index,
												  int i_seq_index,
												  int j_seq_index
												 )
{
	//cerr << "  contextDependence::lt_markov_ratio -> dimensions of lookup table -> ";
	//cerr << env_index << "/" << lookup_table.size() << "   ";
	//cerr << lookup_table.at(env_index).size() << "              j:" << j_seq_index << " i:" << i_seq_index << endl;
	double ratio = lookup_table.at(env_index).at(j_seq_index)->value * lookup_table.at(env_index).at(i_seq_index)->inverse;
	return ratio;
}

double inline contextDependence::markov_ratio (
								 			   vector<short>& i,
								 			   vector<short>& j,
								  	  		   size_t sequence_position
								 	   		  )
{
	double ratio;
	bool profile=false;
	bool seq_end=false;

	if (i.size() < order*2+1)
		if(sequence_position>order) seq_end = true;

	ratio 
	= lookup_context.at(lookup_context_index(j, seq_end))->value 
	* lookup_context.at(lookup_context_index(i, seq_end))->inverse;

	//cerr << "j("; for (vector<short>::iterator it = j.begin(); it != j.end(); ++it) cerr << stateCharacters.at(*it);
	//cerr << ")(val=" << lookup_context.at(lookup_context_index(j, seq_end))->value;
	//cerr << ")  i("; for (vector<short>::iterator it = i.begin(); it != i.end(); ++it) cerr << stateCharacters.at(*it);
	//cerr << ")(inv=" << lookup_context.at(lookup_context_index(i, seq_end))->inverse;
	//cerr << ")" << endl;

	return ratio;
}


#endif