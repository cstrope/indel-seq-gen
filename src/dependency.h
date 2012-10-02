#ifndef _DEPENCENCY_H_
#define _DEPENDENCY_H_

#include "model.h"
#include "tree.h"

extern int changed_site;
extern bool fast_simulation;
extern bool optimize;

using namespace std;

class TNode;
class RateMatrix;

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
	~contextDependence();
	using Counter<contextDependence>::howMany;

	// Accessors
	double return_tuplet_pi(int index, bool inverse = false) { return ( (inverse) ? (1.0/tuplet_pi.at(index)) : tuplet_pi.at(index) ); }
	double return_lookup_table_value(string sequence) { return lookup_table2[sequence]->value; }
	double return_lookup_table_inverse(string sequence) { return lookup_table2[sequence]->inverse; }
	double return_lt_value(int env, int idx) { return lookup_table.at(env).at(idx)->value; }
	double return_lt_inverse(int env, int idx) { return lookup_table.at(env).at(idx)->inverse; }
	double return_order() { return order; }

	void 				set_neutral_lookup_vector(inClade *environment);
	void 				set_sequence_indices(TNode *node, int block_size);
	void 				reset_sequence_indices(TNode *node, int event_site, string subst_event);
	int					getOffset(int env, int codon_pos, short i, short j);
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
	void				get_first_order_indices (int tuplet, int *index1, int *index2);

	// Print functions
	void				report_tuplet_pi();

private:

	// The order of Markov dependence.
	int order;		

	// 
	vector<vector<LookUp*> > lookup_table;
	map<string, LookUp*> 	lookup_table2;

	// Triplet frequencies
	vector<double>			tuplet_pi;
	vector<int>				context_specific_index_offset;	// see allocate_lookup_context_vector();
	//vector<LookUp*>			lookup_context;
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

	Dependency(double dep_order, double dependence_superscript, string& output_files);	// 3oMM FWD
	Dependency(double dep_order, string& dependency_file);	// 3oMM EPC
	Dependency(double dep_order, string &dependency_counts_file, string $neutral_counts_file);	// HDS FWD/EPC?
	Dependency(double dep_order, inClade *global_environment);	// For cases where Neutral model = independent sites model.
	Dependency(double dep_order, int block_size, string& file);
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
	cerr << "Ratio: " << lookup_table.at(env_index).at(j_seq_index)->value << " / " << lookup_table.at(env_index).at(i_seq_index)->value << endl;
	double ratio = lookup_table.at(env_index).at(j_seq_index)->value * lookup_table.at(env_index).at(i_seq_index)->inverse;
	return ratio;
}

/*double inline contextDependence::markov_ratio (
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
*/

#endif