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

#include "inTree.h"

using namespace std;

inTree::inTree(const seqGenOptions *options, int num, inClade *globalEnv)
	   : randomInvariableAssignment( ((options->invariableSites) ? true : false) ),
	     scaleBranches((options->default_branchScale != 1 && options->default_branchScale > 0) ? true : false),
	     scaleTree((options->default_branchScale < 0) ? true : false),
	     ma_columnCollapseMethod('c'),
	     branchScale(options->default_branchScale),
	     proportion_invariable(options->default_proportion_invariable),
	     partitionRate(1.0),
	     ma_numSeqsToUse(-1),
	     rootSeqType(RANDOM),
	     partitionLength(0),
	     treeNum(num),
	     codon_offset(0),
		 label(""),
	     ma_range(""),
	     my_tree(new TTree())
{
	inClade *rootEnv = new inClade("Tree Parameters", options);
	rootEnv->rootEnvSetup(globalEnv);
	my_tree->treeEnv.push_back(rootEnv);
	if(options->invariableSites) my_tree->treeEnv.front()->invariableSites = true;
	else my_tree->treeEnv.front()->invariableSites = false;
	filename.clear();
	invariable_by_partition.clear();
	rootseq_by_partition.clear();
	motif_specs.clear();
	tree.clear();
	input_MA.clear();
	input_MA_motifs.clear();
}

void inTree::perturbTree(double perturb_value)
{
	setPerturb(my_tree->root->branch1, perturb_value);
	setPerturb(my_tree->root->branch2, perturb_value);
}

void inTree::setPerturb(TNode *node, double perturb_value)
{
	node->branch->perturbation = 1.0/perturb_value + (perturb_value - 1.0/perturb_value)*rndu();
	if (node->tipNo == -1) {
		setPerturb(node->branch1, perturb_value);
		setPerturb(node->branch2, perturb_value);
	}
}

bool inTree::CheckPhylogeneticAncestralNodes(string& node_names, int action) 
{
	bool write = true;

	if(action == SET) {
		for (int i = 0; i < my_tree->numTips; i++) node_names += "1";
		node_names += ",";
	}

	CPAN(my_tree->root->branch1, node_names, action, &write);
	if(write) CPAN(my_tree->root->branch2, node_names, action, &write);

	return write;
}

void inTree::SortTaxonTips(vector<string>& otu_names) 
{
	TNode **permute_tips;
	int which_otu;

	permute_tips = (TNode**)malloc(my_tree->numTips * sizeof(TNode*));

	// For the sake of printing out, need all tip lists to be in the same order.
	int i = 0;
	for (vector<TNode*>::iterator itn = my_tree->tips.begin(); itn != my_tree->tips.end(); itn++, i++) {
		permute_tips[i] = (*itn);
	}

	i = 0;
	for (vector<TNode*>::iterator itn = my_tree->tips.begin(); itn != my_tree->tips.end(); itn++, i++) {
		which_otu = 0;
		while ( (my_tree->names.at(i)).compare(otu_names.at(which_otu)) != 0 && which_otu < my_tree->numTips)
			which_otu++;
		if ( (my_tree->names.at(i)).compare(otu_names.at(which_otu)) != 0) {
			cerr << "Taxon " << my_tree->names.at(i) << " not found." << endl;
			exit(EXIT_FAILURE);
		}

		(*itn)=permute_tips[which_otu];
	}
	
	free(permute_tips);
}

void inTree::CPAN(TNode *node, string& node_names, int action, bool *inval) 
{
	
	if(action == SET) {
		for (int i = 0; i < my_tree->numTips; i++) node_names += node->bipartition.at(i);
		node_names += ",";
		if(node->tipNo == -1) {
			CPAN(node->branch1, node_names, action, inval);
			CPAN(node->branch2, node_names, action, inval);
		}
	} else {
		string tmp = "";
		size_t found;
		for (int i = 0; i < my_tree->numTips; i++) tmp += node->bipartition.at(i);		
		found = node_names.find(tmp);
		if(found == string::npos) *inval = false;
		else {
			if(node->tipNo == -1) {
				CPAN(node->branch1, node_names, action, inval);
				CPAN(node->branch2, node_names, action, inval);
			}
		}
	}
}

double inTree::CalculatePathLengths() 
{
	my_tree->root->branch->branch1_max_path = SearchPath(my_tree->root,my_tree->root->branch1,DIST_CALC);
	my_tree->root->branch->branch1_time_relative_length = my_tree->root->branch->length1;

	my_tree->root->branch->branch2_max_path = SearchPath(my_tree->root,my_tree->root->branch2,DIST_CALC);
	my_tree->root->branch->branch2_time_relative_length = my_tree->root->branch->length2;

	return ( (my_tree->root->branch->branch1_max_path > my_tree->root->branch->branch2_max_path) ? 
			  my_tree->root->branch->branch1_max_path : my_tree->root->branch->branch2_max_path );
}

void inTree::Print_Trace_Header(stringstream& header, vector<string> otu_names, seqGenOptions* options) 
{
	int numTaxa = otu_names.size();
	vector<string>::iterator taxon_names;
	taxon_names = otu_names.begin();

	if (options->writeAncestors) {
		header << ++numTaxa << " ";
		for (vector<bool>::iterator it = my_tree->root->bipartition.begin(); it != my_tree->root->bipartition.end(); it++)
			header << (*it);
		header << endl;
	}
	
	Print_Branch(my_tree->root->branch1,header, taxon_names, numTaxa, options);
	Print_Branch(my_tree->root->branch2,header, taxon_names, numTaxa, options);
}

void inTree::Print_Branch(TNode *node, stringstream& header, vector<string>::iterator& taxon_names, int& num, seqGenOptions *options) 
{

	if (node->tipNo == -1) {
		if (options->writeAncestors) {
			header << ++num << " ";
			for (vector<bool>::iterator it = node->bipartition.begin(); it != node->bipartition.end(); it++)
				header << (*it);
			header << endl;
		}
		Print_Branch(node->branch1, header, taxon_names, num, options);
		Print_Branch(node->branch2, header, taxon_names, num, options);
	} else {
		header << (*taxon_names) << " ";
		for (vector<bool>::iterator it = node->bipartition.begin(); it != node->bipartition.end(); it++)
			header << (*it);
		header << endl;
		taxon_names++;
	}

}

double inTree::SearchMinBranch(TNode *node, double curr_min)
{
	double min1, min2;

	if (node->branch->anc_length() < curr_min) curr_min = node->branch->anc_length();
	
	if (node->tipNo==-1) {
		min1 = SearchMinBranch(node->branch1, curr_min);
		min2 = SearchMinBranch(node->branch2, curr_min);

		min1 = SearchMinBranch(node->branch1, curr_min);
		min2 = SearchMinBranch(node->branch2, curr_min);

		if (min1 < curr_min) curr_min = min1;
		if (min2 < curr_min) curr_min = min2;
	} 
	
	return curr_min;
}

double inTree::SearchPath(TNode *anc, TNode *des, int flag) 
{
	if (flag == DIST_CALC)
		des->branch->branch0_time_relative_length = des->branch->anc_length();

    if (des->tipNo==-1) { 
		if(flag == DIST_CALC) 
			des->branch->branch1_time_relative_length = des->branch->length1;
		des->branch->branch1_max_path = SearchPath(des,des->branch1,flag);
		if(flag == DIST_CALC) 
			des->branch->branch2_time_relative_length = des->branch->length2;
		des->branch->branch2_max_path = SearchPath(des,des->branch2,flag);
	} else {
		des->branch->branch1_max_path = 0;
		return des->branch->branch0_time_relative_length;
	}

	return ( (des->branch->branch1_max_path > des->branch->branch2_max_path) ? 
			 (des->branch->branch0_time_relative_length + des->branch->branch1_max_path) : 
			 (des->branch->branch0_time_relative_length + des->branch->branch2_max_path) 
		   );
}

void inTree::NudgeBranches(double step_size, int simulation_type)
{
	my_tree->root->branch->length1 = DoNudge(my_tree->root->branch1, step_size, simulation_type);
	my_tree->root->branch->length2 = DoNudge(my_tree->root->branch2, step_size, simulation_type);
}

double inTree::DoNudge(TNode *node, double step_size, int simulation_type)
{
	double nudge, evenQ;
	int isInt;
	double *P;

	if (simulation_type == DISCRETE_EVOLUTIONARY_STEPS) { 
		P=&node->branch->length0;
	} else if (simulation_type == TIME_RELATIVE_STEPS) {
		P=&node->branch->branch0_time_relative_length;
	} else if (simulation_type == GILLESPIE) {
		P=&node->branch->length0;
		return *P;
	}
	else {
		cerr << "Undefined simulation_type in inTree::DoNudge." << endl;
		exit(EXIT_FAILURE);
	}
	
	isInt = (int)(*P / step_size);
	if (*P - isInt * step_size == 0) nudge = 0;
	else { 
		nudge = step_size - (*P - isInt * step_size); 
		isInt++; 
	}
	
	*P += nudge;
	evenQ = *P / step_size;

	if(node->tipNo == -1) {
		if (simulation_type == DISCRETE_EVOLUTIONARY_STEPS) {
			node->branch1->branch->length0 -= nudge;
			node->branch2->branch->length0 -= nudge;
			node->branch->length1 = DoNudge(node->branch1, step_size, simulation_type);
			node->branch->length2 = DoNudge(node->branch2, step_size, simulation_type);
		} else if (simulation_type == TIME_RELATIVE_STEPS) {
			node->branch1->branch->branch0_time_relative_length -= nudge;
			node->branch2->branch->branch0_time_relative_length -= nudge;
			node->branch->branch1_time_relative_length = DoNudge(node->branch1, step_size, simulation_type);
			node->branch->branch2_time_relative_length = DoNudge(node->branch2, step_size, simulation_type);		
		}
	} else {
		//if(evenQ - (int)evenQ > PREC_EPS) cerr << "Node " << node->tipNo << " not divisible by nudge" << endl;	
	}

	node->numEventsToSimulate = isInt;
	if (simulation_type == DISCRETE_EVOLUTIONARY_STEPS) {
		node->BL_step_size = step_size;
	} else if (simulation_type == TIME_RELATIVE_STEPS) {
		node->BL_step_size = (node->branch->anc_length() / node->branch->branch0_time_relative_length) * step_size;
	}
	
	return *P;
}

void inTree::TimeScale_Tree(double max_path, TNode *node) 
{
	double multiplier;
	
	if(node->tipNo==-1) {
		// Left branch
		multiplier = max_path / node->branch->branch1_max_path;
		node->branch->branch1_time_relative_length *= multiplier;
		BranchScale_Tree(node,node->branch1, multiplier);
		node->branch->branch1_max_path = SearchPath(node,node->branch1,TIME_SCALE);
		// Right branch
		multiplier = max_path / node->branch->branch2_max_path;
		node->branch->branch2_time_relative_length *= multiplier;
		BranchScale_Tree(node,node->branch2, multiplier);
		node->branch->branch2_max_path = SearchPath(node,node->branch2,TIME_SCALE);

		TimeScale_Tree(node->branch->branch1_max_path - node->branch->branch1_time_relative_length, node->branch1);
		TimeScale_Tree(node->branch->branch2_max_path - node->branch->branch2_time_relative_length, node->branch2);
	}
}

void inTree::BranchScale_Tree(TNode *anc, TNode *des, double scale) 
{
	des->branch->branch0_time_relative_length *= scale;
	des->trDistanceFromRoot = anc->trDistanceFromRoot + des->branch->branch0_time_relative_length;
	des->DistanceFromRoot   = anc->DistanceFromRoot   + des->branch->length0;

    if (des->tipNo==-1) { 
		des->branch->branch1_time_relative_length *= scale;
		des->branch->branch2_time_relative_length *= scale;
		BranchScale_Tree(des, des->branch1, scale);
		BranchScale_Tree(des, des->branch2, scale);
	}
}

void inTree::calcDistancesFromRoot()
{
	calcDist(my_tree->root, my_tree->root->branch1);
	calcDist(my_tree->root, my_tree->root->branch2);
}

void inTree::calcDist(TNode *anc, TNode *des) 
{
	des->DistanceFromRoot   = anc->DistanceFromRoot   + des->branch->length0;

    if (des->tipNo==-1) { 
		calcDist(des, des->branch1);
		calcDist(des, des->branch2);
	}
}

//////////
/// 10-7-2010: Branch scaling has been done after timescaling, which is erroneous, since longest
/// path will be incorrectly computed. This function applies branchScaling factors to the branch
/// lengths of the nodes before any calculations take place. Local branch scales take precedence
/// over the global scales, if both are present.
//////////
void inTree::apply_branchScales(double scalar)
{
	my_tree->root->branch->length1 = apply2Subtree(my_tree->root->branch1, scalar);
	my_tree->root->branch->length2 = apply2Subtree(my_tree->root->branch2, scalar);
	my_tree->root->branch->length0 = my_tree->root->branch->length1;
}

double inTree::apply2Subtree(TNode *node, double scalar)
{
	if (node->tipNo == -1) {
		node->branch->length1 = apply2Subtree(node->branch1, scalar);
		node->branch->length2 = apply2Subtree(node->branch2, scalar);
	}
	if (branchScale != 1) node->branch->length0 *= scalar;

	return node->branch->length0;
}

void inTree::Define_Ancestors() 
{
	int node_num = -1;

	my_tree->root->mytipNo = node_num++;
	Define_Anc(my_tree->root, my_tree->root->branch1, &node_num);
	Define_Anc(my_tree->root, my_tree->root->branch2, &node_num);
	if (!my_tree->rooted)
		Define_Anc(my_tree->root, my_tree->root->branch0, &node_num);
}

void inTree::Define_Anc(TNode *anc, TNode *des, int *node_num) 
{
	des->anc = anc;
	des->mytipNo = (*node_num)++;
	if (des->branch1 != NULL) {
		Define_Anc(des, des->branch1, node_num);
		Define_Anc(des, des->branch2, node_num);
	}
}

void inTree::Define_Bipartitions(vector<string>& otu_names) 
{
	vector<bool>::iterator it;

	my_tree->root->tipNo = ROOT_FLAG;
	for (list<TNode*>::iterator itn = my_tree->nodeList.begin(); itn != my_tree->nodeList.end(); itn++) {
		it = (*itn)->bipartition.begin();
		(*itn)->bipartition.insert(it, my_tree->numTips, false);
	}

	// Define ancestors for the Push_Bipartition recursive routine.	
	my_tree->root->bipartition.assign(my_tree->numTips, true);
	my_tree->root->nodeEnv->bipartition.assign(my_tree->numTips, true);

	int which_otu;
	int i = 0;

	for(vector<TNode*>::iterator itn = my_tree->tips.begin(); itn != my_tree->tips.end(); itn++, i++) {
		which_otu=0;
		while ( (my_tree->names.at(i)).compare(otu_names.at(which_otu)) != 0 && which_otu < my_tree->numTips) 
			which_otu++;
		if ( (my_tree->names.at(i)).compare(otu_names.at(which_otu)) != 0) {
			cerr << "Taxon " << my_tree->names.at(i) << " not found." << endl;
			exit(EXIT_FAILURE);
		}
		(*itn)->bipartition.at(which_otu) = true;
		Push_Bipartition((*itn)->anc, which_otu);
	}

	Define_Clade_Bipartitions();
	my_tree->root->tipNo = -1;
}

void inTree::Define_Clade_Bipartitions()
{
	// Define clade bipartitions also.
	for (list<TNode*>::iterator itn = my_tree->nodeList.begin(); itn != my_tree->nodeList.end(); itn++) {
		if ( !(*itn)->clade_label.empty() )
			(*itn)->nodeEnv->bipartition = (*itn)->bipartition;
	}

	for (list<inClade*>::iterator it = my_tree->treeEnv.begin(); it != my_tree->treeEnv.end(); it++) {
		for (list<inMotif*>::iterator it2 = (*it)->my_motifs.begin(); it2 != (*it)->my_motifs.end(); it2++) {
			(*it2)->bipartition = (*it)->bipartition;
		}
	}
}

void inTree::Push_Bipartition(TNode *node, int which_tip)
{
	if(node->tipNo == ROOT_FLAG) ;
	else {	
		node->bipartition.at(which_tip) = true;
		Push_Bipartition(node->anc, which_tip);
	}
}

void inTree::Print_Time_Rel_Tree(ostream& out, bool writeAncestors) 
{
	// Name clades as internal node numbers (ancestors).
	int nodeNo = my_tree->numTips + 2;

	out << "Partition " << treeNum << endl;
	out << "(";
	PrintSubtree(my_tree,my_tree->root->branch1,out, writeAncestors, &nodeNo, true, true);
	out << ",";
	PrintSubtree(my_tree,my_tree->root->branch2,out, writeAncestors, &nodeNo, true, true);
	out << ")" << my_tree->numTips+1 << ";" << endl << endl;
}

void inTree::Print_Tree(ostream& out, bool scale)
{
	int nodeNo = my_tree->numTips + 2;

	if (scale) out << "Partition " << treeNum << endl;
	out << "(";
	PrintSubtree(my_tree,my_tree->root->branch1,out, true, &nodeNo, false, scale);
	out << ",";
	PrintSubtree(my_tree,my_tree->root->branch2,out, true, &nodeNo, false, scale);
	out << ")";
	if (scale) out << my_tree->numTips+1;
	out << ";" << endl << endl;
}

void inTree::PrintSubtree(TTree *tree, TNode *node, ostream& out, bool writeAncestors, int *nodeNo, bool time_rel, bool scale) 
{
	int my_nodeNo;

	if(node->tipNo==-1) {
		my_nodeNo = *nodeNo;
		(*nodeNo)++;
		out << "(";
		PrintSubtree(tree,node->branch1,out, writeAncestors, nodeNo, time_rel, scale);
		out << ",";
		PrintSubtree(tree,node->branch2,out, writeAncestors, nodeNo, time_rel, scale);
		out << ")";
		if (scale) out << my_nodeNo;
		out << ":" << ((time_rel) ? node->branch->branch0_time_relative_length : node->branch->anc_length() );
	} else {
		out << my_tree->names.at(node->tipNo) << ":" << ((time_rel) ? node->branch->branch0_time_relative_length : node->branch->anc_length() );
	}

}

void inTree::perform_checks(string input, seqGenOptions *options)
{
	size_t found;
		
	found = input.find('{');
	if (found != string::npos) {
		my_tree->treeEnv.front()->indelFlag = true;
	}
	
	found = input.find("[:");
	if (found != string::npos) {
		//inputRoot = true;
		proportion_invariable = 0;
	}
	
}

string inTree::tree_only(string input)
{
	size_t found, found2, found_tmp, found_tmp2;
	
	found2 = input.rfind(')');
	if (found2 == string::npos) {
		cerr << "Missing closing parentheses in input guide tree." << endl;
		exit(EXIT_FAILURE);
	}
	found = input.find(']');

	//////////
	/// Remove ""
	//////////
	found_tmp = input.find('"');
	if (found_tmp != string::npos) {
		found_tmp2 = input.find(']');
		if (found_tmp2 != string::npos)
			if (found < found_tmp2) 
				found = found_tmp2;
	}	

	//////////
	/// Remove ##
	//////////
	found_tmp = input.find('#');
	if (found_tmp != string::npos) {
		found_tmp2 = input.find('#', found_tmp+1);
		if (found_tmp2 != string::npos) 
			if (found < found_tmp2) 
				found = found_tmp2;

		//////////
		/// Also need to worry about branch scaling, so this snippet of code gets the appropriate variable/value combo.
		//////////
		string opt_str = input.substr(found_tmp + 1, found_tmp2 - found_tmp - 1);
		list<string> tokens = split(opt_str, ",");
		for (list<string>::iterator tok = tokens.begin(); tok != tokens.end(); tok++) {
			trim(*tok);
			if ((*tok)[0] == 'b') {
				double scale = atof((*tok).substr(1).c_str());
				scaleBranches = true;
				branchScale = scale;
			}
		}
	}

	//////////
	/// Remove {}
	//////////
	found_tmp = input.find('}');
	if (found_tmp != string::npos)
		if (found_tmp > found)
			found = found_tmp;

	return input.substr(found+1, found2 - found + 1);

}

void inTree::parseTree(string input, seqGenOptions *options)
{
	size_t found;
	size_t found2;
	stringstream warning;
	int start_pos = 0;

	if (options->treelength_check) {
		tree = tree_only(input);
		//////////
		/// If branchScale is < 0, then the tree is expected to be tree scaled. However, a
		/// treescaling will return a tree with average path length equal to the branchScale
		/// value, which makes a treelength_check pointless. Thus, unset the branchScale, and 
		/// assume the user wants the unscaled tree treelengths.
		//////////
		if (branchScale < 0) { branchScale = 1; }
		my_tree->ReadTree(tree, &start_pos, my_tree, (my_tree->names).empty());	
			return;
	}
	
	//Find the label
	found = input.find('"');
	if (found != string::npos) {
		found2 = input.find('"', found + 1);
		if (found2 != string::npos) {
			label = input.substr(found + 1, found2 - found - 1);
			input.erase(found, found2 - found + 1);
		} else {
			cerr << "Missing closing \"." << endl;
			exit(EXIT_FAILURE);
		}
	}

	//Find indel information
	parseIndel(input, my_tree->treeEnv.front(), options);

	//Find subseq option overrides
	parsePound(input, my_tree->treeEnv.front(), false, options);

	// Find Root Sequence Options
	found = input.find('[');
	if (found != string::npos) {
		string root_info = "";
		if (input.at(found + 1) == ':') {
			my_tree->treeEnv.front()->invariableSites=true;
			options->inputRoot = true;
			if (proportion_invariable == -1) {
				proportion_invariable = 0;
				my_tree->treeEnv.front()->proportion_invariable=0;
			}
			found2 = input.find(']', found + 2);
			if (found2 != string::npos) {
				root_info = input.substr(found + 2, found2 - (found + 2));
				input.erase(found, found2 - found + 1);
				list<string> maOrRoot = split(root_info, "(");
				if (maOrRoot.size() > 1) {
					rootSeqType = MULTIPLE_ALIGNMENT_ROOT;
					randomInvariableAssignment = false;
					filename = maOrRoot.front();
					maOrRoot.pop_front();
					list<string> ma_tmp = split(*maOrRoot.begin(),")");
					list<string> ma_options = split(ma_tmp.front(), ",");
					ma_tmp.pop_front();
					string chk = trim(ma_tmp.front());
					if(chk != "") {
						partitionRate = atof((*ma_tmp.begin()).substr(1,(*ma_tmp.begin()).length()).c_str());
					}
					if (ma_options.front() != "") {
						ma_range = ma_options.front();
						list<string> check = split(ma_range,":");
						if (check.size() != 2) {
							cerr << "Invalid range specified for \"" << filename << "\"" << endl;
							exit(EXIT_FAILURE);
						}
					}
					ma_options.pop_front();
					if (ma_options.front() != "") {
						ma_numSeqsToUse = atoi(ma_options.front().c_str());
					}
					ma_options.pop_front();
					if (ma_options.front() != "") {
						ma_columnCollapseMethod = ma_options.front()[0];
					}

					ifstream root_input;
					filename = addFilePath(filename, options->filePath);
					root_input.open(filename.c_str(), ifstream::in);
					string readfile, readfile2;
					list<string> tmp;
					// Seems strange, but NEXUS reports conserved sites, but before constructed.
					if (options->fileFormat == NEXUSFormat) {
						if (root_input.good()) {
							getline(root_input, readfile);
							invariable_by_partition = trim(readfile);
							getline(root_input, readfile);
							readfile2 = trim(readfile);
							tmp = split(readfile2, " \t", 1);
							if(tmp.size() == 1) {
								cerr << filename << " format is bad.";
								cerr << "Format is a taxon name and sequence, separated by spaces or tabs." << endl;
								exit(EXIT_FAILURE);
							}
							rootseq_by_partition = trim(tmp.back());
							partitionLength = rootseq_by_partition.size();
						}
					}
					Read_MA(filename, options);
				} else {
					rootSeqType = SINGLE_ROOT;
					randomInvariableAssignment = false;
					ifstream root_input;
					list<string> ma_tmp = split((*maOrRoot.begin()), ",");
					if (ma_tmp.size() == 2) {
						filename = ma_tmp.front();
						partitionRate = atof(ma_tmp.back().c_str());
					} else if (ma_tmp.size() == 1) {
						filename = maOrRoot.front();
					} else {
						cerr << "Error: Unknown root sequence input format.\n\n";
						exit(EXIT_FAILURE);
					}
					filename = addFilePath(filename, options->filePath);
					root_input.open(filename.c_str(), ifstream::in);
					string readfile;
					if (root_input.good()) {
						getline(root_input, readfile);
						partitionLength = /*maxPartitionLength =*/ atoi(readfile.c_str());
						getline(root_input, readfile);
						invariable_by_partition = trim(readfile);
						for (string::iterator it = invariable_by_partition.begin(); it != invariable_by_partition.end(); it++)
							if (!((*it) == '0' || (*it) == '1' || (*it) == '2' || (*it) == '3')) {
								cerr << "Invalid character in root sequence input invariable array: " << (*it) << endl;
								cerr << "NOTE: This error may be caused if the invariable array is not included in root sequence input file." << endl;
								exit(EXIT_FAILURE);
							}
						getline(root_input, readfile);
						rootseq_by_partition = trim(readfile);
						if(rootseq_by_partition.size() != (size_t)partitionLength || invariable_by_partition.size() != (size_t)partitionLength) {
							cerr << filename << " format problem: Stated partitionLength (" << partitionLength
								 << ") not equal to partition lengths of the root sequence or invariable arrays. ("
								 << rootseq_by_partition.size() << " and " << invariable_by_partition.size() << ")" << endl;
							exit(EXIT_FAILURE);
						}
						input_MA.push_back(rootseq_by_partition);	// For motifs, store original.
						Read_MA(filename, options);
					} else {
						cerr << "Cannot open " << maOrRoot.front() << " for reading." << endl;
						exit(EXIT_FAILURE);
					}
					root_input.close();
				}
			} else {
				cerr << "Missing ]." << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			found2 = input.find(']', found + 1);
			if (found2 != string::npos)  {
				rootSeqType = RANDOM;
				root_info = input.substr(found + 1, found2 - (found + 1));
				input.erase(found, found2 - found + 1);
				list<string> ma_tmp = split(root_info, ",");
				if(ma_tmp.size() == 2) {
					partitionLength = atoi(ma_tmp.front().c_str());
					partitionRate = atoi(ma_tmp.back().c_str());
				} else if (ma_tmp.size() == 1) {
					partitionLength = atoi(root_info.c_str());
				} else {
					cerr << "Error: Unknown root sequence input format.\n\n";
					exit(EXIT_FAILURE);
				}
			} else {
				cerr << "Missing ]." << endl;
				exit(EXIT_FAILURE);
			}
		}
	} else {
		rootSeqType = RANDOM;
		partitionLength = 500;
		warning << "Tree " << treeNum << ": No root sequence length or file detected. Assumed a random, length 500 sequence." << endl;
		options->SpoolWarnings(warning.str());
	}
	
	//////////
	/// Find the tree 	
	//////////
	found = input.find('(');
	if (found != string::npos) {
		found2 = input.rfind(')');
		if (found2 != string::npos) {
			tree = input.substr(found, found2 - found + 1);
			input.erase(found, found2 - found + 2);
		}
	}
	
	//my_tree->treeEnv.front()->Print_Environment();
	
	//Check if there is anything left in input
	trim(input);
	if (input != "") {
		cerr << "Unrecognized tree format: " << input << endl;
		exit(EXIT_FAILURE);
	}

	my_tree->ReadTree(tree, &start_pos, my_tree, (my_tree->names).empty());	
	my_tree->root->nodeEnv = my_tree->treeEnv.front();
}

void inTree::parseLineages(seqGenOptions *options) 
{
	fstream file;
    string line;
	size_t found, found1, found2;
	stringstream warning;
	string remove_whitespace = "";
	bool lineage_read_done = false;
	bool motif_read_done = false;
	string all_clades_string = "";
	string all_motifs_string = "";
	string line_append = "";
	bool doMotifs = true, doLineages = true;

	file.open(options->lineageSpecificFile.c_str(), ios::in);
	if (!file.is_open()) {
		cerr << "Unable to open motif file \"" << options->lineageSpecificFile << "\": " << endl;
		exit(EXIT_FAILURE);
	}

	while (!file.eof()) {
    	getline(file, line);
		line_append += line;
	}
	file.close();

	remove_whitespace = line_append;
	remove_whitespace = Remove_Whitespace(remove_whitespace);
	
	found = remove_whitespace.find("MOTIFS");
	found1 = remove_whitespace.find("LINEAGES");
	
	//////////
	/// These two if stmts remove any mention of the MOTIF stuff from remove_whitespace. They will
	/// be taken care of in line_append.
	//////////
	if (found != string::npos && found1 != string::npos) {
		if (found > found1) remove_whitespace.erase(found);
		else remove_whitespace.erase(0,found1);
	}

	if (found1 == string::npos) {
		if (!remove_whitespace.empty()) {
			warning << "Lineage file " << options->lineageSpecificFile.c_str() << " not empty after reading lineages and motifs.";
			warning << "  Content = \"" << remove_whitespace << "\"" << endl;
			options->SpoolWarnings(warning.str());
			warning.flush();
		}
		remove_whitespace.clear();
	}

	found = remove_whitespace.find("LINEAGES={");
	if (found == string::npos) doLineages = false;
  	else {
  		remove_whitespace.erase(found,10);
		if(found != string::npos && remove_whitespace.at(found) != '}') {	// Then we have lineage constraints.
			// Find the end of the lineage file. Marked by ";}" in remove_whitespace.
			found1 = 0;
			do {
				found1 = remove_whitespace.find(';', found1+1);
				if (found1 != string::npos) {
					if (remove_whitespace.at(found1+1) == '}') {	// End of lineages.
						lineage_read_done = true;
						all_clades_string = remove_whitespace.substr(found, found1);
					}
				}
			} while(found1 != string::npos && !lineage_read_done);
			remove_whitespace.erase(found, found1+2-found);
		}
	}

	found = line_append.find("MOTIFS");
	if (found == string::npos) doMotifs = false;
	else {
		line_append.erase(0, found+6);
		found = line_append.find_first_not_of(" \t\f\n\v\r");
		found1 = line_append.find_first_not_of(" \t\v\n\r\f", found+1);
		if (line_append.at(found) != '=' || line_append.at(found1) != '{') {
			cerr << "Error specifying motifs in lineage specification file." << endl;
			cerr << "Format:  MOTIFS = { motifs.... }" << endl;
			exit(EXIT_FAILURE);
		} else line_append.erase(0, found1+1);

		bool done = false;
		while (!done) {
			found1 = line_append.find(';');
			if (found1 == string::npos) {
				cerr << "Cannot specify empty motif file." << endl;
				exit(EXIT_FAILURE);
			}

			found2 = line_append.find_first_not_of(" \t\f\v\n\r", found1+1);
			if (line_append.at(found2) == '}') {		// found motif.
				line_append.erase(found2+1, string::npos);
				done = true;
			}
			if (found2 == string::npos) {
				cerr << "Error parsing motif file. May be missing '}'" << endl;
				exit(EXIT_FAILURE);
			}
			all_motifs_string += line_append.substr(0, found1+1);
			line_append.erase(0,found1+1);
		}
	}

	if (doMotifs) parseMotif(all_motifs_string, options);
	if (doLineages) parseClade(all_clades_string, options);
}

void inTree::parseMotif(string& motif_string, seqGenOptions *options)
{
	size_t found, found1, found_marker, found_name, found_pattern, found_clades;
	string marker, name, pattern, clades;
	inMotif *motif;
	inClade *clade;
	
	list<string> statements = split(motif_string, ";");
	list<string> each_motif;
	string this_motif = "";
	string trimmed;
	string str;

	// PARSE MOTIFS
	for (list<string>::iterator it = statements.begin(); it != statements.end(); it++) {
		trimmed = trim((*it));
		found = trimmed.find(':');
		found1 = trimmed.find('=');
		if (found1 > found) {		// list of clades.
			if (!this_motif.empty()) {
				each_motif.push_back(this_motif);
				this_motif.clear();
			}
		}
		this_motif += trimmed;
	}
	each_motif.push_back(this_motif);

	list<string> each_clade;
	for (list<string>::iterator it = each_motif.begin(); it != each_motif.end(); it++) {
		found_name = (*it).find("NAME");
		statements = split((*it), "=");
		if ( (found_name == string::npos && statements.size() > 3) || statements.size() > 4) {
			cerr << "Unrecognized syntax in motif file: " << endl << (*it) << endl;
			exit(EXIT_FAILURE);
		}
		found_marker = (*it).find("MARKER");
		found_name = (*it).find("NAME");
		found_pattern = (*it).find("PATTERN");

		if (found_marker != string::npos && found_pattern != string::npos) {
			if (found_name != string::npos) {
				// Get marker
				if (found_marker < found_name) {
					if (found_name < found_pattern) marker = (*it).substr(found_marker, found_name-found_marker);
					else marker = (*it).substr(found_marker, found_name-found_marker);
				} else {
					if (found_marker < found_pattern) marker = (*it).substr(found_marker, found_pattern-found_marker);
					else marker = (*it).substr(found_marker);
				}
				// pattern
				if (found_pattern < found_name) {
					if (found_name < found_marker) pattern = (*it).substr(found_pattern, found_name-found_pattern);
					else pattern = (*it).substr(found_pattern, found_name-found_pattern);
				} else {
					if (found_pattern < found_marker) pattern = (*it).substr(found_pattern, found_marker-found_pattern);
					else pattern = (*it).substr(found_pattern);
				}
				// and name
				if (found_name < found_pattern) {
					if (found_pattern < found_marker) name = (*it).substr(found_name, found_pattern-found_name);
					else name = (*it).substr(found_name, found_pattern-found_name);
				} else {
					if (found_name < found_marker) name = (*it).substr(found_name, found_marker-found_name);
					else name = (*it).substr(found_name, found_name-found_pattern);
				}
			} else { 
				if(found_marker < found_pattern) {
					marker = (*it).substr(found_marker, found_pattern-found_marker);
					pattern = (*it).substr(found_pattern);
				} else {
					marker = (*it).substr(found_marker);
					pattern = (*it).substr(found_pattern, found_marker-found_pattern);
				}
			}

			found = (*it).find(marker);
			(*it).erase(found, marker.length());
			marker = Remove_Whitespace(marker);
			motif = new inMotif();
			motif->marker = marker.at(marker.length()-1);

			found = (*it).find(name);
			if (found != string::npos) {
				(*it).erase(found, name.length());
				found = name.find_first_of("=");
				name.erase(0,found+1);
				motif->name = trim(name);
			} else motif->name = "Unnamed";

			char marker_spec = marker.at(marker.length()-1);
			if (isdigit(marker_spec)) motif->sequence_template = true;

			found = (*it).find(pattern);
			(*it).erase(found, pattern.length());
			pattern = Remove_Whitespace(pattern);
			motif->regex = pattern.substr(8);

			found_clades = (*it).find_first_of(':');
			str = (*it).substr(0,found_clades);
			clades = Remove_Whitespace(str);
			each_clade = split(clades,",");
			if (each_clade.size() > 1) {
				for (list<string>::iterator it2 = each_clade.begin(); it2 != each_clade.end(); it2++) {
					if( (*it2).compare("root") == 0) clade = my_tree->treeEnv.front();
					else clade = FindEnvironment(my_tree, (*it2));
					if (clade != NULL) clade->my_motifs.push_back(motif);
				}
			} else {
				if (clades.compare("root") == 0) clade = my_tree->treeEnv.front();
				else clade = FindEnvironment(my_tree, clades);
				if (clade != NULL) clade->my_motifs.push_back(motif);
			}
			(*it).erase(0, found_clades+1);
			(*it) = trim(*it);
			if ( !(*it).empty() ) {
				cerr << "Invalid syntax in motif file. Motif read as: " << endl;
				motif->report();
				cerr << " unread: " << (*it) << endl;
				exit(EXIT_FAILURE);
			}
			motif_specs.push_back(motif);
		} else {
			cerr << "Motifs must have a specified marker and pattern." << endl;
			exit(EXIT_FAILURE);
		}
	} 
}

void inTree::parseClade(string& clade_string, seqGenOptions *options) 
{
	size_t found, found2, found3, found4;
	string clade_str;
	string clade_label;
	inClade *clade;
	stringstream warning;

	// Process each clade.
	list<string> clade_parameters = split(clade_string, ";");
	for(list<string>::iterator it = clade_parameters.begin(); it != clade_parameters.end(); it++) {
		// Gather the label to identify which clade this occurs in:
		found = (*it).find(':');
		if (found == string::npos) {
			cerr << "Error parsing lineage-specific file. May be missing a closed '}'" << endl;
			exit(EXIT_FAILURE);
		} else {
			found2 = (*it).find_first_of("\"{}#");
			if (found2 < found) {
				cerr << "ERROR: You must specify clade membership first in lineages." << endl;
				cerr << "FORMAT:" << endl;
				cerr << "clade_members: \"environment_label\"#substitution#{insertion/deletion};" << endl;
				cerr << "  where clade_members is a comma separated list followed by a colon. Other options do not depend on order." << endl;
				exit(EXIT_FAILURE);
			} else {
				string clade_tmp;
				bool skip;
				list<string> clade_members;
				clade_tmp = (*it).substr(0, found);
				clade_members = split(clade_tmp, ",");
				(*it).erase(0,found+1);
				
				for (list<string>::iterator sit = clade_members.begin(); sit != clade_members.end(); sit++) {
					if ((*sit).find_first_of("(") != string::npos) {
						// Will be a taxon spec.
						string taxon_str;
						taxon_str += (*sit);
						if(sit == clade_members.end()) {
							cerr << "Error processing lineages: Missing ')' in clade specification.";
							exit(EXIT_FAILURE);
						}

						list<string>::iterator sit2 = sit;
						sit2++;
						if (sit2 != clade_members.end()) {
							for (;;) {
								taxon_str += "," + (*sit2);
								found3 = (*sit2).find_first_of(')',0);
								(*sit2).clear();
								sit2++;
								if (found3 != string::npos) {
									break;
								} else if (sit2 == clade_members.end()) {
									cerr << "Error processing lineages: Missing ')' in clade specification.";
									exit(EXIT_FAILURE);
								}
							}
						}
						(*sit).assign(taxon_str);
					} 
					if (sit == clade_members.begin()) {
						clade_label = (*sit);
					}
					else {
						clade_tmp.clear();
						clade_tmp += (*sit) + ":" + (*it);
						clade_parameters.push_back(clade_tmp);
					}
				}
			}
		}

		clade = FindEnvironment(my_tree, clade_label);
		if (clade == NULL) {
			//////////
			/// This little beauty of a code block is for a taxon name only.
			//////////
			found3 = clade_label.find_first_of("(");
			if (found3 != string::npos) {
				string taxon_name;
				string which_trees;

				found4 = clade_label.find_first_of(")");
				taxon_name = clade_label.substr(0,found3);
				which_trees = clade_label.substr(found3+1, found4-(found3+1));

				list<string> tree_nums = split(which_trees,",");
				bool taxon_clade_in_tree = false;
				for (list<string>::iterator sit = tree_nums.begin(); sit != tree_nums.end(); sit++)
					if (this->treeNum == atoi((*sit).c_str())) taxon_clade_in_tree = true;

				if (taxon_clade_in_tree) {
					// Maybe a taxon clade. We'll find out.
					vector<TNode*>::iterator vit = this->my_tree->tips.begin();
					for (vector<string>::iterator sit = this->my_tree->names.begin(); sit != this->my_tree->names.end(); sit++, vit++) {
						if (taxon_name.compare(*sit) == 0) {
							clade = my_tree->AddClade(taxon_name);
							(*vit)->clade_label = taxon_name;
						}
					}
				}		

			//////////
			/// Check to see if it is a valid taxon. If so, it is active for all trees.
			/// If not, it should be a clade name.
			//////////
			} else { 
				vector<TNode*>::iterator vit = this->my_tree->tips.begin();
				bool pushed = false;
				for (vector<string>::iterator sit = this->my_tree->names.begin(); sit != this->my_tree->names.end(); sit++, vit++) {
					if (clade_label.compare(*sit) == 0) {
						clade = my_tree->AddClade(clade_label);
						(*vit)->clade_label = clade_label;
						pushed = true;
					}
				}
				//////////
				/// Not a taxon, so it is a clade. Add this clade to the tree.
				/// Since it was also not found, it is not a motif name, either.
				//////////
				if (!pushed) {
					clade = my_tree->AddClade(clade_label);
				}
			}
		}

		//////////
		/// If it is a taxon, the taxon name will be made into a clade in the above code.
		//////////
		if (clade != NULL) {
			if (clade->processed) {		// Clade specified more than once in lineages.
				warning << "Clade " << clade->clade_name << " specified more than once in lineages.";
				warning << " Using environment of first occurrance in lineages.";
				options->SpoolWarnings(warning.str());
			} else {
				clade->processed = true;	// Mark that clade was found in tree.
				found = (*it).find('"');
				if (found != string::npos) {
					found2 = (*it).find('"', found + 1);
					if (found2 != string::npos) {
						clade->environment_name = (*it).substr(found + 1, found2 - found - 1);
						(*it).erase(found, found2 - found + 1);
					} else {
						cerr << "Missing closing \". in lineage " << clade->clade_name << endl;
						exit(EXIT_FAILURE);
					}
				}

				parseIndel((*it), clade, options);
				parsePound((*it), clade, 1, options);
				// clade->Print_Environment();
			}
		} 
	} // end of clade processing.

	int which_clade;
	if ( (which_clade = CladeUnprocessed(my_tree)) ) {	// Not all clades of tree were processed.		
		warning << "Clade " << CladeName(my_tree, which_clade) << " is an invalid clade for tree " << treeNum << ".";
		warning << " NOTE: Clade names are case-sensitive." << endl;
		options->SpoolWarnings(warning.str());
	}
	// my_tree->report_clades();
}

inMotif *inTree::FindMotif(char marker)
{
	inMotif *thisMotif;
	for (list<inMotif*>::iterator it = motif_specs.begin(); it != motif_specs.end(); it++) {
		if ( (*it)->marker == marker) return (*it);
	}
	return NULL;
}

void inTree::propagateTreePointers()
{
	my_tree->root->nodeEnv = my_tree->treeEnv.front();

	propagateCladePointers(my_tree->root->branch1);
	propagateCladePointers(my_tree->root->branch2);
}

void inTree::propagateCladePointers(TNode *node) 
{
	if ( !(node->nodeEnv = FindEnvironment(my_tree, node->clade_label)))
		node->nodeEnv = node->anc->nodeEnv;

//	node->report(this->my_tree); cerr << " belongs to clade " << node->nodeEnv->environment_name << endl;

	if (node->tipNo == -1) { // Recurse
		propagateCladePointers(node->branch1);
		propagateCladePointers(node->branch2);
	}
}

string inTree::addFilePath(string suffix, string prefix) 
{
	return (prefix + "/" + suffix);
}

void inTree::parseIndel(string& input, inClade *clade, seqGenOptions *options) 
{
	size_t found, found2;
	string insert_lengthDistFile;
	string delete_lengthDistFile;
	string indel_info;

	found = input.find('{');
	if (found != string::npos) {
		clade->indelFlag = true;
		clade->set_in_clade[__indelFlag__] = true;
		found2 = input.find('}', found + 1);
		if (found2 != string::npos) {
			indel_info = input.substr(found + 1, found2 - found - 1);
			input.erase(found, found2 - found + 1);
			list<string> splitcomma = split(indel_info, ",");
			//////////
			/// Make sure it is a number.
			//////////
			//Remove Whitespace.
			trim(splitcomma.front());
			// Verify there is a number and that no non-numbers exist.
			if ( splitcomma.front().find_first_of("0123456789") == string::npos 
				 || splitcomma.front().empty()
				 || splitcomma.front().find_first_not_of("0123456789") != string::npos) {
				cerr << "Incorrect maxIndel specification: \"" << splitcomma.front() << "\"." << endl;
				exit(EXIT_FAILURE);
			}
			//////////
			clade->maxIndel = atoi((*splitcomma.begin()).c_str());
			if(clade->maxIndel == 0) clade->indelFlag = false;	// No indels, simple as that.
			clade->set_in_clade[__maxIndel__] = true;
			splitcomma.pop_front();
			clade->insert_lengthDistribution.assign(clade->maxIndel+1,0);
			clade->delete_lengthDistribution.assign(clade->maxIndel+1,0);
			clade->set_in_clade[__insert_lengthDistribution__] = clade->set_in_clade[__delete_lengthDistribution__] = true;
			clade->set_in_clade[__P_ins___] = clade->set_in_clade[__P_del___] = true;
			found = (*splitcomma.begin()).find('/');
			if (found != string::npos) {
				clade->P_ins_ = atof((*splitcomma.begin()).substr(0, found).c_str());
				clade->P_del_ = atof((*splitcomma.begin()).substr(found + 1).c_str());
			} else clade->P_ins_ = clade->P_del_ = atof((*splitcomma.begin()).c_str());
			splitcomma.pop_front();
			//////////
			/// If splitcomma is empty, then the length distribution files have been left out
			/// and the Chang & Benner length dist files should be used. (4-13-2010)
			//////////
		 	if (!(clade->P_ins_ == 0 && clade->P_del_ == 0) && !splitcomma.empty()) {
				string fileread;
				ifstream dist_input;
//				found = (*splitcomma.begin()).find('/');
				found = (*splitcomma.begin()).find(':');
				if (found != string::npos) {
					insert_lengthDistFile = (*splitcomma.begin()).substr(0, found);
					trim(insert_lengthDistFile);
					insert_lengthDistFile = addFilePath(insert_lengthDistFile, options->filePath);
					dist_input.open(insert_lengthDistFile.c_str(), ifstream::in);
					if (dist_input.good()) {
						getline(dist_input, fileread);
						double sum = 0;
						list<string> temp = split(fileread, ",");
						if ((int)temp.size() < clade->maxIndel) {
							cerr << "Unable to read \"" << insert_lengthDistFile << "\":";
							cerr << " maxIndel is " << clade->maxIndel << " but file supplies only " << temp.size() << " values." << endl;
							exit(EXIT_FAILURE);
						}
						list<string>::iterator it = temp.begin();
						for (int j = 0; j < clade->maxIndel; j++) {
							sum += atof((*it).c_str());
							it++;
						}
						insert_lengthDistFile = "";
						double temp_val;
						stringstream iter_str, temp_val_str;
						it = temp.begin();
						for (int j = 0; j < clade->maxIndel; j++) {
							temp_val = atof((*it).c_str()) / sum;
							iter_str << (j + 1);
							if (temp_val < 0.0001) {
								insert_lengthDistFile += iter_str.str() + ":0." + string(PRECISION - 2, '0');
							} else {
								temp_val_str << temp_val;
								insert_lengthDistFile += iter_str.str() + ":" + temp_val_str.str().substr(0, PRECISION);
							}
							clade->insert_lengthDistribution.at(j+1) = temp_val;
							insert_lengthDistFile += " ";
							iter_str.str("");
							temp_val_str.str("");
							it++;
						}
					} else {
						cerr << "Cannot open " << insert_lengthDistFile << ": Invalid filename." << endl;
						exit(EXIT_FAILURE);
					}
					dist_input.close();
			
					delete_lengthDistFile = (*splitcomma.begin()).substr(found + 1);
					trim(delete_lengthDistFile);
					delete_lengthDistFile = addFilePath(delete_lengthDistFile, options->filePath);
					dist_input.open(delete_lengthDistFile.c_str(), ifstream::in);
					if (dist_input.good()) {
						getline(dist_input, fileread);
						double sum = 0;
						list<string> temp = split(fileread, ",");
						if ((int)temp.size() < clade->maxIndel) {
							cerr << "Unable to read \"" << delete_lengthDistFile << "\":"<< endl;
							cerr << " maxIndel is " << clade->maxIndel << " but file supplies only " << temp.size() << " values." << endl;
							exit(EXIT_FAILURE);
						}
						list<string>::iterator it = temp.begin();
						for (int j = 0; j < clade->maxIndel; j++) {
							sum += atof((*it).c_str());
							it++;
						}
						delete_lengthDistFile = "";
						double temp_val;
						stringstream iter_str, temp_val_str;
						it = temp.begin();
						for (int j = 0; j < clade->maxIndel; j++) {
							temp_val = atof((*it).c_str()) / sum;
							iter_str << (j + 1);
							if (temp_val < 0.0001) {
								delete_lengthDistFile += iter_str.str() + ":0." + string(PRECISION - 2, '0');
							} else {
								temp_val_str << temp_val;
								delete_lengthDistFile += iter_str.str() + ":" + temp_val_str.str().substr(0, PRECISION);
							}
							clade->delete_lengthDistribution.at(j+1) = temp_val;
							delete_lengthDistFile += " ";
							iter_str.str("");
							temp_val_str.str("");
							it++;
						}
					} else {
						cerr << "Cannot open " << delete_lengthDistFile << ": Invalid filename." << endl;
						exit(EXIT_FAILURE);
					}
				} else {
					insert_lengthDistFile = (*splitcomma.begin());
					trim(insert_lengthDistFile);
					insert_lengthDistFile = addFilePath(insert_lengthDistFile, options->filePath);
					dist_input.open(insert_lengthDistFile.c_str(), ifstream::in);
					if (dist_input.good()) {
						getline(dist_input, fileread);
						double sum = 0;
						list<string> temp = split(fileread, ",");
						if ((int)temp.size() < clade->maxIndel) {
							cerr << "Unable to read \"" << insert_lengthDistFile << "\":";
							cerr << " maxIndel is " << clade->maxIndel << " but file supplies only " << temp.size() << " values." << endl;
							exit(EXIT_FAILURE);
						}
						list<string>::iterator it = temp.begin();
						for (int j = 0; j < clade->maxIndel; j++) {
							sum += atof((*it).c_str());
							it++;
						}
						insert_lengthDistFile = "";
						delete_lengthDistFile = "";
						double temp_val;
						stringstream iter_str, temp_val_str;
						it = temp.begin();
						for (int j = 0; j < clade->maxIndel; j++) {
							temp_val = atof((*it).c_str()) / sum;
							iter_str << (j + 1);
							if (temp_val < 0.0001) {
								insert_lengthDistFile += iter_str.str() + ":0." + string(PRECISION - 2, '0');
								delete_lengthDistFile += iter_str.str() + ":0." + string(PRECISION - 2, '0');
							} else {
								temp_val_str << temp_val;
								insert_lengthDistFile += iter_str.str() + ":" + temp_val_str.str().substr(0, PRECISION);
								delete_lengthDistFile += iter_str.str() + ":" + temp_val_str.str().substr(0, PRECISION);
							}
							clade->insert_lengthDistribution.at(j+1) = clade->delete_lengthDistribution.at(j+1) = temp_val;
							insert_lengthDistFile += " ";
							delete_lengthDistFile += " ";
							iter_str.str("");
							temp_val_str.str("");
							it++;
						}
					} else {
						cerr << "Cannot open " << insert_lengthDistFile << ": Invalid filename." << endl;
						exit(EXIT_FAILURE);
					}
				}
				dist_input.close();
			} else {		// using Chang & Benner model
				if ( !splitcomma.empty() ) {
					stringstream warning;
					warning << "Cannot assign a length distribution to the Chang & Benner model of indel creation. Simulation run assuming the Chang and Benner model." << endl;
					options->SpoolWarnings(warning.str());
				}
				double sum = 0.0;
				for(int j=0; j<clade->maxIndel; j++) {
					clade->insert_lengthDistribution.at(j+1) = 2628.0 * pow(j+1,-1.821); 
      				clade->delete_lengthDistribution.at(j+1) = 2628.0 * pow(j+1,-1.821); 
   					sum += clade->insert_lengthDistribution.at(j+1);
				}
				for(int j=0; j<clade->maxIndel; j++) {
					clade->insert_lengthDistribution.at(j+1) /= sum;
					clade->delete_lengthDistribution.at(j+1) /= sum;
				}
			}
		} else {
			cerr << "Missing }." << endl;
			exit(EXIT_FAILURE);
		}
	}
}

void inTree::parsePound(string& input, inClade *clade, bool local, seqGenOptions *options) 
{
	size_t found, found2;

	found = input.find('#');
	if (found != string::npos) {
		found2 = input.find('#', found + 1);
		if (found2 != string::npos) {
			string opt_str = input.substr(found + 1, found2 - found - 1);
			input.erase(found, found2 - found + 1);

			list<string> tokins = split(opt_str, ",");
			for (list<string>::iterator tok = tokins.begin(); tok != tokins.end(); tok++) {
				trim(*tok);
				if ((*tok)[0] == 'a') {		// Gamma rates.
					clade->rateHetero = GammaRates;
					clade->gammaShape = atof((*tok).substr(1).c_str());
					clade->set_in_clade[__rateHetero__] = true;
					if (clade->gammaShape <= 0.0) {
						cerr << "Bad Gamma shape in lineages." << endl;
						exit(EXIT_FAILURE);
					}
				} else if ((*tok)[0] == 'b') {
					if (!local) {
						double scale = atof((*tok).substr(1).c_str());
						//////////
						/// Need to set both bool scaling factors, since this subsequence-specific
						/// parameter setting will override the global setting.
						//////////
						if (scale < 0) {
							scaleTree = true;
							scaleBranches = false;
						} else {
							scaleBranches = true;
							scaleTree = false;
						}
						branchScale = scale;
					} else {
						cerr << "Cannot scale branches in lineages." << (*tok)[0] << (*tok).substr(1).c_str() << "." << endl;
						exit(EXIT_FAILURE);
					}
				} else if ((*tok)[0] == 'c') {
					string codonrates;
					clade->rateHetero = CodonRates;
					clade->numCats = 3;
					clade->set_in_clade[__rateHetero__] = true;
					codonrates = (*tok).substr(1);
					list<string> rates = split(codonrates, ",");
					if (rates.size() != 3) {
						cerr << "You must specify rates for the 3 categories of codon positions in lineages." << endl;
						exit(EXIT_FAILURE);
					}
					int j = 0;
					for (list<string>::iterator it = rates.begin(); it != rates.end(); it++) {
						clade->catRate[j] = atof((*it).c_str());
						if (clade->catRate[j] <= 0) {
							cerr << "Bad Category Rates: " << codonrates << endl;
							exit(EXIT_FAILURE);
						}
						j++;
					}
				} else if ((*tok)[0] == 'g') {
					clade->rateHetero = DiscreteGammaRates;
					clade->numCats = atoi((*tok).substr(1).c_str());
					clade->set_in_clade[__rateHetero__] = true;	
					if (clade->numCats <= 0 || clade->numCats > MAX_RATE_CATS) {
						cerr << "Bad number of gamma categories in lineages." << endl;
					}
				} else if ((*tok)[0] == 'i') {
					proportion_invariable = atof((*tok).substr(1).c_str());
					clade->proportion_invariable = atof((*tok).substr(1).c_str());
					clade->set_in_clade[__proportion_invariable__] = true;
					if (rootSeqType == RANDOM)
						randomInvariableAssignment = true;
					clade->invariableSites = true;
					if (local) clade->set_in_clade[__invariableSites__] = true;
				} else if ((*tok)[0] == 'f') {
					string frequencies = (*tok).substr(1);
					trim(frequencies);
					frequencies = addFilePath(frequencies, options->filePath);
					ifstream freq_input;
					freq_input.open(frequencies.c_str(), fstream::in);
					if (freq_input.good()) {
						getline(freq_input, frequencies);
						list<string> freqs = split(frequencies, ",");
						if(freqs.size() != (size_t)numStates) {
							cerr << "Frequency file " << frequencies << " contains "
								 << freqs.size() << " frequencies. Expected " << numStates << " characters." << endl;
							if (freqs.size() == numStates+1) {
								cerr << "This may be caused by a trailing comma." << endl;
							}
							exit(EXIT_FAILURE);
						}
						double total = 0;
						for(list<string>::iterator it = freqs.begin(); it != freqs.end(); it++) {
							if (atof((*it).c_str()) == 0.0) {
								total += 1.0E-11;
							} else total+=atof((*it).c_str());
						}

						frequencies = "";
						stringstream temp_val_str;
						int i = 0;
						clade->values2Export2Freq.assign(numStates,0);
						clade->set_in_clade[__values2Export2Freq__] = true;
						for(list<string>::iterator it = freqs.begin(); it != freqs.end(); it++, i++) {
							temp_val_str << atof((*it).c_str()) / total;
							clade->values2Export2Freq.at(i) = atof((*it).c_str()) / total;
							frequencies += temp_val_str.str() + (i < (numStates-1) ? "," : "");
							temp_val_str.str("");
						}
					} else {
						cerr << "Cannot open " << frequencies << ": Invalid filename." << endl;
						exit(EXIT_FAILURE);
					}
					freq_input.close();
					if (!isNucModel) CheckInputAAFrequencies(clade->values2Export2Freq);
				} else if ((*tok)[0] == 'm') {
					string str;
					clade->model = -1;
					str = (*tok).substr(1);
					clade->model = FindModel(str);
					clade->set_in_clade[__model__] = true;
					model = FindModel(str);
				} else if ((*tok)[0] == 'n') {
					if (!local) {
						cerr << "Cannot declare a neofunctionalization globally. Try setting parameter in lineages." << endl;
						exit(EXIT_FAILURE);
					}
					if (clade->constraintChange == PSEUDOGENE) {
						cerr << "Cannot declare a partition to be both a pseudogene and a neofunctionalization." << endl;
						exit(EXIT_FAILURE);
					}
					clade->constraintChange = NEOFUNCTIONALIZATION;
					clade->set_in_clade[__constraintChange__] = true;
				} else if ((*tok)[0] == 'p') {
					if (!local) {
						cerr << "Cannot declare a pseudogene globally. Try setting parameter in lineages." << endl;
						exit(EXIT_FAILURE);
					}
					if (clade->constraintChange == NEOFUNCTIONALIZATION) {
						cerr << "Cannot declare a partition to be both a pseudogene and a neofunctionalization." << endl;
						exit(EXIT_FAILURE);
					}
					clade->constraintChange = PSEUDOGENE;
					clade->set_in_clade[__constraintChange__] = true;
				} else if ((*tok)[0] == 'r') {
					if (!local) {
						if (clade->rateHetero == CodonRates)
							clade->rateHetero = NoRates;
						clade->set_in_clade[__rateHetero__] = true;
					} else {
						clade->rateHetero = NoRates;
					}
				} else {
					cerr << "Invalid option argument in #, " << (*tok)[0] << endl;
					exit(EXIT_FAILURE);
				}
			}
		} else {
			cerr << "Missing closing #." << endl;
			exit(EXIT_FAILURE);
		}
	}
}

void inTree::Setup_Tree(seqGenOptions *options) 
{
	size_t found;
	inMotif *thisMotif;
	Reset_inMotif();

	//////////
	/// Set the ancestor of the root to the root itself. Useful for checking to see if the
	/// node under evolution is actually the root node (setting active props, in particular)
	//////////
	my_tree->root->anc = my_tree->root;

	if (rootSeqType == RANDOM) {
		rootseq_by_partition.clear();
		invariable_by_partition.clear();
		partitionLength = getRootSeq_and_Motif(RANDOM, options);
	} else if (rootSeqType == MULTIPLE_ALIGNMENT_ROOT) {
		rootseq_by_partition.clear();
		invariable_by_partition.clear();
		partitionLength = getRootSeq_and_Motif(MULTIPLE_ALIGNMENT_ROOT, options);
		for (vector<string>::iterator it = motifs_by_partition.begin(); it != motifs_by_partition.end(); it++) {
			found = (*it).find_first_not_of("*");
			if (found == string::npos) {
				cerr << "Problem parsing motif in MA input file. Motif consist of all \"*\" characters." << endl;
				exit(EXIT_FAILURE);
			}
			thisMotif = FindMotif((*it).at(found));
			if (thisMotif == NULL) {
				cerr << "Invalid motif marker specified in MA file \"" << filename << "\"" << endl;
				exit(EXIT_FAILURE);
			}
			thisMotif->sitemap = (*it);
		}
	} else {
		rootseq_by_partition.clear();
		invariable_by_partition.clear();
		partitionLength = getRootSeq_and_Motif(SINGLE_ROOT, options);
	}

	setupTreeStructs(partitionLength, options);
	my_tree->root_numSites = partitionLength;
	ConvertRootSequence(my_tree);

	if(my_tree->treeEnv.front()->invariableSites && !randomInvariableAssignment) {
		for(int i = 0; i < partitionLength; i++) {
			my_tree->root->seq_evo.at(i).setInvariableState(atoi(invariable_by_partition.substr(i,1).c_str()));
		}
	}

}

void inTree::Reset_inMotif()
{
	// Because I change the inMotif sitemap with variable sites, need to reset on each run.
	inMotif *thisMotif;
	size_t found;

	for (vector<string>::iterator it = input_MA_motifs.begin(); it != input_MA_motifs.end(); it++) {
		found = (*it).find_first_not_of("*");
		if (found == string::npos) {
			cerr << "Somehow an invalid motif got to this point..." << endl;
			exit(EXIT_FAILURE);
		}
		thisMotif = FindMotif((*it).at(found));
		thisMotif->sitemap = (*it);
	}
}

void inTree::setupTreeStructs(int partitionLength, seqGenOptions *options) 
{
	//////////
	/// Set categories of evolvingSequence sites (if GammaRates or DiscreteGammaRates... perhaps CodonRates, too)
	//////////
	SetCategories(my_tree->root,partitionLength,options);	// Initial categories for root.

	//////////
	/// This will be necessary no matter what, pretty much.
	//////////
	Create_Global_Arrays(my_tree, partitionLength);
}

void inTree::ConvertRootSequence(TTree *tree) 
{
	string stateCharacters;
	if(isNucModel) stateCharacters = "ACGT";
	else stateCharacters = "ARNDCQEGHILKMFPSTWYV";

	//////////
	/// Put the root sequence into the evolving sequence structure.
	//////////
	for(int i = 0; i < partitionLength; i++) {
		if(rootseq_by_partition.at(i) == 'U' && isNucModel) rootseq_by_partition.at(i) = 'T';
		if(stateCharacters.find(rootseq_by_partition.at(i)) > (size_t)numStates) {
			cerr << "Invalid character '" << rootseq_by_partition.at(i) << "' at position " << i 
				 << " in input file:" << endl << rootseq_by_partition << endl;;
			exit(EXIT_FAILURE);
		}
		tree->root->evolvingSequence->evolutionaryAttributes.at(i).setState(stateCharacters.find(rootseq_by_partition.at(i)));
	}
}

void inTree::Read_MA(string ma_filename, seqGenOptions *options)
{
	list<inMotif*> motifs;
	ifstream ma_in;
	string readin, readin_seq;
	size_t found, found2;
	inMotif *thisMotif;
	input_MA.clear();
	input_MA_motifs.clear();
	size_t ma_length;

	motifs.clear();
	motifs_by_partition.clear();
	ma_in.open(ma_filename.c_str());
	if(ma_in.is_open()) {
		getline(ma_in, readin);
		if (rootSeqType == SINGLE_ROOT) getline(ma_in, readin);
		invariable = trim(readin);
		ma_length = invariable.size();
		while(ma_in.good()) {
			getline(ma_in, readin);
			if (!readin.empty()) {
				if (rootSeqType == SINGLE_ROOT) {
					readin_seq = trim(readin);
				} else {
					found = readin.find_first_of(" \t");
					found2 = readin.find_first_not_of(" \t", found);
					readin_seq = readin.substr(found2,readin.npos);
				}
				if((int)readin_seq.size() != ma_length) {
					cerr << "Sequences in multiple alignment file: " << ma_filename << " are of unequal size." << endl;
					exit(EXIT_FAILURE);
				}
				found = readin_seq.find("*");
				string check = "";
				check = readin_seq.substr(0,1);
				check.push_back('.');	// looking for a star. any char or a . is fine.
				found2 = readin_seq.find_first_not_of(check);
				if (found != string::npos || found2 == string::npos) {
					found = readin_seq.find_first_not_of("*.");
					if (found != string::npos) {
						thisMotif = FindMotif(readin_seq.at(found));
						if (thisMotif != NULL) {
							thisMotif->sitemap = readin_seq;
							input_MA_motifs.push_back(readin_seq);
							motifs_by_partition.push_back(readin_seq);		// Building vector to keep MA motifs
						} else {
							stringstream warning;
							warning << "Input alignment \"" << ma_filename << "\" has undefined motif placed: \'" << readin_seq.at(found) << "\'.";
							options->SpoolWarnings(warning.str());
						}
					} else {
						cerr << "Error of unknown cause in list<inMotif*> inTree::Read_MA(string ma_filename)." << endl;
						exit(EXIT_FAILURE);
					}
				} else {
					input_MA.push_back(readin_seq);
				}
			}
		}
	} else {
		cerr << "Could not open input root sequence file: " << ma_filename << endl;
		exit(EXIT_FAILURE);
	}

}

int inTree::getRootSeq_and_Motif(int type, seqGenOptions *options) 
{
	int from_pos, to_pos;
	size_t num_seqs_in_MA = 0;
	size_t num_sites;
	size_t found;
	//////////
	/// no_motifs: flag passed to constructRootMotifSites. Set to true only if a sequence
	/// is RANDOM && -1 flag is not set (proportion of random sequence to set motifs.
	//////////
	bool no_motifs = false;

	if (type == MULTIPLE_ALIGNMENT_ROOT) {
		num_seqs_in_MA = input_MA.size();
		// Initialize sitemap, to be built later.
		for (list<inMotif*>::iterator it = motif_specs.begin(); it != motif_specs.end(); it++) 
			(*it)->sitemap.clear();

		// Parse range
		if(!ma_range.empty()) {
			list<string> splitrange = split(ma_range, ":");		
			from_pos = atoi((*splitrange.begin()).c_str()) - 1;
			splitrange.pop_front();
			to_pos = atoi((*splitrange.begin()).c_str()) - 1;
			if(to_pos > invariable.size()) {
				cerr << "Cannot read in multiple alignment: Position " << to_pos << " exceeds alignment boundary." << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			from_pos = 0;
			to_pos = invariable.size();
		}

		// Set the number of sequences to pick
		vector<string> data_MA;
		if(ma_numSeqsToUse == -1) {
			ma_numSeqsToUse = num_seqs_in_MA;
			for(int i=0; i<ma_numSeqsToUse; i++) 
				data_MA.push_back(input_MA.at(i));
		} else {
			int random_number;
			for(int i=0; i<ma_numSeqsToUse; i++) {
				random_number = rand() % num_seqs_in_MA;
				data_MA.push_back(input_MA.at(random_number));
			}
		}

		string stateChar;
		if(isNucModel) stateChar = "ACGTYRUN-";
		else stateChar = "ARNDCQEGHILKMFPSTWYVBZJX-";
		vector<int> count (stateChar.size());
		int running_total, k, random_number;
		int gap = stateChar.find_first_of('-');
		size_t elem_pos = 0;
		for(int i = 0; i < to_pos - from_pos; i++) {
			count.assign(stateChar.size(),0);
			for(int j = 0; j < ma_numSeqsToUse; j++) {
				elem_pos = stateChar.find_first_of(data_MA.at(j).at(i));
				if(elem_pos > stateChar.size()) {
					cerr << "Invalid character '" << data_MA.at(j).at(i) << "' in input multiple alignment"
						 << " file " << filename << endl;
					exit(EXIT_FAILURE);
				}
				if(isNucModel) {
					if (stateChar.at(elem_pos) == 'Y') 
						count.at(stateChar.find('C'))+=50, count.at(stateChar.find('T'))+=50;
					if (stateChar.at(elem_pos) == 'R') 
						count.at(stateChar.find('A'))+=50, count.at(stateChar.find('G'))+=50;
					else if (stateChar.at(elem_pos) == 'U')
						count.at(stateChar.find('T'))+=100;
					else if (stateChar.at(elem_pos) == 'N') 
						count.at(0)+=25,count.at(1)+=25,count.at(2)+=25,count.at(3)+=25;
					else count.at(elem_pos)+=100;
				} else {
					if (stateChar.at(elem_pos) == 'B')
						count.at(stateChar.find('N'))+=50, count.at(stateChar.find('D'))+=50;
					else if (stateChar.at(elem_pos) == 'Z')
						count.at(stateChar.find('E'))+=50, count.at(stateChar.find('Q'))+=50;
					else if (stateChar.at(elem_pos) == 'J')
						count.at(stateChar.find('L'))+=50, count.at(stateChar.find('I'))+=50;
					else if (stateChar.at(elem_pos) == 'X') 
						for(size_t x = stateChar.find('A'); x <= stateChar.find('V'); x++)
							count.at(x)+=100/20;
					else
						count.at(elem_pos)+=100;
				}
			}
			int max_elem, max_pos;
			max_elem = max_pos = 0;
			// -1 b/c gap is last.
			for(int j = 0; j < (int)stateChar.size()-1; j++) {
				if(count.at(j) > max_elem) {
					max_elem = count.at(j);
					max_pos = j;
				}
			}

			// Translate into a single sequence.
			if(count.at(gap) <= (ma_numSeqsToUse*100/2.)) {
				invariable_by_partition += invariable.at(i);
				if(invariable.at(i) == '1' || invariable.at(i) == '3' || ma_columnCollapseMethod == 'c') {
					rootseq_by_partition += stateChar.at(max_pos);
				} else {
					random_number = rand() % (ma_numSeqsToUse*100 - count.at(gap));
					running_total = 0;
					for(k=0, running_total=0; random_number > running_total; k++) {
						running_total += count.at(k);
					}
					if(k != 0) k--;
					rootseq_by_partition += stateChar.at(k);
				}

				for (vector<string>::iterator it = motifs_by_partition.begin(); it != motifs_by_partition.end(); it++) {
					found = (*it).find_first_not_of("*");
					for (list<inMotif*>::iterator it2 = motif_specs.begin(); it2 != motif_specs.end(); it2++)
						if ((*it2)->marker == (*it).at(found))
							(*it2)->sitemap += (*it).at(i);
				}
			} else {
				bool ismotifpos = false;
				inMotif *which_motif;
				size_t which_regEx_position, j;
				bitset<20> state_bits, tmp_bits;
				state_bits.reset();		// Set all to 0
				string state_set;
				bool non_null = false;
				for (vector<string>::iterator it = motifs_by_partition.begin(); it != motifs_by_partition.end(); it++) {
					if ((*it).at(i) != '*' && (*it).at(i) != '.') {
						tmp_bits.reset();
						which_regEx_position = 1, j = 1;
						for (list<inMotif*>::iterator it2 = motif_specs.begin(); it2 != motif_specs.end(); it2++) {
							if ( (*it).at(i) == (*it2)->marker ) { 
								which_motif = (*it2);
							}
						}

						ismotifpos = true;
						while ((*it).at(i-j) != '*') {
							if( (*it).at(i-j) != '.' ) which_regEx_position++;
							j++;
						}
						tmp_bits |= which_motif->getRegExValues(which_regEx_position);
						if (tmp_bits.count() == stateCharacters.size()) {
							non_null = true;
						} else {
							state_bits |= tmp_bits;
						} 
					}
				}

				if (state_bits.none() && non_null) state_bits.set();
				if (ismotifpos) {
					invariable_by_partition += invariable.at(i);
					if (state_bits.none()) rootseq_by_partition += "-";
					else {
						state_set.clear();
						for (size_t q = 0; q < state_bits.size(); q++) 
							if (state_bits.test(q)) state_set += stateCharacters.at(q);
						random_number = rand() % state_set.size();
						rootseq_by_partition += state_set.at(random_number);
					}

					for (vector<string>::iterator it = motifs_by_partition.begin(); it != motifs_by_partition.end(); it++) {
						found = (*it).find_first_not_of("*");
						for (list<inMotif*>::iterator it2 = motif_specs.begin(); it2 != motif_specs.end(); it2++)
							if ((*it2)->marker == (*it).at(found))
								(*it2)->sitemap += (*it).at(i);
					}
				}
			}
		}
	} else if (type == SINGLE_ROOT) { 
		rootseq_by_partition = input_MA.at(0);	// Reset rootseq (no longer has '-').
		invariable_by_partition = invariable;
		setMultiStateCharacters();
	} else {		// RANDOM ROOT
		char *cstr = new char [partitionLength + 2];
		RandomSequence(cstr, partitionLength, options);
		cstr[partitionLength] = '\0';
		for(int x = 0; x < partitionLength; x++) rootseq_by_partition += stateCharacters[cstr[x]];
		delete[] cstr;
		if (options->random_sequence_proportion_motif)
			motif_specs = my_tree->DrawMotifs(rootseq_by_partition, options->random_sequence_proportion_motif);		
		else no_motifs = true;
		invariable_by_partition.assign(rootseq_by_partition.size(), '0');
		// All are active at root, so set them to the global environment, treeEnv.front()
		for (list<inMotif*>::iterator it = motif_specs.begin(); it != motif_specs.end(); it++)
			my_tree->treeEnv.front()->my_motifs.push_back(*it);
		Define_Clade_Bipartitions();
	}

	if (no_motifs) {
		my_tree->root->constructUnconstrainedSequence(rootseq_by_partition);
	} else {
		my_tree->constructRootMotifSites(motif_specs, rootseq_by_partition, invariable_by_partition);
	}
	
	return rootseq_by_partition.size();
}

void inTree::setMultiStateCharacters()
{
	int three_state_rand;

	for (string::iterator it = rootseq_by_partition.begin(); it != rootseq_by_partition.end(); it++) {
		if (isNucModel) {
			//////////
			/// Set according to DDBJ: http://www.ddbj.nig.ac.jp/sub/ref1-e.html
			//////////
			// Non-critical issue:
			// Make sure that if this is coding sequence, changes do not make a stop codon.
			// Under real-sequence circumstances, this should not be an issue.
			//////////
			switch ((*it)) {
			case 'U':	// T
				(*it) = 'T';
				break;
			case 'N':	// ACGT
				(*it) = stateCharacters[rand() % numStates];
				break;
			case 'M':	// AC -- amino
				((rand() % 2 == 0) ? ((*it) = 'A') : ((*it) = 'C'));
				break;
			case 'R':	// AG -- purine
				((rand() % 2 == 0) ? ((*it) = 'A') : ((*it) = 'G'));
				break;
			case 'W':	// AT
				((rand() % 2 == 0) ? ((*it) = 'A') : ((*it) = 'T'));
				break;
			case 'S':	// CG
				((rand() % 2 == 0) ? ((*it) = 'C') : ((*it) = 'G'));
				break;
			case 'Y':	// CT -- pyramidine
				((rand() % 2 == 0) ? ((*it) = 'C') : ((*it) = 'T'));
				break;
			case 'K':	// GT -- keto
				((rand() % 2 == 0) ? ((*it) = 'T') : ((*it) = 'G'));
				break;
			case 'V':	// ACG -- not T
				three_state_rand = rand() % 3;
				(
				 (three_state_rand % 3 == 0) 
				  ? ((*it) = 'A') 
				  : (
				     (three_state_rand % 3 == 1) 
					  ? ((*it) = 'C')
					  : ((*it) = 'G')
				    )
				);
				break;
			case 'H':	// ACT -- not g
				three_state_rand = rand() % 3;
				(
				 (three_state_rand % 3 == 0) 
				  ? ((*it) = 'A') 
				  : (
				     (three_state_rand % 3 == 1) 
					  ? ((*it) = 'C')
					  : ((*it) = 'T')
				    )
				);
				break;
			case 'D':	// AGT -- not C
				three_state_rand = rand() % 3;
				(
				 (three_state_rand % 3 == 0) 
				  ? ((*it) = 'A') 
				  : (
				     (three_state_rand % 3 == 1) 
					  ? ((*it) = 'T')
					  : ((*it) = 'G')
				    )
				);
				break;
			case 'B':	// CGT -- not A
				three_state_rand = rand() % 3;
				(
				 (three_state_rand % 3 == 0) 
				  ? ((*it) = 'T') 
				  : (
				     (three_state_rand % 3 == 1) 
					  ? ((*it) = 'C')
					  : ((*it) = 'G')
				    )
				);
				break;
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			default:
				break;
			}
		} else {
			switch ((*it)) {
			case 'B':	// ND
				((rand() % 2 == 0) ? ((*it) = 'N') : ((*it) = 'D'));
				break;
			case 'Z':	// EQ
				((rand() % 2 == 0) ? ((*it) = 'E') : ((*it) = 'Q'));
				break;
			case 'J':	// LI
				((rand() % 2 == 0) ? ((*it) = 'L') : ((*it) = 'I'));
				break;
			case 'X':	// ARNDCQEGHILKMFPSTWYV
				(*it) = stateCharacters[rand() % numStates];
				break;
			case 'U':
				(*it) = 'C';
				break;
			default:
				break;
			}
		}
	}
}

list<string> split(string str, string delim, bool mixedDelim, bool merge) 
{
	list<string> tokins;
	size_t found;
	string possibleTokin;
	int p_found, offset;
	if (mixedDelim) {
		offset = 1;
		found = str.find_first_of(delim);
	} else {
		offset = delim.length();
		found = str.find(delim);
	}
	p_found = -offset;
	
	while (found != string::npos) {
		possibleTokin = str.substr(p_found + offset, found - (p_found + offset));
		if ((merge && possibleTokin != "") || !merge) {
			tokins.push_back(possibleTokin);
		}
		p_found = found;
		if (mixedDelim) {
			found = str.find_first_of(delim, p_found + offset);
		} else {
			found = str.find(delim, p_found + offset);
		}
		
	}
	possibleTokin = str.substr(p_found + offset);
	if ((merge && possibleTokin != "") || !merge) {
		tokins.push_back(possibleTokin);
	}
	return tokins;
}

string& trim(string& str) 
{
	char const * delims = " \t\r\n\v\f";
	size_t found;
	//Leading whitespace
	found = str.find_first_not_of(delims);
	str.erase(0, found);
	
	//Trailing whitespace
	found = str.find_last_not_of(delims);
	return str.erase(found + 1);
}

string& Remove_Whitespace(string& str) 
{
	for(int i = 0; i < str.size(); i++) {
		if(isspace(str.at(i)))
			str.erase(i,1);
	}
	return str;
}
