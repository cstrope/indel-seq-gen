#ifndef _MOTIF_H_
#define _MOTIF_H_

//////////
/// Definitions
//////////
#define ACTIVE_MOTIF   4
#define DORMANT_MOTIF  5
#define VARIABLE_MOTIF 6
// For insertions. Left is the position to the left, right for right.
#define LEFT 0
#define RIGHT 1
// siteProperties general positions.
#define FIRST_POSITION 0
#define OTHER_POSITION 1
#define LAST_POSITION 2
#define ISG_INT_MAX 100000000
#define MOTIF 0
#define TEMPLATE 1
//////////

#include <string>
#include <vector>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <set>
#include <memory>
#include <assert.h>
#include "evolve.h"
#include "inTree.h"
#include "model.h"
#include "tree.h"

using namespace std;

class siteProperties;
class activeProperties;
class varSite;
class inTree;
class inMotif;
class inClade;
class motifSite;
class Substitution;
class Insertion;
class Deletion;
class siteConstraints;
class Indel;
class TNode;
class RateMatrix;
class Site;
class siteDependencies;
class Likelihood;
class Dependency;

extern int actual_events, virtual_events;

class inMotif : private Counter<inMotif>
{
public:
	using Counter<inMotif>::howMany;
	string name;
	char marker;
	string regex;
	string sitemap;
	vector<bool> bipartition;
	bool sequence_template;

	inMotif()
		: name (""),
		  marker ('*'),
		  regex (""),
		  sitemap (""),
		  sequence_template (false)
	{ bipartition.clear(); }
	inMotif(string& prosite_name, string& prosite_regex, size_t num_tips) 
		: name (prosite_name),
		  marker ('*'),
		  regex (prosite_regex),
		  sitemap (""),
		  sequence_template (false)
	{ bipartition.assign(true, num_tips); }
	inMotif(inMotif *thisMotif, string& map, char mark) 
		: name (thisMotif->name),
		  marker (mark),
		  regex (thisMotif->regex),
		  sitemap (map),
		  sequence_template (false)
	{ bipartition.assign(thisMotif->bipartition.begin(), thisMotif->bipartition.end()); }

	vector<siteProperties*> enumerateRegEx(TNode *node);
	list<siteRegEx*> parseRegEx();
	void removeLastPosition();
	bitset<20> getRegExValues(size_t which_site);
	bool isTemplate();
	void report();
};

class siteRegEx : private Counter<siteRegEx>
{
public:
	using Counter<siteRegEx>::howMany;
	list<string> allowable_characters;
	string site_regex;
	size_t sites_occupied;
	bool N_term_motif, C_term_motif, last_motif_site_optional;

	siteRegEx(string regex, size_t sites, bool N_term, bool C_term, bool last_site_optional)
			 : site_regex(regex), sites_occupied(sites), N_term_motif(N_term),
			   C_term_motif(C_term), last_motif_site_optional(last_site_optional)
	{ allowable_characters.clear(); }
};

class Substitution : private Counter<Substitution>
{
public:
	using Counter<Substitution>::howMany;
	bitset<20> substitution_bitstring;
		
	Substitution() { substitution_bitstring.set(); }
	Substitution(Substitution* site2copy) 
	{ 
		substitution_bitstring.set(); 
		substitution_bitstring &= site2copy->substitution_bitstring; 
	}
	void copy(Substitution* site2copy);
	void setSiteBits(bitset<20> regex_site);
	string report_bitset();
	void setInvariable(char residue);
	bool siteInvariable();
};

class varSite : private Counter<varSite>
{
public:
	using Counter<varSite>::howMany;
	int min;
	int max;
	set<Deletion*> member_set;
	TNode *my_TNode;
	varSite *descendant_equiv;

	///////////
	/// CONSTRUCTORS
	///////////
	varSite(size_t minl, size_t maxl, TNode *my_node) 
		: min(minl), 
		  max(maxl), 
		  my_TNode(my_node), 
		  descendant_equiv(NULL)
	{
		member_set.clear();
		my_TNode->variable_region_list.push_back(this);
	}
	varSite(varSite *varSite2copy, TNode *my_node)
		: min(varSite2copy->min), 
		  max(varSite2copy->max), 
		  my_TNode(my_node), 
		  descendant_equiv(NULL)
	{
		member_set.clear();
		varSite2copy->descendant_equiv = this;
	}
	
	///////////
	/// MEMBER FUNCTIONS
	///////////
	bool insertion(int size = 1);
	bool deletion(int size = 1);
	bool isUnconstrained();
	void add2Member (Deletion *newMember);
	void removeFromMember (Deletion *thisMember);
	void print_members();
	varSite *set_on_site_ptr(varSite *prev_site);
};

class Deletion : private Counter<Deletion>
{
public:
	using Counter<Deletion>::howMany;
	varSite *my_sequence_template_varSite, *my_motif_varSite;

	//////////
	/// CONSTRUCTORS
	//////////
	Deletion() : my_sequence_template_varSite(NULL), my_motif_varSite(NULL) 
	{ }
	Deletion(varSite *in_template_varSite, varSite *in_motif_varSite) 
		: my_sequence_template_varSite(in_template_varSite), my_motif_varSite(in_motif_varSite)
	{ }
	Deletion(Deletion *site2copy, motifSite *anc, bool isTemplate) 
		: my_sequence_template_varSite
			( ( (anc != NULL) 
				? ( (isTemplate) 
				 	? (site2copy->my_sequence_template_varSite->descendant_equiv)
				 	: NULL
				  )
				: ( (isTemplate)
				 	? site2copy->my_sequence_template_varSite
				 	: NULL
				  )
			  ) 
		    ),
		  my_motif_varSite
		    ( ( (anc != NULL) 
				? ( (isTemplate) 
					? NULL
				 	: (site2copy->my_motif_varSite->descendant_equiv)
				  )
				: ( (isTemplate)
					? NULL
				 	: site2copy->my_motif_varSite
				  )
			  ) 
			)
	{ }

	//////////
	/// MEMBER FUNCTIONS
	//////////
	bool deletionAllowed(int size = 1);
	void copy(Deletion *site2copy, bool isTemplate, bool isDescendant, Deletion *anc_copy);
	void copy(Deletion *site2copy, short type, bool isDescendant, Deletion *anc_copy);
	void setMembership(varSite *mOS, bool isTemplate);
};

class Insertion : private Counter<Insertion>
{
public:
	using Counter<Insertion>::howMany;
	// Each one has a LEFT and RIGHT motif.
	varSite *my_sequence_template_varSite_left, *my_sequence_template_varSite_right;
	varSite *my_motif_varSite_left, *my_motif_varSite_right;
	varSite *on_site_motif_varSite, *on_site_sequence_template_varSite;	
		
	Insertion()
		: my_sequence_template_varSite_left (NULL),
		  my_motif_varSite_left (NULL),
		  my_sequence_template_varSite_right (NULL),
		  my_motif_varSite_right (NULL),
		  on_site_sequence_template_varSite (NULL),
		  on_site_motif_varSite (NULL)
	{ }
	Insertion(varSite *in_template_varSite, varSite *in_motif_varSite)
		: my_sequence_template_varSite_left (in_template_varSite),
		  my_sequence_template_varSite_right (in_template_varSite),
		  my_motif_varSite_left (in_motif_varSite),
		  my_motif_varSite_right (in_motif_varSite),
		  on_site_motif_varSite
		  	(
		  	  (
		  	    (in_motif_varSite != NULL)
		  	    ? (
		  	    	(in_motif_varSite->min == 0)
		  	    	? in_motif_varSite
		  	    	: NULL
		  	      )
		  	    : NULL
		  	  )
		  	),
		  on_site_sequence_template_varSite
		  	(
		  	  (
		  	    (in_template_varSite != NULL)
		  	    ? (
		  	    	(in_template_varSite->min == 0)
		  	    	? in_template_varSite
		  	    	: NULL
		  	      )
		  	    : NULL
		  	  )
		  	)
	{ }
	Insertion(Insertion *site2copy, motifSite *anc, bool isTemplate, Insertion *anc_copy); // complex, w/ return stmts.

	bool insertionAllowed(int size = 1);
	void copy(Insertion *site2copy, bool isTemplate, bool isDescendant, Insertion *anc_copy);
	void copy(Insertion *site2copy);
	void L_ins_copy(Deletion *site2copy, Deletion *prev_site, short type);
	void setMembership(varSite *mOS, bool isTemplate, size_t side);
};

class Indel : private Counter<Indel>
{
public:
	using Counter<Indel>::howMany;
	Insertion    *L_ins_, *R_ins_;	
	Deletion     *del;

	Indel();
	Indel(vector<siteProperties*> prev_site);
	Indel(Indel *site2copy, motifSite *anc, bool isTemplate);
	Indel(varSite *in_template_varSite, varSite *in_motif_varSite);

	void copy(Indel *site2copy, motifSite *prev, bool isTemplate, bool last_site, bool isDescendant, Indel *anc_copy);
	void set(Indel *site2copy, Indel *prev_site, short type, bool isDescendant, Indel *anc_indel);
	void createObjects();
	
	void report();
};

class activeProperties : private Counter<activeProperties>
{
public:
	using Counter<activeProperties>::howMany;
	Substitution *subst;
	Indel		 *indel;
	inMotif		 *fromMotif;
	inMotif		 *fromTemplate;

	activeProperties() 
		: subst(NULL),
		  indel(NULL),
		  fromMotif(NULL),
		  fromTemplate(NULL)
	{ }
	void report();
	void setProperties(bool activeProp, motifSite *ms);
	void set_properties(list<siteProperties*>& site_props, list<siteProperties*>& anc_site_props, TNode *node, activeProperties *prev_site, size_t site_number);
	void setFirstMotifPos(TNode *node);
	void setProps(short type, siteProperties *props, activeProperties *prev_site, short relation, Indel *anc_indel);
	void setLastMotifPos();

};

class motifSite : private Counter<motifSite>
{
public:
	using Counter<motifSite>::howMany;
	list<siteProperties*> site_props;
	TNode *my_node;
	activeProperties active_properties;

	motifSite();
	void setActiveProps(vector<bool>& bipartition_to_test, motifSite *site, motifSite *anc, inClade *clade, bool last_site);
	void Site_deleteMerge(motifSite *prev, motifSite *next, bool edge = false);
	void deleteSiteObjects(bool edge);
	void setNode(TNode *node);
	void copy(motifSite* site2copy);
	void Print_Nulls();
};

class Likelihood : private Counter<Branch>
{
public:
	using Counter<Branch>::howMany;
	vector<double> Li_xi_;			// Vector holding the likelihoods of each character.
	bool touched;

	Likelihood ()
		: touched (false)
	{ Li_xi_.assign(numStates, 0); }
	//////////
	/// Why 2 functions? 
	//////////
	// State likelihoods are calculated once only, and are active throughout the simulation, so they
	// set the Site likelihood variable, thus requiring a void return value:
	void calculateStateLikelihood(
								  TNode *node, 
								  int category, 
								  int sequence_position,
								  double branch_length_scalar
								 );
	//
	// whereas the likelihood for each category must be calculated, which will be used to infer the
	// category that the Site belongs to, and afterwards, will be unnecessary. Thus, it will use
	// only a temporary variable, so the conditional likelihood variable is returned.
	vector<double> calculateCatLikelihood(
										  TNode *node, 
										  int category, 
										  int sequence_position,
										  vector<double> b1_Li,
										  vector<double> b2_Li
										 );
	//
	//////////
};

class Site : private Counter<Site>
{
public:
	vector<double> site_rate_away;
	vector<double> forward_rate_away;	/// Qij^D
	vector<double> e_QijDt;				/// Qij for Pij^\tilde{D}
	vector<double> Pjk0;				/// Numerator
	double Pik0;						/// Denominator
	double Pr_selected_for_substitution;
	int numSubstAway;
	list<siteDependencies*> interactions;
	Likelihood L_i;
	bool calcTauIJ;

protected:
	short 	state;
	short 	value;
	short 	phenotype;
	short 	invariable_state;
	short 	gamma_category;
	double 	gammaRates;
	double 	rij;
	int		lookup_table_sequence_index;
	int		lookup_table_environment_index;
	
public:
	using Counter<Site>::howMany;
	motifSite motif;

	Site() 
		: state(-1), 
		  value(-1), 
		  phenotype(-1), 
		  invariable_state(0), 
		  gamma_category(0), 
		  gammaRates(1),
		  Pik0(0),
		  calcTauIJ(true),
		  Pr_selected_for_substitution(1),
		  lookup_table_sequence_index(-1),
		  lookup_table_environment_index(-1),
		  numSubstAway(0)
	{ }
	void setrij(double r) { rij = r; }
	void  	set_lookup_table_sequence_index(int value) { lookup_table_sequence_index = value; }
	void  	set_lookup_table_environment_index(int value) { lookup_table_environment_index = value; }
	double 	returnrij() { return rij; }
	int  	return_lookup_table_sequence_index() { return lookup_table_sequence_index; }
	int  	return_lookup_table_environment_index() { return lookup_table_environment_index; }
	void setState(short newState) { state = newState; }
	void setCategory(short newCat) { gamma_category = newCat; }
	void setInvariableState(short newState) { invariable_state = newState; }
	void copyInvariableState(Site copy_site) { invariable_state = copy_site.invariable_state; }
	void setGamma(double gamma_rate) { gammaRates = gamma_rate; }
	void InitializeMotif(motifSite *site_type, varSite *template_varSite, varSite *motif_varSite);
	void push_siteProps(siteProperties *site_prop) { motif.site_props.push_back(site_prop); }
	void setSiteRateAway(double max_site_width, double gamma, RateMatrix *rates);
	void printSiteRateAway();
	bool doSubstitution(double value, int step_type, TNode *node, int *codon_position, string& event);
	bool isSet() { return ( (state >= 0 && state < numStates) ? true : false ); }
	bool isDefinedByUser();
	double setSiteRateAway(vector<double>& Qij, RateMatrix *rates);
	double returnGamma() { return gammaRates; }
	short returnCategory() { return gamma_category; }
	short returnState() { return state; }
	short returnInvariableState() { return invariable_state; }
	long double return_epc_ij()
	{ 
		return site_rate_away.at(state) - ( (state == 0) ? 0 : site_rate_away.at(state-1) );
	}
	double forward_rate_away_from_site(Branch *branch);
};

class siteDependencies : private Counter<siteDependencies>
{
public:
	list<Site*> interact_sites;

	siteDependencies() { }

	bool isDependentInteraction();
};

class Sequence : public Site, private Counter<Sequence>
{
public:
	using Counter<Sequence>::howMany;
	vector<Site> 	evolutionaryAttributes;
	TNode 			*my_node;
	double 			Qidot;
	double 			Qidot_k__T__;
	double 			Rij;

	//////////
	/// Constructors
	//////////
	Sequence(
			 Sequence *anc, 
			 TNode *node
			) 
		: evolutionaryAttributes( anc->evolutionaryAttributes ),
		  my_node(node),
		  Qidot(0.0),
		  Qidot_k__T__(0.0),
		  Rij (0.0)
	{ }
	Sequence(
			 TNode *node, 
			 int seqLength
			)
		: evolutionaryAttributes(seqLength, Site()), 
		  my_node(node),
		  Qidot(0.0),
		  Qidot_k__T__(0.0),
		  Rij (0.0)
	{ }

	void init(
			  TNode *des = NULL,
			  string initial_state = "", 
			  inClade *env = NULL, 
			  motifSite *site_type = NULL, 
			  varSite* template_varSite = NULL, 
			  varSite *motif_varSite = NULL
			 );

 	void setRij(double R) { Qidot = Rij = R; }
	double returnRij() { return Rij; }
	void print_sequence();
	void setActiveProps(bool insertion = false);
	void ConvertRootSequence();
	int	 compare_sequence(Sequence *sequence_to_compare);	// Returns the # of differing positions.
	double forward_rate_away_from_sequence(Branch *branch, int event_site, int end_site);
};

class siteProperties : private Counter<siteProperties>
{
public:
	using Counter<siteProperties>::howMany;
	Substitution subst;
	Indel		 indel;
	inMotif		 *fromMotif;
	inMotif		 *fromTemplate;
		
	siteProperties()
		: subst (Substitution()),
		  indel(Indel()), 
		  fromMotif(NULL), 
		  fromTemplate(NULL)
	{ }
	siteProperties(vector<siteProperties*> site, varSite *var, inMotif *infromMotif, bool isTemplate);
	siteProperties(motifSite *site_type, varSite *in_template_varSite, varSite *in_motif_varSite)
		: subst (Substitution()),
		  indel(Indel(in_template_varSite, in_motif_varSite)),
		  fromMotif
		    ( 
		      (in_motif_varSite) 
		  	  ? site_type->active_properties.fromMotif
		  	  : NULL	   
		    ),
		  fromTemplate
		    (
		      (in_template_varSite)
		      ? site_type->active_properties.fromTemplate
		      : NULL
		    )
    { }
	siteProperties(siteProperties* site2copy, motifSite *anc)
		: subst (Substitution(site2copy->subst)),
		  indel
		  	(
		   	  (
		   	    (site2copy->fromTemplate != NULL)
		   	    ? (Indel(&site2copy->indel, anc, true))
		   	    : (
		   	    	(site2copy->fromMotif != NULL)
		   	      	? (Indel(&site2copy->indel, anc, false))
		   	      	: (Indel())
			  	  )
			  )
			),
	  	  fromMotif(site2copy->fromMotif),
	  	  fromTemplate(site2copy->fromTemplate)
	{ }

	void copy(siteProperties* stuff, short which, short type, motifSite *prev, siteProperties *anc, bool last_site);
};

#endif
