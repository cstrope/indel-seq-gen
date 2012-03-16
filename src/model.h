///  Header file for model.cpp ///
#ifndef _MODEL_H_
#define _MODEL_H_

#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <cstring>
#include <math.h>
#include "eigen.h"
#include "evolve.h"
#include "inTree.h"
#include "gamma.h" 
#include "tree.h"

#define NUM_AA 20
#define NUM_AA_REL_RATES 190
#define NUM_NUC 4
#define NUM_NUC_REL_RATES 6

class inClade;
class seqGenOptions;

extern std::string stateCharacters;
extern unsigned int numStates_squared;
extern unsigned int numStates_cubed;

enum { 
	  NONE=-1, 
	  JC69,
	  K80,
	  K81,
	  K81ne,
	  F81, 
	  F84, 
	  HKY, 
	  T92,
	  TN93,
	  TN93eq,
	  TIM,
	  TIMeq,
	  TVM,
	  TVMeq,
	  SYM,
	  GTR, 
	  JTT, 
	  WAG, 
	  PAM, 
	  BLOSUM, 
	  MTREV, 
	  CPREV, 
	  GENERAL, 
	  numModels 
	 };
enum { 
	  NUC_NONE = -1, 
	  NUC_JC69, 
	  NUC_K80,
	  NUC_K81,
	  NUC_K81ne,
	  NUC_F81,
	  NUC_F84, 
	  NUC_HKY, 
	  NUC_T92, 
	  NUC_TN93,
	  NUC_TN93eq,
	  NUC_TIM,
	  NUC_TIMeq,
	  NUC_TVM,
	  NUC_TVMeq,
	  NUC_SYM,
	  NUC_GTR, 
	  AA_JTT, 
	  AA_WAG, 
	  AA_DAYHOFF78, 
	  AA_BLOSUM62, 
	  AA_MTREV24, 
	  AA_CPREV45, 
	  AA_GENERAL, 
	 };
enum { A, C, G, T };
enum { ala, arg, asn, asp, cys, gln, glu, gly, his, ileu, leu, lys, met, phe, pro, ser, thr, trp, tyr, val};
enum { AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, AGA, AGC, AGG, AGT, ATA, ATC, ATG, ATT,
	   CAA, CAC, CAG, CAT, CCA, CCC, CCG, CCT, CGA, CGC, CGG, CGT, CTA, CTC, CTG, CTT,
	   GAA, GAC, GAG, GAT, GCA, GCC, GCG, GCT, GGA, GGC, GGG, GGT, GTA, GTC, GTG, GTT,
	   TAA, TAC, TAG, TAT, TCA, TCC, TCG, TCT, TGA, TGC, TGG, TGT, TTA, TTC, TTG, TTT
	 };
	
extern const std::string nucleotides;
extern const std::string aminoAcids;
extern const std::string codons;

class RateMatrix : private Counter<RateMatrix>
{
public:
	vector<vector<double> > Pij;/// Transition probability matrix for time T, for each of NUM_CATS
	vector<double> Qij;			/// Instantaneous rate matrix, used to calculate delta_t transitions.
	vector<double> kappa;		/// Nucl. subst. matrix ts/tv values
	vector<double> pi;			/// Nucleotide frequencies
	vector<double> Sij;			/// Relative rates, i.e., GTR->a,b,c,d,e,f
	vector<double> Cijk;		/// Matrix whose used to calculate Pij. (NUM_NUC)^3 values.
	double 	Root[NUM_AA];		/// Needed for eigen() & GTR, previously global variable.
	short	num_categories;		/// Num discrete gamma categories
	double	alphaGamma;			/// Gamma Categories alpha value.
	vector<double> catRate;		/// Rates for each of the categories.
	vector<double> freqRate;	/// No idea what this does, yet. But, used in creating catRates.
	int		substitution_model;	/// code for the particular matrix under consideration.

	using Counter<RateMatrix>::howMany;

	RateMatrix()
		: substitution_model (-1),
		  num_categories (1),
		  alphaGamma (-1)
	{ 
		Qij.clear();
		Pij.clear();
		kappa.clear();
		pi.clear();
		Sij.clear();
		catRate.clear();
	}

	void 	PrintTransitionMatrices();
	void	printQij();
	void 	printPij();
	void	printCijk();
	void	printRoot();
	void	printCategories();
	void	fullReport();
	void	InitializeSubstitutionVectors();
	void 	setModelConditions(string& matrix, string& rates_in, string& freq_in);
	int		FindModel(string& theModel);
	void	SetModel(inClade *branch);
	void	SetupMatrix(bool transition_probability_independent_sites_setup);
	void 	SetFrequencies(double *inFrequencies, inClade *branch = NULL);
	void 	CheckInputAAFrequencies(vector<double>& inFreq);
	void	CheckFrequencies();
	void 	SetVector(vector<double>& vec, short state, double len);
	void 	SetNucVector(vector<double>& vec, short state, double len);
	void 	SetGTRVector(vector<double>& vec, short state, double len);
	void	setPij(Site site, double BL, int rateHetero);
	void	SetMatrix(vector<double>& matrix, double len);
	void 	SetRelativeRates(double *inRelativeRate);
	void 	SetAAFreqs(int theModel, inClade *branch);

	static double jttRelativeRates[NUM_AA_REL_RATES];
	static double jttFrequencies[NUM_AA];
	static double wagRelativeRates[NUM_AA_REL_RATES];
	static double wagFrequencies[NUM_AA];
	static double dayhoffRelativeRates[NUM_AA_REL_RATES];
	static double dayhoffFrequencies[NUM_AA];
	static double blosumRelativeRates[NUM_AA_REL_RATES];
	static double blosumFrequencies[NUM_AA];
	static double mtrevRelativeRates[NUM_AA_REL_RATES];
	static double mtrevFrequencies[NUM_AA];
	static double cprevRelativeRates[NUM_AA_REL_RATES];
	static double cprevFrequencies[NUM_AA];

};

extern const char *modelNames[numModels];
extern const char *modelTitles[numModels];

extern int model, isNucModel, numStates;

#endif /* _MODEL_H_ */
