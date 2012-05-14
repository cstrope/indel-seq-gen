#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <string>
#include <cstring>
#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include "tree.h"
#include "propose_path.h"

using namespace std;

class eventTrack;
class Statistics;
class PathProposal;
class PathProbability;

class SampleStatistics : private Counter<SampleStatistics>
{
public:
	using Counter<SampleStatistics>::howMany;
	PathProbability *P_path;
	
	SampleStatistics();
	void   calculateStatistics( list<eventTrack*> *events );
	void   reportSampleStatistics();
	double FWD() { return forward_probability; }
	double EPC() { return epc_probability; }
	double m_i() { return FWD() - EPC(); }
	void   set_w_i( double val ) { w_i = val; }
	double return_w_i() { return w_i; }
	void   set_numEvents( int val ) { numEvents = val; }

private:
	double forward_probability;
	double epc_probability;
	double w_i;
	int    numEvents;
};

class Statistics : private Counter<Statistics>
{
public:
	using Counter<Statistics>::howMany;
	vector<SampleStatistics> sample_stats;
	double importance_sampling_weight;
	
	Statistics(int howMany) 
		: importance_sampling_weight (0), 
		  M (-numeric_limits<double>::max()), 
		  m_i_sum (0),
		  effective_sample_size (0)
	{ sample_stats.assign(howMany, SampleStatistics()); }
	double returnM() { return M; }
	void setM( double m_i ) { M = m_i; }
	void setWeights();
	void reportStatistics();

	void calc_ESS();
	double ESS();
	class ImportanceSampling : private Counter<ImportanceSampling>
	{
	public:
	
		ImportanceSampling() : u(0), r(0), med(0) { }
		// Accessors
		void set_mu( double mu ) { u = mu; }
		void set_sigma( double sigma ) { r = sigma; }
		void set_median( double median ) { med = median; }
		double return_mu() { return u; }
		double return_sigma() { return r; }
		double return_median() { return med; }
		void mu_sigma( vector<double> X );

		class Range : private Counter<Range>
		{
		public:
			Range() : min(numeric_limits<double>::max()), max(-numeric_limits<double>::max()) { }
			void set_min( double minimum ) { min = minimum; }
			void set_max( double maximum ) { max = maximum; }
			double return_min() { return min; }
			double return_max() { return max; }
		private:
			double min;
			double max;
		}range;

	private:
		double u;
		double r;
		double med;
	
	}importance_sampling;


	void MCMC_run(TTree *tree, int num_steps, string out_file, PathProposal *path);
	class MCMC : private Counter<MCMC>
	{
	public:
		PathProposal *current;
		PathProposal *proposed;

		MCMC() { }	
		double J();  					// J(\theta^t | \theta^*) //
		double Pr_data_given_path();	// Pr(X | \theta) //
		double logr(double P, double J, double Pprime, double Jprime);
		bool   accept_proposal(double P, double J, double Pstar, double Jstar);

		// Routines for MCMC step types.
		double	resample_subpath(TTree *tree, TNode *i_0, TNode *k_0, double t_0, double T, double current_cycle_probability, bool *accepted);
		double	rasmus_resample(TTree *tree, TNode *i_0, TNode *k_0, double t_0, double T, double current_cycle_probability);
	}mcmc;


private:
	double M;	/// MAXIMUM IMPORTANCE SAMPLING REPLICATE VALUE.
	double m_i_sum;	// Sum of all m_i for the SampleStatistics.
	double effective_sample_size;
};

/// Function declarations
PathProposal *setPath(PathProposal *path);
void Leakage();
long print_memory(string message, size_t howMany, int object_size);


#endif