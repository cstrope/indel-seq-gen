
#ifndef _pst_h
#define _pst_h

typedef struct pst_struct {
  /* labels are redundant, defined via the traversal */
  char smooth;
  double *p_c;
  struct pst_struct **sonsp; } pst_node;

typedef pst_node *pst_type;

/* public functions */

pst_type make_empty_pst(int absize);

void print_pst (char *AB, int absize, pst_type T, int *nodesp);
/* we DFS the tree - results in infix notation */

void print_leaves_pst(char *AB, int absize, int max_string, 
		      pst_type T, int *nodesp);

void traverse_pst(char *AB, pst_type nodep, char *S, 
		  pst_node **leafpp, int *lrip);
/* output: (*leafpp)= pointer to leaf of traverse
            (*lrip) = index of last read letter of S
	            (if 0 all read, if lenS none read)  */

void add_pst_node (char *AB, int absize, pst_type T, char *S);
/* builds all nodes on path to S if needed */

void add_probs_pst(char *AB, int absize, int max_string, 
		  int m, char **ds_ptr, pst_type T);
/* allocates and sets all gamma_S(c) for all nodes where p-c==NULL*/

void print_probs_pst(char *AB, int absize, int max_string, pst_type T);
/* nice print of all gamma_S(c) */

void print_pst_stat (const char *AB, const int absize, const int max_string, \
		     pst_type T, const int tsize, const int m, char **ds_ptr, \
		     int *nodesp, int *leavesp, FILE* fp);
/* prints main stat features to FILE inc table of most freq/long leaves */

void sum_pst_stat (const char *AB, const int absize, const int max_string, \
		   pst_type T, const int m, char **ds_ptr, int *nodesp, \
		   int *leavesp, int *llp, int *lfp, int *flp, int *ffp);
/* silent major stat params gathering inc single most freq/long leaf */

void smoothen_pst (int absize, double min_prob, pst_type T);
/* p(x) -> p(x) * (1 - absize * min_prob) + min_prob */
/* smoothes to min value min_prob. only those whose smooth==0 */

void smoothP0_pst(const char *AB, const int absize, const int max_string, 
		  const int m, char **ds_ptr, pst_type T, const double b, 
		  const double x0, const double p0[]);
/* smoothes in relation to a given prior - defunct */

void pseudo_counts_smooth(const char *AB, const int absize, 
           const int max_string, const int m, char **ds_ptr, pst_type T, 
	   const double mue, const double q[][23], double Q[]);
/* the new smoother of alg2. we init Q[] within it */

double log10like_on_pst(char *AB, pst_type T, char *S);
/* returns cumulative result */

double log10like_detail(char *AB, pst_type T, char *S, FILE *fp);
/* saves persymbol log10like to fp and returns cumulative result */

void predict_entry (char *AB, pst_type T, char *S, const int lwrite, \
		    FILE *lfp, const int swrite, FILE *sfp, char *D);
/* for predict_fasta_set in file_handle.c */
/* saves long/short format predictions to designated files */

void emit_string (FILE *fp, pst_type T, char *AB, int absize, 
		  int string_len);
/* for emit.c */
/* emits a string to designated file */

void extend_pst(char *AB, int absize, int max_string,
		      pst_type T, int *nodesp, int *new_nodesp);
/* for 2pfa.c */
/* adds dummy nodes with prediction matching their father */
/* so that a PFA can be derived from the PST */
/* see RST-96 appendix B, p.26 */

void print_pfa_form(char *AB, int absize, int max_string, pst_type T,
	       FILE *lfp, FILE *sfp);
/* for 2pfa.c */
/* prints a pfa from the extended_pst */


#endif /*_pst_h */

