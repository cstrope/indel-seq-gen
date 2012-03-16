
#ifndef _pfa_h
#define _pfa_h

#include "pst.h"

typedef struct pfa_node_struct {
  char *label;
  double p_stat;
  double *p_c;
  int *son_c; } pfa_node;

typedef struct pfa_struct {
  int nodes;
  pfa_node *node; } pfa_type;

/* public functions */

pfa_type make_empty_pfa(int absize, int nodes);

void save_pfa (char *outfile, char *AB, int absize, int L_max, pfa_type M);
/* clobbers existing if found */

void load_pfa();

void solve_stat_pfa(pfa_type M, int absize);
/* finds the stationary probabilites and updates the pfa */

void epst2pfa(char *AB, int absize, int max_string, pst_type T, pfa_type M);

#endif /*_pfa_h */

