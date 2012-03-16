
#ifndef _learn_algs_h
#define _learn_algs_h

#include "pst.h"


pst_type grow_tree1 (char *AB, int absize, int m, char **ts_ptr, int *li, \
		    int *pot_nodesp, pst_type T0);
/* classic algorithm */
/* T0 = NULL to grow a new tree, T0 = old tree, to try & extend it */

pst_type grow_tree2 (char *AB, int absize, int m, char **ts_ptr, int *li, \
		    int *pot_nodesp, pst_type T0);
/* bio algorithm - protein oriented */
/* T0 = NULL to grow a new tree, T0 = old tree, to try & extend it */


#endif /*_learn_algs_h */

