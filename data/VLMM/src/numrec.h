
#ifndef _numrec_h
#define _numrec_h

/* functions from Numerical Recipes in C 2nd ed., windows diskette ver 2.08 */

/* converted from float to double */

void ludcmp(double **a, int n, int *indx, double *d);
/* ch.2, p.46: perform LU decomposition. */
 
void lubksb(double **a, int n, int *indx, double b[]);
/* ch.2, p.47: solve LU system */


#endif /*_numrec_h */

