#ifndef _pattern_match_h
#define _pattern_match_h


void count_xSx (const char *S, const int m, char **ds_ptr, 
		int *chi_Sp, int *chi_xSp, int *chi_Sxp);
/* can handle a NULL string */


int count_S   (const char *S, const int m, char **ds_ptr);
int count_Sx  (const char *S, const int m, char **ds_ptr);
double *proball_Sc (const char *AB, int absize, 
                  const int m, char **ds_ptr, char *S);
/* allocates the pst nodes gamma_S(c); conceptually belongs in pst.h*/
/* note: if S==NUL we use S[-1] */

void allchi_Sc (const char *AB, int absize, const int m, char **ds_ptr, 
		char *S, int chi_Sc[]);
     /* chi_Sc should be preallocated, results are strored in it */

int count_diff_S (const char *S, const int m, char **ds_ptr);


#endif /*_pattern_match_h */

