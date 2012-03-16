#ifndef _file_handle_h
#define _file_handle_h

#include "pst.h"

#define STRINGLEN_MAX 7500  /* technical, buffer upper limit */
                            /* mustn't exceed <limits.h> INT_MAX */
                            /* in v1.2 and v2.0 longest string = 5217 */

/* actualy we also assume sum of all string length to fit into an int
   therefore if sizeof(int)=2 expect to fail on long data-set files */


void read_alphabet (int *absizep, char **ABp, const char *filename);
/* revised 11/99 from apfa/c-code with a more robust version */

void read_data_set (char **data_setp, int *mp, char ***ds_ptrp, int **lip,
		    const char *AB, const char *filename);
/* reads all strings. no index. printf's some stat. defunct */

void read_ind_data_set (char **data_setp, int *mp, char ***ds_ptrp, int **lip,
		    const char *AB, const char *filename, const char *indfile);
/* reads strings using an index file. printf's stats. defunct */

void mread_ind_data_set (char **data_setp, int *mp, char ***ds_ptrp, \
			 int **lip, const char *AB, const char *filename, \
			 const char *indfile, int *tez, int *trz, int *trl);
/* silent read of strings using an index. stats are returned to caller */

void read_fasta_set (char **data_setp, int *mp, char ***ds_ptrp, \
		     int **lip, const char *AB, const char *filename, \
		     const char *indfile, int *tez, int *trz, int *trl);
/* reads fasta format files of [comment_line,lines] */
/* if indexed all indices=0 are read, rest ('test') are ignored */
/* otherwise all strings are read */

void predict_fasta_set (const char *long_out, const char *short_out, \
			pst_type T, char *AB, const char *filename, \
			const char *indfile, int *tez, int *trz, int *trl);
/* reads fasta format files of [comment_line,lines] */
/* if indexed all indices=0 are read, rest ('test') are ignored */
/* otherwise all strings are read */
/* the file is loaded to memory single string at a time */ 
/* a string is loaded, predicted, overwritten */

void save_pst (char *outfile, char *AB, int absize, int L_max, pst_type T);
/* clobbers existing if found. Saves absize, AB, L-max and the tree itself */
/* does not save the smooth indicator */
/* saving L-max is useful. When you load it you get a depth bound*/

void load_pst (char *pstfile, char **ABp, int *absizep, int *L_maxp, pst_type *Tp);
/* restores the above mentioned parameters along with the tree */
/* marks smooth=1 for all nodes */

void process_test_set (char *filename, char *AB, pst_type T);
/* reads all strings. no index. printf's results. defunct */

void fileproc_test_set (char *filename, char *resfile, char *AB, pst_type T);
/* silent read of all strings. no index. appends results to file. defunct */

void pred_onall_save (char *filename, char *resfile, char *AB, pst_type T, char detail);
/* silent read of all strings. no index. writes results to file. */

void prediction_stat (char *AB, pst_type T, const int n, const char *sfile, 
       const char *ifile, const char *s0file, const char *i0file, const int s1,
       const int s2, int s3, double *non1p, double *non2p, double *non12p, 
       double *iso1p, double *iso2p, double *iso12p);
/* silent read strings by indices. returns stats */


void pred_stat_save (char *AB, pst_type T, const int n, const char *sfile, 
       const char *ifile, const char *s0file, const char *i0file, int s1,
       int s2, int s3, double *non1p, double *non2p, double *non12p, 
       double *iso1p, double *iso2p, double *iso12p, const char *wfile,
       const char *w0file, const int tsize);
/* silent read strings by indices. saves results for matlab/draw.c
   and returns stats  - PFAM damaged!!! */

#endif /*_file_handle_h */



