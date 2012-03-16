
/* add includes */
/* update Makefile */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pfa.h"
#include "numrec.h"

pfa_type make_empty_pfa(int absize, int nodes)
{
  pfa_type M;
  int i;
  
  M.nodes = nodes;
  
  if ((M.node = (pfa_node *) malloc (nodes*sizeof(pfa_node))) == NULL)
    {
      fprintf (stderr,"failed to allocate pfa array\n");
      exit(1);
    }
  
  for (i=0; i<nodes; i++)
    {
      M.node[i].label = NULL; /* marks a vacant slot */
      M.node[i].p_c = NULL;   /* in 2pfa will point to the tree vectors */
      if ((M.node[i].son_c = (int *) malloc (absize*sizeof(int))) == NULL)
	{
	  fprintf (stderr,"failed to allocate a pfa node sons array\n");
	  exit(1);
	}
    }
  return(M);
}


void save_pfa (char *outfile, char *AB, int absize, int L_max, pfa_type M)
{
  FILE *fp;
  int i,j;

  if ((fp = fopen(outfile, "w")) == NULL)
      {
	perror("failed to open <pfa-file> for write");
	exit(1);
      }

  fprintf(fp, "absize\t%d\nab\t%s\nL_max\t%d\nnodes\t%d\n", absize, AB, L_max, M.nodes);

  for (i=0; i<M.nodes; i++)
    {
      fprintf(fp,"\n[%d] label: <%s> p_stat: %1.8f\n", i, M.node[i].label, M.node[i].p_stat);
      for (j=0; j<absize; j++)
	fprintf(fp," %1.8f <%s> --%c--> <%s> [%d]\n", M.node[i].p_c[j], M.node[i].label, AB[j], 
		M.node[M.node[i].son_c[j]].label, M.node[i].son_c[j]);
    }

#ifdef REMOUT
  for (i=0; i<M.nodes; i++)
    {
      fprintf(fp,"%s %1.8f ", M.node[i].label, M.node[i].p_stat);
      for (j=0; j<absize; j++)
	fprintf(fp,"%1.8f %d ", M.node[i].p_c[j], M.node[i].son_c[j]);
    }
#endif

  fclose(fp);
  return;
}



int pfa_index(pfa_type M, pst_node *Tnode, const char *label)
/* scanning left to right if a node with this label exists - return its index
   otherwise insert it at the first (leftmost) free node and return its index */
{
  int i;

  for(i=0; i<M.nodes; i++)
    {
      if (M.node[i].label == NULL)
	{
	  if ((M.node[i].label = (char *) malloc ((strlen(label)+1)*sizeof(char))) == NULL)
	    {
	      fprintf (stderr,"failed to allocate a pfa node label\n");
	      exit(1);
	    }
	  strcpy(M.node[i].label, label);
	  M.node[i].p_c = Tnode->p_c;
	  return(i);
	}
      if (strcmp(M.node[i].label, label) == 0)
	return(i);
    }
  fprintf(stderr,"c-code bug: impossible state in pfa_index\n");
  exit(1);
}


void epst2pfa_aux(char *AB, int absize, int max_string, pst_type T0, pst_node *T,
		   char *S, int lenS, pfa_type M)
/* T is the current node while T0 is the tree root */
{
  const char *ROOT = "-"; /* root node label (mustn't collide with alphabet!)*/
  int i,j;
  int lri;
  pst_node *nodep;

  j = pfa_index(M, T, (lenS==0 ? ROOT : S));

  for(i=0; i<absize; i++)
    {
      S[lenS] = AB[i];
      traverse_pst(AB, T0, S, &nodep, &lri);
      M.node[j].son_c[i] = pfa_index(M, nodep, (lri==lenS+1 ? ROOT : S+lri));
      S[lenS] = '\0';
    }

  for(i=0; i<absize; i++)
    if ((T->sonsp)[i] != NULL)
      {
	S[-1]=AB[i];
	epst2pfa_aux(AB,absize,max_string,T0,(T->sonsp)[i],S-1,lenS+1,M);
      }
}


void epst2pfa(char *AB, int absize, int max_string, pst_type T, pfa_type M)
{
  char *S;
  
  /* +1 for the terminating '\0' and +1 for the next symbol */
  if ((S = (char *) malloc (max_string+2)) == NULL)
    {
      fprintf (stderr,"failed to allocate epst2pfa print buffer\n");
      exit(1);
    }

  S[max_string+1] = '\0';
  S[max_string] = '\0';
  epst2pfa_aux (AB, absize, max_string, T, T, S+max_string, 0, M);
  free(S);
}


void solve_stat_pfa(pfa_type M, int absize)
/* uses LU decomposition - NUMREC book ch. 2 routines */
{
  double **a, *b, d; /* we match names with NUMREC */
  int *indx;         /* NUMREC desires: a[1.n][1.n] b[1.n] indx[1.n] */
  int i,j;
 
  if ((a = (double **) malloc ((M.nodes+1)*sizeof(double *))) == NULL)
    {
      fprintf (stderr,"failed to allocate transition matrix backbone\n");
      exit(1);
    }

  for (i=1; i<=M.nodes; i++) /* we dont bother with a[0] */
    {
      if ((a[i] = (double *) malloc ((M.nodes+1)*sizeof(double))) == NULL)
	{
	  fprintf (stderr,"failed to allocate transition matrix row\n");
	  exit(1);
	}
      for (j=1; j<=M.nodes; j++) /* we dont care about a[?][0] */
	a[i][j] = 0.0; /* calloc not trusted on double */
    }

  if ((b = (double *) malloc ((M.nodes+1)*sizeof(double))) == NULL)
    {    
      fprintf (stderr,"failed to allocate stationary prob. vector\n");
      exit(1);
    }
  
  if ((indx = (int *) malloc ((M.nodes+1)*sizeof(int))) == NULL)
    {    
      fprintf (stderr,"failed to allocate row permutation vector\n");
      exit(1);
    }
  
  for (j=0; j<M.nodes; j++)
    for (i=0; i<absize; i++)
      a[M.node[j].son_c[i]+1][j+1] = M.node[j].p_c[i];  /* a_ij=pr(j->i) */
                                               /* but the matrix is offset by 1 */

  /* convert transition matrix into a linear eq system */
  for (i=1; i<=M.nodes; i++)
    {
      a[i][i] -= 1.0;    /* ax=x => (a-I)x=0 */
      a[M.nodes][i] = 1.0;      /* the next two lines override the last eq */
      b[i] = 1.0*(i==M.nodes);  /* (which is dependant on the others) to   */
    }                             /* the nomralization rule sum(x_i) = 1     */

#ifdef REMARK_OUT
  for (i=1; i<=M.nodes; i++)
    {
      for (j=1; j<=M.nodes; j++)
	printf("%6.3f ",a[i][j]);
      printf("   x_%d   %6.3f\n", i, b[i]);
    }
#endif

  ludcmp(a, M.nodes, indx, &d);  /* NUMREC: LU decomp a onto itself */
  lubksb(a, M.nodes, indx, b);   /* NUMREC: solve LU system into b  */

  for (i=0;i<M.nodes;i++)
    M.node[i].p_stat = b[i+1];

  free(indx);
  free(b);
  for (i=1; i<=M.nodes; i++)
    free(a[i]);
  free(a);

  return;

}

