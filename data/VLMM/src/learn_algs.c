 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "learn_params.h"  /* learning params come as GLOBALS */
#include "learn_algs.h"
#include "priority_queue.h"
#include "pattern_match.h"

/* local */

double exp_chi_S (char *AB, char *S, int N[])
/* DEFUNCT */
{
  int i,s;
  double e;

  s = strlen(S);
  e = N[s];
  for (i=0; i<s; i++)
    e *= p0[strchr(AB,S[i])-AB];
  return(e);
}


void stock_pqueue1 (int absize, pst_type T0, char *AB, char *S, int m,
		  char **ts_ptr, pQ_type *Qp, int *N, int *pot_nodesp)
{
  int chi_S, chi_xS, chi_Sx, i;
  double p_S;
  
  (*pot_nodesp) += 1;
  S[0] = AB[0]; /* so it won't contain garbage */
  if (strlen(S)>L_max) /* we won't consider the sons of an L deep leave */
    return;
  for (i=0; i<absize; i++)
    {
      S[0] = AB[i];
      if ((T0->sonsp)[i] != NULL)
	stock_pqueue1 (absize, (T0->sonsp)[i], AB, S-1, m, ts_ptr, Qp, N, pot_nodesp);
      else /* a candidate for rejuvination */
	{
	  count_xSx (S, m, ts_ptr, &chi_S, &chi_xS, &chi_Sx);
	  p_S = ((double) chi_S)/N[strlen(S)];      
	  if (p_S >= p_min)
	    {
	      stack_pqueue (Qp, S, p_S, chi_Sx, chi_xS, 
			      count_Sx(S+1,m,ts_ptr));
	    }
	}
    }
  return;
}


void stock_pqueue2 (int absize, pst_type T0, char *AB, char *S, int m,
		  char **ts_ptr, pQ_type *Qp, int *N, int *pot_nodesp)
{
  int chi_S, chi_xS, chi_Sx, i;
  
  (*pot_nodesp) += 1;
  S[0] = AB[0]; /* so it won't contain garbage */
  if (strlen(S)>L_max) /* we won't consider the sons of an L deep leave */
    return;
  for (i=0; i<absize; i++)
    {
      S[0] = AB[i];
      if ((T0->sonsp)[i] != NULL)
	stock_pqueue2 (absize, (T0->sonsp)[i], AB, S-1, m, ts_ptr, Qp, N, pot_nodesp);
      else /* a candidate for rejuvination */
	{
	  count_xSx (S, m, ts_ptr, &chi_S, &chi_xS, &chi_Sx);
	  if (1.0*count_diff_S(S,m,ts_ptr)/m >= shows) 
	    {
	      stack_pqueue (Qp, S, 0.0, chi_Sx, chi_xS, 
			      count_Sx(S+1,m,ts_ptr));
	    }
	}
    }
  return;
}



/* public */


pst_type grow_tree1(char *AB, int absize, int m, char **ts_ptr, int *li, \
		   int *pot_nodesp, pst_type T0)

{

  pst_type T;
  pQ_type Q; /* used as a LIFO stack */
  int i;
  char *S;
  int chi_S, chi_xS, chi_Sx, chi_sufSx;
  int chi_Sc, chi_sufSc;
  int chi_cS, chi_xcS, chi_cSx;
  double p_S, p_cS;
  int insert;
  int lenS;
  double lr;
  int *N;

  if ((N = (int *) malloc ((L_max+1)*sizeof(int))) == NULL)
    {
      perror("failed to allocate N\n");
      exit(1);
    }
  N[0] = 0;  /* unused */
  for (lenS=1; lenS<=L_max; lenS++) {
    N[lenS] = 0;
    for (i=1; i<=m; i++) 
      N[lenS] += (li[i]-(lenS-1))*(li[i]>=lenS);
  }  /* accounts correctly for strings who are shorter than the suffix */


  Q = make_empty_pqueue();
  if ((S = (char *) malloc (L_max+3)) == NULL)
    {  /* one upfront, one at the end, and one for \0 */
      perror("failed to allocate S");
      exit(1);
    }
  S = S+1; /* so we can add one letter upfront */

  if (T0 == NULL)
    {
      T = make_empty_pst(absize);
      (*pot_nodesp) = 1; /* count the root */
      S[1] = '\0';
      for (i=0; i<absize; i++)
	{
	  S[0] = AB[i];
	  count_xSx (S, m, ts_ptr, &chi_S, &chi_xS, &chi_Sx);
	  p_S = ((double) chi_S)/N[1];      
	  if (p_S >= p_min)
	    {
	      stack_pqueue (&Q, S, p_S, chi_Sx, chi_xS, N[1]);
	    }
	}
    }
  else
    {
      (*pot_nodesp) = 0;
      S[L_max+1] = '\0';
      stock_pqueue1 (absize, T0, AB, S+L_max, m, ts_ptr, &Q, N, pot_nodesp);
      T = T0;
    }

  while (not_empty_pqueue(Q))
    {
      dequeue_pqueue(&Q, S, &p_S, &chi_Sx, &chi_xS, &chi_sufSx);      
      (*pot_nodesp) += 1;
      insert = 0;
      lenS = strlen(S);
      S[lenS+1]='\0'; /* in case we want to add a suffix letter */
      for (i=0; i<absize; i++)
	{
	  if (insert==0)
	    {
	      S[lenS]=AB[i]; /* turn S to Sc */
 	      chi_Sc = count_S (S, m, ts_ptr);
	      if (((double) chi_Sc)/chi_Sx >= (1.0+alpha)*gamma_min)
		{
		  chi_sufSc = count_S (S+1, m, ts_ptr);
		  lr = (((double) chi_Sc)/chi_Sx) / 
		       (((double) chi_sufSc)/chi_sufSx);
		  insert = ((lr>=p_ratio) || ((1.0/lr)>=p_ratio));
		}
	    }
	  if (lenS<L_max)
	    {
	      S[lenS]='\0'; /* turn back to S */
	      S[-1]=AB[i];  /* now S-1 is cS */
	      count_xSx (S-1, m, ts_ptr, &chi_cS, &chi_xcS, &chi_cSx);
	      p_cS = ((double) chi_cS)/N[lenS+1];
	      if (p_cS >= p_min)
		{
		  stack_pqueue (&Q, S-1, p_cS, chi_cSx, chi_xcS, chi_Sx); 
		}
	    }
	}
      if (insert==1)
	{
	  S[lenS]='\0'; /* turn back to S */
	  add_pst_node(AB, absize, T, S);  /* adds midway nodes if needed */
	}

    }
  
  return(T);
}




pst_type grow_tree2(char *AB, int absize, int m, char **ts_ptr, int *li, \
		   int *pot_nodesp, pst_type T0)

{

  pst_type T;
  pQ_type Q;  /* used as a LIFO stack */
  int i;
  char *S;
  int chi_S, chi_xS, chi_Sx, chi_sufSx;
  int chi_Sc, chi_sufSc;
  int chi_cS, chi_xcS, chi_cSx;
  double dummy;
  int insert;
  int lenS;
  double lr;
  int *N;

  if ((N = (int *) malloc ((L_max+1)*sizeof(int))) == NULL)
    {
      perror("failed to allocate N\n");
      exit(1);
    }
  N[0] = 0;  /* unused */
  for (lenS=1; lenS<=L_max; lenS++) {
    N[lenS] = 0;
    for (i=1; i<=m; i++) 
      N[lenS] += (li[i]-(lenS-1))*(li[i]>=lenS);
  }  /* accounts correctly for strings who are shorter than the suffix */


  Q = make_empty_pqueue();
  if ((S = (char *) malloc (L_max+3)) == NULL)
    {  /* one upfront, one at the end, and one for \0 */
      perror("failed to allocate S");
      exit(1);
    }
  S = S+1; /* so we can add one letter upfront */

  
  if (T0 == NULL)
    {
      (*pot_nodesp) = 1; /* count the root */
      T = make_empty_pst(absize);
      S[1] = '\0';
      for (i=0; i<absize; i++)
	{
	  S[0] = AB[i];
	  count_xSx (S, m, ts_ptr, &chi_S, &chi_xS, &chi_Sx);
	  if (1.0*count_diff_S(S,m,ts_ptr)/m >= shows) 
	    {
	      stack_pqueue (&Q, S, 0.0, chi_Sx, chi_xS, N[1]);
	    }
	}
    }
  else
    {
      (*pot_nodesp) = 0;
      S[L_max+1] = '\0';
      stock_pqueue2 (absize, T0, AB, S+L_max, m, ts_ptr, &Q, N, pot_nodesp);
      T = T0;
    }

  while (not_empty_pqueue(Q))
    {
      dequeue_pqueue(&Q, S, &dummy, &chi_Sx, &chi_xS, &chi_sufSx);      
      (*pot_nodesp) += 1;
      insert = 0;
      lenS = strlen(S);
      S[lenS+1]='\0'; /* in case we want to add a suffix letter */
      for (i=0; i<absize; i++)
	{
	  if (insert==0)
	    {
	      S[lenS]=AB[i]; /* turn S to Sc */
 	      chi_Sc = count_S (S, m, ts_ptr);
	      if (1.0*chi_Sc/chi_Sx >= gama)
		{
		  chi_sufSc = count_S (S+1, m, ts_ptr);
		  lr = (((double) chi_Sc)/chi_Sx) / 
		       (((double) chi_sufSc)/chi_sufSx);
		  insert = ((lr>=p_ratio) || ((1.0/lr)>=p_ratio));
		}
	    }
	  if (lenS<L_max)
	    {
	      S[lenS]='\0'; /* turn back to S */
	      S[-1]=AB[i];  /* now S-1 is cS */
	      count_xSx (S-1, m, ts_ptr, &chi_cS, &chi_xcS, &chi_cSx);
	      if (1.0*count_diff_S(S-1,m,ts_ptr)/m >= shows) 
		{
		  stack_pqueue (&Q, S-1, 0.0, chi_cSx, chi_xcS, chi_Sx); 
		}
	    }
	}
      if (insert==1)
	{
	  S[lenS]='\0'; /* turn back to S */
	  add_pst_node(AB, absize, T, S);  /* adds midway nodes if needed */
	}
    }
  
  return(T);
}


