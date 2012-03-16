
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pst.h"
#include "pattern_match.h"
#include "priority_queue.h"

/* local */

int ABind (char *AB, char c)
{
  return(strchr(AB,c)-AB);
}

pst_node *make_pst_node (int absize)
{
  pst_node *nodep;
  int i;

  if ((nodep = (pst_node *) malloc (sizeof(pst_node))) == NULL)
    {
      fprintf (stderr,"failed to allocate pst node\n");
      exit(1);
    }
  nodep->p_c = NULL; /* in add_probs_pst we init and fill */
  if ((nodep->sonsp = (pst_node **) malloc (absize*sizeof(pst_node *))) == NULL)
    {
      fprintf (stderr,"failed to allocate pst node sons pointers\n");
      exit(1);
    }
  for (i=0; i<absize; i++)
    (nodep->sonsp)[i] = NULL;

  nodep->smooth = 0;

  return(nodep);
}


pst_node *add_pst_leaf(char *AB,int absize, pst_node *leafp, char c)
/* returns a pointer to the new leaf */

{
  pst_node *newleafp;

  newleafp = make_pst_node(absize);
  (leafp->sonsp)[ABind(AB,c)] = newleafp;
  return(newleafp);
}



void add_probs_aux (char *AB, int absize, int m, char **ds_ptr,
		       pst_type T, char *S)
/* S must be prorerly inited by add_probs_pst
   it points to the string\0 of this level */

{
  int i;

  if (T->p_c == NULL)
    T->p_c = proball_Sc (AB, absize, m, ds_ptr, S); /* func in pattern_match */
  
  for (i=0; i<absize; i++)
    {
      if ((T->sonsp)[i] != NULL)
	{
	  S[-1]=AB[i];
	  add_probs_aux(AB,absize,m,ds_ptr,(T->sonsp)[i],S-1);
	}
    }
  
}


void print_probs_aux (char *AB, int absize, pst_type T, char *S)
/* S must be prorerly inited by print_probs_pst 
   it points to the string\0 of this level */

{
  int i;

  printf("<%s>",S);
  for (i=0; i<absize; i++)
    {
      if (i % 7 == 0)
	printf("\n\t");
      if (T->p_c[i] == 0.0)
	printf("    --    ");
      else
	printf("%c %6.4f  ", AB[i], T->p_c[i]);
    }
  printf("\n\n");
  
  for (i=0; i<absize; i++)
    {
      if ((T->sonsp)[i] != NULL)
	{
	  S[-1]=AB[i];
	  print_probs_aux(AB,absize,(T->sonsp)[i],S-1);
	}
    }
  
}



void print_leaves_aux (char *AB, int absize, pst_type T, char *S, int *nodesp)
/* S must be prorerly inited by print_leaves_pst */
{
  int i;
  int leaf;

  leaf=1;
  for (i=0; i<absize; i++)
    {
      if ((T->sonsp)[i] != NULL)
	{
	  leaf=0;
	  S[-1]=AB[i];
	  print_leaves_aux(AB,absize,(T->sonsp)[i],S-1, nodesp);
	}
    }
  
  if (leaf)
    {
      printf("%s  ",S);
      (*nodesp)++;
    }
}


void print_stat_aux (const char *AB, const int absize, pst_type T, char *S, \
		     pQ_type *Qfp, pQ_type *Qlp, const int m, char **ds_ptr, \
		     int *nodesp, int *leavesp)
{
  int i, leaf, freq, len;

  (*nodesp)++;
  leaf=1;
  for (i=0; i<absize; i++)
    {
      if ((T->sonsp)[i] != NULL)
	{
	  leaf=0;
	  S[-1]=AB[i];
	  print_stat_aux(AB,absize,(T->sonsp)[i],S-1,Qfp,Qlp,m,ds_ptr, \
			 nodesp,leavesp);
	}
    }
  
  if (leaf)
    {
      (*leavesp)++;
      len = strlen(S);
      freq = count_S(S,m,ds_ptr);
      enqueue_pqueue(Qfp,S,1.0*freq+((1e-10)*len),len,freq,0);
      enqueue_pqueue(Qlp,S,1.0*len+((1e-10)*freq),len,freq,0);
    }
  return;
}





void sum_stat_aux (const char *AB, const int absize, pst_type T, char *S, \
		     const int m, char **ds_ptr, int *nodesp, int *leavesp, \
		     int *llp, int *lfp, int *flp, int *ffp)
{
  int i, leaf, freq, len;

  (*nodesp)++;
  leaf=1;
  for (i=0; i<absize; i++)
    {
      if ((T->sonsp)[i] != NULL)
	{
	  leaf=0;
	  S[-1]=AB[i];
	  sum_stat_aux(AB, absize, (T->sonsp)[i], S-1, m, ds_ptr, \
			 nodesp, leavesp, llp, lfp, flp, ffp);
	}
    }
  
  if (leaf)
    {
      (*leavesp)++;
      len = strlen(S);
      freq = count_S(S,m,ds_ptr);
      if ((len > *llp) || ((len == *llp) && (freq>*lfp)))
	{
	  *llp = len;
	  *lfp = freq;
	}
      if ((freq > *ffp) || ((freq == *ffp) && (len>*flp)))
	{
	  *flp = len;
	  *ffp = freq;
	}
    }
  return;
}



/* public */

pst_type make_empty_pst(int absize)
{
  return(make_pst_node(absize));
}


void print_pst (char *AB, int absize, pst_type T, int *nodesp)
/* we DFS the tree */
{
  int i;
  int leaf;

  (*nodesp)++;
  leaf=1;
  for (i=0; i<absize; i++)
    {
      if ((T->sonsp)[i] != NULL)
	{
	  if (leaf)
	    {
	      leaf=0;
	      printf("[");
	    }
	  printf(" %c", AB[i]);
	  print_pst(AB,absize,(T->sonsp)[i], nodesp);
	}
    }
  if (!leaf)
    printf(" ]");
}


void traverse_pst(char *AB, pst_type nodep, char *S, 
		  pst_node **leafpp, int *lrip)

/* output: (*leafpp)= pointer to leaf of traverse
            (*lrip) = index of last read letter of S
	            (if 0 all read, if lenS none read)  */
/* handles empty string correctly.
   reads string right to left, suffix-wise! */
{
  pst_node *nextp;

  (*lrip) = strlen(S);
  if ((*lrip) > 0)
    nextp = (nodep->sonsp)[ABind(AB,S[(*lrip)-1])];
  else
    nextp = NULL;
  while (nextp != NULL) /* the tree continues */
    {
      nodep = nextp;
      (*lrip)--;
      if (*lrip == 0) /* the string ends */
	break;
      nextp = (nodep->sonsp)[ABind(AB,S[(*lrip)-1])];
    }
  (*leafpp) = nodep;
}


void add_pst_node (char *AB, int absize, pst_type T, char *S)
/* builds all nodes on path to S if needed */

{
  int lri; /* index of last read char in S (reading right to left) */
  pst_node *leafp;

  traverse_pst(AB, T, S, &leafp, &lri);

  while (lri>0)
    {
      leafp = add_pst_leaf(AB, absize, leafp, S[lri-1]);
      lri--;
    }
  return;
}

void print_leaves_pst(char *AB, int absize, int max_string, 
		      pst_type T, int *nodesp)

{
  char *S;

  if ((S = (char *) malloc (max_string+1)) == NULL)
    {
      fprintf (stderr,"failed to allocate leaves print buffer\n");
      exit(1);
    }

  S[max_string]='\0';
  print_leaves_aux (AB, absize, T, S+max_string, nodesp);
  free(S);
}


void add_probs_pst(char *AB, int absize, int max_string, 
		  int m, char **ds_ptr, pst_type T)

{
  char *S;

  /* deep down inside, proball_Sc will need the extra char */
  if ((S = (char *) malloc (max_string+2)) == NULL)
    {
      fprintf (stderr,"failed to allocate add pst probs buffer\n");
      exit(1);
    }

  S[max_string+1]='\0';
  add_probs_aux (AB, absize, m, ds_ptr, T, S+max_string+1);
  free(S);
}




void print_probs_pst(char *AB, int absize, int max_string, pst_type T)

{
  char *S;

  if ((S = (char *) malloc (max_string+1)) == NULL)
    {
      fprintf (stderr,"failed to allocate print pst probs buffer\n");
      exit(1);
    }

  S[max_string]='\0';
  print_probs_aux (AB, absize, T, S+max_string);
  free(S);
}


void print_pst_stat (const char *AB, const int absize, const int max_string, \
		     pst_type T, const int tsize, const int m, char **ds_ptr, \
		     int *nodesp, int *leavesp, FILE* fp)
{
  pQ_type Qf;
  pQ_type Ql;
  char *S;
  double dummy1;
  int i, dummy2, len, freq, prevl, prevf;

  if ((S = (char *) malloc (max_string+1)) == NULL)
    {
      fprintf (stderr,"failed to allocate stat print buffer\n");
      exit(1);
    }

  S[max_string]='\0';
  (*nodesp) = 0;
  (*leavesp) = 0;
  Qf = make_empty_pqueue();
  Ql = make_empty_pqueue();
  print_stat_aux(AB,absize,T,S+max_string,&Qf,&Ql,m,ds_ptr,nodesp,leavesp);

  fprintf(fp, "\nthe %d most frequent leaves:", tsize);
  fprintf(fp, "%*sthe %d longest leaves:", max_string-9+6, " ", tsize);
  fprintf(fp, "\n\n");
  fprintf(fp, "      len   freq   %*s", max_string, "suffix");
  fprintf(fp, "            len   freq   %*s\n", max_string, "suffix");
  fprintf(fp, "      ---   ----   %*s", max_string, "------");
  fprintf(fp, "            ---   ----   %*s\n", max_string, "------");

  prevf = 0;
  prevl = 0;
  for (i=1; (i<=tsize) && (not_empty_pqueue(Qf) || not_empty_pqueue(Ql)); i++)
    {
      if (not_empty_pqueue(Qf))
	{
	  dequeue_pqueue(&Qf,S,&dummy1,&len,&freq,&dummy2);
	  if (freq == prevf) 
	    fprintf(fp, "( -) ");
	  else
	    fprintf(fp, "(%2d) ", i);
	  fprintf(fp, " %3d  %5d   %*s", len, freq, max_string, S);
	  fprintf(fp, "  ||  ");
	  prevf = freq;
	}

      if (not_empty_pqueue(Ql))
	{
	  dequeue_pqueue(&Ql,S,&dummy1,&len,&freq,&dummy2);
	  if (len == prevl) 
	    fprintf(fp, "( -) ");
	  else
	    fprintf(fp, "(%2d) ", i);
	  fprintf(fp, " %3d  %5d   %*s\n", len, freq, max_string, S);
	  prevl = len;
	}
    }
  
  fprintf(fp, "\n");
  free_pqueue(Qf);
  free_pqueue(Ql);
  free(S);
  return;

}



void sum_pst_stat (const char *AB, const int absize, const int max_string, \
		   pst_type T, const int m, char **ds_ptr, int *nodesp, \
		   int *leavesp, int *llp, int *lfp, int *flp, int *ffp)
{
  char *S;

  if ((S = (char *) malloc (max_string+1)) == NULL)
    {
      fprintf (stderr,"failed to allocate stat print buffer\n");
      exit(1);
    }

  S[max_string]='\0';

  *nodesp = 0;
  *leavesp = 0;
  *llp = 0;
  *lfp = 0;
  *flp = 0;
  *ffp = 0;

  sum_stat_aux(AB,absize,T,S+max_string,m,ds_ptr,nodesp,leavesp, \
	       llp,lfp,flp,ffp);

  free(S);
  return;

}



void smoothen_pst (int absize, double min_prob, pst_type T)
/* p(x) -> p(x) * (1 - absize * min_prob) + min_prob */
{
  int i;

  if (T->smooth == 0)
    {
      for(i=0; i<absize; i++)
	T->p_c[i] = T->p_c[i] * (1 - absize * min_prob) + min_prob;
      T->smooth = 1;
    }
  
  for (i=0; i<absize; i++)
    {
      if ((T->sonsp)[i] != NULL)
	  smoothen_pst(absize, min_prob, (T->sonsp)[i]);
    }
}


double lambda_fun (double b, double x0, double x)

{
  return((tanh(b*(x-x0))-tanh(-b*x0))/2.0);
}


void smoothP0_aux (const char *AB, const int absize, pst_type T, char *S, 
		   const int m, char **ds_ptr, const double b, 
		   const double x0, const double p0[])

{
  int i;
  double lambda;

  lambda = lambda_fun(b, x0, count_Sx(S,m,ds_ptr));
  for(i=0; i<absize; i++)
    {
      T->p_c[i] = lambda * T->p_c[i] + (1-lambda) * p0[i];
    }
  for (i=0; i<absize; i++)
    {
      if ((T->sonsp)[i] != NULL)
	{
	  S[-1] = AB[i];
	  smoothP0_aux(AB, absize, (T->sonsp)[i], S-1, m, ds_ptr, b, x0 ,p0);
	}
    }
}


void smoothP0_pst(const char *AB, const int absize, const int max_string, 
		  const int m, char **ds_ptr, pst_type T, const double b, 
		  const double x0, const double p0[])
{
  char *S;

  if ((S = (char *) malloc (max_string+1)) == NULL)
    {
      fprintf (stderr,"failed to allocate S buffer for smoothing\n");
      exit(1);
    }

  S[max_string]='\0';
  smoothP0_aux (AB, absize, T, S+max_string, m, ds_ptr, b, x0, p0);
  free(S);
}



void pcounts_smooth_aux(const char *AB, const int absize, pst_type T, 
           char *S, const int m, char **ds_ptr, const double mue, 
	   const double q[][23], const double Q[], int chi_Sc[])

{


  int i,c, chi_Sx;
  double B_S, b_c;

  allchi_Sc(AB, absize, m, ds_ptr, S, chi_Sc);

  if (T->smooth == 0) 
    {
      chi_Sx = 0;
      B_S = 0.0;
      for (i=0; i<absize; i++)
	{
	  chi_Sx += chi_Sc[i];
	  B_S += (chi_Sc[i] > 0);
	}
      B_S *= mue;
      
      for (c=0; c<absize; c++)
	{
	  b_c = 0.0;
	  for (i=0; i<absize; i++)
	    b_c += (chi_Sc[i]*q[i][c])/(chi_Sx*Q[i]);
	  b_c *= B_S;
	  T->p_c[c] = (chi_Sc[c]+b_c)/(chi_Sx+B_S);
	}
      T->smooth = 1;
    }


  for (i=0; i<absize; i++)
    {
      if ((T->sonsp)[i] != NULL)
	{
	  S[-1] = AB[i];
	  pcounts_smooth_aux(AB, absize, (T->sonsp)[i], S-1, m, ds_ptr, 
			     mue, q, Q, chi_Sc );
	}
    }
  
  return;
}



void pseudo_counts_smooth(const char *AB, const int absize, 
           const int max_string, const int m, char **ds_ptr, pst_type T, 
	   const double mue, const double q[][23], double Q[])

{
  int i,j, *chi_Sc;
  char *S;

  if ((S = (char *) malloc (max_string+1)) == NULL)
    {
      fprintf (stderr,"failed to allocate S buffer for smoothing\n");
      exit(1);
    }

  S[max_string]='\0';

  if ((chi_Sc = (int *) malloc (absize*sizeof(int))) == NULL)
    {
      fprintf (stderr,"failed to allocate letter count buffer for smoothing\n");
      exit(1);
    }

  for (i=0; i<absize; i++)
    {
      Q[i] = 0.0;
      for (j=0; j<absize; j++)
	Q[i] += q[i][j];
    }

  pcounts_smooth_aux (AB, absize, T, S+max_string, m, ds_ptr, mue, q, Q, chi_Sc);

  free(S);
  free(chi_Sc);

  return;
}



double log10like_on_pst(char *AB, pst_type T, char *S)

{
  char c;
  int i;
  int lri;
  int lenS;
  pst_node *leafp;
  double loglike;

  loglike = 0.0;
  lenS = strlen(S);
  for (i=0; i<lenS; i++)
    {
      c = S[i];
      S[i] = '\0';
      traverse_pst(AB, T, S, &leafp, &lri);
      loglike += log10(leafp->p_c[ABind(AB,c)]);
      S[i] = c;
    }
  return(loglike);
}


double log10like_detail(char *AB, pst_type T, char *S, FILE *fp)

{
  char c;
  int i;
  int lri;
  int lenS;
  pst_node *leafp;
  double loglike, persymbol;

  loglike = 0.0;
  lenS = strlen(S);
  for (i=0; i<lenS; i++)
    {
      c = S[i];
      S[i] = '\0';
      traverse_pst(AB, T, S, &leafp, &lri);
      persymbol = log10(leafp->p_c[ABind(AB,c)]);
      fprintf(fp, "%1.8f\n", persymbol);
      loglike += persymbol;
      S[i] = c;
    }
  return(loglike);
}


void predict_entry (char *AB, pst_type T, char *S, const int lwrite, 
		    FILE *lfp, const int swrite, FILE *sfp, char *D)
{
  char c;
  int i;
  int lri;
  int lenS, sumdepth;
  pst_node *leafp;
  double loglike, persymbol;

  sumdepth = 0;
  loglike = 0.0;
  lenS = strlen(S);
  if (lwrite)
    fprintf(lfp, "%7d\t%s\n", lenS, D);
  /* to read it in MATLAB: len = fscanf(fd, '%d', 1); */
  /*                       fgetl(fd); */
  /*                       d = fscanf(fd, '%d %f', [2 len]); */
  for (i=0; i<lenS; i++)
    {
      c = S[i];
      S[i] = '\0';
      traverse_pst(AB, T, S, &leafp, &lri);
      persymbol = leafp->p_c[ABind(AB,c)]; /* the prob */
      if (lwrite)
	fprintf(lfp, "%5d\t%1.8f\n", i-lri, persymbol);
      sumdepth += (i - lri);
      loglike += log10(persymbol);
      S[i] = c;
    }
  if (swrite)
    fprintf(sfp, "%7d\t%6d\t%1.8f\t%-40s<\n", lenS, sumdepth, loglike, D);
  /* to read it in MATLAB: [v] = fscanf(fd, '%d\t%d\t%f\t%*41c\n', [3 inf]); */
  return;
}



void emit_string (FILE *fp, pst_type T, char *AB, int absize, 
		  int string_len)
{
  int i,c;
  int lri;
  pst_node *leafp;
  double cumprob,rnd;
  char *S;
  
  if ((S = (char *) malloc (string_len+1)) == NULL)
    {
      fprintf (stderr,"failed to allocate new string buffer\n");
      exit(1);
    }

  for (i=0; i<string_len; i++)
    {
      S[i] = '\0';
      traverse_pst(AB, T, S, &leafp, &lri);
      rnd = 1.0*rand()/RAND_MAX;
      cumprob = 0.0;
      for (c=0; c<absize; c++)
	{
	  cumprob += leafp->p_c[c];
	  if (rnd <= cumprob)
	    break;
	}
      S[i] = AB[c];
    }

  fprintf(fp, "%s\n", S);

  free(S);
  return;
}


pst_node *add_dummy_pst_leaf(char *AB,int absize, pst_node *leafp, char c)
/* returns a pointer to the new leaf */

{
  pst_node *newleafp;

  newleafp = make_pst_node(absize); /* with smooth=0 */
  newleafp->p_c = leafp->p_c; /* physicially sharing his father's pred. vector */
  (leafp->sonsp)[ABind(AB,c)] = newleafp;
  return(newleafp);
}


void add_dummy_pst_node (char *AB, int absize, pst_type T, char *S, 
			 int *new_nodesp)
/* builds all nodes on path to S if needed */

{
  int lri; /* index of last read char in S (reading right to left) */
  pst_node *leafp;

  
  traverse_pst(AB, T, S, &leafp, &lri);

  while (lri>0)
    {
      leafp = add_dummy_pst_leaf(AB, absize, leafp, S[lri-1]);
      lri--;
      (*new_nodesp)++;
    }
  return;
}


void extend_pst_aux (char *AB, int absize, pst_type T0, pst_node *T, char *S, 
		     int *nodesp, int *new_nodesp)
/* S must be prorerly inited by print_leaves_pst */
/* T0 will maintain the root and T the current node */
{
  char c;
  int i;
  int leaf;

  (*nodesp)++;
  leaf=1;
  for (i=0; i<absize; i++)
    {
      /* original leaves have smooth=1 while dummy leaves get smooth=0 */
      /* we ignore the latter and find the leaves of the original tree */
      if (((T->sonsp)[i] != NULL) && (((T->sonsp)[i])->smooth==1))
	{
	  leaf=0;
	  S[-1]=AB[i];
	  extend_pst_aux(AB,absize,T0,(T->sonsp)[i],S-1, nodesp, new_nodesp);
	}
    }
  
  if (leaf)
    {
      /* add all prefixes of the leaf */
      for (i=1; i<strlen(S);i++)
	{
	  c = S[i];
	  S[i] = '\0';
	  add_dummy_pst_node(AB,absize,T0,S,new_nodesp);
	  S[i] = c;
	}
    }
}


void extend_pst(char *AB, int absize, int max_string,
		      pst_type T, int *nodesp, int *new_nodesp)

{
  char *S;

  if ((S = (char *) malloc (max_string+1)) == NULL)
    {
      fprintf (stderr,"failed to allocate leaves print buffer\n");
      exit(1);
    }

  S[max_string]='\0';
  extend_pst_aux (AB, absize, T, T, S+max_string, nodesp, new_nodesp);
  free(S);
}


void print_pfa_form_aux(char *AB, int absize, int max_string, pst_type T0, 
			pst_node *T, char *S, int lenS, FILE *lfp, FILE *sfp)
/* T is the current node while T0 is the tree root */
{
  const char *ROOT = "-"; /* root node label */
  int i;
  int lri;
  pst_node *nodep;

  fprintf(lfp,"%s\n", (lenS==0 ? ROOT : S));

  for(i=0; i<absize; i++)
    {
      fprintf(sfp, "%-*s    %c    ", max_string, (lenS==0 ? ROOT : S), AB[i]);
      S[lenS] = AB[i];
      traverse_pst(AB, T0, S, &nodep, &lri);
      fprintf(sfp, "%-*s     %1.8f\n", max_string, (lri==lenS+1 ? ROOT : S+lri), 
	      (T->p_c)[i]);
      S[lenS] = '\0';
    }

  for(i=0; i<absize; i++)
    if ((T->sonsp)[i] != NULL)
      {
	S[-1]=AB[i];
	print_pfa_form_aux(AB,absize,max_string,T0,(T->sonsp)[i],S-1,lenS+1,lfp,sfp);
      }
}


void print_pfa_form(char *AB, int absize, int max_string, pst_type T,
	       FILE *lfp, FILE *sfp)
/* lfp - labels, sfp - sparse matrix */
{
  char *S;
  
  /* +1 for the terminating '\0' and +1 for the next symbol */
  if ((S = (char *) malloc (max_string+2)) == NULL)
    {
      fprintf (stderr,"failed to allocate leaves print buffer\n");
      exit(1);
    }

  S[max_string+1] = '\0';
  S[max_string] = '\0';
  print_pfa_form_aux (AB, absize, max_string, T, T, S+max_string, 0, lfp, sfp);
  free(S);
}


  
