
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "pattern_match.h"


/* local */

int ABindex (const char *AB, char c)
{
  return(strchr(AB,c)-AB);
}


/* public */

void count_xSx (const char *S, const int m, char **ds_ptr, 
		int *chi_Sp, int *chi_xSp, int *chi_Sxp)

{
  int lenS;
  int i;
  char *p;
  
  *chi_Sp = 0;
  *chi_xSp = 0;
  *chi_Sxp = 0;
  lenS = strlen(S);
  
  if (lenS == 0)
    for(i=1; i<=m; i++)
      {
	(*chi_xSp) += strlen(ds_ptr[i]);
	(*chi_Sxp) += strlen(ds_ptr[i]);
      }
  else
    {
      for (i=1; i<=m; i++)
	{
	  p = ds_ptr[i];
	  while ((p=strstr(p,S))!=NULL)
	    {
	      (*chi_Sp)++;
	      (*chi_xSp) += (*(p-1)!='\0');
	      (*chi_Sxp) += (*(p+lenS) != '\0');
	      p++; /* so we don't find the same substring again */
	    }
	}
    }
  return;
}


int count_S (const char *S, const int m, char **ds_ptr)

{
  int chi_S;
  int dummy_xS;
  int dummy_Sx;
  count_xSx (S, m, ds_ptr, &chi_S, &dummy_xS, &dummy_Sx);
  return(chi_S);
}

int count_Sx (const char *S, const int m, char **ds_ptr)

{
  int dummy_S;
  int dummy_xS;
  int chi_Sx;
  count_xSx (S, m, ds_ptr, &dummy_S, &dummy_xS, &chi_Sx);
  return(chi_Sx);
}

double *proball_Sc (const char *AB, int absize, 
                  const int m, char **ds_ptr, char *S)
/* if S==NUL we use S[-1] */

{
  int lenS;
  int i;
  char *p;
  double *p_Sc;
  double count_sum;

  if ((p_Sc = (double *) malloc (absize*sizeof(double))) == NULL)
    {
      perror("failed to allocate a tree node probability buffer\n");
      exit(1);
    }

  lenS = strlen(S);
  if (lenS > 0)
    {
      for (i=0; i<absize; i++)
	p_Sc[i] = 0.0;
      for (i=1; i<=m; i++)
	{
	  p = ds_ptr[i];
	  while ((p=strstr(p,S)) != NULL)
	    {
	      if (p[lenS] != '\0')
		p_Sc[ABindex(AB,p[lenS])] += 1.0;
	      p++; /* so we don't find the same substring again */
	    }
	}
    }  
  else
    {
      for (i=0; i<absize; i++)
	{
	  S[-1] = AB[i];
	  p_Sc[i] = (double) count_S(S-1, m, ds_ptr);
	}
    }

  count_sum = 0;
  for (i=0; i<absize; i++)
    count_sum += p_Sc[i];
  for (i=0; i<absize; i++)
    p_Sc[i] = p_Sc[i]/count_sum;
  /* this is where we normalize. no smoothing (at least here) */
  return(p_Sc);
}



int count_diff_S (const char *S, const int m, char **ds_ptr)

{
  int i,c;
  
  if (strlen(S) == 0)
    return(0);

  c=0;
  for(i=1; i<=m; i++)
    c += (strstr(ds_ptr[i],S) != NULL);

  return(c);
}


void allchi_Sc (const char *AB, int absize, const int m, char **ds_ptr, 
		char *S, int chi_Sc[])
/* if S==NUL we use S[-1] */

{
  int lenS;
  int i;
  char *p;

  lenS = strlen(S);
  if (lenS > 0)
    {
      for (i=0; i<absize; i++)
	chi_Sc[i] = 0;
      for (i=1; i<=m; i++)
	{
	  p = ds_ptr[i];
	  while ((p=strstr(p,S)) != NULL)
	    {
	      if (p[lenS] != '\0')
		chi_Sc[ABindex(AB,p[lenS])] += 1;
	      p++; /* so we don't find the same substring again */
	    }
	}
    }  
  else
    {
      for (i=0; i<absize; i++)
	{
	  S[-1] = AB[i];
	  chi_Sc[i] = count_S(S-1, m, ds_ptr);
	}
    }

  return;
}


