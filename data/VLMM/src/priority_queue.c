
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "priority_queue.h"


/* local */

pQ_node *make_pqueue_node (const char *fromS, const double p_S,
		     const int chi_Sx, const int chi_xS, const int chi_sufSx)

{
  pQ_node *pQ_np;

  if ((pQ_np = (pQ_node *) malloc (sizeof(pQ_node))) == NULL)
    {
      fprintf (stderr,"failed to allocate priority queue node\n");
      exit(1);
    }
  if ((fromS != NULL) && ((pQ_np->S = (char *) malloc(strlen(fromS)+1)) == NULL))
    {
      fprintf (stderr,"failed to allocate priority queue node string\n");
      exit(1);
    }
  if (fromS == NULL)
    pQ_np->S = NULL;
  else
    strcpy(pQ_np->S, fromS);
  pQ_np->p_S = p_S;
  pQ_np->chi_Sx = chi_Sx;
  pQ_np->chi_xS = chi_xS;
  pQ_np->chi_sufSx = chi_sufSx;
  pQ_np->nextp = NULL;

  return(pQ_np);
}

void free_pqueue_node (pQ_node *pQ_np)

{     
  if (pQ_np->S != NULL)
    free(pQ_np->S);
  free(pQ_np);
  return;
}


/* public */

pQ_type make_empty_pqueue ()

{
  return (NULL);
}


int not_empty_pqueue (pQ_type pQ)
/* returns (0) if queue is empty or (1) otherwise */

{
  return (pQ != NULL);
}


void dequeue_pqueue (pQ_type *pQp, char *S, double *p_Sp,
		     int *chi_Sxp, int *chi_xSp, int *chi_sufSxp)
/* dequeues first element of (*pQp). Assumes non empty queue! */
/* S is assumed pre-allocated and sufficient in size or NULL*/

{
  pQ_type newpQ;

  if (((*pQp)->S != NULL) && (S != NULL))
    strcpy(S,(*pQp)->S);
  (*p_Sp) = (*pQp)->p_S;
  (*chi_Sxp) = (*pQp)->chi_Sx;
  (*chi_xSp) = (*pQp)->chi_xS;
  (*chi_sufSxp) = (*pQp)->chi_sufSx;

  newpQ = (*pQp)->nextp;
  free_pqueue_node(*pQp);
  (*pQp) = newpQ;
  return;
}

void enqueue_pqueue (pQ_type *pQp, const char *fromS, const double p_S,
		     const int chi_Sx, const int chi_xS, const int chi_sufSx)

{
  pQ_node *pQnp;
  pQ_node *prevp;

  pQnp = make_pqueue_node (fromS, p_S, chi_Sx, chi_xS, chi_sufSx);

  if ((*pQp == NULL) || (p_S > (*pQp)->p_S)) 
    { /* the queue is empty or we insert as first */
      pQnp->nextp = (*pQp);
      (*pQp) = pQnp;
      return;
    }

  for (prevp = (*pQp);
       (prevp->nextp != NULL) && (p_S <= prevp->nextp->p_S);
       prevp = prevp->nextp);

  /* we insert as last or in the middle */
  pQnp->nextp = prevp->nextp; 
  prevp->nextp = pQnp;

  return;
}


void stack_pqueue (pQ_type *pQp, const char *fromS, const double p_S,
		     const int chi_Sx, const int chi_xS, const int chi_sufSx)

{
  pQ_node *pQnp;

  pQnp = make_pqueue_node (fromS, p_S, chi_Sx, chi_xS, chi_sufSx);

  /* we insert as first */
  pQnp->nextp = (*pQp);
  (*pQp) = pQnp;

  return;
}


void free_pqueue (pQ_type pQ)
{
  pQ_type newpQ;
  while (not_empty_pqueue(pQ))
    {
      newpQ = pQ->nextp;
      free_pqueue_node(pQ);
      pQ = newpQ;
    }
  return;
}
