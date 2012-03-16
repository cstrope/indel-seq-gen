#ifndef _priority_queue_h
#define _priority_queue_h

typedef struct pQ_struct {
  char *S;
  double p_S;
  int chi_Sx;
  int chi_xS;
  int chi_sufSx;
  struct pQ_struct *nextp; } pQ_node;

typedef pQ_node *pQ_type;

/* p_S is the insertion key. the biggest p_S is the 
   first element of the queue, the first to be dequeued */


/* public functions */

pQ_type make_empty_pqueue ();

int not_empty_pqueue (pQ_type pQ);

void enqueue_pqueue (pQ_type *pQp, const char *fromS, const double p_S,
		     const int chi_Sx, const int chi_xS, const int chi_sufSx);
/* using p_S for key. the biggest p_s is first to pop */

void stack_pqueue (pQ_type *pQp, const char *fromS, const double p_S,
		     const int chi_Sx, const int chi_xS, const int chi_sufSx);
/* stacks upfront regardless of p_S. turns the queue into a LIFO stack */

void dequeue_pqueue (pQ_type *pQp, char *S, double *p_Sp,
		     int *chi_Sxp, int *chi_xSp, int *chi_sufSxp);
/* dequeues first element of (*pQp). Assumes non empty queue! */
/* S is assumed to be pre-allocated adequately or NULL*/

void free_pqueue (pQ_type pQ);
/* may be empty */

#endif /*_priority_queue_h */

