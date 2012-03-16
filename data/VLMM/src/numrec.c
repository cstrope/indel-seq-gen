
#include <math.h> 
#include <stdio.h>
#include <stdlib.h>

/* pulled out of nrutil.c */

#define NR_END 1 
#define FREE_ARG char* 

void nrerror(char error_text[]) 
/* Numerical Recipes standard error handler */ 
{ 
  fprintf(stderr,"Numerical Recipes run-time error...\n"); 
  fprintf(stderr,"%s\n",error_text); 
  fprintf(stderr,"...now exiting to system...\n"); 
  exit(1); 
} 
      
double *dvector(long nl, long nh) 
     /* allocate a double vector with subscript range v[nl..nh] */ 
{ 
  double *v; 
  
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double))); 
  if (!v) nrerror("allocation failure in dvector()"); 
  return v-nl+NR_END; 
} 

void free_dvector(double *v, long nl, long nh) 
/* free a double vector allocated with dvector() */ 
{ 
  free((FREE_ARG) (v+nl-NR_END)); 
} 



/* NUMREC ch.2 routines */

void ludcmp(double **a, int n, int *indx, double *d) 
{ 
  int i,imax,j,k; 
  double big,dum,sum,temp; 
  double *vv; 
  const double TINY = 1.0e-20;
  
  vv=dvector(1,n);  /* jill: changed from vector to dvector */
  *d=1.0; 
  for (i=1;i<=n;i++) { 
    big=0.0; 
    for (j=1;j<=n;j++) 
      if ((temp=fabs(a[i][j])) > big) big=temp; 
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp"); 
    vv[i]=1.0/big; 
  } 
  for (j=1;j<=n;j++) { 
    for (i=1;i<j;i++) { 
      sum=a[i][j]; 
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j]; 
      a[i][j]=sum; 
    } 
    big=0.0; 
    for (i=j;i<=n;i++) { 
      sum=a[i][j]; 
      for (k=1;k<j;k++) 
	sum -= a[i][k]*a[k][j]; 
      
      a[i][j]=sum; 
      if ( (dum=vv[i]*fabs(sum)) >= big) { 
	big=dum; 
	imax=i; 
      } 
    } 
    if (j != imax) { 
      for (k=1;k<=n;k++) { 
	dum=a[imax][k]; 
	a[imax][k]=a[j][k]; 
	a[j][k]=dum; 
      } 
      *d = -(*d); 
      vv[imax]=vv[j]; 
    } 
    indx[j]=imax; 
    if (a[j][j] == 0.0) a[j][j]=TINY; 
    if (j != n) { 
      dum=1.0/(a[j][j]); 
      for (i=j+1;i<=n;i++) a[i][j] *= dum; 
    } 
  } 
  free_dvector(vv,1,n); /* jill: changed from vector to dvector */
} 


void lubksb(double **a, int n, int *indx, double b[]) 
{ 
  int i,ii=0,ip,j; 
  double sum; 
  
  for (i=1;i<=n;i++) { 
    ip=indx[i]; 
    sum=b[ip]; 
    b[ip]=b[i]; 
    if (ii) 
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j]; 
    else if (sum) ii=i; 
    b[i]=sum; 
  } 
  for (i=n;i>=1;i--) { 
    sum=b[i]; 
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j]; 
    b[i]=sum/a[i][i]; 
  } 
}
