
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "file_handle.h"
#include "priority_queue.h"

/* local */

int file_length (const char *filename)
{
  struct stat stbuf;

  if (stat(filename,&stbuf) == -1)
    {
      perror("failed to 'stat' for file length");
      exit(1);
    }

  return(stbuf.st_size);
}


void rec_save_pst (FILE *fp, int absize, pst_type T)
/* hard to follow by eye. easy to follow by load recursion */
{
  int i;

  for(i=0; i<absize; i++)
    {
      fprintf(fp,"%1.8f %d ",T->p_c[i], (T->sonsp)[i] != NULL);
      if ((T->sonsp)[i] != NULL)
	rec_save_pst(fp, absize, (T->sonsp)[i]);
    }
}


void rec_load_pst (FILE *fp, int absize, pst_type *Tp)

{
  int i;
  int son;

  (*Tp) = make_empty_pst(absize);

  (*Tp)->smooth = 1;

  if (((*Tp)->p_c = (double *) malloc (absize*sizeof(double))) == NULL)
    {
      perror("failed to allocate a tree node probability buffer\n");
      exit(1);
    }

  for(i=0; i<absize; i++)
    {
      fscanf(fp,"%lf %d ",(*Tp)->p_c+i, &son);
      if (son)
	rec_load_pst(fp, absize, (*Tp)->sonsp+i);
    }
}



/* public */



void read_alphabet (int *absizep, char **ABp, const char *filename)
{
  FILE *fp;

  if ((fp = fopen(filename, "r")) == NULL)
      {
	perror("failed to open <ab-file> for read");
	exit(1);
      }
  if ((*ABp = (char *) malloc (file_length(filename))) == NULL)
    {
      perror("failed to allocate AB\n");
      exit(1);
    }
  fscanf(fp,"%s",*ABp);
  *absizep = strlen(*ABp);
  fclose(fp);
  return;
}


void read_data_set (char **data_setp, int *mp, char ***ds_ptrp, int **lip,
		    const char *AB, const char *filename)
{
  FILE *fp;
  char cbuf[STRINGLEN_MAX+1];
  int cbuflen;
  int i;
  int loc;

  if ((fp = fopen(filename, "r")) == NULL)
      {
	perror("failed to open <train-set> file for read");
	exit(1);
      }

  printf("train-set file name: %s\n\n", filename);
  printf("reject for: ");
  i = 0;
  *mp = 0;
  loc = 0;  /* calcs li[0] = total string length */
  while (fscanf(fp,"%s",cbuf) != EOF)
    {
      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "single string size exceeds STRINGLEN_MAX = %d (may update and re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}
      if (strspn(cbuf,AB)!=cbuflen)
	{
	  i = 1;
	  printf ("%c", cbuf[strspn(cbuf,AB)]);
	}
      else
	{
	  (*mp)++;
	  loc += cbuflen;
	}
    }
  if (i==0)
    printf ("<none rejected>");
  printf("\n\nvalid strings = %d\ntotal lengths = %d\n", *mp, loc);
  rewind(fp);
  if ((*data_setp = (char *) malloc (loc+(*mp)+1)) == NULL)
    {
      perror("failed to allocate data set storage char array\n");
      exit(1);
    }
  if (((*ds_ptrp = (char **) malloc (((*mp)+1)*sizeof(char *))) == NULL) ||
      ((*lip     = (int *)   malloc (((*mp)+1)*sizeof(int   ))) == NULL))
    {
      perror("failed to allocate data set storage related variables\n");
      exit(1);
    }
  (*lip)[0] = loc;
  (*ds_ptrp)[0] = NULL;    /* we don't use it */
  (*data_setp)[0] = '\0';  /* we need it for efficient string count */

  i = 0;
  loc = 1; /* now calcs relative pos within *data_setp, [0] contains '\0' */
  while (fscanf(fp,"%s", cbuf) != EOF)
    {
      if (strspn(cbuf,AB)==(cbuflen=strlen(cbuf)))
	{
	  (*lip)[++i] = cbuflen;
	  strcpy((*data_setp)+loc,cbuf); /* trailing \0 is copied */
	  (*ds_ptrp)[i] = (*data_setp)+loc;
	  loc += cbuflen+1;
	}
    }

  fclose(fp);
  return;

}


void read_ind_data_set (char **data_setp, int *mp, char ***ds_ptrp, int **lip,
		    const char *AB, const char *filename, const char *indfile)
{
  FILE *fp;
  FILE *ifp;
  char cbuf[STRINGLEN_MAX+1];
  int cbuflen;
  int i;
  int ind;
  int loc;
  int tn;

  if ((fp = fopen(filename, "r")) == NULL)
      {
	perror("failed to open <strings_file> for read");
	exit(1);
      }

  if ((ifp = fopen(indfile, "r")) == NULL)
      {
	perror("failed to open <index_file> for read");
	exit(1);
      }

  printf("train-set strings file name: %s\n\n", filename);
  printf("train-set index   file name: %s\n\n", indfile);
  printf("reject for: ");
  tn = 0;
  i = 0;
  *mp = 0;
  loc = 0;  /* calcs li[0] = total string length */
  while (fscanf(fp,"%s",cbuf) != EOF)
    {
      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "single string size exceeds STRINGLEN_MAX = %d (may update and re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}

      if ((fscanf(ifp,"%d",&ind)) == EOF)
	{
	  fprintf(stderr,"<index_file> is shorter then <string_file>\n");
	  exit(1);
	}
      

      if (strspn(cbuf,AB)!=cbuflen)
	{
	  i = 1;
	  printf ("%c", cbuf[strspn(cbuf,AB)]);
	}
      else
	{
	  if (ind == 0) {   /* train */
	    (*mp)++;
	    loc += cbuflen;
	  }
	  else
	    tn++;
	}
    }
  if (i==0)
    printf ("<none rejected>");
  printf("\n\nfound %d valid train strings of total length %d\n", *mp, loc);
  printf("train string avg. length %0.1f\n", 1.0*loc/(*mp));
  printf("(plus %d valid unused test strings, %0.1f%% of the cluster)\n", \
	 tn, 100.0*tn/(tn+(*mp)));

  rewind(fp);
  rewind(ifp);

  if ((*data_setp = (char *) malloc (loc+(*mp)+1)) == NULL)
    {
      perror("failed to allocate data set storage char array\n");
      exit(1);
    }
  if (((*ds_ptrp = (char **) malloc (((*mp)+1)*sizeof(char *))) == NULL) ||
      ((*lip     = (int *)   malloc (((*mp)+1)*sizeof(int   ))) == NULL))
    {
      perror("failed to allocate data set storage related variables\n");
      exit(1);
    }
  (*lip)[0] = loc;
  (*ds_ptrp)[0] = NULL;    /* we don't use it */
  (*data_setp)[0] = '\0';  /* we need it for efficient string count */

  i = 0;
  loc = 1; /* now calcs relative pos within *data_setp, [0] contains '\0' */
  while (fscanf(fp,"%s", cbuf) != EOF)
    {
      fscanf(ifp,"%d", &ind);
      if ((strspn(cbuf,AB)==(cbuflen=strlen(cbuf))) && (ind==0))
	{
	  (*lip)[++i] = cbuflen;
	  strcpy((*data_setp)+loc,cbuf); /* trailing \0 is copied */
	  (*ds_ptrp)[i] = (*data_setp)+loc;
	  loc += cbuflen+1;
	}
    }

  fclose(fp);
  fclose(ifp);

  return;

}



void mread_ind_data_set (char **data_setp, int *mp, char ***ds_ptrp, \
			 int **lip, const char *AB, const char *filename, \
			 const char *indfile, int *tez, int *trz, int *trl)
{
  FILE *fp;
  FILE *ifp;
  char cbuf[STRINGLEN_MAX+1];
  int cbuflen;
  int i;
  int ind;
  int loc;
  int tn;

  if ((fp = fopen(filename, "r")) == NULL)
      {
	perror("failed to open <strings_file> for read");
	exit(1);
      }

  if ((ifp = fopen(indfile, "r")) == NULL)
      {
	perror("failed to open <index_file> for read");
	exit(1);
      }

  tn = 0;
  i = 0;
  *mp = 0;
  loc = 0;  /* calcs li[0] = total string length */
  while (fscanf(fp,"%s",cbuf) != EOF)
    {
      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "single string size exceeds STRINGLEN_MAX = %d (may update and re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}

      if (strspn(cbuf,AB)!=cbuflen)
	{
	  fprintf(stderr,"found a rejected string\n");
	  exit(1);
	}

      if ((fscanf(ifp,"%d",&ind)) == EOF)
	{
	  fprintf(stderr,"<index_file> is shorter then <string_file>\n");
	  exit(1);
	}
      

      if (ind == 0) {   /* train */
	(*mp)++;
	loc += cbuflen;
      }
      else
	tn++;
      
    }

  *tez = tn;
  *trz = *mp;
  *trl = loc;
    
  rewind(fp);
  rewind(ifp);

  if ((*data_setp = (char *) malloc (loc+(*mp)+1)) == NULL)
    { /* format: \0 string \0 ... \0 string \0 */
      perror("failed to allocate data set storage char array\n");
      exit(1);
    }
  if (((*ds_ptrp = (char **) malloc (((*mp)+1)*sizeof(char *))) == NULL) ||
      ((*lip     = (int *)   malloc (((*mp)+1)*sizeof(int   ))) == NULL))
    {
      perror("failed to allocate data set storage related variables\n");
      exit(1);
    }
  (*lip)[0] = loc;
  (*ds_ptrp)[0] = NULL;    /* we don't use it */
  (*data_setp)[0] = '\0';  /* we need it for efficient string count */

  i = 0;
  loc = 1; /* now calcs relative pos within *data_setp, [0] contains '\0' */
  while (fscanf(fp,"%s", cbuf) != EOF)
    {
      fscanf(ifp,"%d", &ind);
      if ((strspn(cbuf,AB)==(cbuflen=strlen(cbuf))) && (ind==0))
	{
	  (*lip)[++i] = cbuflen;
	  strcpy((*data_setp)+loc,cbuf); /* trailing \0 is copied */
	  (*ds_ptrp)[i] = (*data_setp)+loc;
	  loc += cbuflen+1;
	}
    }

  fclose(fp);
  fclose(ifp);

  return;

}



void read_fasta_set (char **data_setp, int *mp, char ***ds_ptrp, \
		     int **lip, const char *AB, const char *filename, \
		     const char *indfile, int *tez, int *trz, int *trl)
{
  FILE *fp;
  FILE *ifp;
  char cbuf[STRINGLEN_MAX+1];
  int cbuflen;
  int i;
  int ind;
  int loc;
  int indexed, flag, state, incoming;

  if ((fp = fopen(filename, "r")) == NULL)
      {
	perror("failed to open <strings_file> for read");
	exit(1);
      }

  indexed = (*indfile != '@');

  if ((indexed) && ((ifp = fopen(indfile, "r")) == NULL))
      {
	perror("failed to open <index_file> for read");
	exit(1);
      }

  /* 1st pass - count size to allocate, verify content */

  *tez = 0;
  *mp = 0;
  loc = 0;  /* calcs li[0] = total string length */
  flag = 1;
  state = -1;

  while (flag)
    {

      do
	{
	  if (fgets(cbuf, STRINGLEN_MAX, fp) == NULL)
	    {
	      incoming = 0; /* EOF or error :-( */
	      cbuflen = -1;
	    }
	  else
	    {
	      if ((cbuflen = strlen(cbuf)) >= STRINGLEN_MAX-1) /* max poss */
		{
		  fprintf(stderr, "single LINE size seems to exceeds STRINGLEN_MAX = %d\n (you may update and re-compile)\n", STRINGLEN_MAX);
		  exit(1);
		}
	      if (cbuf[cbuflen-1] == '\n')
		{
		  cbuf[--cbuflen] = '\0';
		}
	      if (cbuflen > 0)
		if (*cbuf == '>')
		  incoming = 1; /* comment line */
		else
		  {
		    incoming = 2; /* data line */
		    if (strspn(cbuf,AB)!=cbuflen)
		      {
			fprintf(stderr,"found a rejected data string:\n%s\n", cbuf);
			exit(1);
		      }
		  }
	    }
	}
      while (cbuflen == 0);
	
      switch (state)
	{
	case -1: /* start of file */
	  switch(incoming)
	    {
	    case 0:
	      fprintf(stderr,"Read fasta: empty string_file\n");
	      exit(1);
	      break;
	    case 1:
	      state = 0;
	      break;
	    case 2:
	      fprintf(stderr,"Read fasta: data start without comment line\n");
	      exit(1);
	      break;
	    default:
	      fprintf(stderr,"c-code bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }
	  break;

	case  0: /* after comment */
	  switch(incoming)
	    {
	    case 0:
	      fprintf(stderr,"Read fasta: comment line without data\n");
	      exit(1);
	      break;
	    case 1:
	      fprintf(stderr,"Read fasta: two consecutive comment lines\n");
	      exit(1);
	      break;
	    case 2:
	      if ((indexed) && ((fscanf(ifp,"%d",&ind)) == EOF))
		{
		  fprintf(stderr,"index_file is shorter then string_file\n");
		  exit(1);
		}
	      if ((indexed) && (ind))
		{
		  (*tez)++;
		  state = 2;
		}
	      else /* non-indexed || ind == 0 */
		{
		  (*mp)++;
		  loc += cbuflen;
		  state = 1;
		}
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;
	  
	case  1: /* in goodata */
	  switch(incoming)
	    {
	    case 0:
	      flag = 0;
	      break;
	    case 1:
	      state = 0;
	      break;
	    case 2:
	      loc += cbuflen;
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;

	case  2: /* in badata */
	  switch(incoming)
	    {
	    case 0:
	      flag = 0;
	      break;
	    case 1:
	      state = 0;
	      break;
	    case 2:
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;
	  
	default:
	  fprintf(stderr,"Program bug: impossible state in read fasta\n");
	  exit(1);
	  break;
	}
    }

  *trz = *mp;
  *trl = loc;
    
  rewind(fp);
  if (indexed)
    rewind(ifp);

  if ((*data_setp = (char *) malloc (loc+(*mp)+1)) == NULL)
    { /* format: \0 string \0 ... \0 string \0 */
      perror("failed to allocate data set storage char array\n");
      exit(1);
    }
  if (((*ds_ptrp = (char **) malloc (((*mp)+1)*sizeof(char *))) == NULL) ||
      ((*lip     = (int *)   malloc (((*mp)+1)*sizeof(int   ))) == NULL))
    {
      perror("failed to allocate data set storage related variables\n");
      exit(1);
    }
  (*lip)[0] = loc;
  (*ds_ptrp)[0] = NULL;    /* we don't use it */
  (*data_setp)[0] = '\0';  /* we need it for efficient string count */
  
  /* 2nd pass - read into memory */

  i = 0; /* string number, like *mp above */
  loc = 1; /* now calcs relative pos within *data_setp, [0] contains '\0' */
  flag = 1;
  state = -1;

  while (flag)
    {

      do
	{
	  if (fgets(cbuf, STRINGLEN_MAX, fp) == NULL)
	    {
	      incoming = 0; /* EOF or error :-( */
	      cbuflen = -1;
	    }
	  else
	    {
	      if ((cbuflen = strlen(cbuf)) >= STRINGLEN_MAX-1) /* max poss */
		{
		  fprintf(stderr, "single string size roughly exceeds STRINGLEN_MAX = %d (may update and re-compile)\n", STRINGLEN_MAX);
		  exit(1);
		}
	      if (cbuf[cbuflen-1] == '\n')
		{
		  cbuf[--cbuflen] = '\0';
		}
	      if (cbuflen > 0)
		if (*cbuf == '>')
		  incoming = 1; /* comment line */
		else
		  {
		    incoming = 2; /* data line */
		    if (strspn(cbuf,AB)!=cbuflen)
		      {
			fprintf(stderr,"found a rejected data string:\n%s\n", cbuf);
			exit(1);
		      }
		  }
	    }
	}
      while (cbuflen == 0);
      
      switch (state)
	{
	case -1: /* start of file */
	  switch(incoming)
	    {
	    case 0:
	      fprintf(stderr,"Read fasta: empty string_file\n");
	      exit(1);
	      break;
	    case 1:
	      state = 0;
	      break;
	    case 2:
	      fprintf(stderr,"Read fasta: data start without comment line\n");
	      exit(1);
	      break;
	    default:
	      fprintf(stderr,"c-code bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }
	  break;

	case  0: /* after comment */
	  switch(incoming)
	    {
	    case 0:
	      fprintf(stderr,"Read fasta: comment line without data\n");
	      exit(1);
	      break;
	    case 1:
	      fprintf(stderr,"Read fasta: two consecutive comment lines\n");
	      exit(1);
	      break;
	    case 2:
	      if ((indexed) && ((fscanf(ifp,"%d",&ind)) == EOF))
		{
		  fprintf(stderr,"index_file is shorter then string_file\n");
		  exit(1);
		}
	      if ((indexed) && (ind))
		state = 2;
	      else /* non-indexed || ind == 0 */
		{
		  (*lip)[++i] = cbuflen;
		  strcpy((*data_setp)+loc,cbuf); /* trailing \0 is copied */
		  (*ds_ptrp)[i] = (*data_setp)+loc;
		  loc += cbuflen; /* trailing \0 to be trampled unless last */
		  state = 1;
		}
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;
	  
	case  1: /* in goodata */
	  switch(incoming)
	    {
	    case 0:
	      loc++; /* last trailing \0 */
	      flag = 0;
	      break;
	    case 1:
	      loc++; /* string trailing \0 */
	      state = 0;
	      break;
	    case 2:
	      (*lip)[i] += cbuflen;
	      strcpy((*data_setp)+loc,cbuf); /* trailing \0 is copied */
	      loc += cbuflen; /* trailing \0 to be trampled unless last */
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;

	case  2: /* in badata */
	  switch(incoming)
	    {
	    case 0:
	      flag = 0;
	      break;
	    case 1:
	      state = 0;
	      break;
	    case 2:
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;
	  
	default:
	  fprintf(stderr,"Program bug: impossible state in read fasta\n");
	  exit(1);
	  break;
	}
    }


  fclose(fp);
  if (indexed)
    fclose(ifp);

  return;

}




void predict_fasta_set (const char *long_out, const char *short_out, \
			pst_type T, char *AB, const char *filename, \
			const char *indfile, int *tez, int *trz, int *trl)
{
  FILE *fp;
  FILE *ifp;
  FILE *lfp;
  FILE *sfp;
  char cbuf[STRINGLEN_MAX+1];
  int cbuflen;
  int ind;
  int loc;
  int indexed, flag, state, incoming, lwrite, swrite;
  char *entryp;
  int maxloc;
  char D[41]; /* first 40 chars of description line are echoed */

  if ((fp = fopen(filename, "r")) == NULL)
      {
	perror("failed to open <strings_file> for read");
	exit(1);
      }

  indexed = (*indfile != '@');

  lwrite = (*long_out != '@');
  swrite = (*short_out != '@');
  if (!(lwrite || swrite))
      {
	fprintf(stderr,"either long/short save format must be specified.\n");
	exit(1);
      }

  if ((indexed) && ((ifp = fopen(indfile, "r")) == NULL))
      {
	perror("failed to open <index_file> for read");
	exit(1);
      }

  /* 1st pass - count size to allocate, verify content */

  *trl = 0;    /* total predicting length */
  *tez = 0;    /* ignored size */
  *trz = 0;    /* predicting size */
  maxloc = 0;  /* calcs max entry length */
  flag = 1;
  state = -1;

  while (flag)
    {

      do
	{
	  if (fgets(cbuf, STRINGLEN_MAX, fp) == NULL)
	    {
	      incoming = 0; /* EOF or error :-( */
	      cbuflen = -1;
	    }
	  else
	    {
	      if ((cbuflen = strlen(cbuf)) >= STRINGLEN_MAX-1) /* max poss */
		{
		  fprintf(stderr, "single LINE size seems to exceeds STRINGLEN_MAX = %d\n (you may update and re-compile)\n", STRINGLEN_MAX);
		  exit(1);
		}
	      if (cbuf[cbuflen-1] == '\n')
		{
		  cbuf[--cbuflen] = '\0';
		}
	      if (cbuflen > 0)
		if (*cbuf == '>')
		  incoming = 1; /* comment line */
		else
		  {
		    incoming = 2; /* data line */
		    if (strspn(cbuf,AB)!=cbuflen)
		      {
			fprintf(stderr,"found a rejected data string:\n%s\n", cbuf);
			exit(1);
		      }
		  }
	    }
	}
      while (cbuflen == 0);
	
      switch (state)
	{
	case -1: /* start of file */
	  switch(incoming)
	    {
	    case 0:
	      fprintf(stderr,"Read fasta: empty string_file\n");
	      exit(1);
	      break;
	    case 1:
	      state = 0;
	      break;
	    case 2:
	      fprintf(stderr,"Read fasta: data start without comment line\n");
	      exit(1);
	      break;
	    default:
	      fprintf(stderr,"c-code bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }
	  break;

	case  0: /* after comment */
	  switch(incoming)
	    {
	    case 0:
	      fprintf(stderr,"Read fasta: comment line without data\n");
	      exit(1);
	      break;
	    case 1:
	      fprintf(stderr,"Read fasta: two consecutive comment lines\n");
	      exit(1);
	      break;
	    case 2:
	      if ((indexed) && ((fscanf(ifp,"%d",&ind)) == EOF))
		{
		  fprintf(stderr,"index_file is shorter then string_file\n");
		  exit(1);
		}
	      if ((indexed) && (ind))
		{
		  (*tez)++;
		  state = 2;
		}
	      else /* non-indexed || ind == 0 */
		{
		  (*trz)++;
		  (*trl) += cbuflen;
		  loc = cbuflen; /* restarted */
		  state = 1;
		}
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;
	  
	case  1: /* in goodata */
	  switch(incoming)
	    {
	    case 0:
	      if (loc > maxloc)
		maxloc = loc;
	      flag = 0;
	      break;
	    case 1:
	      if (loc > maxloc)
		maxloc = loc;
	      state = 0;
	      break;
	    case 2:
	      (*trl) += cbuflen;
	      loc += cbuflen;
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;

	case  2: /* in badata */
	  switch(incoming)
	    {
	    case 0:
	      flag = 0;
	      break;
	    case 1:
	      state = 0;
	      break;
	    case 2:
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;
	  
	default:
	  fprintf(stderr,"Program bug: impossible state in read fasta\n");
	  exit(1);
	  break;
	}
    }

  rewind(fp);
  if (indexed)
    rewind(ifp);

  if ((entryp = (char *) malloc (maxloc+1+1)) == NULL)
    { /* format: \0 string \0 */
      perror("failed to allocate max entry storage char array\n");
      exit(1);
    }
  entryp[0] = '\0';  /* we need it for efficient string count */

  /* 2nd pass - read entry into memory, predict, repeat */

  if ((lwrite) && ((lfp = fopen(long_out, "w")) == NULL))
      {
	perror("failed to open <long_out> for write");
	exit(1);
      }

  if ((swrite) && ((sfp = fopen(short_out, "w")) == NULL))
      {
	perror("failed to open <short_out> for write");
	exit(1);
      }

  flag = 1;
  state = -1;
  D[40] = '\0';

  while (flag)
    {

      do
	{
	  if (fgets(cbuf, STRINGLEN_MAX, fp) == NULL)
	    {
	      incoming = 0; /* EOF or error :-( */
	      cbuflen = -1;
	    }
	  else
	    {
	      if ((cbuflen = strlen(cbuf)) >= STRINGLEN_MAX-1) /* max poss */
		{
		  fprintf(stderr, "single string size roughly exceeds STRINGLEN_MAX = %d (may update and re-compile)\n", STRINGLEN_MAX);
		  exit(1);
		}
	      if (cbuf[cbuflen-1] == '\n')
		{
		  cbuf[--cbuflen] = '\0';
		}
	      if (cbuflen > 0)
		if (*cbuf == '>')
		  incoming = 1; /* comment line */
		else
		  {
		    incoming = 2; /* data line */
		    if (strspn(cbuf,AB)!=cbuflen)
		      {
			fprintf(stderr,"found a rejected data string:\n%s\n", cbuf);
			exit(1);
		      }
		  }
	    }
	}
      while (cbuflen == 0);
      
      switch (state)
	{
	case -1: /* start of file */
	  switch(incoming)
	    {
	    case 0:
	      fprintf(stderr,"Read fasta: empty string_file\n");
	      exit(1);
	      break;
	    case 1:
	      strncpy(D, cbuf, 40);
	      state = 0;
	      break;
	    case 2:
	      fprintf(stderr,"Read fasta: data start without comment line\n");
	      exit(1);
	      break;
	    default:
	      fprintf(stderr,"c-code bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }
	  break;

	case  0: /* after comment */
	  switch(incoming)
	    {
	    case 0:
	      fprintf(stderr,"Read fasta: comment line without data\n");
	      exit(1);
	      break;
	    case 1:
	      fprintf(stderr,"Read fasta: two consecutive comment lines\n");
	      exit(1);
	      break;
	    case 2:
	      if ((indexed) && ((fscanf(ifp,"%d",&ind)) == EOF))
		{
		  fprintf(stderr,"index_file is shorter then string_file\n");
		  exit(1);
		}
	      if ((indexed) && (ind))
		state = 2;
	      else /* non-indexed || ind == 0 */
		{
		  loc = 1; 
         /* restarted, calcs relative pos within entryp, [0] contains '\0' */
		  strcpy(entryp+loc,cbuf); /* trailing \0 is copied */
		  loc += cbuflen; 
         /* restarted, trailing \0 to be trampled unless last */
		  state = 1;
		}
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;
	  
	case  1: /* in goodata */
	  switch(incoming)
	    {
	    case 0:
	      predict_entry (AB, T, entryp+1, lwrite, lfp, swrite, sfp, D);
	      flag = 0;
	      break;
	    case 1:
	      predict_entry (AB, T, entryp+1, lwrite, lfp, swrite, sfp, D);
	      strncpy(D, cbuf, 40);
	      state = 0;
	      break;
	    case 2:
	      strcpy(entryp+loc,cbuf); /* trailing \0 is copied */
	      loc += cbuflen; /* trailing \0 to be trampled unless last */
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;

	case  2: /* in badata */
	  switch(incoming)
	    {
	    case 0:
	      flag = 0;
	      break;
	    case 1:
	      strncpy(D, cbuf, 40);
	      state = 0;
	      break;
	    case 2:
	      break;
	    default:
	      fprintf(stderr,"Program bug: impossible state in read fasta\n");
	      exit(1);
	      break;
	    }

	  break;
	  
	default:
	  fprintf(stderr,"Program bug: impossible state in read fasta\n");
	  exit(1);
	  break;
	}
    }


  fclose(fp);
  if (indexed)
    fclose(ifp);
  if (lwrite)
    fclose(lfp);
  if (swrite)
    fclose(sfp);

  return;

}




void save_pst (char *outfile, char *AB, int absize, int L_max, pst_type T)
/* clobbers existing if found */
{
  FILE *fp;

  if ((fp = fopen(outfile, "w")) == NULL)
      {
	perror("failed to open <tree-file> for write");
	exit(1);
      }
  fprintf(fp, "%d %s %d ", absize, AB, L_max);
  rec_save_pst (fp, absize, T);
  fclose(fp);
  return;
}


void load_pst (char *pstfile, char **ABp, int *absizep, int *L_maxp, 
	       pst_type *Tp)
{
  FILE *fp;

  if ((fp = fopen(pstfile, "r")) == NULL)
      {
	perror("failed to open <tree-file> for read");
	exit(1);
      }
  fscanf(fp, "%d ", absizep);
  if ((*ABp = (char *) malloc ((*absizep)+1)) == NULL)
    {
      perror("failed to allocate AB buffer\n");
      exit(1);
    }
  fscanf(fp,"%s %d ",*ABp, L_maxp);
  rec_load_pst (fp, *absizep, Tp);
  fclose(fp);
  return;
}




void process_test_set (char *filename, char *AB, pst_type T)
{
  FILE *fp;
  char cbuf[STRINGLEN_MAX+1];
  int cbuflen;
  int i;

  if ((fp = fopen(filename, "r")) == NULL)
      {
	perror("failed to open <test-set> file for read");
	exit(1);
      }

  printf("#no.\tlength\tlog10likeli\n");
  i = 0;
  while (fscanf(fp,"%s",cbuf) != EOF)
    {
      printf("(%d)\t", ++i);
      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "string size exceeds STRINGLEN_MAX = %d (may update & re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}
      printf("%d\t", strlen(cbuf));
      if (strspn(cbuf,AB)!=cbuflen)
	{
	  printf ("\t\treject for: %c\n", cbuf[strspn(cbuf,AB)]);
	}
      else
	{
	  printf("%g\n", log10like_on_pst(AB, T, cbuf));
	}
    }

  printf("\n");
  fclose(fp);
  return;

}



void fileproc_test_set (char *filename, char *resfile, char *AB, pst_type T)
{
  FILE *fp;
  FILE *rfp;
  char cbuf[STRINGLEN_MAX+1];
  int cbuflen;


  if ((fp = fopen(filename, "r")) == NULL)
      {
	perror("failed to open <strings_file> for read");
	exit(1);
      }

   if ((rfp = fopen(resfile, "a")) == NULL)
      {
	perror("failed to open <results_file> for append");
	exit(1);
      }

  while (fscanf(fp,"%s",cbuf) != EOF)
    {
      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "string size exceeds STRINGLEN_MAX = %d (may update & re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}
      if (strspn(cbuf,AB)!=cbuflen)
	{
	  fprintf (rfp,"%1.8f ", 0.0);
	}
      else
	{
	  fprintf(rfp, "%1.8f ", log10like_on_pst(AB, T, cbuf));
	}
    }


  fclose(fp);
  fclose(rfp);
  return;

}


void pred_onall_save (char *filename, char *resfile, char *AB, pst_type T, char detail)
{
  FILE *fp;
  FILE *rfp;
  char cbuf[STRINGLEN_MAX+1];
  int cbuflen;


  if ((fp = fopen(filename, "r")) == NULL)
      {
	perror("failed to open <strings_file> for read");
	exit(1);
      }

   if ((rfp = fopen(resfile, "w")) == NULL)
      {
	perror("failed to open <results_file> for write");
	exit(1);
      }

  while (fscanf(fp,"%s",cbuf) != EOF)
    {

      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "string size exceeds STRINGLEN_MAX = %d (may update & re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}
      if (strspn(cbuf,AB)!=cbuflen)
	{
	  fprintf(stderr,"found a rejected string\n");
	  exit(1);
	}

      if (detail)
	{
	  log10like_detail(AB, T, cbuf, rfp);
	  fprintf(rfp, "%1.8f\n", 1.0); /* separator */
	}
      else
	fprintf(rfp, "%1.8f\n", log10like_on_pst(AB, T, cbuf));

    }


  fclose(fp);
  fclose(rfp);
  return;

}

void prediction_stat (char *AB, pst_type T, const int n, const char *sfile, 
       const char *ifile, const char *s0file, const char *i0file, const int s1,
       const int s2, int s3, double *non1p, double *non2p, double *non12p, 
       double *iso1p, double *iso2p, double *iso12p)

{
  int e1, e2, e3, i, cbuflen, dum1, dum2;
  char cbuf[STRINGLEN_MAX+1], f1, f2, f12;
  double a, amin, iso3;
  pQ_type Q;
  FILE *fp, *ifp;

  amin = -DBL_MAX;
  Q = make_empty_pqueue();
  e1 = 0;      /* train  */
  e2 = 0;      /* test   */
  e3 = 0;      /* others */

  if ((fp = fopen(s0file, "r")) == NULL) 
    {
      perror("failed to open others strings file for read");
      exit(1);
    }
  if ((ifp = fopen(i0file, "r")) == NULL)
    {
      perror("failed to open others names file for read");
      exit(1);
    }

  dum1 = 0;
  while ((dum1<s3) && (fscanf(fp,"%s",cbuf) != EOF))
    {
      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "single string size exceeds STRINGLEN_MAX = %d (may update and re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}
      
      if (strspn(cbuf,AB)!=cbuflen)
	{
	  fprintf(stderr,"found a rejected string\n");
	  exit(1);
	}
      
      if ((fscanf(ifp,"%*s %*s %*d %d",&i)) == EOF)
	             /* id acc len cln */
	{
	  fprintf(stderr,"names file is shorter then others string file\n");
	  exit(1);
	}
      
      if (i != n)
	{
	  dum1++;
	  a = log10like_on_pst(AB, T, cbuf)/cbuflen;
	  if (a>amin)
	    amin = a;
	  enqueue_pqueue(&Q, NULL, a, 3, 0,0);
	}
    }
  fclose(fp);
  fclose(ifp);
  

/* we disable so that PFAM could be run...
  if (dum1 < s3)
    {
      fprintf(stderr,"not enough strings in 0.strings\n");
      exit(1);
    }
*/    
  
  if ((fp = fopen(sfile, "r")) == NULL) 
    {
      perror("failed to open cluster <strings_file> for read");
      exit(1);
    }
  if ((ifp = fopen(ifile, "r")) == NULL)
    {
      perror("failed to open cluster <index_file> for read");
      exit(1);
    }
  
  while (fscanf(fp,"%s",cbuf) != EOF)
    {
      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "single string size exceeds STRINGLEN_MAX = %d (may update and re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}

      if (strspn(cbuf,AB)!=cbuflen)
	{
	  fprintf(stderr,"found a rejected string\n");
	  exit(1);
	}

      if ((fscanf(ifp,"%d",&i)) == EOF)  /* 0=train and 1=test */
	{
	  fprintf(stderr,"<index_file> is shorter then <string_file>\n");
	  exit(1);
	}
      
      a = log10like_on_pst(AB, T, cbuf)/cbuflen;
      if (a<=amin)
	{
	  enqueue_pqueue(&Q, NULL, a, i+1, 0,0);
	  ((i == 0) ? e1++ : e2++);
	}
    }

  fclose(fp);
  fclose(ifp);

  *non1p = 1.0*e1; /* /s1; */
  *non2p = 1.0*e2; /* /s2; */
  *non12p = 1.0*(e1+e2); /* /(s1+s2); */
  *iso1p = *non1p;
  *iso2p = *non2p;
  *iso12p = *non12p;
  iso3 = 0.0;

  f1 = (*iso1p > iso3);
  f2 = (*iso2p > iso3);
  f12 = (*iso12p > iso3);

  while (f1+f2+f12 > 0)
    {
      dequeue_pqueue(&Q, NULL, &a, &i, &dum1, &dum2);
      switch(i)
	{
	  case 1:
	    e1--;
	    if (f1)
	      *iso1p = 1.0*e1; /* /s1; */
	    if (f12)
	      *iso12p = 1.0*(e1+e2); /* /(s1+s2); */
	    break;
	  case 2:
	    e2--;
	    if (f2)
	      *iso2p = 1.0*e2; /* /s2; */
	    if (f12)
	      *iso12p = 1.0*(e1+e2); /* /(s1+s2); */
	    break;
	  case 3:  
	    e3++;
	    iso3 = 1.0*e3; /* /s3; */
	    break;
	}
      
      f1 = (*iso1p > iso3);
      f2 = (*iso2p > iso3);
      f12 = (*iso12p > iso3);
    }


  free_pqueue(Q);
  free(cbuf);

  return;
}



void pred_stat_save (char *AB, pst_type T, const int n, const char *sfile, 
       const char *ifile, const char *s0file, const char *i0file, int s1,
       int s2, int s3, double *non1p, double *non2p, double *non12p, 
       double *iso1p, double *iso2p, double *iso12p, const char *wfile,
       const char *w0file, const int tsize)

{
  int e1, e2, e3, i, cbuflen, dum1, dum2;
  char cbuf[STRINGLEN_MAX+1], f1, f2, f12, nam[20], acc[20]; /* PFAM disabled comb[40]; */
  double a, amin, iso3, y;
  pQ_type Q; /* PFAM disabled Qh; */
  FILE *fp, *ifp, *wfp;

  amin = -DBL_MAX;
  Q = make_empty_pqueue();
  e1 = 0;      /* train  */
  e2 = 0;      /* test   */
  e3 = 0;      /* others */

  if ((wfp = fopen(w0file, "w")) == NULL) 
    {
      perror("failed to open others output file for write");
      exit(1);
    }
  if ((fp = fopen(s0file, "r")) == NULL) 
    {
      perror("failed to open others strings file for read");
      exit(1);
    }
  if ((ifp = fopen(i0file, "r")) == NULL)
    {
      perror("failed to open others names file for read");
      exit(1);
    }
  dum1 = 0;
  while (fscanf(fp,"%s",cbuf) != EOF)
    {
      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "single string size exceeds STRINGLEN_MAX = %d (may update and re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}
      
      if (strspn(cbuf,AB)!=cbuflen)
	{
	  fprintf(stderr,"found a rejected string\n");
	  exit(1);
	}
      
      if ((fscanf(ifp,"%s %s %d %d",nam, acc, &dum2, &i)) == EOF)
	{
	  fprintf(stderr,"names file is shorter then others string file\n");
	  exit(1);
	}
      
      if ((i != n) && (dum1<s3))
	{
	  dum1++;
	  y = log10like_on_pst(AB, T, cbuf);
	  fprintf(wfp, "%f\n", y);
/*	  a = y/cbuflen;
	  if (a>amin)
	    amin = a;
	  sprintf(comb, "%15s %10s", nam, acc);
	  enqueue_pqueue(&Q, comb, a, 3, i, dum2); */
	}             /* id_acc kind cluster len */
      else
	fprintf(wfp, "%f\n", 1.0);
    }
  fclose(fp);
  fclose(ifp);
  fclose(wfp);

/*  
  printf("\nthe %d most related others strings:\n \t slope\t  len\tcluster\tname\n", tsize/2);
  Qh = Q;
  for (i=1; (i<=tsize/2) && (Qh!=NULL); i++)
    {
      printf("(%d)\t%5.3f\t%5d\t%5d\t%s\n", i, -Qh->p_S, Qh->chi_sufSx, 
	     Qh->chi_xS,Qh->S);
      Qh = Qh->nextp;
    }

  if (dum1 < s3)
    {
      fprintf(stderr,"not enough strings in 0.strings\n");
      exit(1);
    }
*/    
  
  if ((wfp = fopen(wfile, "w")) == NULL) 
    {
      perror("failed to open cluster output file for write");
      exit(1);
    }
  if ((fp = fopen(sfile, "r")) == NULL) 
    {
      perror("failed to open cluster <strings_file> for read");
      exit(1);
    }
  if ((ifp = fopen(ifile, "r")) == NULL)
    {
      perror("failed to open cluster <index_file> for read");
      exit(1);
    }
  
  while (fscanf(fp,"%s",cbuf) != EOF)
    {
      if ((cbuflen = strlen(cbuf)) > STRINGLEN_MAX)
	{
	  fprintf(stderr, "single string size exceeds STRINGLEN_MAX = %d (may update and re-compile)\n", STRINGLEN_MAX);
	  exit(1);
	}

      if (strspn(cbuf,AB)!=cbuflen)
	{
	  fprintf(stderr,"found a rejected string\n");
	  exit(1);
	}

      if ((fscanf(ifp,"%d",&i)) == EOF)  /* 0=train and 1=test */
	{
	  fprintf(stderr,"<index_file> is shorter then <string_file>\n");
	  exit(1);
	}
      
      y = log10like_on_pst(AB, T, cbuf);
      fprintf (wfp, "%f\n", y);
/*      a = y/cbuflen;
      if (a<=amin)
	{
	  enqueue_pqueue(&Q, NULL, a, i+1, 0,0);
	  ((i == 0) ? e1++ : e2++);
	} */
    }

  fclose(fp);
  fclose(ifp);
  fclose(wfp);

  return; /* PFAM!!! */

  *non1p = 1.0*e1; /* /s1; */
  *non2p = 1.0*e2; /* /s2; */
  *non12p = 1.0*(e1+e2); /* /(s1+s2); */
  *iso1p = *non1p;
  *iso2p = *non2p;
  *iso12p = *non12p;
  iso3 = 0.0;

  f1 = (*iso1p > iso3);
  f2 = (*iso2p > iso3);
  f12 = (*iso12p > iso3);

  while (f1+f2+f12 > 0)
    {
      dequeue_pqueue(&Q, NULL, &a, &i, &dum1, &dum2);
      switch(i)
	{
	  case 1:
	    e1--;
	    if (f1)
	      *iso1p = 1.0*e1; /* /s1; */
	    if (f12)
	      *iso12p = 1.0*(e1+e2); /* /(s1+s2); */
	    break;
	  case 2:
	    e2--;
	    if (f2)
	      *iso2p = 1.0*e2; /* /s2; */
	    if (f12)
	      *iso12p = 1.0*(e1+e2); /* /(s1+s2); */
	    break;
	  case 3:  
	    e3++;
	    iso3 = 1.0*e3; /* /s3; */
	    break;
	}
      
      f1 = (*iso1p > iso3);
      f2 = (*iso2p > iso3);
      f12 = (*iso12p > iso3);
    }


  free_pqueue(Q);
  free(cbuf);

  return;
}


