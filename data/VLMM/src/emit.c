
/* Emitting a string from an existing PST */
/* this is a GENERAL USE version */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "file_handle.h"
#include "learn_algs.h"


int main (int argc, char **argv)
{


  int absize, string_len;
  int max_L, seed;
  char *AB, *pst_file, *out_file;
  pst_type T;
  FILE *fp;


  switch (argc)
    {
    case 1:  /* short reminder */
      system("clear");
      printf(" (clear screen)\n");
      fprintf(stderr,"GENERAL call format:\n\nString Emission from a given PST:\n\n");
      fprintf(stderr," %s <pst_file> <out_file/@> <string_len> [<random_seed>]\n\n\n",argv[0]);
      fprintf(stderr,"<out_file>    will be written in FASTA format.\n");
      fprintf(stderr,"              if starts with '#' (truncated) the string is appended to out_file\n");
      fprintf(stderr,"              if starts with '@' stdout is used\n");
      fprintf(stderr,"<string_len>  the length of the emitted string\n");
      fprintf(stderr,"<random_seed> if unspecified time is used to init srand\n\n");
      exit(1);

    case 4: case 5:
      pst_file = argv[1];
      out_file = argv[2];
      string_len = atoi(argv[3]);
      if (argc == 5)
	seed = atoi(argv[4]);
      else
	seed = time(NULL);
      break;
      
    default:
      fprintf(stderr,"Syntax Error: type %s for a short reminder\n",argv[0]);
      exit(1);

    }
  

  /* emit */


  load_pst (pst_file, &AB, &absize, &max_L, &T);

  if (*out_file == '@')
    fp = stdout;
  else 
    if (*out_file == '#')
      {
	if ((fp = fopen(out_file+1, "a")) == NULL)
	  {
	    perror("failed to open <out_file> for append");
	    exit(1);
	  }
      }
    else
      if ((fp = fopen(out_file, "w")) == NULL)
	{
	  perror("failed to open <out_file> for write");
	  exit(1);
	}
  
  srand(seed);

  fprintf(fp, ">pst_file: %s   emitted string length: %d   srand seed: %d\n", pst_file, string_len, seed);

  emit_string(fp, T, AB, absize, string_len);

  fclose(fp);

  return(0);
  
}    

