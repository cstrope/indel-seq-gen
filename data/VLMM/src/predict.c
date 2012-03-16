
/* Prdeicting from an existing PST */
/* this is a GENERAL USE version */

#include <stdlib.h>
#include <stdio.h>
#include "file_handle.h"
#include "learn_algs.h"


int main (int argc, char **argv)
{


  int absize, test_size;
  int train_size, train_len, max_L;
  char *AB, *pst_file, *strings_file, *index_file;
  char *long_out, *short_out, *log_file;
  pst_type T;
  FILE *fp;


  switch (argc)
    {
    case 1:  /* short reminder */
      system("clear");
      printf(" (clear screen)\n");
      fprintf(stderr,"GENERAL call format:\n\nPrediction using a given PST:\n\n");
      fprintf(stderr," %s <pst_file> <strings_file> <index_file/@> ...\n",argv[0]);
      fprintf(stderr,"(outputs): ... <long_out/@> <short_out/@> <log_file/@>\n\n\n");
      fprintf(stderr,"<strings_file> should be in FASTA format.\n");
      fprintf(stderr,"<index_file>   entries: 0 = read, other = ignore (no index => read all)\n");
      fprintf(stderr,"<.../@>        if file name starts with '@' it is ignored.\n");
      fprintf(stderr,"<long_out>     {string_len comment {depth like} x string_len} x no_of_strings\n");
      fprintf(stderr,"<short_out>    {string_len sum_depth sum_log10like comment} x no_of_strings\n\n");
      exit(1);

    case 7:
      pst_file = argv[1];
      strings_file = argv[2];
      index_file = argv[3];
      long_out = argv[4];
      short_out = argv[5];
      log_file = argv[6];
      break;
      
    default:
      fprintf(stderr,"Syntax Error: type %s for a short reminder\n",argv[0]);
      exit(1);

    }
  

  /* predict */

  load_pst (pst_file, &AB, &absize, &max_L, &T);

  predict_fasta_set(long_out, short_out, T, AB, strings_file, index_file, 
		    &test_size, &train_size, &train_len);

  if (*log_file != '@')
    {
      if ((fp = fopen(log_file, "w")) == NULL)
        {
          perror("failed to open <log_file> for write");
          exit(1);
        }
      fprintf(fp, "\n");
      fprintf(fp, "strings_file: %s\n", strings_file);
      fprintf(fp, " index_file : %s\n\n", index_file);
      fprintf(fp, "  pst_file  : %s\n", pst_file);
      fprintf(fp, "  long_out  : %s\n", long_out);
      fprintf(fp, " short_out  : %s\n\n", short_out);
      fprintf(fp, "predict_size predict_len    avg_predict_len   ignored_size\n");
      fprintf(fp, "------------ -----------    ---------------   ------------\n");
      fprintf(fp, "  %6d    %9d           %7.1f          %4d\n",
              train_size, train_len, (1.0*train_len)/train_size, test_size);
      fprintf(fp, "\n\n");
      fclose(fp);
    }


  return(0);
  
}    

