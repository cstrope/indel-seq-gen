
/* Training a PST */
/* alg-1 is a GENERAL USE version */
/* alg-2 is a biological version */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "file_handle.h"
#include "learn_algs.h"

extern double p_min;
extern double alpha;
extern double gamma_min;
extern double p_ratio;
extern int L_max;

extern double shows;
extern double gama;
extern double mue;
extern double Q[];
extern const double q[23][23];

int main (int argc, char **argv)
{


  double smooth_p, input_shows=-1;  /*init to avoid warning */
  int absize, m, *li, nodes, leaves, pot_nodes, test_size, alg;
  int train_size, train_len, longl, longf, freql, freqf, dummy;
  char *AB, *train_set, *ab_file, *strings_file, *index_file;
  char *stats_file, *old_pst, *new_pst, *stats_entry;
  int stats_size;
  char *stats_sep;
  char **ts_ptr;
  pst_type T;
  FILE *fp;


  switch (argc)
    {
    case 1:  /* short reminder */
      system("clear");
      printf(" (clear screen)\n");
      fprintf(stderr,"GENERAL call format:\n\nGrowing/extending trees:\n\n");
      fprintf(stderr," %s <ab_file> <strings_file> <index_file/@> ...\n",argv[0]);
      fprintf(stderr,"\t... <old_pst/@> <new_pst/@> <stats_file[:size]/@> <[#]stats_entry/@> ...\n\n");
      fprintf(stderr,"(alg-1): ... 1 <p_min> <alpha> <gamma_min> <p_ratio> <L_max>\n");
      fprintf(stderr,"(alg-2): ... 2 <q_min/N_min> <gama> <p_ratio> <L_max>\n");
      fprintf(stderr,"\n\n");
      fprintf(stderr,"<strings_file> should be in FASTA format.\n");
      fprintf(stderr,"<index_file>   entries: 0 = read, other = ignore (no index => read all)\n");
      fprintf(stderr,"<.../@>        if file name starts with '@' it is ignored.\n");
      fprintf(stderr,"<stats_file>   add ':size' to control size of best leaf table.\n");
      fprintf(stderr,"<stats_entry>  if name begins with '#' (truncated) table file is created.\n");
      fprintf(stderr,"               otherwise a one liner is appended to the named file.\n\n");
      fprintf(stderr,"<q_min/N_min>  alg-2: if param<1 then q_min(%%) else N_min(#).\n");
      fprintf(stderr,"<ab_file>      alg-2: use ONLY with ab.PROTEINS_BZX !\n\n");
      exit(1);

    case 13: /* alg-2 */
      ab_file = argv[1];
      strings_file = argv[2];
      index_file = argv[3];
      old_pst = argv[4];
      new_pst = argv[5];
      stats_file = argv[6];
      stats_entry = argv[7];
      alg = atoi(argv[8]);
      /* the split point for other algs */
      if (alg != 2)
	{
	  fprintf(stderr,"Alg 2 syntax but not Alg 2.\n");
	  exit(1);
	}
      input_shows = atof(argv[9]); /* shows in calced out of it */
      gama = atof(argv[10]);
      p_ratio = atof(argv[11]);
      L_max = atoi(argv[12]);
      break;
      
    case 14: /* alg-1 */
      ab_file = argv[1];
      strings_file = argv[2];
      index_file = argv[3];
      old_pst = argv[4];
      new_pst = argv[5];
      stats_file = argv[6];
      stats_entry = argv[7];
      alg = atoi(argv[8]);
      /* the split point for other algs */
      if (alg != 1)
	{
	  fprintf(stderr,"Alg 1 syntax but not Alg 1.\n");
	  exit(1);
	}
      p_min = atof(argv[9]);
      alpha = atof(argv[10]);
      gamma_min = atof(argv[11]);
      p_ratio = atof(argv[12]);
      L_max = atoi(argv[13]);
      break;
      
    default:
      fprintf(stderr,"Syntax Error: type %s for a short reminder\n",argv[0]);
      exit(1);

    }
  

  /* train */

  T = NULL;
  if (*old_pst != '@')
    load_pst (old_pst, &AB, &absize, &dummy, &T);
  else
    read_alphabet(&absize,&AB,ab_file);

  /* now we can check whether alg-1 smoothing is possible */
  if (alg == 1) {
    smooth_p = 1 - absize*gamma_min;
    if (smooth_p <= 0)
      {
	fprintf(stderr, "Alg-1: cannot smooth. Decrease gamma_min.\n\n");
	exit(1);
      }
  }

  /* load strings into memory */
  read_fasta_set(&train_set, &m, &ts_ptr, &li, AB, strings_file, \
		 index_file, &test_size, &train_size, &train_len);

  if (alg == 2)
    shows = (input_shows>=1) ? input_shows/m : input_shows;
  /* when we divide show/m we effectively turn it into N_min */
  /* because alg-2 compares if (diff_chi_S/m >= shows) ...   */
  if (alg == 1)
    T = grow_tree1(AB, absize, m, ts_ptr, li, &pot_nodes, T);
  else
    T = grow_tree2(AB, absize, m, ts_ptr, li, &pot_nodes, T);

  if (*stats_file != '@' || *stats_entry != '@')
    sum_pst_stat(AB, absize, L_max, T, m, ts_ptr, &nodes, &leaves, \
		 &longl, &longf, &freql, &freqf);

  if (*stats_entry == '#') /* create table */
    {
      if ((fp = fopen(stats_entry+1, "w")) == NULL)
	{
	  perror("failed to open <stats_entry> file for create");
	  exit(1);
	}
      fprintf(fp, "\n");
      fprintf(fp, "strings_file: %s\n", strings_file);
      fprintf(fp, " index_file : %s\n\n", index_file);
      fprintf(fp, "train_size   train_len    avg_train_len   test_size\n");
      fprintf(fp, "----------   ---------    -------------   ---------\n");
      fprintf(fp, "  %4d      %7d           %6.1f       %4d\n",
	      train_size, train_len, (1.0*train_len)/train_size, test_size);
      fprintf(fp, "\n\n");
      fprintf(fp, "alg           parameters             pot_nodes  nodes  leaves nodes/leaf  ");
      fprintf(fp, "long_len long_freq  freq_len freq_freq     train_len\n");
      fprintf(fp, "---  -----------------------------   ---------  -----  ------ ----------  -------- ---------  -------- ---------     ---------\n");
    }
  else /* or open for append to table */
    if ((*stats_entry != '@') && ((fp = fopen(stats_entry, "a")) == NULL))
      {
	perror("failed to open <stats_entry> file for append");
	exit(1);
      }

  if (*stats_entry != '@') /* write entry and close */
    { 
      if (alg == 1)
	fprintf(fp, "%2d   %8.6f %2g %6.4f %6.3f %3d  ", alg, p_min, alpha, gamma_min, p_ratio, L_max);
      else
	fprintf(fp, "%2d   %10.6f  %6.4f %6.3f %3d  ", alg, input_shows, gama, p_ratio, L_max);
      fprintf(fp, "  %6d    %5d  %5d    %5.2f",pot_nodes, nodes, leaves,
	      (1.0*nodes)/leaves);
      fprintf(fp, "       %3d    %5d        %3d    %5d         %7d\n", longl, longf, freql, freqf, train_len);
      fclose(fp);
    }

  if (*stats_file != '@')
    {
       stats_sep = strchr(stats_file,':');
       if (stats_sep == NULL)
	 stats_size = 30;    /* default size */
       else
	 {
	    *stats_sep = '\0';
	    stats_size = atoi(stats_sep+1);
	 }
       if ((fp = fopen(stats_file, "w")) == NULL)
	{
	  perror("failed to open <stats_file> for write");
	  exit(1);
	}
      fprintf(fp, "\n");
      fprintf(fp, "strings_file: %s\n", strings_file);
      fprintf(fp, " index_file : %s\n", index_file);
      fprintf(fp, "  pst_file  : %s\n\n", new_pst);
      fprintf(fp, "train_size   train_len    avg_train_len   test_size\n");
      fprintf(fp, "----------   ---------    -------------   ---------\n");
      fprintf(fp, "  %4d      %7d           %6.1f       %4d\n",
	      train_size, train_len, (1.0*train_len)/train_size, test_size);
      fprintf(fp, "\n\n");
      fprintf(fp, "alg           parameters             pot_nodes  nodes  leaves nodes/leaf\n");
      fprintf(fp, "---  -----------------------------   ---------  -----  ------ ----------\n");
      if (alg == 1)
	fprintf(fp, "%2d   %8.6f %2g %6.4f %6.3f %3d  ", alg, p_min, alpha, gamma_min, p_ratio, L_max);
      else
	fprintf(fp, "%2d   %10.6f  %6.4f %6.3f %3d  ", alg, input_shows, gama, p_ratio, L_max);
      fprintf(fp, "  %6d    %5d  %5d    %5.2f\n\n\n",pot_nodes, nodes, leaves,
	      (1.0*nodes)/leaves);
      print_pst_stat (AB, absize, L_max, T, stats_size, m, ts_ptr, &nodes, &leaves, fp); 
      fclose(fp);
    }

  add_probs_pst(AB, absize, L_max, m, ts_ptr, T);

  if (alg == 1)
    smoothen_pst(absize, gamma_min, T);
  else
    pseudo_counts_smooth(AB, absize, L_max, m, ts_ptr, T, mue, q, Q);


  /* save PST */
  if (*new_pst != '@')
    save_pst(new_pst, AB, absize, L_max, T);
  
  return(0);
  
}    

