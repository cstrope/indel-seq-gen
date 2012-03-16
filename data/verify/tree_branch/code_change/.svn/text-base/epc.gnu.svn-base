set term postscript eps color

set xlabel "MCMC step number"
set ylabel "EPC probability"

plot "one_round_fwdepc-0.5-0.0.results" usi 1:3 w points ti "0.0", \
     "one_round_fwdepc-0.5-0.5.results" usi 1:3 w points ti "0.5", \
     "one_round_fwdepc-0.5-1.0.results" usi 1:3 w points ti "1.0"
