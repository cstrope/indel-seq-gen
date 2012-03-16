set term postscript eps color

set xlabel "MCMC cycle number"
set ylabel "Forward probability"

#set yrange [-650:-540]
#set xrange [0:10000]

set output "5-1.0_F594.218.eps"
plot "one_round_fwdepc-5-1.0.results" usi 1:7 w points ti "Forward Probability per MCMC cycle"

#set yrange [210:250]
#set ylabel "Number of Substitutions"
#set output "5-1.0_S230.eps"
#plot "one_round_fwdepc-5-1.0.results" usi 1:9 w points ti "Substitutions per MCMC cycle"
