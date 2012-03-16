set term postscript eps color
set title "Events per Branch Length Simulated, 100 replicates"
set xlabel "Branch length simulated"
set ylabel "Number of substitution events"
set xrange [0:2.1]
set style line 1 lt 1 lw 3
set style line 2 lt 3 lw 3
set output "run_14.56.41.eps"
plot "run_14.56.41.avg_hits" usi 1:2:3 w yerrorbars ti "EPC" ls 1, \
 "run_14.56.41.avg_hits" usi 1:4:5 w yerrorbars ti "FWD" ls 2
