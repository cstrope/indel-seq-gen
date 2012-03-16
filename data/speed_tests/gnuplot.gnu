set term postscript eps color

set xrange [0:0.1]

plot "run_10.45.52.bins" usi 1:2 w linespoints, \
 "run_10.45.52.bins" usi 1:3 w linespoints
