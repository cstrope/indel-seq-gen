set term postscript eps color

set log y

set xrange [0:1]

plot \
"run_10.59.39-0.38.dat" usi 1:2 w points ti "0.38", \
"run_10.59.39-1.dat" usi 1:2 w points ti "1.0", \
"run_10.59.39-2.dat" usi 1:2 w points ti "2.0", \
"run_10.59.39-100.dat" usi 1:2 w points ti "100", \
"run_10.59.39-1e-06.dat" usi 1:2 w points ti "0.000001", \
"run_10.59.39-0.38.dat" usi 1:3 w points ti "0.38", \
"run_10.59.39-1.dat" usi 1:3 w points ti "1.0", \
"run_10.59.39-2.dat" usi 1:3 w points ti "2.0", \
"run_10.59.39-100.dat" usi 1:3 w points ti "100", \
"run_10.59.39-1e-06.dat" usi 1:3 w points ti "0.000001"
