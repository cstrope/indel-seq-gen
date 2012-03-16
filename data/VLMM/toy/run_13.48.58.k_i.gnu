set term postscript eps color
set xrange [0:1]
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.001
set output "run_13.48.58.k_i.eps"
plot "run_13.48.58.k_i.dat" w boxes
