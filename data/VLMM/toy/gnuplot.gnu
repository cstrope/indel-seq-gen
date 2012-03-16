set term postscript eps color
set auto x
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.83105
set output "run_16.2.5.eps"
plot "run_16.2.5" w boxes
