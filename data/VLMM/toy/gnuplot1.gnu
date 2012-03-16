set term postscript eps color
set auto x
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 8.882
set output "run_10.53.46.1.eps"
plot "run_10.53.46.1" w boxes
