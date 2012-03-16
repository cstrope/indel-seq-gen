set term postscript eps color
set title "Percent accepted path lengths, 2-1"
set xlabel "Subpath length (t_E - t_B)"
set ylabel "Percent M-H acceptance"
set output "2-1_accepted_histo.eps"
plot "2-1_accepted_histo.dat" w points ti "2-1"
