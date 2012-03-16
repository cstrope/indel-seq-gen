set term postscript eps color
set title "Change in Probability versus path length, 2-1"
set xlabel "Subpath length (t_z - t_{z-1})"
set ylabel "Forward probability change"
set output "2-1_pvl_plot.eps"
plot "2-1_pvl_plot.dat" w points ti "2-1"
