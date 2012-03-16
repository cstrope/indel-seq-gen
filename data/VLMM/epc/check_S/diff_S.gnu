set term postscript eps color

set xrange [0:1]

set output "diff_S.eps"
plot "run_13.29.16.0.1.dat" w linespoints ti "0.1", \
"run_13.29.16.0.2.dat" w linespoints ti "0.2", \
"run_13.29.16.0.3.dat" w linespoints ti "0.3", \
"run_13.29.16.0.4.dat" w linespoints ti "0.4", \
"run_13.29.16.0.5.dat" w linespoints ti "0.5", \
"run_13.29.16.0.6.dat" w linespoints ti "0.6", \
"run_13.29.16.0.7.dat" w linespoints ti "0.7", \
"run_13.29.16.0.8.dat" w linespoints ti "0.8", \
"run_13.29.16.0.9.dat" w linespoints ti "0.9", \
"run_13.29.16.1.1.dat" w linespoints ti "1.1", \
"run_13.29.16.1.2.dat" w linespoints ti "1.2", \
"run_13.29.16.1.3.dat" w linespoints ti "1.3", \
"run_13.29.16.1.4.dat" w linespoints ti "1.4", \
"run_13.29.16.1.5.dat" w linespoints ti "1.5", \
"run_13.29.16.1.6.dat" w linespoints ti "1.6", \
"run_13.29.16.1.7.dat" w linespoints ti "1.7", \
"run_13.29.16.1.8.dat" w linespoints ti "1.8", \
"run_13.29.16.1.9.dat" w linespoints ti "1.9", \
"run_13.29.16.1.dat" w linespoints ti "1.0"
