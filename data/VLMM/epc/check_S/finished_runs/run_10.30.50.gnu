set term postscript eps color

set xrange [0:1]
set log y

set output "run_10.30.50.eps"
plot "run_10.30.50.0.1.dat" w linespoints ti "0.1", \
"run_10.30.50.0.2.dat" w linespoints ti "0.2", \
"run_10.30.50.0.3.dat" w linespoints ti "0.3", \
"run_10.30.50.0.4.dat" w linespoints ti "0.4", \
"run_10.30.50.0.5.dat" w linespoints ti "0.5", \
"run_10.30.50.0.6.dat" w linespoints ti "0.6", \
"run_10.30.50.0.7.dat" w linespoints ti "0.7", \
"run_10.30.50.0.8.dat" w linespoints ti "0.8", \
"run_10.30.50.0.9.dat" w linespoints ti "0.9", \
"run_10.30.50.1.1.dat" w linespoints ti "1.1", \
"run_10.30.50.1.2.dat" w linespoints ti "1.2", \
"run_10.30.50.1.3.dat" w linespoints ti "1.3", \
"run_10.30.50.1.4.dat" w linespoints ti "1.4", \
"run_10.30.50.1.5.dat" w linespoints ti "1.5", \
"run_10.30.50.1.6.dat" w linespoints ti "1.6", \
"run_10.30.50.1.7.dat" w linespoints ti "1.7", \
"run_10.30.50.1.8.dat" w linespoints ti "1.8", \
"run_10.30.50.1.9.dat" w linespoints ti "1.9", \
"run_10.30.50.1.dat" w linespoints ti "1.0"
