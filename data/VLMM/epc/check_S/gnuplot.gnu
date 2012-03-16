set term postscript eps color

set ylabel "Forward log Probability"
set y2label "Number of Substitutions"
set y2range [300:600]
set y2tic 300,50
set xlabel "Transition rates Scale factor (S)"

set output "gnuplot.eps"
plot "run_16.8.12.dat" using 1:2:3 axes x1y1 noti with yerrorbars, \
     "run_16.8.12.dat" using 1:4:5 axes x1y2 noti with yerrorbars
