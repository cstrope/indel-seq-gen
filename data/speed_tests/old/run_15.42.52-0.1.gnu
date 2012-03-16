set term postscript eps color

set output "run_15.42.52-0.1.eps"

plot "run_15.42.52-0.1.bins" usi 1:2 w points ti "EPC",\
"run_15.42.52-0.1.bins" usi 1:3 w points ti "FWD"
