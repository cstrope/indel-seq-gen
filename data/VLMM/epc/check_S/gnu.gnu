set term postscript eps color

set log y
set xrange [2.75:3]

plot "qidotvqidot_k_" usi 1:2 w points ti "qi.t", "qidotvqidot_k_" usi 1:3 w points ti "qi.t+dt"
