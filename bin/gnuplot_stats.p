set terminal qt persist size 700,500
set xlabel "Prob. of choosing nearest neighb."
show xlabel
set ylabel "Cost"
show ylabel
set xrange [0:100]
show xrange
set yrange [0:75000]
show yrange
plot "../out/stats_GRASP.dat" with points