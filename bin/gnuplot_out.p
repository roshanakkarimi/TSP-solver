set terminal qt persist size 500,500
plot "../out/h_greedy.dat" using 2:3:1 with labels, \
"" skip 51 with linespoints
