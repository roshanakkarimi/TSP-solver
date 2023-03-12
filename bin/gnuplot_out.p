set terminal qt persist size 500,500
plot "../data/../data/att48.tsp" using 2:3 skip 6 with points
set title "data example"
