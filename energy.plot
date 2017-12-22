set grid
set title "Energy vs Time"
set xlabel "Time"
set ylabel "Energy"
plot "argon_energy.dat" using 1:2 with line title "Potential U", \
     "argon_energy.dat" using 1:3 with line title "Knetics K", \
     "argon_energy.dat" using 1:4 with line title "Total (K+U) "
