set datafile separator ","
set title "Altura vs Tiempo"
set xlabel "Tiempo (s)"
set ylabel "Altura (m)"
plot "trajectory.csv" using 1:3 with lines lw 2 title "Altura"
pause -1

