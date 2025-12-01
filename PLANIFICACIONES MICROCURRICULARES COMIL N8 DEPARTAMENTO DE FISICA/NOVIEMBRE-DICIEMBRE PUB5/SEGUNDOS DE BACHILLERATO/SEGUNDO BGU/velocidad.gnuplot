# velocidad
set datafile separator ","
set title "Velocidad del cohete"
set xlabel "Tiempo (s)"
set ylabel "Velocidad (m/s)"
plot "trajectory.csv" using 1:(sqrt($4*$4 + $5*$5)) with lines lw 2 title "Velocidad"
pause -1

