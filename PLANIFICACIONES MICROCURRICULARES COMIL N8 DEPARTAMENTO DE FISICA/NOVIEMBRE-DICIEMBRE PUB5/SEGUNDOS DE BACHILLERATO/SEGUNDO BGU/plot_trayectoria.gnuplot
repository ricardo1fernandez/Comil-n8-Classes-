# plot_trayectoria.gnuplot
set datafile separator ","
set title "Trayectoria del Cohete de Agua (1L, 1/3 agua, 3 bar)"
set xlabel "Distancia horizontal (m)"
set ylabel "Altura (m)"
set grid
set key left top
plot "trajectory.csv" using 2:3 with lines lw 2 title "Trayectoria"
pause -1

