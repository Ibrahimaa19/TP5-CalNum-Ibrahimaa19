set terminal pngcairo size 1000,500 font "Arial,12"
set output 'perf_TP_loglog.png'

set title "Performances des méthodes directes (échelle log-log)"
set xlabel "Taille du système n"
set ylabel "Temps d'exécution (secondes)"

set logscale xy
set grid
set format x "10^{%L}"
set format y "10^{%L}"

set key top left

# Trace les courbes avec lignes épaisses et points
plot "perf_TP.txt" using 1:2 with linespoints lw 3 pt 7 title "DGBTRF + DGBTRS", \
     "perf_TP.txt" using 1:3 with linespoints lw 3 pt 5 title "LU tridiagonale", \
     "perf_TP.txt" using 1:4 with linespoints lw 3 pt 9 title "DGBSV"
