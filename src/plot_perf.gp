# Fichier simple pour tracer les performances
set terminal png size 800,400
set output 'graph_simple.png'

set title "Performances des méthodes directes"
set xlabel "Taille n"
set ylabel "Temps (secondes)"
set grid

# Trace les 3 courbes
plot "mesures_lapack.txt" with linespoints title "LAPACK", \
     "mesures_manuel.txt" with linespoints title "Mon LU", \
     "mesures_dgbsv.txt" with linespoints title "dgbsv"