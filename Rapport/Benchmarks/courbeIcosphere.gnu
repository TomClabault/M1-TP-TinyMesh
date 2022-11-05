# Définition du format de sortie
set terminal pngcairo

# Définition du nom du fichier de sortie
set output 'courbeIcosphere.png'

# Définition des titres de la figure
set title 'Temps de génération \& Nb triangles / Nombre subdivisions' # titre
set xlabel 'Nombre de subdivisions (icosphère de subdivisions)'     # nom de l'axe des abscisses
set ylabel 'Temps de génération (ms)'   # nom de l'axe des ordonnées
set y2label 'Nombre de triangle (*10^3)'   # nom de l'axe des ordonnées
set y2tics nomirror

plot 'benchmarkDurationIcosphere.dat' axis x1y1 with linespoints title "Temps génération", 'benchmarkTriangleIcosphere.dat' axis x1y2 with linespoints title "Nb triangles"
