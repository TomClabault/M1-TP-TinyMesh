# Définition du format de sortie
set terminal pngcairo

# Définition du nom du fichier de sortie
set output 'courbeCylinder.png'

# Définition des titres de la figure
set title 'Temps de génération \& Nb triangles / Nombre subdivisions' # titre
set xlabel 'Nombre de subdivisions cylindre'      # nom de l'axe des abscisses
set ylabel 'Temps de génération (µs)'             # nom de l'axe des ordonnées
set y2label 'Nombre de triangles'             # nom de l'axe des ordonnées
set y2tics nomirror

plot 'benchmarkDurationCylinder.dat' axis x1y1 with lines title "Temps génération", 'benchmarkTriangleCylinder.dat' axis x1y2 with lines title "Nb triangles"
