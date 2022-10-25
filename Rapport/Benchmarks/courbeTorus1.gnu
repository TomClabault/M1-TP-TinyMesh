# Définition du format de sortie
set terminal pngcairo

# Définition du nom du fichier de sortie
set output 'courbeTorus1.png'

# Définition des titres de la figure
set title 'Temps de génération / Nombre subdivisions' # titre
set xlabel 'Nombre de ring subdivisions'      # nom de l'axe des abscisses
set ylabel 'Temps de génération (ms)'         # nom de l'axe des ordonnées
set y2label 'Nombre de triangles'
set y2tics nomirror

plot 'benchmarkDurationTorus1.dat' axis x1y1 with lines title "Temps de génération", 'benchmarkTriangleTorus1.dat' axis x1y2 with lines title "Nb triangles"
