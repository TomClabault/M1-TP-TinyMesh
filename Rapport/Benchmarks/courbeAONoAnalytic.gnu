# Définition du format de sortie
set terminal pngcairo

# Définition du nom du fichier de sortie
set output 'courbeAONoAnalytic.png'

# Définition des titres de la figure
set title 'Temps de calcul AO (mesh intersection) \& Nb triangles / Nombre subdivisions' # titre
set xlabel 'Nombre de subdivisions icosphere'      # nom de l'axe des abscisses
set ylabel 'Temps de calcul AO (ms)'             # nom de l'axe des ordonnées
set y2label 'Nombre de triangles'             # nom de l'axe des ordonnées
set y2tics nomirror

plot 'benchmarkAOnoAnalyticIntersection.dat' axis x1y1 with lines title "Temps de calcul AO", 'benchmarkAOnoAnalyticIntersection.datTriangles' axis x1y2 with lines title "Nb triangles"
