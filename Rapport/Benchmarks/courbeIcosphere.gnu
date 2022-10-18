# Définition du format de sortie
set terminal pngcairo

# Définition du nom du fichier de sortie
set output 'courbeIcosphere.png'

# Définition des titres de la figure
set title 'Temps de génération / Nombre subdivisions' # titre
set xlabel 'Nombre de subdivisions'                    # nom de l'axe des abscisses
set ylabel 'Temps de génération log4(ms)'                    # nom de l'axe des ordonnées
set logscale y 4

plot 'benchmarkIcosphere.dat' with linespoints title "Icosphere"
