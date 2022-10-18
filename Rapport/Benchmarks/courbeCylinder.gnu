# Définition du format de sortie
set terminal pngcairo

# Définition du nom du fichier de sortie
set output 'courbeCylinder.png'

# Définition des titres de la figure
set title 'Temps de génération / Nombre subdivisions' # titre
set xlabel 'Nombre de subdivisions cylindre'      # nom de l'axe des abscisses
set ylabel 'Temps de génération (ms)'             # nom de l'axe des ordonnées

plot 'benchmarkCylinder.dat' with lines title "Cylinder"
