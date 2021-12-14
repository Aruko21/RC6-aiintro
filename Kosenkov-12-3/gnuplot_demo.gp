reset
set term wxt

color_set1 = '#b54f1e'
color_set2 = '#565af5'

color_line = '#3e3c3d'

set border linewidth 1.5
set pointsize 1
set style line 1 lc rgb color_line pt 7

set palette defined (0 color_set1, 1 color_set2)

set grid

unset key

unset colorbox

set tics scale 0.1
set xtics 1
set ytics 1
set yrange[-1:9]
set xrange[-4:8]
set xlabel 'X_1'
set ylabel 'X_2'

plot 'learn_data.dat' using 1:2:3 w p pt 5 lc palette z

pause -1