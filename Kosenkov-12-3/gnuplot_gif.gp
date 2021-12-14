reset
set term gif animate delay 50
set output 'result.gif'

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

N = system("wc -l weights.dat")

do for [i=1:N] {
    w_0 = system("sed -n " . i . "p weights.dat | awk '{print $1}'")
    w_1 = system("sed -n " . i . "p weights.dat | awk '{print $2}'")
    w_2 = system("sed -n " . i . "p weights.dat | awk '{print $3}'")

    LABEL = "i = " . i
    set obj 10 rect at graph 0.1,0.1 size char strlen(LABEL), char 1
    set obj 10 fillstyle empty border -1 front
    set label 10 at graph 0.1,0.1 LABEL front center

    f(x) = -1/(w_2) * (w_1 * x + w_0)

    plot f(x) w l ls 1, 'learn_data.dat' using 1:2:3 w p pt 5 lc palette z, 'test_data.dat' using 1:2:3 w p pt 3 lc palette z
}
