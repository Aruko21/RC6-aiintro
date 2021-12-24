reset
set term wxt

color_set1 = '#b54f1e'
color_set2 = '#565af5'

color_line = '#3e3c3d'

set border linewidth 1.5
set pointsize 1

# стиль для тестовых данных с классом 1 (класс 0 отображается дефолтным стилем)
set style line 1 lc rgb '#808080' pt 7   # circle gray

# стиль для весов
set style line 2 lc rgb '#fc0362' pt 9   # triangle red

set hidden3d
set dgrid3d 20,20 qnorm 2

unset key

unset colorbox

set tics scale 0.1
set xtics 1
set ytics 1
set yrange[-4:4]
set xrange[-4:4]
set zrange[-2:2]
set xlabel 'x'
set ylabel 'y'

DOTS = system("wc -l init_plot.dat")
N = system("wc -l plots.dat") / (DOTS + 1)

do for [i=0:N:100] {
    LABEL = "i = " . (i+1)
    set obj 10 rect at graph 0.1,0.1 size char strlen(LABEL), char 1
    set obj 10 fillstyle empty border -1 front
    set label 10 at graph 0.1,0.1 LABEL front center

    splot 'init_plot.dat' w l ls 1, 'plots.dat' every :::i::i w l ls 2
    pause 0.01
}
pause mouse
