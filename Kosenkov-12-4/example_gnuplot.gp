reset
set term x11
set size square

set border linewidth 1.5
set pointsize 1.7

# стиль для тестовых данных с классом 1 (класс 0 отображается дефолтным стилем)
set style line 1 lc rgb '#808080' pt 7   # circle gray

# стиль для весов
set style line 2 lc rgb '#fc0362' pt 9   # triangle red

unset key

set tics scale 0.1
set xtics 1
set ytics 1
set yrange[-3:4]
set xrange[-3:3]
set xlabel 'X_1'
set ylabel 'X_2'

set grid
set macro

set style increment user

N = system("wc -l coordinates.dat") * 5

do for [i=0:N] {
    # номер итерации обучения
    LABEL = "i = " . (i+1)
    set obj 10 rect at -1.5,-1.5 size char strlen(LABEL), char 1
    set obj 10 fillstyle empty border -1 front
    set label 10 at -1.5,-1.5 LABEL front center


    plot 'coordinates.dat' u 1:2 w p ls 1, 'weights.dat' every :::i::i using 2:3 w p ls 2, 'weights.dat' every :::i::i using 2:3:1 with labels notitle

    pause 0.1
}
pause mouse