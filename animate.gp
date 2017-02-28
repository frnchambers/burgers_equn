#!/usr/bin/gnuplot

reset

datafile='data/burgers.dat'


# png
set terminal pngcairo size 400,400 enhanced font 'Verdana,10'

system('mkdir png')

set yrange[-1:1]

stats datafile skip 1

n=0
do for [i=1:50] {
    n=n+1
    set output sprintf('png/burgers%03.0f.png',n)
    plot datafile index i u 1:2 w linespoints title columnheader(1)
}

system('convert -delay 10 -loop 0 png/*.png burgers.gif')
system('rm -rf png/')
