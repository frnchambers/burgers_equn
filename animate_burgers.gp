#!/usr/bin/gnuplot

reset

datafile='data/solution.dat'


# png
set terminal pngcairo size 400,400 enhanced font 'Verdana,10'

system('mkdir png')

#set yrange[-1:1]

stats datafile skip 1 nooutput
set xrange [STATS_min_x:STATS_max_x]

n=0
do for [i=1:(STATS_blocks-2)] {
    n=i+1
    set output sprintf('png/burgers%03.0f.png',n)
    set key autotitle columnheader
    plot datafile index i u 1:2 w linespoints title ''
}

system('convert -delay 10 -loop 0 png/*.png burgers.gif')
system('rm -rf png/')
system('ristretto burgers.gif')

