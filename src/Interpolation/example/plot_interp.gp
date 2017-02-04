#!/usr/bin/gnuplot

reset

set terminal qt size 900,900 enhanced persist

# set terminal postscript eps size 24cm,18cm enhanced color font 'Ubuntu,15' linewidth 1
# set output 'noburst_T-1x8.eps'

# set terminal dumb


# set lmargin 10
# set rmargin 10

# set bmargin 4
# set tmargin 1

set grid
set key right top

method=INTERP
mesh=MESH

pointsfile='data/points_'.mesh.'.dat'
linefile='data/line_'.mesh.'_'.method.'.dat'
derivs1='data/derivs1_'.mesh.'_'.method.'.dat'
derivs2='data/derivs2_'.mesh.'_'.method.'.dat'

set multiplot layout 3,2 rowsfirst title 'Analysis of '.method.' Interpolation'
set title ''

#set xrange[-1.1:1.1]
set xlabel 'x'


## ----------- First plot ----------- ##
unset log y
set format y '%.1f'
set ylabel 'f(x)'
plot linefile   u 1:2 w lines title 'interpolation', \
     pointsfile u 1:2 w points title 'data points'

set log y
set format y '10^{%L}'
set ylabel '{/Symbol D} f(x)'
plot linefile u 1:(abs($4)) w lines title ''


# ## ---------- Second plot ---------- ##
# unset log y
# set format y '%.1f'
# set ylabel 'df(x)/dx'
# plot derivs1 u 1:2 w points title 'interpolation points', \
#      derivs1 u 1:3 w lines title 'true poines'

# set log y
# set format y '10^{%L}'
# set ylabel '{/Symbol D} df(x)/dx'
# plot derivs1 u 1:(abs($4)) w lines title ''


# ## ---------- Third plot ---------- ##
# unset log y
# set format y '%.1f'
# set ylabel 'd^2f(x)/dx^2'
# plot derivs2 u 1:2 w points title 'interpolation points', \
#      derivs2 u 1:3 w lines title 'true poines'

# set log y
# set format y '10^{%L}'
# set ylabel '{/Symbol D} d^2f(x)/dx^2'
# plot derivs2 u 1:(abs($4)) w lines title ''



