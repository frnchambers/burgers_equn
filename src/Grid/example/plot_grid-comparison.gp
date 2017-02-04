#!/usr/bin/gnuplot

reset

set grid
set key left top

#### choose parameters
## good choice here
m1=0.0
m2=0.7
s=0.0
# m1=0.0
# m2=1.0
# s=0.0 # cannot be <0 otherwise multivalued

#### solved linear system for these
a=15.0/8.0  - m1*7.0/8.0  - m2     + s/8.0
b=-21.0/4.0 + m1*9.0/4.0  + m2*3.0 - s/4.0
c=35.0/8.0  - m1*11.0/8.0 - m2*3.0 + s/8.0
d=m2

T7z(x) = (a*(x**7)) + (b*(x**5)) + (c*(x**3)) + (d*x)
shiny(x) = sinh(x*2)/sinh(2)
cheb_1(x) = -cos( pi / 2.0 * ( x+1.0 ) )

zm=1.8

f(t) = (1.0-zm)*t*t + t


set terminal qt 0 size 400,400 enhanced persist


set xrange[-1:1]
set yrange[-1:1]
# plot T7z(x)


set terminal qt 1 size 1000,400 enhanced persist
set multiplot layout 1,3

set xlabel 'z'
set ylabel 'x=g(z)'

set xrange[-1:1]
set yrange[-1:1]
set title 'Full Range'
plot T7z(x), shiny(x), cheb_1(x) title 'cheb_1(x)', x

unset key
set ylabel ''

set xrange[-1:-0.5]
set yrange[-1:-0.5]
set title 'LHS Boundary'
plot T7z(x), shiny(x), cheb_1(x) title 'cheb_1(x)', x

set xrange[-0.25:0.25]
set yrange[-0.25:0.25]
set title 'Center'
plot T7z(x), shiny(x), cheb_1(x) title 'cheb_1(x)', x

