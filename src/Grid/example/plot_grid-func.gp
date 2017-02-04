#!/usr/bin/gnuplot

reset

set terminal qt size 600,1000 enhanced persist

set grid
set key right top

meshtype=MESH
datafile='data/gridpoints_'.meshtype.'.dat'

set multiplot layout 2,1 title meshtype

set xlabel ''
set ylabel 'x=g(z)'
plot datafile u 1:2 w linespoints title '', \
     datafile u (0):2 w points title ''

set xlabel 'i'
set ylabel 'w(i)'
plot datafile u 1:3 w linespoints title ''



# #### choose parameters
# m1=0.5
# m2=0.7
# s=1.0

# #### solved for these
# a=15.0/8.0  - m1*7.0/8.0  - m2     + s/8.0
# b=-21.0/4.0 + m1*9.0/4.0  + m2*3.0 - s/4.0
# c=35.0/8.0  - m1*11.0/8.0 - m2*3.0 + s/8.0
# d=m2

# g(x) = (a*(x**7)) + (b*(x**5)) + (c*(x**3)) + (d*x)
# gp(x) = 7*a*(x**6) + 5*b*(x**4) + 3*c*(x**2) + d

# plot g(x)





