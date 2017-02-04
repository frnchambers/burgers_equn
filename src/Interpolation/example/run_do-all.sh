#!/bin/bash

set -e


## ---- compile the things ---- ##
make run_lagrange run_spline run_GSL-spline

## ---- run the things ---- ##
echo "lagrange:"
./run_lagrange $1
echo "spline:"
./run_spline $1
echo "GSL-spline:"
./run_GSL-spline $1

## ---- plot the things ---- ##
gnuplot -e "METHOD='lagrange'" plot_interp.gp
gnuplot -e "METHOD='spline'" plot_interp.gp
gnuplot -e "METHOD='GSL-spline'" plot_interp.gp


