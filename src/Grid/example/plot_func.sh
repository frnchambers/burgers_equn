#!/bin/bash

if [[ $# -eq 0 ]] ; then
    gnuplot -e "MESH='Sinh'" plot_grid-func.gp
    gnuplot -e "MESH='Cheb_1'" plot_grid-func.gp
    gnuplot -e "MESH='Cheb_2'" plot_grid-func.gp
else
    gnuplot -e "MESH='$1'" plot_grid-func.gp
fi


