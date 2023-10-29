#!/usr/bin/gnuplot
set terminal png
set key autotitle columnhead
set datafile separator ","

set output "plot_2A.png"
set xlabel "p"
set logscale x 10
set grid
plot "data_min.csv" using 1:2 with lines lw 1

set output "plot_2B.png"
set xlabel "p"
set logscale x 10
set grid
plot "data_min.csv" using 1:3 with lines lw 1