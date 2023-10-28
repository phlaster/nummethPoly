#!/usr/bin/gnuplot
set terminal png
set key autotitle columnhead
set datafile separator ","

set output "plot_1A.png"
set xlabel "x"
set grid
plot "data_baza_a.csv" using 1:2 with lines lw 1,\
    "data_baza_b.csv" using 1:2 with lines lw 1,\
    "data_baza_b.csv" using 1:3 with lines lw 1 lc "black",

set output "plot_1B.png"
set xlabel "x"
set grid
plot "data_baza_a.csv" using 1:2 with lines lw 1,\
    "data_baza_a.csv" using 1:3 with lines lw 1

set output "plot_1C.png"
set xlabel "x"
set grid
plot "data_baza_b.csv" using 1:2 with lines lw 1,\
    "data_baza_b.csv" using 1:3 with lines lw 1

set output "plot_1D.png"
set xlabel "x"
set grid
plot "data_baza_a.csv" using 1:4 with lines lw 1,\
    "data_baza_b.csv" using 1:4 with lines lw 1,\

set output "plot_1E.png"
set xlabel "h"
set logscale x 10
set logscale y 10
set grid
plot "data_baza_c.csv" using 1:3 with lines lw 1,\
    "data_baza_c.csv" using 1:2 with lines lw 1,\