#!/bin/bash

gnuplot << EOF
set terminal pngcairo color enhanced font 'Verdana,12' size 1440, 1024
set output 'selection-logscale.png'
set title 'Comparison of selection algorithms.'
set border 3
set key top left
set xtics nomirror
set ytics nomirror
set grid ytics
set xlabel 'selectivity of the predicate (in %)'
set ylabel 'overall time (in ms)'
set xrange [0.01:*]
set yrange [0:*]
set logscale x 1.2
set format x '%5.2f'
plot 'selection.dat' using 1:2 with lines title 'branching', \
     ''              using 1:3 with lines title 'branch free (mult)', \
     ''              using 1:4 with lines title 'branch free (select)', \
     ''              using 1:5 with lines title 'branch free (predicated)'
EOF

gnuplot << EOF
set terminal pngcairo color enhanced font 'Verdana,12' size 1440, 1024
set output 'selection-linear.png'
set title 'Comparison of selection algorithms.'
set border 3
set key top left
set xtics nomirror
set ytics nomirror
set grid ytics
set xlabel 'selectivity of the predicate (in %)'
set ylabel 'overall time (in ms)'
set xrange [0.01:*]
set yrange [0:*]
set format x '%5.2f'
plot 'selection.dat' using 1:2 with lines title 'branching', \
     ''              using 1:3 with lines title 'branch free (mult)', \
     ''              using 1:4 with lines title 'branch free (select)', \
     ''              using 1:5 with lines title 'branch free (predicated)'
EOF
