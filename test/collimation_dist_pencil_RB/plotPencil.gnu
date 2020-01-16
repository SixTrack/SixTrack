set term png size 1200,800 notransparent enhanced
set output 'pencilbeam_distr_RB.png'
set grid

set multiplot layout 2,2 rowsfirst

set xlabel 'x [mm]'
set ylabel 'xp [mm]'
plot \
     'pencilbeam_distr_RB.dat' index 0 using ($1+$2*0.6):2 with points notitle

set xlabel 'sigmv [mm]'
set ylabel 'E [MeV]'
plot \
     'pencilbeam_distr_RB.dat' index 0 using 5:6 with points notitle

set xlabel 'y [mm]'
set ylabel 'yp [mm]'
set xrange [-1.01:-0.99]
plot \
     'pencilbeam_distr_RB.dat' index 0 using ($3+$4*0.6):4 with points notitle

set xlabel 'y [mm]'
set ylabel 'yp [mm]'
set xrange [0.99:1.01]
plot \
     'pencilbeam_distr_RB.dat' index 0 using ($3+$4*0.6):4 with points notitle

unset multiplot