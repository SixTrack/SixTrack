set term png size 1200,800 notransparent enhanced
set output 'pencilbeam_distr_type3.png'
set grid

set multiplot layout 2,2 rowsfirst

set xlabel 'x [mm]'
set ylabel 'xp [mm]'
plot \
     'pencilbeam_distr_type3.dat' index 0 using ($7+$9*0.6):9 with points notitle

set xlabel 'sigmv [mm]'
set ylabel 'E [MeV]'
plot \
     'pencilbeam_distr_type3.dat' index 0 using 11:14 with points notitle

set xlabel 'y [mm]'
set ylabel 'yp [mm]'
set xrange [-1.01:-0.99]
plot \
     'pencilbeam_distr_type3.dat' index 0 using ($8+$10*0.6):10 with points notitle

set xlabel 'y [mm]'
set ylabel 'yp [mm]'
set xrange [0.99:1.01]
plot \
     'pencilbeam_distr_type3.dat' index 0 using ($8+$10*0.6):10 with points notitle

unset multiplot