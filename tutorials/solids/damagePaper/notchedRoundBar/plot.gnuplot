set term pngcairo dashed size 1024,768 font "Arial,12"
set datafile separator " "

set style line 1 linecolor rgb 'black' linetype 6 linewidth 1 ps 1
set style line 2 linecolor rgb 'black' linetype 4 linewidth 1 ps 1
set style line 3 linecolor rgb 'red' linetype 6 linewidth 1 ps 1

set grid
#set xrange [0:50]
set yrange [:36]
set xtics
set ytics
#set ytics 0.002
set output "forceVsDisp.png"
set xlabel "Displacement (in mm)"
set ylabel "Force (in kN)"
set key right bottom;

# Wedge is 5 degrees so multiply force by (360/5) = 72
plot "postProcessing/0/solidForcesup.dat" u ($1*2):($3*72/1e3) w lp ls 1 t "This work"
