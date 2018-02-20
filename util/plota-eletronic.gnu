
set key font "Helvetica,15"

set xrange [-1.2:1.2]
set xtics font "Helvetica,15"

set yrange [-1.2:1.2]
set ytics font "Helvetica,15"

set style line 1 lt 1 lc rgb "black" lw 5

plot "yTrajectory.csv" using 1:2 title "angle = 20" with linespoints ls 1