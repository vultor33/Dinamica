
set key font "Helvetica,15"

set xrange [0:305]
set xtics font "Helvetica,15"

set yrange [0.7:3]
set ytics font "Helvetica,15"

set style line 1 lt 1 lc rgb "black" lw 5

plot "rppGraph.csv" using 1:2 title "angle = 5" with linespoints ls 1