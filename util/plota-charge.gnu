
set xtics font "Helvetica,6"

set boxwidth 0.5
set style fill solid
set grid
set style line 1 lc rgb "blue"
plot "electronPlot.csv" using 2:xtic(1) with boxes ls 1 notitle

