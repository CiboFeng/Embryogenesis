set terminal postscript enhanced color font "Arial-bold, 25" solid
set output "TAD1.eps"

set xlabel "Genomic Distance (Mb)"
set ylabel "TAD Signal"

set border lw 2

set xrange [0.0:40.0]
set yrange [-2:2]

set ytics 2

set mxtics 2
set mytics 2

set size 1.0,0.6

p "TAD.dat" u ($1/10):2 w l lt -1 lw 3 notitle, 0 w l lw 3 lc rgb "red" notitle

!ps2pdf TAD1.eps && rm TAD1.eps

set output "TAD2.eps"
set xrange [40.0:85.7]

p "TAD.dat" u ($1/10):2 w l lt -1 lw 3 notitle, 0 w l lw 3 lc rgb "red" notitle

!ps2pdf TAD2.eps && rm TAD2.eps
