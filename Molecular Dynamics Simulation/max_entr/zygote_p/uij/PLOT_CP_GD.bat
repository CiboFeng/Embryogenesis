set terminal postscript enhanced color font "Arial-bold, 25"
set output "CP_GD.eps"

set xlabel "Genomic Distance (Mb)"
set ylabel "Contact Probability"

set logscale

set border lw 2

set xrange [0.1:]
set yrange [0.0001:]

set format x "10^{%T}"
set format y "10^{%T}"

set label "s^{-1.0}" at graph 0.35,0.55

set size 0.9,1.0

p "CP_GD.dat" u 1:2 w l lt -1 lw 6 notitle,\
  [0.5:7] 1/x**1.0*0.08 w l lt 18 lw 8 notitle

!ps2pdf CP_GD.eps && rm CP_GD.eps
