set terminal postscript enhanced color font "Arial-bold, 25" solid
set output "CP_GD.eps"

set xlabel "Genomic Distance (Mb)"
set ylabel "Contact Probability"

load "/hpc2hdd/home/chu-amat/cbfengphy/gnuplot-palettes/parula.pal"
set logscale

set border lw 2

set xrange [0.1:57.7]
set yrange [0.0001:1]
set format x "10^{%T}"
set format y "10^{%T}"

#set label "s^{-1.0}" at graph 0.35,0.55

set size 0.9,1.0

set key bottom left
set key samplen -1
set key spacing 1.5
set key textcolor variable

p "CP_GD_HiC.dat"    u 1:2 w l lt -1 lc rgb "black"   lw 4 title "Hi-C",\
  "CP_GD_Sim.dat"    u 1:2 w l ls 15 lw 6 title "Sim",\

#  [0.5:7] 1/x**1.0*0.08 w l lt 18 lw 8 notitle

!ps2pdf CP_GD.eps && rm CP_GD.eps

set output "TAD1.eps"
unset logscale
set ytics format "%.0f"
set xtics format "%.0f"
set yrange [-3.5:1.5]
set ytics 1
set size 1.0,0.5
set ylabel "Insulation Score"
set xrange [0:37.8]
set xtics 10

p "InsulationScore_HiC.dat"    u ($1/10):2 w l lt -1 lc rgb "black" lw 4 notitle,\
  "InsulationScore_Sim.dat"    u ($1/10):2 w l ls 15 lw 6 notitle

!ps2pdf TAD1.eps && rm TAD1.eps

set output "TAD2.eps"
set yrange [-3.5:1.5]
set ytics 1
set size 1.0,0.5
set xrange [20:57.8]
set xtics 10

p "InsulationScore_HiC.dat"    u ($1/10):2 w l lt -1 lc rgb "black" lw 4 notitle,\
  "InsulationScore_Sim.dat"    u ($1/10):2 w l ls 15 lw 6 notitle

!ps2pdf TAD2.eps && rm TAD2.eps

set output "CM.eps"
set yrange [-0.15:0.15]
set ytics format "%.1f"
set ytics 0.1
set size 1.0,0.5
set ylabel "Compartment Profiles"
set xrange [1:57.7]
set xtics 20

p "<cat -n CM_HiC.dat" u ($1/1):2 w l lt -1 lc rgb "black" lw 4 notitle,\
  "<cat -n CM_Sim.dat" u ($1/1):2 w l ls 15 lw 6 notitle

!ps2pdf CM.eps && rm CM.eps
