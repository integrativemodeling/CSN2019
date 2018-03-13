reset
set terminal pdfcairo enhanced color font "Arial, 48" size 20,10
set border lw 10 lc rgb "#484848"
set output "Example.pdf"
unset key 


set bmargin at screen 0.10
set tmargin at screen 0.98

RSIDE=0.98
DX=0.45

set tics scale 0.0
set multiplot

set offset 0.0, 0.0, graph 0.05, graph 0.05

set ylabel "Percentage Satisfied (%)" tc rgb "#484848"
set xlabel "Crosslink Pairs" tc rgb "#484848" 

set yr[0:101]
set xr[-1:142]
set format y "%.1f"
set xtics font "Arial, 9" rotate by 90 offset 0.0, 0.0 right

set rmargin at screen RSIDE-DX 
set lmargin at screen RSIDE-2*DX + 0.025
plot "Inter_XLS.txt" usi 1:3:xtic(2) w p pt 7 ps 2 lc rgb "blue"


unset ylabel
unset xlabel
unset ytic
set yr[0:101]
set xr[-1:172]
set rmargin at screen RSIDE
set lmargin at screen RSIDE-DX  + 0.0
set xtics font "Arial, 9" rotate by 90 offset 0.0, 0.0 right

plot "Intra_XLS.txt" usi 1:3:xtic(2) w p pt 7 ps 2 lc rgb "red"

unset multiplot 

