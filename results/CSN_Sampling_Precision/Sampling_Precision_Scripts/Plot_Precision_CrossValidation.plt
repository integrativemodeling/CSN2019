#############################################
######## Written by Ilan E. Chemmama ########
########   Andrej Sali Laboratory    ########
########    UC - San Francisco       ########
#############################################

reset

set terminal pngcairo transparent enhanced color font "ArialNarrow-Bold, 39" size 1028,1028
set border lw 5 lc rgb "#484848"

set encoding iso_8859_1
set palette defined (0.0 "#B30000", 0.25 "#E34A33",  0.50 "#FC8D59", 0.75 "#FDCC8A", 1.0 "#FFFFFF")

#set palette rgb -15,-5,-7

set cbrange [10.0:60.0]
#set zrange [0.0:60.0]

set xrange [0:8]
set yrange [13:40]

set xlabel "" tc rgb "#484848" 
set ylabel "Precision (\305)" tc rgb "#484848" 

unset xlabel
unset ylabel
unset xtics
unset ytics 

unset key
#unset colorbox

set tics scale 0

set size 1.0, 0.8
set output (sprintf("%s", ARG2))
plot sprintf("%s", ARG1) usi 1:3:3 w p palette pt 7 ps 2 notitle
set output




