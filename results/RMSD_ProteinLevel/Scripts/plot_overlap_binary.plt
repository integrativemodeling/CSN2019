reset

set terminal pngcairo transparent enhanced color font "ArialNarrow-Bold, 39" size 1028,1028
set border lw 5 lc rgb "#484848"

set encoding iso_8859_1
#set palette defined (0.2 "#FEF0D9", 0.4 "#FDCC8A", 0.6 "#FC8D59", 0.8 "#E34A33", 1.0 "#B30000")
set palette defined (0.0 "#000000", 1.0 "#ff0000")
#set palette rgb -15,-5,-7

#set cbrange [0.0:3.0]
#set zrange [0.0:3.0]

set xrange [-0.5:7.5]
set yrange [-0.5:7.5]

set xlabel " " tc rgb "#484848" 
set ylabel "Subunit" tc rgb "#484848" 

#unset xlabel
unset ylabel

set xtics ("CSN1" 0, "CSN2" 1, "CSN3" 2, "CSN4" 3, "CSN5" 4, "CSN6" 5, "CSN7" 6, "CSN8" 7) rotate by 90 center offset 0, -1.0
set ytics ("CSN1" 0, "CSN2" 1, "CSN3" 2, "CSN4" 3, "CSN5" 4, "CSN6" 5, "CSN7" 6, "CSN8" 7) center offset -1.80, 0

unset key
unset colorbox

set output (sprintf("%s", ARG2))
plot sprintf("%s", ARG1) matrix w image
set output




