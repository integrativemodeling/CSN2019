reset
set terminal pngcairo transparent enhanced color font "Arial-Bold, 48" size 1028, 1028
set border lw 5 lc rgb "#484848"

set encoding iso_8859_1

set xlabel "Crosslink distance ({\305})" tc rgb "#484848"
set ylabel "Number of Crosslinks" tc rgb "#484848" 

stats sprintf("%s", ARG1) usi 2 prefix "dss"
stats sprintf("%s", ARG2) usi 2 prefix "dhs"
stats sprintf("%s", ARG3) usi 2 prefix "bms"

set xr[0:100]
set yr[0:1.05]
set format y "%.2f"
set format x "%.1f"

set xtics 0, 25, 100 nomirror
unset ytics

set linetype 5 dashtype 4 lw 5
set linetype 6 dashtype 2 lw 5

set arrow nohead from 30.0,0 to 30.0,1.05 lw 5 lt 5 lc rgb "#2828BE" back filled
set arrow nohead from 30.0,0 to 30.0,1.05 lw 5 lt 6 lc rgb "#116325" back filled
set arrow nohead from 45.0,0 to 45.0,1.05 lw 5 lt 6 lc rgb "#2828BE" back filled

set output "Together.png"
set key font "Arial-Bold, 30" tc rgb "#484848"


plot sprintf("%s", ARG1) usi 1:($2 / dss_max) w histeps lc rgb "#2828BE" lw 5 title "DSSO", \
     sprintf("%s", ARG2) usi 1:($2 / dhs_max) w histeps lc rgb "#116325" lw 5 title "DHSO", \
     sprintf("%s", ARG3) usi 1:($2 / bms_max) w histeps lc rgb "#C63F05" lw 5 title "BMSO"

set output

