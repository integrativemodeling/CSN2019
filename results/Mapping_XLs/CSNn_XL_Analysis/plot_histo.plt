#############################################
########   Plot score distribution   ########
######## Written by Ilan E. Chemmama ########
########   Andrej Sali Laboratory    ########
########    UC - San Francisco       ########
#############################################


### To run the script
### gnuplot -c plot_histo.plt <input file> <output file with png extension>


reset
set terminal pngcairo transparent enhanced color font "Arial-Bold, 48" size 1028, 1028
set border lw 5 lc rgb "#484848"

set encoding iso_8859_1

set xlabel "Crosslink distance ({\305})" tc rgb "#484848"
set ylabel "Number of Crosslinks" tc rgb "#484848" 

stats sprintf("%s", ARG1) usi 2 prefix "u"
set xr[0:100]
set yr[0:1.05]

set format y "%.2f"
set format x "%.1f"

max = u_max
set xtics 0, 25, 100 nomirror
unset ytics

set linetype 5 dashtype 2 lw 5
set arrow nohead from sprintf("%s", ARG2),0 to sprintf("%s", ARG2),1.05 lw 5 lt 5 lc rgb sprintf("#%s", ARG4) back filled

set output sprintf("%s", ARG3)
plot sprintf("%s", ARG1) usi 1:($2/max) w histeps lw 5 lc rgb sprintf("#%s", ARG4) notitle
set output

