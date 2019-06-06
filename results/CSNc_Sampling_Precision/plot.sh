#!/usr/bin/bash

#############################################
########  Return sampling conv test  ########
######## Written by Ilan E. Chemmama ########
########   Andrej Sali Laboratory    ########
########    UC - San Francisco       ########
#############################################

for dir in BMS DHS DSS BMS_DHS DSS_BMS DSS_DHS DSS_DHS_BMS; 
do
    cd ${dir}
    rm *.pdf *.png 
    gnuplot -c ../Plot_Cluster_Population.plt Cop9_$dir
    gnuplot -c ../Plot_Convergence_NM.plt  Cop9_$dir
    gnuplot -c ../Plot_Convergence_SD.plt Cop9_$dir
    gnuplot -c ../Plot_Convergence_TS.plt Cop9_$dir
    cd ..
done

