#!/usr/bin/bash

gnuplot -c Plot_Cluster_Population.plt $1
gnuplot -c Plot_Convergence_NM.plt  $1 
gnuplot -c Plot_Convergence_SD.plt $1 
gnuplot -c Plot_Convergence_TS.plt $1

