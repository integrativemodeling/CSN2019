#!/usr/bin/bash

for xlset in DSS DHS BMS; 
do 
    cat ../run_analysis/XLs_satisfaction_${xlset}_Intra_cluster_1.csv | sed 's/|/,/g' |  awk 'BEGIN{FS=OFS=","}{print $4, $5, $6, $7,$11, $12, $13, $14, $15}' > Cop9.${xlset}.Intra.csv

    cat ../run_analysis/XLs_satisfaction_${xlset}_Inter_cluster_1.csv | sed 's/|/,/g' |  awk 'BEGIN{FS=OFS=","}{print $4, $5, $6, $7,$11, $12, $13, $14, $15}' > Cop9.${xlset}.Inter.csv

done 