#!/usr/bin/bash


for file in `ls ../run_analysis/XLs_distances_*_Int*_cluster_1.csv`;
do
    echo $file
    grep -v "Cross" $file > tmp
    awk 'BEGIN{FS=","}{for(i=2; i<NF; i++) print $i}' tmp > $(echo $file | sed 's/XLs_distances_/ /g' | awk '{print $2}' )
    rm tmp
done

cat BMS_* > BMS.txt
python get_histo.py BMS.txt BMS_H.txt

cat DHS_* > DHS.txt
python get_histo.py DHS.txt DHS_H.txt 

cat DSS_* > DSS.txt
python get_histo.py DSS.txt DSS_H.txt 


gnuplot -c plot_histo.plt DSS_H.txt 30.0 DSS_Dist.png 2828BE
gnuplot -c plot_histo.plt DHS_H.txt 30.0 DHS_Dist.png 116325
gnuplot -c plot_histo.plt BMS_H.txt 45.0 BMS_Dist.png C63F05

gnuplot -c plot_histo_together.plt DSS_H.txt DHS_H.txt BMS_H.txt