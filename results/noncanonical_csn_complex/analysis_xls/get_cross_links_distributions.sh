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
