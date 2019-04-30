#!/usr/bin/bash

#for f in `ls -1v ???_*.csv`; 
#do 
#    python get_histo.py $f $(echo $f | sed 's/_cluster_1.csv/_H.txt/g'); 
#done

gnuplot -c plot_histo.plt BMS_Intra_H.txt 45.0 BMS_Intra_Dist.png C63F05
gnuplot -c plot_histo.plt BMS_Inter_H.txt 45.0 BMS_Inter_Dist.png C63F05

gnuplot -c plot_histo.plt DHS_Intra_H.txt 30.0 DHS_Intra_Dist.png 116325
gnuplot -c plot_histo.plt DHS_Inter_H.txt 30.0 DHS_Inter_Dist.png 116325

gnuplot -c plot_histo.plt DSS_Intra_H.txt 30.0 DSS_Intra_Dist.png 2828BE
gnuplot -c plot_histo.plt DSS_Inter_H.txt 30.0 DSS_Inter_Dist.png 2828BE

montage DSS_Intra_Dist.png DHS_Intra_Dist.png BMS_Intra_Dist.png DSS_Inter_Dist.png DHS_Inter_Dist.png BMS_Inter_Dist.png -tile 3x2 -geometry +0+0 -background none separated_csn1-9.png