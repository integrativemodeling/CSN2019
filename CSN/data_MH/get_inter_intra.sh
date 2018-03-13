#!/usr/bin/bash

awk 'BEGIN{FS=OFS=","}{if(NR==1) print $0; else{if($1!=$2) print $0;}}' Cop9.BMS.csv > Cop9.BMS.Inter.csv
awk 'BEGIN{FS=OFS=","}{if(NR==1) print $0; else{if($1==$2) print $0;}}' Cop9.BMS.csv > Cop9.BMS.Intra.csv
awk 'BEGIN{FS=OFS=","}{if(NR==1) print $0; else{if($1==$2) print $0;}}' Cop9.DSS.csv > Cop9.DSS.Intra.csv
awk 'BEGIN{FS=OFS=","}{if(NR==1) print $0; else{if($1!=$2) print $0;}}' Cop9.DSS.csv > Cop9.DSS.Inter.csv
awk 'BEGIN{FS=OFS=","}{if(NR==1) print $0; else{if($1==$2) print $0;}}' Cop9.DHS.csv > Cop9.DHS.Intra.csv
awk 'BEGIN{FS=OFS=","}{if(NR==1) print $0; else{if($1!=$2) print $0;}}' Cop9.DHS.csv > Cop9.DHS.Inter.csv
