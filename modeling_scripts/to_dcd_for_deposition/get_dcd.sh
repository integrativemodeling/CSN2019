#!/bin/bash

for dir in DSS DHS BMS DSS_DHS DSS_BMS BMS_DHS;
do
    python to_dcd.py /scratch/Cop9/Cop9_${dir}/prefilter/Cluster_1/sample_A/ A_${dir}.txt A_${dir}.dcd
    python to_dcd.py /scratch/Cop9/Cop9_${dir}/prefilter/Cluster_1/sample_B/ B_${dir}.txt B_${dir}.dcd

done
