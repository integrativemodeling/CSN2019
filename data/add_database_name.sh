#!/usr/bin/bash

for file in *.csv
do
    awk 'BEGIN{print "prot1,prot2,res1,res2,id"}{print}' $file > temp0
    mv temp0 $file
done
