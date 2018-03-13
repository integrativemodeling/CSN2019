#!/bin/bash
#$ -l arch=linux-x64 
#$ -l h_rt=300:0:0
#$ -l mem_free=2G
#$ -S /bin/bash 
#$ -N pref
#$ -cwd 
#$ -j y 
#$ -R y
#$ -t 1-43

source ~/.bash_profile
i=$(($SGE_TASK_ID - 1))

while read p;
do

    dir=$(echo $p | awk '{print $1}')
    rep=$(echo $p | awk '{print $2}')
    frm=$(echo $p | awk '{print $4}')
    frn=$(echo $p | awk '{print $4 + 1}')
    
    process_output.py -f ../modeling${dir}/output/stat.${rep}.out -n ${frn} > XLS/${dir}_${rep}_${frm}.txt

done < XLS/GSM_${i}
