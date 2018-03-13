#!/bin/bash
#$ -l arch=linux-x64 
#$ -l h_rt=20:0:0
#$ -l mem_free=2G
#$ -S /bin/bash 
#$ -N pref
#$ -cwd 
#$ -j y 
#$ -R y
#$ -t 1-12

source ~/.bash_profile
i=$(($SGE_TASK_ID - 1))

while read p;
do
    
    tp=$(echo $p | awk 'BEGIN{FS="|"}{print $2}')
    
    id=$(echo $p | awk 'BEGIN{FS="|"}{print $3}' | awk 'BEGIN{FS="."}{print $1}')

    p1=$(echo $p | awk 'BEGIN{FS="|"}{print $4}')
    r1=$(echo $p | awk 'BEGIN{FS="|"}{print $5}')

    p2=$(echo $p | awk 'BEGIN{FS="|"}{print $6}')
    r2=$(echo $p | awk 'BEGIN{FS="|"}{print $7}')

    grep ${p} XLS/Stats/*.txt > XLS/${tp}_${p1}_${r1}_${p2}_${r2}_${id}.txt

done < XLS/List${i}
