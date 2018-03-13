#!/bin/bash
#$ -l arch=linux-x64 
#$ -l h_rt=4:0:0
#$ -l mem_free=2G
#$ -S /bin/bash 
#$ -N Slice
#$ -cwd 
#$ -j y 
#$ -R y
#$ -t 1-26

source ~/.bash_profile
i=$(($SGE_TASK_ID - 1))

module load imp/2.8.0

path=/netapp/sali/ilan/Cop9/Cop9_MH_Inter_Intra

while read p;
do
    rmf_direc=$(echo $p | awk '{print $1}')
    rmf_repli=$(echo $p | awk '{print $2}')
    rmf_frame=$(echo $p | awk '{print $4}')

    rmf_slice ${path}/modeling${rmf_direc}/output/rmfs/${rmf_repli}.rmf3 ${rmf_direc}_${rmf_repli}_${rmf_frame}.rmf3 -f ${rmf_frame}

done < tmp_${i}

