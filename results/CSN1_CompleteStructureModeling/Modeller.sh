#!/bin/bash
#$ -l arch=linux-x64 
#$ -l h_rt=18:0:0
#$ -l mem_free=2G
#$ -S /bin/bash 
#$ -N csn
#$ -cwd 
#$ -j y 
#$ -R y
#$ -t 1-100

module load modeller
i=$(($SGE_TASK_ID - 1))

mkdir Run_${i}
cp CSN1_Modeller.py Run_${i}
cp CSN.pdb Run_${i}
cp CSN1_CSN1.ali Run_${i}
cd Run_${i}

python CSN1_Modeller.py $i 

cd ..