#!/bin/bash
#$ -l arch=linux-x64 
#$ -l h_rt=300:0:0
#$ -l mem_free=2G
#$ -S /bin/bash 
#$ -N pref
#$ -cwd 
#$ -j y 
#$ -R y
#$ -t 1-60

source ~/.bash_profile
i=$(($SGE_TASK_ID - 1))
DIR0=/netapp/sali/ilan/Cop9/Cop9_MH_Inter_Intra/modeling${i}/output/

for rep in `seq 0 7`; 
do 
    nlines=$(wc -l ${DIR0}/stat.${rep}.out | awk '{print $1}')
    for f in `seq 1 $((nlines - 1))`;
    do

	process_output.py -f ${DIR0}/stat.${rep}.out -n ${f} | grep -E "Total_Score|ConnectivityRestraint|Distance|rmf|ExcludedVolumeSphere_CSN|CrossLinkingMassSpectrometryRestraint_Data_Score|MonteCarlo_Nframe|^CrossLinkingMassSpectrometryRestraint_Psi_PSI" > tmp_${i}

	rmffile=$(grep rmf_file tmp_${i} | awk '{print $2}')
	frameid=$(grep rmf_frame_index tmp_${i} |  awk '{print $2}')
	mc_steps=$(grep MonteCarlo_Nframe tmp_${i} | awk '{print $2}')

	totsco=$(grep Total_Score tmp_${i} | awk '{print $2}')
	evscor=$(grep ExcludedVolumeSphere_CSN tmp_${i} | awk '{print $2}')
	connect=$(grep ConnectivityRestraint tmp_${i} | awk '{sum+=$2}END{print sum}')

	dssa=$(grep CrossLinkingMassSpectrometryRestraint_Data_Score_DSS_Inter tmp_${i} | awk '{print $2}')
	dhsa=$(grep CrossLinkingMassSpectrometryRestraint_Data_Score_DHS_Inter tmp_${i} | awk '{print $2}')
	bmsa=$(grep CrossLinkingMassSpectrometryRestraint_Data_Score_BMS_Inter tmp_${i} | awk '{print $2}')

	dssb=$(grep CrossLinkingMassSpectrometryRestraint_Data_Score_DSS_Intra tmp_${i} | awk '{print $2}')
        dhsb=$(grep CrossLinkingMassSpectrometryRestraint_Data_Score_DHS_Intra tmp_${i} | awk '{print $2}')
        bmsb=$(grep CrossLinkingMassSpectrometryRestraint_Data_Score_BMS_Intra tmp_${i} | awk '{print $2}')

	sdssa=$(grep CrossLinkingMassSpectrometryRestraint_Distance_\|DSS_Inter tmp_${i} | awk 'BEGIN{sum=0}{if($2<35) sum++}END{print sum}')
	sdhsa=$(grep CrossLinkingMassSpectrometryRestraint_Distance_\|DHS_Inter tmp_${i} | awk 'BEGIN{sum=0}{if($2<35) sum++}END{print sum}')
	sbmsa=$(grep CrossLinkingMassSpectrometryRestraint_Distance_\|BMS_Inter tmp_${i} | awk 'BEGIN{sum=0}{if($2<45) sum++}END{print sum}')
	
	sdssb=$(grep CrossLinkingMassSpectrometryRestraint_Distance_\|DSS_Intra tmp_${i} | awk 'BEGIN{sum=0}{if($2<35) sum++}END{print sum}')
        sdhsb=$(grep CrossLinkingMassSpectrometryRestraint_Distance_\|DHS_Intra tmp_${i} | awk 'BEGIN{sum=0}{if($2<35) sum++}END{print sum}')
        sbmsb=$(grep CrossLinkingMassSpectrometryRestraint_Distance_\|BMS_Intra tmp_${i} | awk 'BEGIN{sum=0}{if($2<45) sum++}END{print sum}')


	pdssa=$(grep CrossLinkingMassSpectrometryRestraint_Psi_PSI_DSS_Inter tmp_${i} | awk '{print $2}')
	pdhsa=$(grep CrossLinkingMassSpectrometryRestraint_Psi_PSI_DHS_Inter tmp_${i} | awk '{print $2}')
	pbmsa=$(grep CrossLinkingMassSpectrometryRestraint_Psi_PSI_BMS_Inter tmp_${i} | awk '{print $2}')
		
	pdssb=$(grep CrossLinkingMassSpectrometryRestraint_Psi_PSI_DSS_Intra tmp_${i} | awk '{print $2}')
        pdhsb=$(grep CrossLinkingMassSpectrometryRestraint_Psi_PSI_DHS_Intra tmp_${i} | awk '{print $2}')
        pbmsb=$(grep CrossLinkingMassSpectrometryRestraint_Psi_PSI_BMS_Intra tmp_${i} | awk '{print $2}')



	echo $i $rep $rmffile $frameid $mc_steps " | " $totsco $evscor $connect " | " $dssa $dhsa $bmsa $sdssa $sdhsa $sbmsa $pdssa $pdhsa $pbmsa " | " $dssb $dhsb $bmsb $sdssb $sdhsb $sbmsb $pdssb $pdhsb $pbmsb >> Scores_R${i}.txt	
	rm tmp_${i}
    done
done

