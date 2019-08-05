#!/usr/bin/bash

dir0=DSSO_BMSO_DHSO
for dir in BMSO DHSO DSSO BMSO_DHSO DSSO_BMSO DSSO_DHSO DSSO_BMSO_DHSO_P;
do
    for i in `seq 0 7`;
    do
	for j in `seq $i 7`;
	do
	    
	    mean1=$(cat RMSD_WRT_Centroid_${dir}/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1;b++}END{print sum/b}');
	    devi1=$(cat RMSD_WRT_Centroid_${dir}/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1; sum2+=$1*$1; b++}END{print sqrt(sum2/b - sum/b*sum/b)}');
	    
	    mean2=$(cat RMSD_WRT_Centroid_${dir0}/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1;b++}END{print sum/b}');
	    devi2=$(cat RMSD_WRT_Centroid_${dir0}/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1; sum2+=$1*$1; b++}END{print sqrt(sum2/b - sum/b*sum/b)}');

	    cent1=$(cat structures/RMSD_CSN_vs_$(echo $dir | sed 's/_/./g')/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1;b++}END{print sum/b}');
	    
	    echo ${dir} RB$i RB$j $mean1 $devi1 $mean2 $devi2 $cent1;
	done;

    done > Results.${dir}.vs.CSN.txt;
    
    cat Results.${dir}.vs.CSN.txt | \
	sed 's/RB7/CSN8/g' | \
	sed 's/RB6/CSN3/g' | \
	sed 's/RB5/CSN2/g' | \
	sed 's/RB4/CSN1/g' | \
	sed 's/RB3/CSN7/g' | \
	sed 's/RB2/CSN6/g' | \
	sed 's/RB1/CSN5/g' | \
	sed 's/RB0/CSN4/g' | \
	sort -n -k1 -k2 > Results.${dir}.vs.CSN.Rbs.txt
    rm Results.${dir}.vs.CSN.txt;
done

dir0=DSSO_BMSO_DHSO
for dir in Crys CRL1A CRL4A;
do
    for i in `seq 0 7`;
    do
	for j in `seq $i 7`;
	do
	
	    mean2=$(cat RMSD_WRT_Centroid_${dir0}/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1;b++}END{print sum/b}');
	    devi2=$(cat RMSD_WRT_Centroid_${dir0}/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1; sum2+=$1*$1; b++}END{print sqrt(sum2/b - sum/b*sum/b)}');
	    cent1=$(cat structures/RMSD_CSN_vs_$(echo $dir | sed 's/_/./g')/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1;b++}END{print sum/b}');

	    echo ${dir} RB$i RB$j $mean2 $devi2 $cent1;
	done;
    
    done > Results.${dir}.vs.CSN.txt;
done
cat Results.${dir}.vs.CSN.txt | \
    sed 's/RB7/CSN8/g' | \
    sed 's/RB6/CSN3/g' | \
    sed 's/RB5/CSN2/g' | \
    sed 's/RB4/CSN1/g' | \
    sed 's/RB3/CSN7/g' | \
    sed 's/RB2/CSN6/g' | \
    sed 's/RB1/CSN5/g' | \
    sed 's/RB0/CSN4/g' | \
    sort -n -k1 -k2 > Results.${dir}.vs.CSN.Rbs.txt
rm Results.${dir}.vs.CSN.txt;

#!/usr/bin/bash

dir0=DSSO_BMSO_DHSO_P
for dir in CRL1A CRL4A;
do
    for i in `seq 0 7`;
    do
	for j in `seq $i 7`;
	do
	    
	        mean2=$(cat RMSD_WRT_Centroid_${dir0}/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1;b++}END{print sum/b}');
		    devi2=$(cat RMSD_WRT_Centroid_${dir0}/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1; sum2+=$1*$1; b++}END{print sqrt(sum2/b - sum/b*sum/b)}');
		        cent1=$(cat structures/RMSD_CSN_vs_$(echo $dir | sed 's/_/./g')/RMSD_pairs_all_versus_centroid_${i}_${j}.csv | awk '{sum+=$1;b++}END{print sum/b}');

			    echo ${dir} RB$i RB$j $mean2 $devi2 $cent1;
			    done;
    
    done > Results.${dir}.vs.CSNn.txt;
    
    cat Results.${dir}.vs.CSNn.txt | \
	sed 's/RB7/CSN8/g' | \
	sed 's/RB6/CSN3/g' | \
	sed 's/RB5/CSN2/g' | \
	sed 's/RB4/CSN1/g' | \
	sed 's/RB3/CSN7/g' | \
	sed 's/RB2/CSN6/g' | \
	sed 's/RB1/CSN5/g' | \
	sed 's/RB0/CSN4/g' | \
	sort -n -k1 -k2 > Results.${dir}.vs.CSNn.Rbs.txt
    rm Results.${dir}.vs.CSNn.txt;
done
