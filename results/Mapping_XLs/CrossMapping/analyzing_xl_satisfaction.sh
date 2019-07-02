#!/usr/bin/bash

for fi in `ls *CSN?.txt`;
do
    ndssre=$(grep DSSO $fi | awk '{if($2!=$4) print $1}' | wc -l)
    ndhsre=$(grep DHSO $fi | awk '{if($2!=$4) print $1}' | wc -l)
    nbmsre=$(grep BMSO $fi | awk '{if($2!=$4) print $1}' | wc -l)

    ndssra=$(grep DSSO $fi | awk '{if($2==$4) print $1}' | wc -l)
    ndhsra=$(grep DHSO $fi | awk '{if($2==$4) print $1}' | wc -l)
    nbmsra=$(grep BMSO $fi | awk '{if($2==$4) print $1}' | wc -l)

    dssra=$(grep DSSO $fi | awk '{if($2==$4) print $0}' | grep ', 0,' | wc -l)
    dhsra=$(grep DHSO $fi | awk '{if($2==$4) print $0}' | grep ', 0,' | wc -l)
    bmsra=$(grep BMSO $fi | awk '{if($2==$4) print $0}' | grep ', 0,' | wc -l)
    
    dssre=$(grep DSSO $fi | awk '{if($2!=$4) print $0}' | grep ', 0,' | wc -l)
    dhsre=$(grep DHSO $fi | awk '{if($2!=$4) print $0}' | grep ', 0,' | wc -l)
    bmsre=$(grep BMSO $fi | awk '{if($2!=$4) print $0}' | grep ', 0,' | wc -l)
    
    root=$(echo $fi | sed 's/.txt/.S.txt/g')
    
    echo "DSSO" $dssre $dssra $ndssre $ndssra $(echo "scale=2;1-$dssra/$ndssra" | bc) $(echo "scale=2;1-$dssre/$ndssre" | bc) > $root
    echo "DHSO" $dhsre $dhsra $ndhsre $ndhsra $(echo "scale=2;1-$dhsra/$ndhsra" | bc) $(echo "scale=2;1-$dhsre/$ndhsre" | bc) >> $root
    echo "BMSO" $bmsre $bmsra $nbmsre $nbmsra $(echo "scale=2;1-$bmsra/$nbmsra" | bc) $(echo "scale=2;1-$bmsre/$nbmsre" | bc) >> $root

done