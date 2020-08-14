#!/bin/bash

wm='Test1D_5km'
#'Test1D_10km'

#grid resolution
r='5000'
#'10000'

#shape function
whichsf='gimpm' # smpm'

#initial particles per grid cell
ppe='9' # 4 16'


for i in $whichsf	 
do
    for j in $ppe     
    do
  		echo $i
		echo $j
		echo $wm
		echo $r
		
		sed  "s/<sf>/$i/; s/<ppe>/$j/; s/<m>/$wm/; s/<gr>/$r/g" \
		     frontprop1d.sif > fp_5k.sif
		
		ElmerSolver fp_5k.sif
		#rm fp_5k.sif
    done
done

