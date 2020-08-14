#!/bin/bash

wm='Test1D_5km'
#'Test1D_10km'

#grid resolution
r='5000'
#'10000'

#shape function
whichsf='gimpm' #smpm'

#initial particles per grid cell
ppe='9' #9 16'

#use particle reweighting?
reweight='false' # true'

for i in $whichsf	 
do
    for j in $ppe   
    do
	for k in $reweight	 
	do

		echo $wm
		echo $r
  		echo $i
		echo $j		
		echo $k
		
		sed  "s/<sf>/$i/; s/<ppe>/$j/; s/<m>/$wm/; s/<gr>/$r/; s/<rw>/$k/g" \
		     steady1d.sif > steady_5k.sif
		
		ElmerSolver steady_5k.sif
		#rm steady_5k.sif
	done
    done
done
