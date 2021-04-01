#!/bin/bash

wm=(Test1D_10km Test1D_5km Test1D_2.5km Test1D_1.25km Test1D_0.625km)

#grid resolution
res=(10000 5000 2500 1250 625)

#shape function
whichsf='gimpm smpm'

#initial particles per grid cell
ppe='4 9 16'

#use particle reweighting?
reweight='false true'

#for m in $wm
for ((i = 0; i < ${#wm[@]}; ++i)); do
    m=(${wm[$i]})
    r=(${res[$i]})
    for s in $whichsf
    do
	for k in $reweight
	do
	    if [ $s == 'gimpm' ] && [ $k == 'true' ]
	    then
		echo "Skipping reweight==true for gimpm"
		continue
	    fi
	    for j in $ppe
	    do

		echo $m
		echo $r
  		echo $s
		echo $k
		echo $j

		sed  "s/<sf>/$s/; s/<ppe>/$j/; s/<m>/$m/; s/<gr>/$r/; s/<rw>/$k/g" \
		     steady1d.sif > steady_1d.sif

		ElmerSolver steady_1d.sif
	    done
	done
    done
done
