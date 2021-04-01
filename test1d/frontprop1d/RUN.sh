#!/bin/bash

wm=(Test1D_5km Test1D_2.5km Test1D_1.25km Test1D_0.625km)

#grid resolution
res=(5000 2500 1250 625)

#shape function
whichsf='gimpm smpm'

#initial particles per grid cell
ppe='4 9 16'



for ((i = 0; i < ${#wm[@]}; ++i)); do
    m=(${wm[$i]})
    r=(${res[$i]})
    for i in $whichsf
    do
	for j in $ppe
	do
  	    echo $i
	    echo $j
	    echo $m
	    echo $r

	    sed  "s/<sf>/$i/; s/<ppe>/$j/; s/<m>/$m/; s/<gr>/$r/g" \
		 frontprop1d.sif > fp.sif

	    ElmerSolver fp.sif
	    #rm fp_5k.sif
	done
    done
done
