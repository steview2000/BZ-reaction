#!/bin/bash

DATASET=1302113
for i in {0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.58,0.60,0.61,0.62,0.63,0.64,0.65,0.70,0.72,0.74,0.800}

do
	DATASET=$((DATASET+1))
	echo "${DATASET}S| $i|">>/home/stwe/Data/HomoForcing/HomoForcingSim.dat
done

