#!/bin/bash

trap "exit" INT

DATASET=1308034
for i in {0.2,0.25,0.3,0.35}
do
	DATASET=$((DATASET+1))
	echo $DATASET
	cat include/header0.h >header.h
#	echo "#define I0 0.001">>header.h
#	echo "#define f_f $i">>header.h
	echo "#define EPSILON $i" >>header.h
#	echo "#define F $i">>header.h
	echo "#define DATASET \"${DATASET}S\" ">>header.h
	mv header.h include/
	make
#	cp a-Start.dat a-1.dat
#	cp c-Start.dat c-1.dat
	nice -n 19 ./BZ-single
	cp a-1.dat ~/Data/${DATASET}S/
	cp c-1.dat ~/Data/${DATASET}S/

done

echo "Done!! (BZ-single)">/home/stwe/.notify
