#!/bin/bash
# Skript to couple two (different) spirals both having the same spiral frequencies 

trap "exit" INT

# I0=0
#cp c-1Start.dat c-1.dat
#cp a-1Start.dat a-1.dat
#cp c-2Start.dat c-2.dat
#cp a-2Start.dat a-2.dat
#cp c-3Start.dat c-3.dat
#cp a-3Start.dat a-3.dat
#
cat include/header0.h > include/header.h
echo '#define DATASET "1402240S"'>> include/header.h
echo '#define I0 0.01'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_triple


# I0=0.001
#cp c-1Start.dat c-1.dat
#cp a-1Start.dat a-1.dat
#cp c-2Start.dat c-2.dat
#cp a-2Start.dat a-2.dat
#cp c-3Start.dat c-3.dat
#cp a-3Start.dat a-3.dat
#
#cat include/header0.h > include/header.h
#echo '#define DATASET "1402241S"'>> include/header.h
#echo '#define I0 0.02'>>include/header.h
#echo '#endif' >>include/header.h
#make
#nice -n 19 ./BZ_triple
##
##
### I0=0.002
##cp c-1Start.dat c-1.dat
##cp a-1Start.dat a-1.dat
##cp c-2Start.dat c-2.dat
##cp a-2Start.dat a-2.dat
##cp c-3Start.dat c-3.dat
##cp a-3Start.dat a-3.dat
##
##
#cat include/header0.h > include/header.h
#echo '#define DATASET "1402242S"'>> include/header.h
#echo '#define I0 0.03'>>include/header.h
#echo '#endif' >>include/header.h
#make
#nice -n 19 ./BZ_triple
##
##
### I0=0.005
##cp c-1Start.dat c-1.dat
##cp a-1Start.dat a-1.dat
##cp c-2Start.dat c-2.dat
##cp a-2Start.dat a-2.dat
##cp c-3Start.dat c-3.dat
##cp a-3Start.dat a-3.dat
##
##
##
#cat include/header0.h > include/header.h
#echo '#define DATASET "1402243S"'>> include/header.h
#echo '#define I0 0.04'>>include/header.h
#echo '#endif' >>include/header.h
#make
#nice -n 19 ./BZ_triple
##
##
### I0=0.008
##cp c-1Start.dat c-1.dat
##cp a-1Start.dat a-1.dat
##cp c-2Start.dat c-2.dat
##cp a-2Start.dat a-2.dat
##cp c-3Start.dat c-3.dat
##cp a-3Start.dat a-3.dat
##
#cat include/header0.h > include/header.h
#echo '#define DATASET "1402244S"'>> include/header.h
#echo '#define I0 0.05'>>include/header.h
#echo '#endif' >>include/header.h
#make
#nice -n 19 ./BZ_triple
##
##
## I0=0.01
##cp c-1Start.dat c-1.dat
##cp a-1Start.dat a-1.dat
##cp c-2Start.dat c-2.dat
##cp a-2Start.dat a-2.dat
##cp c-3Start.dat c-3.dat
##cp a-3Start.dat a-3.dat
#
#cat include/header0.h > include/header.h
#echo '#define DATASET "1402245S"'>> include/header.h
#echo '#define I0 0.06'>>include/header.h
#echo '#endif' >>include/header.h
#make
#nice -n 19 ./BZ_triple
##
#
## I0=0.02
##cp c-1Start.dat c-1.dat
##cp a-1Start.dat a-1.dat
##cp c-2Start.dat c-2.dat
##cp a-2Start.dat a-2.dat
##cp c-3Start.dat c-3.dat
##cp a-3Start.dat a-3.dat
#
#cat include/header0.h > include/header.h
#echo '#define DATASET "1402246S"'>> include/header.h
#echo '#define I0 0.07'>>include/header.h
#echo '#endif' >>include/header.h
#make
#nice -n 19 ./BZ_triple
#
#
## I0=0.05
##cp c-1Start.dat c-1.dat
##cp a-1Start.dat a-1.dat
##cp c-2Start.dat c-2.dat
##cp a-2Start.dat a-2.dat
##cp c-3Start.dat c-3.dat
##cp a-3Start.dat a-3.dat
#
#cat include/header0.h > include/header.h
#echo '#define DATASET "1402247S"'>> include/header.h
#echo '#define I0 0.08'>>include/header.h
#echo '#endif' >>include/header.h
#make
#nice -n 19 ./BZ_triple


# I0=0.1
#cp c-1Start.dat c-1.dat
#cp a-1Start.dat a-1.dat
#cp c-2Start.dat c-2.dat
#cp a-2Start.dat a-2.dat
#cp c-3Start.dat c-3.dat
#cp a-3Start.dat a-3.dat

#cat include/header0.h > include/header.h
#echo '#define DATASET "1309028-2S"'>> include/header.h
#echo '#define I0 0.1'>>include/header.h
#echo '#endif' >>include/header.h
#make
#nice -n 19 ./BZ_triple

# I0=0.2
#cp c-1Start.dat c-1.dat
#cp a-1Start.dat a-1.dat
#cp c-2Start.dat c-2.dat
#cp a-2Start.dat a-2.dat
#cp c-3Start.dat c-3.dat
#cp a-3Start.dat a-3.dat

#cat include/header0.h > include/header.h
#echo '#define DATASET "1309029-2S"'>> include/header.h
#echo '#define I0 0.2'>>include/header.h
#echo '#endif' >>include/header.h
#make
#nice -n 19 ./BZ_triple
#

