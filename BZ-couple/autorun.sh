#!/bin/bash
# Skript to couple two (different) spirals both having the same spiral frequencies 

trap "exit" INT

./CreateCoresAt -50 50 0.4
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502100S"'>> include/header.h
echo '#define I0 0.0'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 0.8
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502101S"'>> include/header.h
echo '#define I0 0.01'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple
#./CreateCoresAt -50 50 1.2
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502102S"'>> include/header.h
echo '#define I0 0.02'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple


#./CreateCoresAt -50 50 1.6
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502103S"'>> include/header.h
echo '#define I0 0.03'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 2.0
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502104S"'>> include/header.h
echo '#define I0 0.04'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple


#./CreateCoresAt -50 50 2.4
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502105S"'>> include/header.h
echo '#define I0 0.05'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 2.8
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502106S"'>> include/header.h
echo '#define I0 0.06'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 3.2
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502107S"'>> include/header.h
echo '#define I0 0.07'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 3.6
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502108S"'>> include/header.h
echo '#define I0 0.08'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 4.0
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502109S"'>> include/header.h
echo '#define I0 0.09'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple



#./CreateCoresAt -50 50 4.4
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502111S"'>> include/header.h
echo '#define I0 0.10'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 4.8
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502112S"'>> include/header.h
echo '#define I0 0.11'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 5.2
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502113S"'>> include/header.h
echo '#define I0 0.12'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 5.4
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502114S"'>> include/header.h
echo '#define I0 0.13'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

#./CreateCoresAt -50 50 5.8
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502115S"'>> include/header.h
echo '#define I0 0.14'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple

cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502116S"'>> include/header.h
echo '#define I0 0.15'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502117S"'>> include/header.h
echo '#define I0 0.16'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502118S"'>> include/header.h
echo '#define I0 0.17'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple
#./CreateCoresAt -50 50 6.2
cp ~/Data/999Core3/a-1.dat a-1.dat
cp ~/Data/999Core3/c-1.dat c-1.dat
cp ~/Data/999Core3/a-2.dat a-2.dat
cp ~/Data/999Core3/c-2.dat c-2.dat

cat include/header0.h > include/header.h
echo '#define DATASET "1502119S"'>> include/header.h
echo '#define I0 0.2'>>include/header.h
echo '#endif' >>include/header.h
make
nice -n 19 ./BZ_couple
