O =		obj
S =		src
LIB = 	include

all	: avTiff flatten background 

avPower	: $O/avPower.o $O/tiffFunc.o
	cc -o avPower $O/avPower.o $O/tiffFunc.o -ltiff


avTiff	: $O/avTiff.o $O/tiffFunc.o
	cc -o avTiff $O/avTiff.o $O/tiffFunc.o -ltiff


flatten : $O/tiffFunc.o $O/flatten.o
	cc -o flatten $O/flatten.o $O/tiffFunc.o -ltiff

background : $O/tiffFunc.o $O/background.o
	cc -o background $O/background.o $O/tiffFunc.o -ltiff

tiffTest : $O/tiffFunc.o $O/tiffTest.o
	cc -o tiffTest $O/tiffTest.o $O/tiffFunc.o -ltiff

$O/avTiff.o : $S/avTiff.c $(LIB)/tiffFunc.h
	cc -c -g $S/avTiff.c
	mv avTiff.o $O/	

$O/tiffFunc.o : $S/tiffFunc.c $(LIB)/tiffFunc.h
	cc -c -g $S/tiffFunc.c
	mv tiffFunc.o $O/

$O/flatten.o :	$S/flatten.c $(LIB)/tiffFunc.h
	cc -c -g $S/flatten.c
	mv flatten.o $O/

$O/background.o :	$S/background.c $(LIB)/tiffFunc.h
	cc -c -g $S/background.c
	mv background.o $O/
