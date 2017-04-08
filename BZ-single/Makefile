all: BZ-single

OBJ = ./obj
SRC = ./src
LIB = ./include
FLAG = /home/sweiss/analysis/obj/tiffFunc.o -lm -ltiff

BZ-single:	${OBJ}/BZ-single.o ${OBJ}/gnuplot_i.o ${OBJ}/numFunc.o
	cc -O3 -o BZ-single ${OBJ}/BZ-single.o ${OBJ}/gnuplot_i.o ${OBJ}/numFunc.o ${FLAG} 

${OBJ}/BZ-single.o: ${SRC}/BZ-single.c ${LIB}/header.h
	cc -c -O3 ${SRC}/BZ-single.c
	mv BZ-single.o ${OBJ}/


${OBJ}/numFunc.o: ${SRC}/numFunc.c ${LIB}/header.h
	cc -c -O3 ${SRC}/numFunc.c
	mv numFunc.o ${OBJ}
