all: BZ_triple 

OBJ = ./obj
SRC = ./src
LIB = ./include
FLAG = /home/stwe/analysis/obj/tiffFunc.o -lm -ltiff

BZ_triple:	${OBJ}/BZ_triple.o ${OBJ}/numFunc.o
	cc -O3 -o BZ_triple ${OBJ}/BZ_triple.o ${OBJ}/numFunc.o ${FLAG} 


${OBJ}/BZ_triple.o: ${SRC}/BZ_triple.c ${LIB}/header.h
	cc -c -O3 ${SRC}/BZ_triple.c
	mv BZ_triple.o ${OBJ}/

${OBJ}/numFunc.o: ${SRC}/numFunc.c ${LIB}/header.h
	cc -c -O3 ${SRC}/numFunc.c
	mv numFunc.o ${OBJ}
