all: BZ_couple CreateCoresAt

OBJ = ./obj
SRC = ./src
LIB = ./include
FLAG = ../analysisPic/obj/tiffFunc.o -lm -ltiff

BZ_couple:	${OBJ}/BZ_couple.o ${OBJ}/numFunc.o
	cc -O3 -o BZ_couple ${OBJ}/BZ_couple.o ${OBJ}/numFunc.o ${FLAG} 

CreateCoresAt:	${OBJ}/CreateCoresAt.o ${OBJ}/numFunc.o
	cc -O3 -o CreateCoresAt ${OBJ}/CreateCoresAt.o ${OBJ}/numFunc.o ${FLAG} 


${OBJ}/CreateCoresAt.o: ${SRC}/CreateCoresAt.c ${LIB}/header.h
	cc -c -O3 ${SRC}/CreateCoresAt.c
	mv CreateCoresAt.o ${OBJ}/


${OBJ}/BZ_couple.o: ${SRC}/BZ_couple.c ${LIB}/header.h
	cc -c -O3 ${SRC}/BZ_couple.c
	mv BZ_couple.o ${OBJ}/

${OBJ}/numFunc.o: ${SRC}/numFunc.c ${LIB}/header.h
	cc -c -O3 ${SRC}/numFunc.c
	mv numFunc.o ${OBJ}
