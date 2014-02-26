HEADERS = dataIO.h params.h
OPTIONS = -g -Wall
OUTPUT  = ../bin/calc_auc

all: ${OUTPUT}

calc_auc.o: calc_auc.c $(HEADERS)
	gcc -c $(OPTIONS) calc_auc.c	

dataIO.o: dataIO.c $(HEADERS)
	gcc -c $(OPTIONS) dataIO.c
	
${OUTPUT}: calc_auc.o dataIO.o
	gcc $(OPTIONS) -o ${OUTPUT} calc_auc.o dataIO.o -lm

clean:
	rm -f *.o
	rm -f ${OUTPUT}



