OPTIONS = -O2
#OPTIONS = -fopenmp -O2
OUTPUT_DIR = ../bin/

all: MKDIR ../bin/hmnb_train ../bin/hmnb_pred ../bin/mnb_pred \
		../bin/mnb_train ../bin/nb_train ../bin/nb_pred \
		../bin/nbc_train ../bin/nbc_pred

MKDIR:
	mkdir -p OUTPUT_DIR

clean:
	rm -f *.o
	rm -f ${OUTPUT_DIR}hmnb_train
	rm -f ${OUTPUT_DIR}mnb_train
	rm -f ${OUTPUT_DIR}nb_train
	rm -f ${OUTPUT_DIR}nbc_train
	rm -f ${OUTPUT_DIR}hmnb_pred
	rm -f ${OUTPUT_DIR}mnb_pred
	rm -f ${OUTPUT_DIR}nb_pred
	rm -f ${OUTPUT_DIR}nbc_pred

../bin/hmnb_pred: hmnb_pred.o 
	g++ ${OPTIONS} -o ${OUTPUT_DIR}hmnb_pred hmnb_pred.o 

../bin/mnb_pred: mnb_pred.o 
	g++ ${OPTIONS} -o ${OUTPUT_DIR}mnb_pred mnb_pred.o 

../bin/nb_pred: nb_pred.o 
	g++ ${OPTIONS} -o ${OUTPUT_DIR}nb_pred nb_pred.o 

../bin/nbc_pred: nbc_pred.o 
	g++ ${OPTIONS} -o ${OUTPUT_DIR}nbc_pred nbc_pred.o 

../bin/hmnb_train: hmnb_train.o initTheta.o loadData.o digamma_function.o \
		sumlogexp.o asa121.o gamma_function.o ln_gamma_function.o
	g++ ${OPTIONS} -o ${OUTPUT_DIR}hmnb_train hmnb_train.o initTheta.o loadData.o digamma_function.o \
		sumlogexp.o asa121.o gamma_function.o ln_gamma_function.o \

../bin/mnb_train: mnb_train.o initTheta.o loadData.o digamma_function.o \
		sumlogexp.o asa121.o gamma_function.o ln_gamma_function.o
	g++ ${OPTIONS} -o ${OUTPUT_DIR}mnb_train mnb_train.o initTheta.o loadData.o digamma_function.o \
		sumlogexp.o asa121.o gamma_function.o ln_gamma_function.o \

../bin/nb_train: nb_train.o initTheta.o loadData.o digamma_function.o \
		sumlogexp.o asa121.o gamma_function.o ln_gamma_function.o
	g++ ${OPTIONS} -o ${OUTPUT_DIR}nb_train nb_train.o initTheta.o loadData.o digamma_function.o \
		sumlogexp.o asa121.o gamma_function.o ln_gamma_function.o \

../bin/nbc_train: nbc_train.o initTheta.o loadData.o digamma_function.o \
		sumlogexp.o asa121.o gamma_function.o ln_gamma_function.o
	g++ ${OPTIONS} -o ${OUTPUT_DIR}nbc_train nbc_train.o initTheta.o loadData.o digamma_function.o \
		sumlogexp.o asa121.o gamma_function.o ln_gamma_function.o \

hmnb_train.o: hmnb_train.cpp model.h
	g++ ${OPTIONS} -c -o hmnb_train.o hmnb_train.cpp

mnb_train.o: mnb_train.cpp
	g++ ${OPTIONS} -c -o mnb_train.o mnb_train.cpp

nb_train.o: nb_train.cpp
	g++ ${OPTIONS} -c -o nb_train.o nb_train.cpp

nbc_train.o: nbc_train.cpp
	g++ ${OPTIONS} -c -o nbc_train.o nbc_train.cpp

hmnb_pred.o: hmnb_pred.cpp model.h
	g++ ${OPTIONS} -c -o hmnb_pred.o hmnb_pred.cpp

mnb_pred.o: mnb_pred.cpp 
	g++ ${OPTIONS} -c -o mnb_pred.o mnb_pred.cpp

nb_pred.o: nb_pred.cpp
	g++ ${OPTIONS} -c -o nb_pred.o nb_pred.cpp

nbc_pred.o: nbc_pred.cpp
	g++ ${OPTIONS} -c -o nbc_pred.o nbc_pred.cpp

loadData.o: loadData.cpp model.h
	g++ ${OPTIONS} -c -o loadData.o loadData.cpp

initTheta.o: initTheta.cpp model.h
	g++ ${OPTIONS} -c -o initTheta.o initTheta.cpp

digamma_function.o: digamma_function.c
	g++ ${OPTIONS} -c -o digamma_function.o digamma_function.c

sumlogexp.o: sumlogexp.cpp
	g++ ${OPTIONS} -c -o sumlogexp.o sumlogexp.cpp

asa121.o: asa121.c
	g++ ${OPTIONS} -c -o asa121.o asa121.c

gamma_function.o: gamma_function.c
	g++ ${OPTIONS} -c -o gamma_function.o gamma_function.c

ln_gamma_function.o: ln_gamma_function.c
	g++ ${OPTIONS} -c -o ln_gamma_function.o ln_gamma_function.c
