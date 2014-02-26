#!/bin/bash
K=5
FOLD=2
FIX_HYP=0
PRED_MODE=2

cd outputs
echo "split data"
../bin/split ../data/demo.csv demo ${FOLD}

echo "change format"
for ((i=0; i < ${FOLD}; i++))
do
	../bin/changeFormat ${i}.a train_${i}.dat 0
	../bin/changeFormat ${i}.b test_${i}.dat 1
	rm ${i}.a ${i}.b
done

mkdir -p model
echo "train mnb"
for ((i=0; i<${FOLD}; i++))
do
	F=${K}HMNB${i}
	../bin/hmnb_train train_${i}.dat model/${F} ${FIX_HYP} ${K}
done

echo "test hmnb"
mkdir -p result
for ((i=0; i<${FOLD}; i++))
do
	F=${K}HMNB${i}
	OF=result/${F}
	mkdir -p ${OF}
	../bin/hmnb_pred model/${F}/ test_${i}.dat ${OF}/test_${i}.outputs ${PRED_MODE}
	../bin/hmnb_pred model/${F}/ ../data/bad.dat ${OF}/bad_${i}.outputs ${PRED_MODE}
	cat ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs > ${OF}/outputs_${i}
	../bin/combine ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs ${OF}/targets_${i}
	rm ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs 
	../bin/calc_auc ${OF}/targets_${i} ${OF}/outputs_${i} > ${OF}/auc_${i}
	cat ${OF}/auc_${i}
done


## buggy
#echo "train mnb"
#for ((i=0; i<${FOLD}; i++))
#do
#	F=${K}MNB${i}
#	../bin/mnb_train train_${i}.dat model/${F} ${K}
#done

#echo "test mnb"
#mkdir -p result
#for ((i=0; i<${FOLD}; i++))
#do
#	F=${K}MNB${i}
#	OF=result/${F}
#	mkdir -p ${OF}
#	../bin/mnb_pred model/${F}/ test_${i}.dat ${OF}/test_${i}.outputs ${PRED_MODE}
#	../bin/mnb_pred model/${F}/ ../data/bad.dat ${OF}/bad_${i}.outputs ${PRED_MODE}
#	cat ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs > ${OF}/outputs_${i}
#	../bin/combine ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs ${OF}/targets_${i}
#	rm ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs 
#	../bin/calc_auc ${OF}/targets_${i} ${OF}/outputs_${i} > ${OF}/auc_${i}
#	cat ${OF}/auc_${i}
#done

echo "train nb"
for ((i=0; i<${FOLD}; i++))
do
	F=${K}NB${i}
	../bin/nb_train train_${i}.dat model/${F} ${K}
done

echo "test nb"
mkdir -p result
for ((i=0; i<${FOLD}; i++))
do
	F=${K}NB${i}
	OF=result/${F}
	mkdir -p ${OF}
	../bin/nb_pred model/${F}/ test_${i}.dat ${OF}/test_${i}.outputs
	../bin/nb_pred model/${F}/ ../data/bad.dat ${OF}/bad_${i}.outputs
	cat ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs > ${OF}/outputs_${i}
	../bin/combine ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs ${OF}/targets_${i}
	rm ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs 
	../bin/calc_auc ${OF}/targets_${i} ${OF}/outputs_${i} > ${OF}/auc_${i}
	cat ${OF}/auc_${i}
done

echo "train nbc"
for ((i=0; i<${FOLD}; i++))
do
	F=NBC${i}
	../bin/nbc_train train_${i}.dat model/${F}
done

echo "test nbc"
mkdir -p result
for ((i=0; i<${FOLD}; i++))
do
	F=NBC${i}
	OF=result/${F}
	mkdir -p ${OF}
	../bin/nbc_pred model/${F}/ test_${i}.dat ${OF}/test_${i}.outputs
	../bin/nbc_pred model/${F}/ ../data/bad.dat ${OF}/bad_${i}.outputs
	cat ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs > ${OF}/outputs_${i}
	../bin/combine ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs ${OF}/targets_${i}
	rm ${OF}/test_${i}.outputs ${OF}/bad_${i}.outputs 
	../bin/calc_auc ${OF}/targets_${i} ${OF}/outputs_${i} > ${OF}/auc_${i}
	cat ${OF}/auc_${i}
done
