all:
	make -C src all
	make -C calc_auc all 
	make -C preprocess all

clean: reset
	make -C src clean
	make -C calc_auc clean
	make -C preprocess clean
	

reset:
	rm -rf outputs/*
