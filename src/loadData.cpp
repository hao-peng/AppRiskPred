/**********
Function for loading the model
Author: Hao Peng (penghATpurdue.edu)
Last updated time: Feb 25, 2014
**********/
#include "model.h"


int loadData(char *input) {
	FILE *fp;
	fp = fopen(input, "r");
	if (fp == NULL) {
		printf("cannot open file: %s\n", input);
		exit(-1);
	}
	fscanf(fp, "%d %d", &M, &L);
	N = (int *) malloc(sizeof(int)*M);
	R = (int ***) malloc(sizeof(int**)*M);
	for (int m = 0; m < M; m++) {
		fscanf(fp, "%d", N+m);
		R[m] = (int**) malloc(sizeof(int*)*N[m]);
		for (int n = 0; n < N[m]; n++) {
			R[m][n] = (int*) malloc(sizeof(int)*L);
			for (int l = 0; l < L; l++) {
				fscanf(fp, "%d", &R[m][n][l]);
			}
		}
	}
	fclose(fp);
	return 1;
}
