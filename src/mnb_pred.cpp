/**********
Prediction for MNB model
Author: Hao Peng (penghATpurdue.edu)
Last updated time: Feb 25, 2014
**********/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "asa121.h"

#define OUT_LOOP 8000
#define EPS 1e-10
#define EQUAL(a,b) ((a-b)<EPS && a-b>-EPS)

#define UPDATE_ALPHA 0
#define FIX_ALPHA 1

//#define DEBUG
#define PRIOR
//#define NEW_PRIOR


#define MAX_ESTIMATE 0
#define MEAN_ESTIMATE 1
#define SUM_ESTIMATE 2


int K = 0;
int M = 0;
int L = 0;
int *N = NULL;
int ***R = NULL;

double **mytheta = NULL;
double **r = NULL;
double log_likelihood = 0;
int method = MAX_ESTIMATE;

void load(const char *input) {
	char filename[100];
	FILE *fp;


	strcpy(filename, input);
	strcat(filename, "/r");
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fscanf(fp, "%d %d\n", &M, &K);
	r = (double **) malloc(sizeof(double*)*M);
	N = (int *) malloc(sizeof(int)*M);
	for (int m = 0; m < M; m++) {
		r[m] = (double *) malloc(sizeof(double)*K);
		for (int k = 0; k < K; k++) {
			fscanf(fp, "%lf", &r[m][k]);
		}
	}
	fclose(fp);

	strcpy(filename, input);
	strcat(filename, "/theta");
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fscanf(fp, "%d %d\n", &K, &L);
	mytheta = (double **) malloc(sizeof(double *)*K);
	for (int k = 0; k < K; k++) {
		mytheta[k] = (double *) malloc(sizeof(double)*L);
		for (int l = 0; l < L; l++) {
			fscanf(fp, "%lf", &mytheta[k][l]);
		}
	}
	fclose(fp);
}

void estimate(char *data, char *output) {
	FILE *fp, *ofp;
	fp = fopen(data, "r");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}

	ofp = fopen(output, "w");
	if (ofp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}

	int NN;

	fscanf(fp, "%d %d", &NN, &L);
	int *newR = (int *) malloc(sizeof(int)*L);
	double *p = (double *) malloc(sizeof(double)*M);
	double maxP = 0;

	for (int i = 0; i < NN; i++){
		for (int j = 0; j < L; j++) {
			fscanf(fp, "%d", &newR[j]);
		}
#pragma omp parallel shared(M,N,K,L,p,newR,r)
		{
#pragma omp for schedule(dynamic,3)
			for (int m = 0; m < M; m++) {
				p[m] = 0;
				for (int k = 0; k < K; k++) {
					double q = r[m][k];
					for (int l = 0; l < L; l++) {
#ifdef NEW_PRIOR
						if (l == 56) continue;
#endif
						if (newR[l]) {
							q *= mytheta[k][l];
						} else {
							q *= (1 - mytheta[k][l]);
						}
					}
					p[m] += q;
				}
			}
		}
		maxP = 0;
		if (method == 0) {
			for (int m = 0; m < M; m++) {
				if (maxP < p[m]) {
					maxP = p[m];
				}
			}	
		} else if (method == 1) {
			for (int m = 0; m < M; m++ ){
				maxP += p[m];
			}
		}
		fprintf(ofp, "%e\n", maxP);
	}

	free(p);
	free(newR);
	fclose(fp);
}

int main(int argc, char **argv) {
	if (argc != 4 && argc != 5) {
		printf("Usage: %s input data output [method]\nmethod 0: max\nmode 1: sum", argv[0]);
		printf("method: different evaluation stretages for multiple groups\n");
		return 0;
	}
	if (argc == 5) {
		method = atoi(argv[4]);
	}

	printf("load\n");
	load(argv[1]);
	estimate(argv[2], argv[3]);
	return 0;
}
