/**********
Prediction for NB model
Author: Hao Peng (penghATpurdue.edu)
Last updated time: Sep 17, 2015
**********/
#ifdef USE_OMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "asa121.h"

//#define IN_LOOP 5000
#define OUT_LOOP 8000
#define EPS 1e-10
#define EQUAL(a,b) ((a-b)<EPS && a-b>-EPS)

#define UPDATE_ALPHA 0
#define FIX_ALPHA 1

//#define DEBUG
//#define PRIOR
#define NEW_PRIOR



#define MAX_ESTIMATE 0
#define MEAN_ESTIMATE 1
#define SUM_ESTIMATE 2


int K = 0;
int M = 0;
int L = 0;
int *N = NULL;
int ***R = NULL;

double **mytheta = NULL;
double *r = NULL;
double log_likelihood = 0;
int method = MAX_ESTIMATE;

void load(const char *input) {
	char filename[100];
	FILE *fp;


	strcpy(filename, input);
	strcat(filename, "/r");
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("cannot open file r\n");
		exit(-1);
	}
	fscanf(fp, "%d\n", &K);
	r = (double *) malloc(sizeof(double)*K);
	for (int k = 0; k < K; k++) {
		fscanf(fp, "%lf", &r[k]);
	}
	fclose(fp);

	strcpy(filename, input);
	strcat(filename, "/theta");
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("cannot open file theta\n");
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
		printf("cannot open file %s\n", data);
		exit(-1);
	}

	ofp = fopen(output, "w");
	if (ofp == NULL) {
		printf("cannot open file %s\n", output);
		exit(-1);
	}

	int NN;

	fscanf(fp, "%d %d", &NN, &L);
	int *newR = (int *) malloc(sizeof(int)*L);
	double *p = (double *) malloc(sizeof(double)*M);

	for (int i = 0; i < NN; i++){
		double maxP = 0;
		for (int j = 0; j < L; j++) {
			fscanf(fp, "%d", &newR[j]);
		}
		for (int k = 0; k < K; k++) {
			double q = r[k];
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
			maxP += q;
		}
		fprintf(ofp, "%e\n", maxP);
	}

	free(p);
	free(newR);
	fclose(fp);
}

int main(int argc, char **argv) {
	if (argc != 4) {
		printf("Usage: %s input data output\n", argv[0]);
		return 0;
	}

	//printf("load\n");
	load(argv[1]);
	estimate(argv[2], argv[3]);
	return 0;
}
