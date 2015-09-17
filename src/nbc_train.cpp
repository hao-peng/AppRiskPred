/**********
Training for NBC model
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

int loadData(char *input);

void initTheta();

//double loglikelihood();

double DiGamma_Function( double x );

double Ln_Gamma_Function(double x);

double logsumexp(double nums[], int ct);

int K = 0;
int M = 0;
int L = 0;
int *N = NULL;
int ***R = NULL;

#ifdef PRIOR
double A = 1.0;
double B = 1.0;
#endif
#ifdef NEW_PRIOR
double *A = 0;
double *B = 0;
#endif

double **mytheta = NULL;

void initTrain() {
	mytheta = (double **) malloc(sizeof(double*)*M);
	for (int m = 0; m < M; m++ ){
		mytheta[m] = (double *) malloc(sizeof(double)*L);
		for (int l = 0; l < L; l++) {
			mytheta[m][l] = rand()/(double) RAND_MAX;
		}
	}

#ifdef NEW_PRIOR
	A = (double *) malloc(sizeof(double)*L);
	B = (double *) malloc(sizeof(double)*L);
	for (int l = 0; l < L; l++ ){
		A[l] = 1;
		B[l] = 1;
	}
	A[56] = B[56] = 0; //Internet
	B[1] = B[2] = B[21] = B[67] = B[26] = B[69] = B[113] = B[76] = B[87] = B[54] = 3;
	B[68] = B[71] = B[74] = B[81] = B[82] = B[84] = B[83] = B[73] = B[64] = B[112] = B[116] = B[120] = B[116] = B[120] = B[116] = B[114] = B[65] = B[47] = B[19] = B[20] = 2;
#endif
}

void train() {
	/* initialize output */
	printf("init train\n");
	initTrain();

	/* update theta */
#pragma omp parallel shared(N,L,M,mytheta,R)
	{
#pragma omp for schedule(dynamic,1)
		for (int m = 0; m < M; m++) {
			for (int l = 0; l < L; l++) {
				double rR = 0;
				double sum_r = 0;
#ifdef PRIOR
				rR += A;
				sum_r += A + B;
#endif
#ifdef NEW_PRIOR
				if (A[l] == 0) continue;
				rR += A[l];
				sum_r += A[l] + B[l];
#endif
				for (int n = 0; n < N[m]; n++) {
					rR += R[m][n][l];
					sum_r += 1;
				}
				mytheta[m][l] = rR / sum_r;
				if (EQUAL(rR,0.0)) {
					mytheta[m][l] = 0;
				}
				if (mytheta[m][l] < 0 || mytheta[m][l] > 1 || mytheta[m][l] != mytheta[m][l]) {
					printf("error %lf %lf\n", rR, sum_r);
				}
			}
		}
	}
#ifdef DEBUG
		printf("theta:\n");
		for (int m = 0; m < M; m++) {
			for (int l = 0; l < L; l++) {
				printf("%lf ", mytheta[m][l]);
			}
			printf("\n");
		}
#endif
}


void save(const char *output) {
	char filename[100];
	FILE *fp;

	strcpy(filename, output);
	strcat(filename, "/theta");
	fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fprintf(fp, "%d %d\n", M, L);
	for (int m = 0; m < M; m++) {
		for (int l = 0; l < L; l++) {
			fprintf(fp, "%e ", mytheta[m][l]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

#ifdef PRIOR
	strcpy(filename, output);
	strcat(filename, "AB");
	fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fprintf(fp, "%lf %lf\n", A, B);
	fclose(fp);
#endif
	/*
	strcpy(filename, output);
	strcat(filename, "/loglikelihood");
	fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fprintf(fp, "%e\n", loglikelihood());
	fclose(fp);
	*/
}
int main(int argc, char **argv) {
	if (argc != 3) {
		printf("Usage: %s input output\n", argv[0]);
		return 0;
	}

	/* craete output folder */
	char command[1000];
	strcpy(command, "mkdir ");
	strcat(command, argv[2]);
	system(command);
	/* set seed */
	srand(0);
	//printf("load data\n");
	loadData(argv[1]);
	train();
	save(argv[2]);
	return 0;
}
