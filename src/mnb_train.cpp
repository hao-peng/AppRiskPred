/**********
Training for MNB model
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

//#define IN_LOOP 5000
#define OUT_LOOP 8000
#define EPS 1e-10
#define EQUAL(a,b) ((a-b)<EPS && a-b>-EPS)

#define UPDATE_ALPHA 0
#define FIX_ALPHA 1

//#define DEBUG
#define PRIOR
//#define NEW_PRIOR

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
double **r = NULL;

void initTrain() {
	initTheta();
	r = (double **) malloc(sizeof(double*)*M);
	for (int m = 0; m < M; m++) {
		r[m] = (double*) malloc(sizeof(double)*K);
		for (int k = 0; k < K; k++) {
			r[m][k] = 0;
		}
	}
#ifdef NEW_PRIOR
	int total = 0;
	for (int m = 0; m < M; m++) {
		total += N[m];
	}
	A = (double *) malloc(sizeof(double)*L);
	B = (double *) malloc(sizeof(double)*L);
	for (int l = 0; l < L; l++ ){
		A[l] = 1;
		B[l] = 1;
	}
	A[56] = B[56] = 0; //Internet
	B[1] = B[2] = B[21] = B[67] = B[26] = B[69] = B[113] = B[76] = B[87] = B[54] = 2*total;
	B[68] = B[71] = B[74] = B[81] = B[82] = B[84] = B[83] = B[73] = B[64] = B[112] = B[116] = B[120] = B[116] = B[120] = B[116] = B[114] = B[65] = B[47] = B[19] = B[20] = total;
#endif

}

void train() {
	/* initialize output */
	printf("init train\n");
	initTrain();

	/* initialize temp variables */
	double maxDiff = 0;
	double **log_myrho = (double **) malloc(sizeof(double*)*M);
	for (int m = 0; m < M; m++) {
		log_myrho[m] = (double *) malloc(sizeof(double)*K);
	}
	double **old_mytheta = (double **) malloc(sizeof(double*)*K);
	double **log_mytheta = (double **) malloc(sizeof(double*)*K);
	double **log_inv_mytheta = (double **) malloc(sizeof(double*)*K);
	for (int k = 0; k < K; k++) {
		old_mytheta[k] = (double *) malloc(sizeof(double)*L);
		for (int l = 0; l < L; l++) old_mytheta[k][l] = 0;
		log_mytheta[k] = (double *) malloc(sizeof(double)*L);
		log_inv_mytheta[k] = (double *) malloc(sizeof(double)*L);
	}
	double *g = (double *) malloc(sizeof(double)*K);
	double *q = (double *) malloc(sizeof(double)*K);

	for (int out_iter = 0; out_iter < OUT_LOOP; out_iter++) {
		if (out_iter % 100 == 0) printf("Iter: %d\n", out_iter);
		for (int k = 0; k < K; k++) {
			for (int l = 0; l < L; l++) {
#ifdef NEW_PRIOR
				if (A[l] == 0) continue;
#endif
				log_mytheta[k][l] = log(mytheta[k][l]);
				//printf("%lf ", log(mytheta[k][l]));
				log_inv_mytheta[k][l] = log(1-mytheta[k][l]);
			}
			//printf("\n");
		}
		/* e-step */
		//printf("in iter: %d\n", in_iter);
#pragma omp parallel shared(M,N,K,L,log_myrho,log_mytheta,log_inv_mytheta,r)
		{
#pragma omp for schedule(dynamic,3)
			for (int m = 0; m < M; m++) {
				/* computer r */
				for (int k = 0; k < K; k++) {
					log_myrho[m][k] = 0;
				}
				for (int n = 0; n < N[m]; n++) {
					for (int k = 0; k < K; k++) {
						for (int l = 0; l < L; l++) {
#ifdef NEW_PRIOR
							if (A[l] == 0) continue;
#endif
							if (R[m][n][l]) {
								log_myrho[m][k] += log_mytheta[k][l];
							} else {
								log_myrho[m][k] += log_inv_mytheta[k][l];
							}
						}
					}
				}
				double log_sum_rho;
				for (int k = 0; k < K; k++) {
					log_sum_rho = logsumexp(log_myrho[m], K);
				}
				for (int k = 0; k < K; k++) {
					r[m][k] = exp(log_myrho[m][k] - log_sum_rho);
				}
			}
		}
		/* m-step */
		/* update theta */
#pragma omp parallel shared(K,N,L,M,mytheta,r,R)
		{
#pragma omp for schedule(dynamic,1)
			for (int k = 0; k < K; k++) {
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

					for (int m = 0; m < M; m++) {
						for (int n = 0; n < N[m]; n++) {
							rR += r[m][k]*R[m][n][l];
							sum_r += r[m][k];
						}
					}
					mytheta[k][l] = rR / sum_r;
					if (EQUAL(rR,0.0)) {
						mytheta[k][l] = 0;
					}
					if (mytheta[k][l] < 0 || mytheta[k][l] > 1 || mytheta[k][l] != mytheta[k][l]) {
						printf("error %lf %lf\n", rR, sum_r);
					}
				}
			}
		}
		maxDiff = 0;
		for (int k = 0; k < K; k++ ){
			for (int l = 0; l < L; l++) {
#ifdef NEW_PRIOR
				if (A[l] == 0) continue;
#endif
				double diff = old_mytheta[k][l] - mytheta[k][l];
				if (diff > maxDiff) maxDiff = diff;
			   if	(-diff > maxDiff) maxDiff = -diff;
				old_mytheta[k][l] = mytheta[k][l];
			}
		}
		if (maxDiff < 1e-6) {
			printf("finished.\n");
			break;
		}

#ifdef DEBUG
		printf("theta:\n");
		for (int k = 0; k < K; k++) {
			for (int l = 0; l < L; l++) {
				printf("%lf ", mytheta[k][l]);
			}
			printf("\n");
		}
#endif
	}

	if (maxDiff >= 1e-3) {
		printf("not converg\n");
	}

	/* free temp variables */
	free(g);
	free(q);
	for (int k = 0; k < K; k++) {
		free(log_inv_mytheta[k]);
		free(log_mytheta[k]);
		free(old_mytheta[k]);
	}
	free(old_mytheta);
	free(log_inv_mytheta);
	free(log_mytheta);
	for (int m = 0; m < M; m++) {
		free(log_myrho[m]);
	}
	free(log_myrho);
}


void save(const char *output) {
	char filename[100];
	FILE *fp;

	strcpy(filename, output);
	strcat(filename, "/r");
	fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fprintf(fp, "%d %d\n", M, K);
	for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {
			fprintf(fp, "%e ", r[m][k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	strcpy(filename, output);
	strcat(filename, "/theta");
	fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fprintf(fp, "%d %d\n", K, L);
	for (int k = 0; k < K; k++) {
		for (int l = 0; l < L; l++) {
			fprintf(fp, "%e ", mytheta[k][l]);
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
	if (argc == 3) {
		K = 2;
	} else if (argc == 4) {
		K = atoi(argv[3]);
	} else {
		printf("Usage: %s input output [K]\n", argv[0]);
		return 0;
	}

	/* craete output folder */
	char command[1000];
	strcpy(command, "mkdir ");
	strcat(command, argv[2]);
	system(command);
	/* set seed */
	srand(0);
	printf("load data\n");
	loadData(argv[1]);
	train();
	save(argv[2]);
	return 0;
}
