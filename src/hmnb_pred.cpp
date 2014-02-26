/**********
Prediction for HMNB model
Author: Hao Peng (penghATpurdue.edu)
Last updated time: Feb 25, 2014
**********/

#include "model.h"

#define MAX_ESTIMATE 0
#define MEAN_ESTIMATE 1
#define SUM_ESTIMATE 2


int K = 0;
int M = 0;
int L = 0;
int *N = NULL;
int ***R = NULL;

double *myalpha = NULL;
double **mytheta = NULL;
double **mybeta = NULL;
double ***r = NULL;
double log_likelihood = 0;
int method = MAX_ESTIMATE;

void load(const char *input) {
	char filename[100];
	FILE *fp;

	strcpy(filename, input);
	strcat(filename, "/beta");
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fscanf(fp, "%d %d\n", &M, &K);
	mybeta = (double **) malloc(sizeof(double*)*M);
	for (int m = 0; m < M; m++) {
		mybeta[m] = (double *) malloc(sizeof(double)*K);
		for (int k = 0; k < K; k++) {
			fscanf(fp, "%lf", &mybeta[m][k]);
		}
	}
	fclose(fp);

	strcpy(filename, input);
	strcat(filename, "/r");
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fscanf(fp, "%d %d\n", &M, &K);
	r = (double ***) malloc(sizeof(double**)*M);
	N = (int *) malloc(sizeof(int)*M);
	for (int m = 0; m < M; m++) {
		fscanf(fp, "%d\n", &N[m]);
		r[m] = (double **) malloc(sizeof(double*)*N[m]);
		for (int n = 0; n < N[m]; n++) {
			r[m][n] = (double *) malloc(sizeof(double)*K);
			for (int k = 0; k < K; k++) {
				fscanf(fp, "%lf", &r[m][n][k]);
			}
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

	strcpy(filename, input);
	strcat(filename, "/alpha");
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fscanf(fp, "%d", &K);
	myalpha = (double *) malloc(sizeof(double)*K);
	for (int k = 0; k < K; k++) {
		fscanf(fp, "%lf", &myalpha[k]);
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
		if (method == MAX_ESTIMATE || method == SUM_ESTIMATE) {
#pragma omp parallel shared(M,N,K,L,p,mybeta,newR)
			{
#pragma omp for schedule(dynamic,1)
				for (int m = 0; m < M; m++) {
					double sum_beta = 0;
					for (int k = 0; k < K; k++) {
						sum_beta += mybeta[m][k];
					}
					p[m] = 0;
					for (int k = 0; k < K; k++) {
						double q = 1;
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
						q *= mybeta[m][k]/sum_beta;
						p[m] += q;
					}
				}
			}
			if (method == MAX_ESTIMATE) {
				maxP = 0;
				for (int m = 0; m < M; m++) {
					if (p[m] > maxP) {
						maxP = p[m];
					}
				}
			} else {
				maxP = 0;
				for (int m = 0; m < M; m++) {
					maxP += p[m];
				}
			}
		} else {
			/* method two */
			double sum_alpha = 0;
			for (int k = 0; k < K; k++) {
				sum_alpha += myalpha[k];
			}
			maxP = 0;
			for (int k = 0; k < K; k++) {
				double q = 1;
				for (int l = 0; l < L; l++) {
					if (newR[l]) {
						q *= mytheta[k][l];
					} else {
						q *= (1 - mytheta[k][l]);
					}
				}
				q *= myalpha[k]/sum_alpha;
				maxP += q;
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
		printf("Usage: %s input data output [method]\nmethod 0: max\nmethod 1: mean\nmode 2: sum", argv[0]);
		printf("method: different evaluation stretages for multiple groups\n");
		return 0;
	}
	if (argc == 5) {
		method = atoi(argv[4]);
	}

	//printf("load\n");
	load(argv[1]);
	estimate(argv[2], argv[3]);
	return 0;
}
