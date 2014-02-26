#include "train.h"

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

void modelprob(char *data, char *output) {
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

	for (int i = 0; i < NN; i++){
		int m;
		fscanf(fp, "%d", &m); // category
		int sumR = 0;
		for (int j = 0; j < L; j++) {
			fscanf(fp, "%d", &newR[j]);
			sumR +=  newR[j];
		}
		double sum_beta = 0;
		for (int k = 0; k < K; k++) {
			sum_beta += mybeta[m][k];
		}
		double p = 0;
		for (int k = 0; k < K; k++) {
			double q = 1;
			for (int l = 0; l < L; l++) {
				if (newR[l]) {
					q *= mytheta[k][l];
				} else {
					q *= (1 - mytheta[k][l]);
				}
			}
			q *= mybeta[m][k]/sum_beta;
			p += q;
		}
		fprintf(ofp, "%e %d\n", p, sumR);
	}

	free(newR);
	fclose(fp);
}

int main(int argc, char **argv) {
	if (argc != 4) {
		printf("Usage: estimate input data output\n");
		return 0;
	}
	printf("load\n");
	load(argv[1]);
	modelprob(argv[2], argv[3]);
	return 0;
}
