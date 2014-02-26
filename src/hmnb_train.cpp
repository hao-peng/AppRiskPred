/**********
Training for HMNB model
Author: Hao Peng (penghATpurdue.edu)
Last updated time: Feb 25, 2014
**********/

#include "model.h"

int K = 0;
int M = 0;
int L = 0;
int *N = NULL;
int ***R = NULL;

#ifdef NEW_PRIOR
// this is for heuristic initialization for the priors
double *A = 0;
double *B = 0;
#else
double A = 1.0;
double B = 1.0;
#endif

double *myalpha = NULL;
double **mytheta = NULL;
double **mybeta = NULL;
double ***r = NULL;

int mode = UPDATE_ALPHA;

void initTrain() {
	initTheta();
	mybeta = (double **) malloc(sizeof(double*)*M);
	for (int m = 0; m < M; m++) {
		mybeta[m] = (double*) malloc(sizeof(double)*K);
		for (int k = 0; k < K; k++) {
			mybeta[m][k] = 1;
		}
	}
	myalpha = (double *) malloc(sizeof(double)*K);
	for (int k = 0; k < K; k++) {
		myalpha[k] = 1.0/M;
	}
	r = (double ***) malloc(sizeof(double**)*M);
	for (int m = 0; m < M; m++) {
		r[m] = (double**) malloc(sizeof(double*)*N[m]);
		for (int n = 0; n < N[m]; n++) {
			r[m][n] = (double*) malloc(sizeof(double)*K);
			for (int k = 0; k < K; k++) {
				r[m][n][k] = 0;
			}
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
	A[56] = B[56] = 0; //Internet do not use this
	B[1] = B[2] = B[21] = B[67] = B[26] = B[69] = B[113] = B[76] = B[87] = B[54] = 2*total; // high prior for this
	B[68] = B[71] = B[74] = B[81] = B[82] = B[84] = B[83] = B[73] = B[64] = B[112] = B[116] = B[120] = B[116] = B[120] = B[116] = B[114] = B[65] = B[47] = B[19] = B[20] = total; // low prior for this
#endif

}

void train() {
	/* initialize output */
	printf("init train\n");
	initTrain();

	/* initialize temp variables */
	double *myalpha_new = (double *) malloc(sizeof(double)*K);
	double *psi_sum_beta = (double *) malloc(sizeof(double)*M);
	double *psi_myalpha = (double *) malloc(sizeof(double)*K);
	double **log_myrho = (double **) malloc(sizeof(double*)*M);
	double **psi_mybeta = (double **) malloc(sizeof(double*)*M);
	for (int m = 0; m < M; m++) {
		log_myrho[m] = (double *) malloc(sizeof(double)*K);
		psi_mybeta[m] = (double *) malloc(sizeof(double)*K);
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

	double maxDiff = 0;

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
//		for (int in_iter = 0; in_iter < IN_LOOP; in_iter++) {
			//printf("in iter: %d\n", in_iter);
#pragma omp parallel shared(M,N,K,L,mybeta,psi_mybeta,log_myrho,log_mytheta,log_inv_mytheta,r)
			{
#pragma omp for schedule(dynamic,1)
				for (int m = 0; m < M; m++) {
					/* computer r */
					double sum_beta = 0;
					for (int k = 0; k < K; k++) {
						sum_beta += mybeta[m][k];
						psi_mybeta[m][k] = DiGamma_Function(mybeta[m][k]);
					}
					psi_sum_beta[m] = DiGamma_Function(sum_beta);
					for (int n = 0; n < N[m]; n++) {
						for (int k = 0; k < K; k++) {
							log_myrho[m][k] = psi_mybeta[m][k]-psi_sum_beta[m];
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
						double log_sum_rho = logsumexp(log_myrho[m], K);
						for (int k = 0; k < K; k++) {
							r[m][n][k] = exp(log_myrho[m][k] - log_sum_rho);
						}
					}

					/* compute mybeta */
					for (int k = 0; k < K; k++) {
						mybeta[m][k] = myalpha[k];
						for (int n = 0; n < N[m]; n++) {
							mybeta[m][k] = mybeta[m][k] + r[m][n][k];
						}
					}
				}
			}
			/*
			printf("beta:\n");
			for (int m = 0; m < M; m++) {
				for (int k = 0; k < K; k++) {
					printf("%lf ", mybeta[m][k]);
				}
				printf("\n");
			}
			*/
//		}
#ifdef DEBUG
		printf("beta:\n");
		for (int m = 0; m < M; m++) {
			for (int k = 0; k < K; k++) {
				printf("%lf ", mybeta[m][k]);
			}
			printf("\n");
		}
#endif
		/* m-step */
		if (out_iter != OUT_LOOP - 1) {
			/* update alpha */
			if (mode == UPDATE_ALPHA) {
				for (int m = 0; m < M; m++) {
					double sum_beta = 0;
					for (int k = 0; k < K; k++) {
						sum_beta += mybeta[m][k];
						psi_mybeta[m][k] = DiGamma_Function(mybeta[m][k]);
					}
					psi_sum_beta[m] = DiGamma_Function(sum_beta);
				}
				int converge = 0;
				for (int iter = 0; iter < 1000; iter++) {
					double sum_alpha = 0;
					for (int k = 0; k < K; k++) {
						sum_alpha += myalpha[k];
						psi_myalpha[k] = DiGamma_Function(myalpha[k]);
					}
					double psi_sum_alpha = DiGamma_Function(sum_alpha);
					int fault;
					for (int k = 0; k < K; k++) {
						g[k] = M * (psi_sum_alpha - psi_myalpha[k]);
						for (int m = 0; m < M; m++) {
							g[k] += psi_mybeta[m][k] - psi_sum_beta[m];
						}
						q[k] = -M * trigamma(myalpha[k], &fault);
					}
					double z = M * trigamma(sum_alpha, &fault);
					double gq = 0;
					double rq = 0;
					for (int k = 0; k < K; k++) {
						gq = gq + g[k] / q[k];
						rq = rq + 1 / q[k];
					}
					double b = gq / (1 / z + rq);
					for (int k = 0; k < K; k++) {
						myalpha_new[k] = myalpha[k] - (g[k] - b) / q[k];
						if (myalpha_new[k] < 0) {
							printf("warning alpha small than zero\n");
						}
					}
#ifdef DEBUG
					printf("alpha:\n");
					for (int k = 0; k < K; k++) {
						printf("%lf ", myalpha[k]);
					}
					printf("\n");
#endif

					converge = 1; 
					for (int k = 0; k < K; k++) {
						double diff = myalpha_new[k] - myalpha[k];
						if (diff > 1e-6 || diff < -1e-6) {
							converge = 0;
							break;
						}
					}
					if (converge) {
						break;
					}

					double *tmpalpha = myalpha;
					myalpha = myalpha_new;
					myalpha_new = tmpalpha;
				}
				if (!converge) {
					printf("warning: not converge\n");
				}
			}

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
								rR += r[m][n][k]*R[m][n][l];
								sum_r += r[m][n][k];
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
				printf("Finished.\n");
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
		free(psi_mybeta[m]);
		free(log_myrho[m]);
	}
	free(psi_mybeta);
	free(log_myrho);
	free(psi_sum_beta);
	free(psi_myalpha);
	free(myalpha_new);
}


void save(const char *output) {
	char filename[100];
	FILE *fp;

	strcpy(filename, output);
	strcat(filename, "/beta");
	fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fprintf(fp, "%d %d\n", M, K);
	for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {
			fprintf(fp, "%e ", mybeta[m][k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	strcpy(filename, output);
	strcat(filename, "/r");
	fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fprintf(fp, "%d %d\n", M, K);
	for (int m = 0; m < M; m++) {
		fprintf(fp, "%d\n", N[m]);
		for (int n = 0; n < N[m]; n++) {
			for (int k = 0; k < K; k++) {
				fprintf(fp, "%e ", r[m][n][k]);
			}
			fprintf(fp, "\n");
		}
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

	strcpy(filename, output);
	strcat(filename, "/alpha");
	fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("cannot open file\n");
		exit(-1);
	}
	fprintf(fp, "%d\n", K);
	for (int k = 0; k < K; k++) {
		fprintf(fp, "%e ", myalpha[k]);
	}
	fprintf(fp, "\n");
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
	if (argc == 4) {
		K = 2;
	} else if (argc == 5) {
		K = atoi(argv[4]);
	} else {
		printf("Usage: %s input output mode [K]\n mode 0: update alpha\n mode 1: fix alpha\n", argv[0]);
		printf("update alpha: Hierachical Mixture of Naive Bayes (HMNB)\n");
		printf("fix alpha: Mixture of Naive Bayes (MNB)\n");
		return 0;
	}
	mode = atoi(argv[3]);

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
