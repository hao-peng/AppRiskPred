/**********
Compute loglikelihood
Author: Hao Peng (penghATpurdue.edu)
Last updated time: Feb 25, 2014i
**********/
#include "model.h"

double loglikelihood() {
	double *log_lls = (double *) malloc(sizeof(double) * M);

	double sum_alpha = 0;
	for (int k = 0; k < K; k++) {
		sum_alpha += myalpha[k];
	}

	/*
	for (int k = 0; k < K; k++ ) {
		for (int l = 0; l < L; l++) {
			printf("%lf ", log(mytheta[k][l]));
		}
		printf("\n");
	}
	*/

#pragma omp parallel shared(M,N,K,L,log_lls,myalpha,mybeta,r,R,sum_alpha)
	{
#pragma omp for schedule(dynamic,1)
		for (int m = 0; m < M; m++) {
			double sum_beta = 0;;
			for (int k = 0; k < K; k++) {
				sum_beta += mybeta[m][k];
			}
			log_lls[m] = Ln_Gamma_Function(sum_alpha);
			for (int k = 0; k < K; k++) {
				log_lls[m] -= Ln_Gamma_Function(myalpha[k]) + (myalpha[k] - 1) * (DiGamma_Function(mybeta[m][k]) - DiGamma_Function(sum_beta));
			}
			for (int n = 0; n < N[m]; n++) {
				for (int k = 0; k < K; k++) {
					log_lls[m] += r[m][n][k] * (DiGamma_Function(mybeta[m][k]) - DiGamma_Function(sum_beta));
					for (int l = 0; l < L; l++) {
						if (R[m][n][l]) {
							//if (mytheta[k][l] < EPS) {
							//	log_lls[m] += r[m][n][k] * EPS;
							//} else {
								log_lls[m] += r[m][n][k] * log(mytheta[k][l]);
							//}
						} else {
							//if (1 - mytheta[k][l] < EPS) {
							//	log_lls[m] += r[m][n][k] * EPS;
							//} else {
								log_lls[m] += r[m][n][k] * log(1 - mytheta[k][l]);
							//}
						}
					}
				}
			}
			log_lls[m] -= Ln_Gamma_Function(sum_beta);
			for (int k = 0; k < K; k++) {
				log_lls[m] += Ln_Gamma_Function(mybeta[m][k]) - (mybeta[m][k] - 1) * (DiGamma_Function(mybeta[m][k]) - DiGamma_Function(sum_beta));
			}
			for (int n = 0; n < N[m]; n++) {
				for (int k = 0; k < K; k++) {
					log_lls[m] -= r[m][n][k] * log(r[m][n][k]);
				}
			}
		}
	}
	double total = 0;
	for (int m = 0; m < M; m++) {
		total += log_lls[m];
	}
	return total;
}
