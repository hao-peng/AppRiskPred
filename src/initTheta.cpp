/**********
function to initial theta with kmeans
Author: Hao Peng (penghATpurdue.edu)
Last updated time: Feb 25, 2014
**********/

#include "model.h"

void initTheta() {
	mytheta = (double **) malloc(sizeof(double*)*K);
	for (int k = 0; k < K; k++) {
		mytheta[k] = (double*) malloc(sizeof(double)*L);
	}

	/* k mean */
	double **means = (double **) malloc(sizeof(double*)*K);
	for (int k = 0; k < K; k++) {
		means[k] = (double *) malloc(sizeof(double)*L);
	}
	for (int k = 0; k < K; k++) {
		int m = rand()%M;
		int n = rand()%N[m];
		for (int l = 0; l < L; l++) {
			means[k][l] = R[m][n][l];
		}
	}

	int **cluster = (int **) malloc(sizeof(int*)*M);
	for (int m = 0; m < M; m++) {
		cluster[m] = (int *) malloc(sizeof(int)*N[m]);
	}

	double *dists = (double *) malloc(sizeof(double)*M);

	int *count = (int *) malloc(sizeof(int)*K);
	double currTotal = 0, lastTotal = 0;
	for (int iter = 0; iter < 1000; iter++) {
		currTotal = 0;
#pragma omp parallel shared(cluster, means, dists) 
		{
#pragma omp for schedule(dynamic,1)
		for (int m = 0; m < M; m++) {
			dists[m] = 0;
			for (int n = 0; n < N[m]; n++) {
				cluster[m][n] = 0;
				double minDist = DBL_MAX;
				for (int k = 0; k < K; k++) {
					double dist = 0;
					for (int l = 0; l < L; l++) {
						double tmp = means[k][l] - R[m][n][l];
						dist += tmp * tmp;
					}
					if (dist < minDist)  {
						minDist = dist;
						cluster[m][n] = k;
					}
				}
				dists[m] += minDist;
			}
		}
		}
		for (int m = 0; m < M; m++) {
			currTotal += dists[m];
		}

		for (int k = 0; k < K; k++) {
			for (int l = 0; l < L; l++) {
				means[k][l] = 0;
			}
			count[k] = 0;
		}
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < N[m]; n++) {
				int k = cluster[m][n];
				for (int l = 0; l < L; l++) {
					means[k][l] += R[m][n][l];
				}
				count[k]++;
			}
		}
		for (int k = 0; k < K; k++) {
			for (int l = 0; l < L; l++) {
				means[k][l] /= count[k];
			}
		}

		//printf("%lf\n", currTotal);
		/*
		if (EQUAL(currTotal, lastTotal)) {
			printf("equal\n");
			break;
		}
		*/
		lastTotal = currTotal;
	}

	for (int k = 0; k < K; k++) {
		for (int l = 0; l < L; l++) {
			mytheta[k][l] = 1; //normalize
		}
		count[k]++;
	}
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N[m]; n++) {
			int k = cluster[m][n];
			for (int l = 0; l < L; l++) {
				mytheta[k][l] += R[m][n][l];
			}
		}
	}
	for (int k = 0; k < K; k++) {
		for (int l = 0; l < L; l++) {
			mytheta[k][l] /= count[k];
		}
	}

	/*
	printf("theta:\n");
	for (int k = 0; k < K; k++) {
		for (int l = 0; l < L; l++) {
			printf("%lf ", mytheta[k][l]);
		}
		printf("\n");
	}
	*/

	free(count);
	for (int m = 0; m < M; m++) {
		free(cluster[m]);
	}
	free(cluster);
	for (int k = 0; k < K; k++) {
		free(means[k]);
	}
	free(means);

}
