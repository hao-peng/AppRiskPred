/**********
header files for HMNB model
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

#ifdef PRIOR
extern double A;
extern double B;
#endif

#ifdef NEW_PRIOR
extern double *A;
extern double *B;
#endif

extern int K;
extern int M;
extern int L;
extern int *N;
extern int ***R;

extern double *myalpha;
extern double **mytheta;
extern double **mybeta;
extern double ***r;


int loadData(char *input);

void initTheta();

//double loglikelihood();

double DiGamma_Function( double x );

double Ln_Gamma_Function(double x);

double logsumexp(double nums[], int ct);
