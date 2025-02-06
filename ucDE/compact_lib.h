
#include <math.h>
#include <stdio.h>
#include "erf.h"
#include "mt64.h"


double randFrom(double num_min, double num_max);
double sampleSolution(double mu, double sigma);
double denorm(double x, double low_lim, double up_lim);
void updatePv(double *muV, double *stdV, double *winner, double *looser, double Np, int len);
void updatePvB(double *muV, double *stdV, double *winner, double *looser, double Np, int len);
void randSample(int *arr_smp, int n_pop, int n_smp);
