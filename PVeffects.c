#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <stdlib.h>
#include "cDE.h"
#define UsedDimension 10
#define SVPLOT  0

void cec17_test_func(double *, double *,int,int,int);
void init_problem(int nx, int func_num);
double fast_pow(double x, int p);
double feval(double *x, int len);
double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;
int GNVars;
int MaxFEval = 0;
double fopt[1];
int sel_fun;


double sphere(double *xx, int d) {
    double sum = 0.0f;

    for (int i = 0; i < d; i++) {
        double xi = (100.0F)*xx[i];
        sum += xi * xi;
    }

    return sum;
}

double feval(double *x, int len){
    double f, x_min, x_max, xr, x_scaled[UsedDimension];
    int i=0;
    //Set bound constrains to evaluate objective function
    x_min   = -100;
    x_max   = 100;
    xr  = x_max - x_min;
    for(i=0; i<len; i++){
        x_scaled[i]=((x[i] + 1.0)*xr)/2.0 + x_min;
    }
    cec17_test_func(x_scaled, &f, len, 1, sel_fun);
    return f;
}


int main(void){
    double *xsol, result;
    int dim =UsedDimension;
    
    char flie_name[50];
    uint64_t seed = 0xF1D02A87BCB12D29;

    
    set_seed(seed);

    xsol= (double*)calloc(dim, sizeof(double));
    sel_fun = 4;
    printf("Problem_%d, ", sel_fun);
    init_problem(dim, sel_fun);
    double np   = 1.0/300.0;//1.0/82.0;//
    uint64_t cr = MAX_RAND*(uint64_t)floor(pow(0.5,1/(0.25*((double)dim))));
    result  = cDE(xsol, 0.5, cr, 3*dim*1e4, np, dim, feval);
    printf("%.4e\n", result);
    free(xsol);
    xsol=NULL;
    return 0;
}