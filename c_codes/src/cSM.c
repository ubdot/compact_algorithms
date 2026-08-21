#include "cSM.h"
double cSM(double *elite, uint64_t budget, double Np, double B, uint64_t dim, double (*feval)(double *, int)){
    uint64_t i;
    double *mu_v, *std_v, *lbest, *xoff;
    double fittele, fittlb, fittxoff;
    uint64_t budget1, totFE=0, totFE_b=0;
    //Allocate memory
    mu_v    = (double*)calloc(dim, sizeof(double));
    if(mu_v==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    std_v   = (double*)calloc(dim, sizeof(double));
    if(std_v==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    lbest   = (double*)calloc(dim, sizeof(double));
    if(lbest==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    xoff    = (double*)calloc(dim, sizeof(double));
    if(xoff==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    //Calculate number of function evaluations according to the budget
    budget1 = (uint64_t)(0.2 * (double)budget);
    //Sample initial solutions elite and local best
    sampleSolution(0.0, 10.0, elite, -3, -1, dim);
    sampleSolution(0.0, 10.0, lbest, -3, -1, dim);
    //init the PV for the compact module
    for(i=0; i<dim; i++){
        std_v[i]    = 10.0;
        xoff[i]     = elite[i];
    }

    fittele = feval(elite, dim);
    fittlb  = feval(lbest, dim);
    totFE   += 2;

    //Single solution approach
    while (totFE < budget1)
    {
        //Snum operator over the solution
        i           = SNUM(xoff, budget1, totFE, B, dim);
        fittxoff    = feval(xoff, dim);
        totFE++;
        //update solution
        if(fittxoff <= fittele){
            fittele     = fittxoff;
            elite[i]    = xoff[i];
        }
    }
    //Pass information to the compact module
    for(i=0; i<dim; i++){
        mu_v[i]     = -elite[i];
    }
    //Value used into the MNUM operator
    budget1 = budget - budget1;
    //Until the maximun number of function evalutaions is reached execute compact module
    while (totFE < budget)
    {
        for(i=0; i<dim; i++){
            xoff[i] = singleSampleSolution(mu_v[i], std_v[i]);
        }
        MNUM(xoff, budget1, totFE_b, B, dim); //Mnum operator on the best solution
        fittxoff    = feval(xoff, dim);
        totFE_b++;
        totFE++;
        //Update PV model according the results
        if(fittxoff <= fittlb){
            fittlb     = fittxoff;
            for(i=0; i<dim; i++){
                lbest[i]    = xoff[i];
            }
            if(totFE < budget){
                updatePvB(mu_v, std_v, xoff, lbest, Np, dim);
            }
        }else{
            if(totFE < budget){
                updatePvB(mu_v, std_v, lbest, xoff, Np, dim);
            }
        }
    }
    if(fittlb <= fittele){
        fittele = fittlb;
        for(i=0; i<dim; i++){
            elite[i]    = lbest[i];
        }
    }

    cleanup:
    free(mu_v);
    free(std_v);
    free(xoff);
    free(lbest);
    return fittele;
}
