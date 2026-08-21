#include "cPSO.h"
double cPSO(double *xgb, double *c, int budget, double np, int dim, double (*feval)(double *, int)){
    int i, j, k, m;
    double fxgb, fxt, fxlb, r1, r2;
    double *mu_v, *std_v, *xlb, *xt, *vt;
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
    xlb     = (double*)calloc(dim, sizeof(double));
    if(xlb==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    xt      = (double*)calloc(dim, sizeof(double));
    if(xt==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    vt      = (double*)calloc(dim, sizeof(double));
    if(vt==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }

    //Initialize PV
    for(i=0; i< dim; i++){
        mu_v[i] = 0.0f;
        std_v[i]= 10.0f;
    }

    //Initialize solution, global best and velocity
    sampleSolution(mu_v[0], std_v[0], xt,  -3, -1, dim);
    sampleSolution(mu_v[0], std_v[0], vt,  -3, -1, dim);
    sampleSolution(mu_v[0], std_v[0], xgb, -3, -1, dim);
    
    fxgb =  feval(xgb, dim);
    i    = 1;
    while(i < budget){
        for(j=0; j<dim; j++){
            //Sample local best
            xlb[j]  = singleSampleSolution(mu_v[j],std_v[j]);
            //Update velocity
            r1      = doubleRand();
            r2      = doubleRand();
            vt[j]   = c[0]*vt[j] + c[1]*r2*(xlb[j]-xt[j]) + c[2]*r1*(xgb[j]-xt[j]); 
            if(vt[j] > 1.0f){
                vt[j] = 1.0f;
            }else if(vt[j] < -1.0f){
                vt[j] = -1.0f;
            }
            //update position
            xt[j]   = xt[j] + vt[j];
            //repair boundaries
            while(xt[j] > 1.0f){
                xt[j] -= 2.0f;
            }
            while(xt[j] < -1.0f){
                xt[j] += 2.0f;
            }
        }
        fxt =  feval(xt, dim);
        i++;
        //update and competition
        if(i < budget){
            fxlb=  feval(xlb, dim);
            i++;
            if(fxt < fxlb){
                updatePv(mu_v, std_v, xt, xlb, np, dim);
            }else{
                updatePv(mu_v, std_v, xlb, xt, np, dim);
            }
        }

        if(fxt < fxgb){
            fxgb = fxt;
            for(j=0; j<dim; j++){
                xgb[j] = xt[j];
            }
        }
    }

    //clenup memory
    cleanup:
    free(xlb);
    free(xt);
    free(vt);
    free(mu_v);
    free(std_v);
    return fxgb;
}