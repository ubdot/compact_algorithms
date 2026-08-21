#include "cDE_light.h"
double cDE_light(double *elite, double f, double cr, int budget, double np, int dim, double (*feval)(double *, int)){
    int i, j, k, m;
    double *mu_v, *std_v, *trial, pop;
    double ftrial, felite;
    double logCr, scaled_F;
    int xoverL;
    
    //allocate memory
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
    trial   = (double*)calloc(dim, sizeof(double));
    if(trial==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }

    //Adjust used values, scaling factor and crossover
    scaled_F    = 1.0 + f*f;
    logCr       = logf(cr);
    //Init the PV
    for(i=0; i<dim; i++){
        mu_v[i]     = 0.0;
        std_v[i]    = 10.0;
    }
    //Sample initial solution
    sampleSolution(mu_v[0], std_v[0], elite, -3, -1, dim);
    felite = feval(elite, dim);
    //compact differential evolution light
    for(i=0; i<budget-1; i++){
        //Calculate elements to be replaced on trial solution
        xoverL  = (int)(round(log(doubleRand())/logCr));
        j   = (int)(next() % (uint64_t)dim);
        k   = j;
        m   = 0;
        //Exponential crossover with mutant vector, created from the PV
        do{
            //mutant element
            trial[k]    = singleSampleSolution(mu_v[k], scaled_F*std_v[k]);
            //repair bounds
            while(trial[k]>1.0){
                trial[k]-=2.0;
            }
            while(trial[k]<-1.0){
                trial[k]+=2.0;
            }
            k++;
            m++;
            if(k >= dim){
                k   = 0;
            }
            
        }while(m <= xoverL && k != j);
        m   = k;
        //Copy the rest of elements on the trial solution
        while(k != j){
            trial[k]    = elite[k];
            k++;
            if(k>=dim){
                k=0;
            }
        }
        
          
        //Competition
        ftrial  = feval(trial, dim);
        if(ftrial < felite){
            //update model
            updatePvB(mu_v, std_v, trial, elite, np, dim);
            //copy the modified solution
            while(j != m){
                elite[j]  = trial[j];
                j++;
                if(j>=dim){
                    j=0;
                }
            }
            felite = ftrial;
            
        }else{
            //update the model
            updatePvB(mu_v, std_v, elite, trial, np, dim);
        }
        

    }
    //clean up memory
    cleanup:
    free(mu_v);
    free(std_v);
    free(trial);
    return felite;
}