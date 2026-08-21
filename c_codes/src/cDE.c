#include "cDE.h"
double cDE(double *elite, double f, uint64_t cr, int budget, double np, int dim, double (*feval)(double *, int)){
    int i, j, k, m;
    double *mu_v, *std_v, *trial, *pop;
    double ftrial, felite;

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
    pop     = (double*)calloc(3, sizeof(double));
    if(pop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }

    //Initialize memory of the PV
    for(i=0; i<dim; i++){
        mu_v[i]     = 0.0f;
        std_v[i]    = 10.0f;
    }

    //sample initial solution
    sampleSolution(mu_v[0], std_v[0], elite, -3, -1, dim);
    felite = feval(elite, dim);

    //Differential evolution on the elite
    for(i=0; i<budget-1; i++){
        j   = (int)(next() % (uint64_t)dim);
        k   = j;
        //Exponential crossover with mutant vector, created from the PV
        do{
            //mutant element
            sampleSolution(mu_v[k], std_v[k], pop, -3, 0, 1);
            trial[k]    = pop[0] + f*(pop[1] - pop[2]);
            //repair bounds
            while(trial[k]>1){
                trial[k]-=2;
            }
            while(trial[k]<-1){
                trial[k]+=2;
            }
            //if dim is reached, continues on the first element 
            k++;
            if(k >= dim){
                k   = 0;
            }
        }while(next() < cr && k != j);
        m   = k;
        //Copy the rest of the elements on the trial solution
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
            //copy the modified solution if improved
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
    free(pop);
    free(mu_v);
    free(std_v);
    free(trial);
    return felite;
}