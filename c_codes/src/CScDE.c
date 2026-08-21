#include "CScDE.h"

double CScDE(double *elite, int budget, double np, uint64_t nWaves, double freq, uint64_t dim, double (*feval)(double *, int)){
    int i, j, k, m;
    double *mu_v, *std_v, *trial, *pop;
    double ftrial, felite, F_it, CR_it, v_sin;
    
    //Allocating memory for temporary variables
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

    //initialize PV
    for(i=0; i<dim; i++){
        mu_v[i]     = 0.0f;
        std_v[i]    = 10.0f;
    }
    //evaluate elite solution
    sampleSolution(mu_v[0], std_v[0], elite, -3, -1, dim);
    felite = feval(elite, dim);

    for(i=0; i<budget-1; i++){
        //Calculate crossover probability and scaling factor according the paper
        F_it    = 0.0f;
        CR_it   = 0.0f;
        for(j=1; j<=nWaves; j++){
            v_sin   = (0.5f*sinf(PI_2*freq*(double)i/(double)j) + 1.0f);
            F_it    = F_it + v_sin;
            CR_it   = CR_it + v_sin;
        }
        CR_it   = 0.7f + 0.1f*(1.0f/(double)nWaves)*CR_it;
        F_it    = (1.0f/nWaves)*F_it;

        j   = (int)(next() % (uint64_t)dim);
        k   = j;
        do{
            //mutant element
            sampleSolution(mu_v[k], std_v[k], pop, -3, 0, 1);
            trial[k]    = pop[0] + F_it*(pop[1] - pop[2]);
            //fix bounds
            while(trial[k]>1){
                trial[k]-=2.0f;
            }
            while(trial[k]<-1){
                trial[k]+=2.0f;
            }
            k++;
            if(k >= dim){
                k   = 0;
            }
        }while(doubleRand() < CR_it && k != j);
        m   = k;
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
            updatePv(mu_v, std_v, trial, elite, np, dim);
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
            updatePv(mu_v, std_v, elite, trial, np, dim);
        }
    }
    cleanup:
    free(pop);
    free(mu_v);
    free(std_v);
    free(trial);
    pop     = NULL;
    mu_v    = NULL;
    std_v   = NULL;
    trial   = NULL;
    return felite;
}