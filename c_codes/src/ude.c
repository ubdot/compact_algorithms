#include "ude.h"

double ude_mod(double *solution, double f, uint64_t cr, uint64_t budget, uint64_t num_exe, int dim, double (*feval)(double *, int)){
    int i, j, l, k, n_eval=0, elite=0;
    double *trial, *pop, *fpop;
    double ftrial;
    uint64_t i_smp[4];
    trial   = (double*)calloc(dim, sizeof(double));
    if(trial==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    fpop    = (double*)calloc(NPOP, sizeof(double));
    if(fpop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    pop     = (double*)calloc((NPOP*dim), sizeof(double));
    if(pop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }

    for(i=0; i<NPOP; i++){
        for(j=0; j<dim; j++){
            pop[i*dim+j] = 2.0f*doubleRand() - 1.0f; 
        }
        fpop[i] = feval(&pop[i*dim], dim);
        n_eval++;
        if(i!=elite && fpop[i] < fpop[elite]){
            elite = i;
        }
    }

    while(n_eval < budget){
        for(i=0; i<num_exe; i++){
            for(l=0; l<NPOP; l++){
                i_smp[0]    = l;
                randSample(i_smp, NPOP, 3);
                j   = (int)(next() % (uint64_t)dim);
                for(k=0; k<dim; k++){
                    //create mutant element if needed only
                    if(next() < cr || k==j){
                        trial[k]    = pop[i_smp[1]*dim + k] + f*(pop[i_smp[2]*dim + k] - pop[i_smp[3]*dim + k]);
                        //fix bounds
                        while(trial[k]>1.0){
                            trial[k]-=2.0;
                        }
                        while(trial[k]<-1.0){
                            trial[k]+=2.0;
                        }
                    }else{
                        trial[k]    = pop[l*dim + k];
                    }
                }
                ftrial  = feval(trial, dim);
                n_eval++;
                if(ftrial < fpop[l]){
                    for(k=0; k<dim; k++){
                        if(pop[l*dim+k]  != trial[k]){
                            pop[l*dim+k]  = trial[k];
                        }
                    }
                    fpop[l] = ftrial;
                    if(l != elite && fpop[l]<fpop[elite]){
                        elite = l;
                    }
                }
                if(n_eval >= budget){
                    break;
                }
            }
            if(n_eval >= budget){
                break;
            }
        }
        if(n_eval < budget){
            k = elite;
            for(i=0; i<NPOP; i++){
                if(i != elite){
                    for(j=0; j<dim; j++){
                        pop[i*dim+j] = 2.0*doubleRand() - 1.0; 
                    }
                    fpop[i] = feval(&pop[i*dim], dim);
                    n_eval++;
                    if(i!=k && fpop[i] < fpop[k]){
                        k = i;
                    }
                }
                if(n_eval >= budget){
                    break;
                }
            }
            elite = k;
        }
    }
    
    ftrial  = fpop[elite];
    for(i=0; i<dim; i++){
        solution[i] = pop[elite*dim + i];
    }
    
    cleanup:
    free(pop);
    free(fpop);
    free(trial);
    return ftrial;
}

double udeerm(double *solution, double f, double cr, int budget, int GR, int dim, double (*feval)(double *, int)){
    int i, j, k, m, nfevals = 0, G, elite=0, w1, w2;
    uint64_t i_smp[4];
    uint64_t uCR=(uint64_t)(((double)MAX_RAND) * cr);
    double *trial, *pop, *fpop, ftrial;
    
    pop     = (double*)calloc(dim*pop_size, sizeof(double));
    if(pop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    fpop    = (double*)calloc(pop_size, sizeof(double));
    if(fpop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    trial   = (double*)calloc(dim, sizeof(double));
    if(trial==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }

    for(i=0; i< pop_size; i++){
        for(j=0; j<dim; j++){
            pop[i*dim+j]    = 2.0*doubleRand()-1.0;
        }
        fpop[i] = feval(&pop[i*dim], dim); 
        nfevals++;
        if(i>0 && fpop[i] < fpop[elite]){
            elite=i;
        }
    }
    while (nfevals < budget){
        for(G=0; G<GR; G++){
            for(i=0; i< pop_size; i++){
                j   = (int)(next() % (uint64_t)dim);
                i_smp[0]    = i;
                randSample(i_smp, pop_size, 3);
                for(k=0; k<dim; k++){
                    if(j==k || next()<uCR){
                        trial[k] =  pop[i_smp[1]*dim + k] + f*(pop[i_smp[2]*dim + k] - pop[i_smp[3]*dim + k]);
                    }else{
                        trial[k] = pop[i*dim+k];
                    }
                    if(trial[k]<-1.0){
                        trial[k]=-1.0;
                    }else if(trial[k]>1.0){
                        trial[k]=1.0;
                    }
                }
                ftrial  = feval(trial, dim);
                nfevals++;
                if(ftrial < fpop[i]){
                    fpop[i] = ftrial;
                    for(k=0; k<dim; k++){
                        if(pop[i*dim+k]!=trial[k]){
                            pop[i*dim+k]=trial[k];
                        }
                    }
                    if(i!=elite && fpop[i] < fpop[elite]){
                        elite=i;
                    }
                }
                if(nfevals >= budget){
                    break;
                }
            }
            if(nfevals >= budget){
                break;
            }
        }
        if(nfevals >= budget){
            break;
        }
        w1=0;
        for(i=1; i<pop_size; i++){
            if(fpop[i]>fpop[w1]){
                w1=i;
            }
        }
        if(w1 != 0){
            w2=0;
        }else{
            w2=1;
        }
        for(i=0; i<pop_size; i++){
            if(i==w2 || i==w1) continue;
            if(fpop[i]>fpop[w2]){
                w2=i;
            }
        }
        for(k=0; k<dim; k++){
            trial[0] = pop[elite*dim+k] + f*(pop[w1*dim+k] - pop[w2*dim+k]);
            if(trial[0]<-1.0){
                trial[0]=-1.0;
            }else if(trial[0]>1.0){
                trial[0]=1.0;
            }
            trial[1] = pop[elite*dim+k] + f*(pop[w2*dim+k] - pop[w1*dim+k]);
            if(trial[1]<-1.0){
                trial[1]=-1.0;
            }else if(trial[1]>1.0){
                trial[1]=1.0;
            }
            
            pop[w1*dim+k]   = trial[0];
            pop[w2*dim+k]   = trial[1];
        }
        fpop[w1] = feval(&pop[w1*dim], dim);
        nfevals++;
        if(fpop[w1] < fpop[elite]){
            elite=w1;
        }            
        if(nfevals >= budget){
            break;
        }
        fpop[w2] = feval(&pop[w2*dim], dim);
        nfevals++;
        if(fpop[w2] < fpop[elite]){
            elite=w2;
        }            
        if(nfevals >= budget){
            break;
        }
    }
    ftrial = fpop[elite];
    for(i=0; i<dim; i++){
        solution[i] = pop[elite*dim + i];
    }

    cleanup:
    free(pop);
    free(trial);
    free(fpop);
    return ftrial;
}

double ude_fix(double *solution, double f, double cr, int budget, int GR, int dim, double (*feval)(double *, int)){
    int i, j, k, m, nfevals = 0, G, elite=0, w1, w2;
    uint64_t i_smp[4];
    uint64_t uCR=(uint64_t)(((double)MAX_RAND) * cr);
    double *trial, *pop, *fpop, ftrial;
    
    pop     = (double*)calloc(dim*pop_size, sizeof(double));
    if(pop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    fpop    = (double*)calloc(pop_size, sizeof(double));
    if(fpop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    trial   = (double*)calloc(dim, sizeof(double));
    if(trial==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }

    for(i=0; i< pop_size; i++){
        for(j=0; j<dim; j++){
            pop[i*dim+j]    = 2.0*doubleRand()-1.0;
        }
        fpop[i] = feval(&pop[i*dim], dim); 
        nfevals++;
        if(i>0 && fpop[i] < fpop[elite]){
            elite=i;
        }
    }
    while (nfevals < budget){
        for(G=0; G<GR; G++){
            for(i=0; i< pop_size; i++){
                j   = (int)(next() % (uint64_t)dim);
                i_smp[0]    = i;
                randSample(i_smp, pop_size, 3);
                for(k=0; k<dim; k++){
                    if(j==k || next()<uCR){
                        trial[k] =  pop[i_smp[1]*dim + k] + f*(pop[i_smp[2]*dim + k] - pop[i_smp[3]*dim + k]);
                    }else{
                        trial[k] = pop[i*dim+k];
                    }
                    if(trial[k]<-1.0){
                        trial[k]=-1.0;
                    }else if(trial[k]>1.0){
                        trial[k]=1.0;
                    }
                }
                ftrial  = feval(trial, dim);
                nfevals++;
                if(ftrial < fpop[i]){
                    fpop[i] = ftrial;
                    for(k=0; k<dim; k++){
                        if(pop[i*dim+k]!=trial[k]){
                            pop[i*dim+k]=trial[k];
                        }
                    }
                    if(i!=elite && fpop[i] < fpop[elite]){
                        elite=i;
                    }
                }
                if(nfevals >= budget){
                    break;
                }
            }
            if(nfevals >= budget){
                break;
            }
        }
        if(nfevals >= budget){
            break;
        }
        w1=0;
        for(i=1; i<pop_size; i++){
            if(fpop[i]>fpop[w1]){
                w1=i;
            }
        }
        if(w1 != 0){
            w2=0;
        }else{
            w2=1;
        }
        for(i=0; i<pop_size; i++){
            if(i==w2 || i==w1) continue;
            if(fpop[i]>fpop[w2]){
                w2=i;
            }
        }
        for(k=0; k<dim; k++){
            pop[w1*dim+k]   = 2.0*doubleRand()-1.0;
            pop[w2*dim+k]   = 2.0*doubleRand()-1.0;
        }
        fpop[w1] = feval(&pop[w1*dim], dim);
        nfevals++;
        if(fpop[w1] < fpop[elite]){
            elite=w1;
        }            
        if(nfevals >= budget){
            break;
        }
        fpop[w2] = feval(&pop[w2*dim], dim);
        nfevals++;
        if(fpop[w2] < fpop[elite]){
            elite=w2;
        }            
        if(nfevals >= budget){
            break;
        }
    }
    ftrial = fpop[elite];
    for(i=0; i<dim; i++){
        solution[i] = pop[elite*dim + i];
    }
    
    cleanup:
    free(pop);
    free(trial);
    free(fpop);
    pop    = NULL;
    trial   = NULL;
    fpop   = NULL;
    return ftrial;
}