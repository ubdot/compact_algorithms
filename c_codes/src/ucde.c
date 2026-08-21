#include "ucde.h"
#include <stdio.h>

double ucde(double *solution, double f, uint64_t cr_bin, int num_exe, double Dubudget, double np, uint64_t dim, int maxFE, double (*feval)(double *, int)){
    double *mu_v, *std_v, *pop, *fpop, return_var;
    size_t i, j, l, k, m;
    const size_t trial = 5*dim; 
    size_t i_smp[4];
    uint64_t tot_gen_upop, elite;
    uint64_t nfeval = 0;
    uint64_t ubudget = (int)ceil((Dubudget * (double)maxFE-1)/(5.0*(double)num_exe + 4.0));
    uint64_t cr_exp = (uint64_t)(powf(0.5, 1/((double)dim*0.25)) * MAX_RAND);
    
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
    fpop    = (double*)calloc((NPOP+1), sizeof(double));
    if(fpop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }
    pop     = (double*)calloc(((NPOP+1)*dim), sizeof(double));
    if(pop==NULL){
        printf("Error allocating memory");
        goto cleanup;
    }

    elite       = 0;
    //initialize PV
    for(i=0; i<dim; i++){
        mu_v[i] = 0.0;
        std_v[i]= 100.0;
        sampleSolution(0.0, 10.0, pop, 5, i, dim);
    }
    //initialize pop
    for(i=0; i<NPOP; i++){
        fpop[i] = feval(&pop[i*dim], dim);
        nfeval++;
        if(i!=elite && fpop[i] < fpop[elite]){
            elite = i;
        }
    }

    for (tot_gen_upop = 0; tot_gen_upop < ubudget; tot_gen_upop++){
        //micro population module
        for(l=0; l<num_exe; l++){
            for(i=0; i<NPOP; i++){
                //Bin crossover with the mutant vector
                i_smp[0]    = i;
                randSample(i_smp, NPOP, 3);
                i_smp[0]    = i*dim;
                i_smp[1]    = i_smp[1]*dim;
                i_smp[2]    = i_smp[2]*dim;
                i_smp[3]    = i_smp[3]*dim;
                j   = (size_t)(next() % dim);
                for(k=0; k<dim; k++){
                    double temp;
                    //create mutant element if needed only
                    if(next() < cr_bin || k==j){
                        temp    = pop[i_smp[1] + k] + f*(pop[i_smp[2] + k] - pop[i_smp[3] + k]);
                        //fix bounds
                        if(temp>1.0){
                            temp-=2.0;
                        }
                        else if(temp<-1.0){
                            temp+=2.0;
                        }
                    }else{
                        temp  = pop[i_smp[0] + k];
                    }
                    pop[trial +k] = temp;
                }
                
                //Competition and update of the model
                fpop[5]  = feval(&pop[trial], dim);
                nfeval++;
                if(fpop[5] < fpop[i]){
                    //update model
                    //updatePv(mu_v, std_v, trial, &pop[i*dim], np, dim);
                    for(k=0; k<dim; k++){
                        double win, los, mu;
                        if(pop[i_smp[0]+k]  != pop[trial + k]){
                            win = pop[trial + k];
                            los = pop[i_smp[0]+k];
                            mu  = mu_v[k];
                            /*if((mu > win) && (mu - win > 1.0)){
                                win = win + 2.0;
                            }else if((mu < win) && (win - mu > 1.0)){
                                win = win - 2.0;
                            }
                            if((mu > los) && (mu - los > 1.0)){
                                los = los + 2.0;
                            }else if((mu < los) && (los - mu > 1.0)){
                                los = los - 2.0;
                            }*/
                            mu_v[k] += np*(win-los);
                            std_v[k] = fabs(std_v[k] + mu*mu - mu_v[k]*mu_v[k] + np*(win*win - los*los));

                            if(mu_v[k] < -1.0){
                                mu_v[k] = mu_v[k]+2.0;
                            }
                            else if (mu_v[k] > 1.0){
                                mu_v[k] = mu_v[k]-2.0;
                            }
                            if(std_v[k] > 100.0){
                                std_v[k] = 100.0;
                            }

                            pop[i_smp[0]+k]  = pop[trial + k];
                        }
                    }
                    fpop[i] = fpop[5];
                    //update elite position
                    if(i != elite && fpop[i]<fpop[elite]){
                        elite = i;
                    }
                }else{
                    //update model
                    //updatePv(mu_v, std_v, &pop[i_smp[0]], trial, np, dim);
                    
                    for(k=0; k<dim; k++){
                        double win, los, mu;
                        if(pop[i_smp[0]+k]  != pop[trial + k]){
                            win = pop[i_smp[0]+k];
                            los = pop[trial + k];
                            mu  = mu_v[k];
                            /*if((mu > win) && (mu - win > 1.0)){
                                win = win + 2.0;
                            }else if((mu < win) && (win - mu > 1.0)){
                                win = win - 2.0;
                            }
                            if((mu > los) && (mu - los > 1.0)){
                                los = los + 2.0;
                            }else if((mu < los) && (los - mu > 1.0)){
                                los = los - 2.0;
                            }*/
                            
                            mu_v[k] += np*(win-los);
                            std_v[k] = fabs(std_v[k] + mu*mu - mu_v[k]*mu_v[k] + np*(win*win - los*los));
                            if(mu_v[k] < -1.0){
                                mu_v[k] = mu_v[k]+2.0;
                            }
                            else if (mu_v[k] > 1.0){
                                mu_v[k] = mu_v[k]-2.0;
                            }
                            if(std_v[k] > 100.0){
                                std_v[k] = 100.0;
                            }
                        }
                    }
                }
            }
        }

        //restart the worst elements in pop 
        if(tot_gen_upop < ubudget-1){
            for(i=0; i<dim; i++){
                sampleSolution(mu_v[i], sqrt(std_v[i]), pop, elite, i, dim);
            }
            
            for(i=0; i<NPOP; i++){
                if(elite == i) {continue;}
                fpop[i] = feval(&pop[i*dim], dim);
                nfeval++;
                if(i!=elite && fpop[i] < fpop[elite]){
                    elite = i;
                }
            }
        }
    }
    
    for(i=0; i<3; i++){
        if(i >= elite)
            i_smp[i]    = i+1;
        else
            i_smp[i]    = i;
    }
    for(i=0; i<dim;i++){
        std_v[i] = sqrt(std_v[i]);
    }
    while(nfeval < maxFE){
        j   = (size_t)(next() % dim);
        k   = j;
        //Exponential crossover with mutant vector, created from the PV
        do{
            //mutant element
            sampleSolution(mu_v[k], std_v[k], pop, -elite, k, dim);
            pop[trial + k]    = pop[i_smp[0]*dim + k] + f*(pop[i_smp[1]*dim + k] - pop[i_smp[2]*dim + k]);
            //fix bounds
            if(pop[trial + k]>1.0){
                pop[trial + k]-=2.0;
            }
            else if(pop[trial + k]<-1.0){
                pop[trial + k]+=2.0;
            }
            k++;
            if(k >= dim){
                k   = 0;
            }
        }while(next() < cr_exp && k != j);
        m   = k;
        while(k != j){
            pop[trial + k]    = pop[elite*dim + k];
            k++;
            if(k>=dim){
                k=0;
            }
        }
        
          
        //Competition
        fpop[5]  = feval(&pop[trial], dim);
        nfeval++;
        if(fpop[5] < fpop[elite]){
            //update model
            updatePv(mu_v, std_v, &pop[trial], &pop[elite*dim], np, dim);
            //copy the modified solution
            while(j != m){
                pop[elite*dim+j]  = pop[trial + j];
                j++;
                if(j>=dim){
                    j=0;
                }
            }
            fpop[elite] = fpop[5];
        }else{
            //update the model
            updatePv(mu_v, std_v, &pop[elite*dim], &pop[trial], np, dim);
        }
        

    }
    return_var  = fpop[elite];

    for(i=0;i<dim;i++){
        solution[i] = pop[elite*dim+i];
    }

    cleanup:
    free(pop);
    free(fpop);
    free(mu_v);
    free(std_v);
    return return_var;
}
