#include <stdio.h>
#include "erf.h"
#include "compact_lib.h"
#include <time.h>
#include <malloc.h>
#include <stdlib.h>

#define SVPLOT  0

void cec17_test_func(double *, double *,int,int,int);
double fast_pow(double x, int p);
double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;
int GNVars;
int MaxFEval = 0;
double fopt[1];
int sel_fun;

double feval(double *x, int len){
    double f, *x1, x_min, x_max, xr;
    int i=0;
    //Set bound constrains to evaluate objective function
    switch(sel_fun){
    case 32://rosenbrock
        x_min   = -5;
        x_max   = 10;
        break;
    case 33://ackley
        x_min   = -32.768;
        x_max   = 32.768;
        break;
    case 34://rastrigin
        x_min   = -5.12;
        x_max   = 5.12;
        break;
    case 35://schwefel
        x_min   = -500;
        x_max   = 500;
        break;
    default://sphere
        x_min   = -100;
        x_max   = 100;
        break;
    }
    
    xr  = x_max - x_min;
    for(i=0; i<len; i++){
        x[i]=((x[i] + 1.0)*xr)/2.0 + x_min;
    }
    cec17_test_func(x, &f, len, 1, sel_fun);
    for(i=0; i<len; i++){
        x[i]=((x[i] - x_min)*2.0)/xr - 1.0;
    }
    return f;
}

int main(int argc, char* argv[]){
    int D, elite, n_pop, max_eval, tot_eval, upop_cross, num_exe;
    int i, j, k, l, m;
    unsigned long long seed_num;
    int i_smp[4], ubudget;
    double *mu_v, *std_v, *trial, *pop, *fpop;
    double Np, CR, F, ftrial, xr, xs, xt, CR_b;
    double rand_aux, upop_budget;
    
    //Load values to local variables
    sscanf(argv[1], "%d", &sel_fun);    //sel_fun = 6;//
    sscanf(argv[2], "%lf", &F);         //F       = 0.5;//
    sscanf(argv[3], "%lf", &CR_b);      //CR_b    = 0.9;//
    sscanf(argv[4], "%lf", &CR);        //CR      = 0.25;//
    sscanf(argv[5], "%lf", &upop_budget);//upop_budget = 0.8;//
    sscanf(argv[6], "%d", &num_exe);    //;num_exe  = 2;//
    sscanf(argv[7], "%lf", &Np);        //Np=10;//
    sscanf(argv[8], "%ull", &seed_num); //seed_num=100;//
    sscanf(argv[9], "%d", &D);          //D       = 10;                       //

#if defined(SVPLOT) && SVPLOT == 1
    int plt = 10*D;
#endif

    //set seed  for the random generator
    init_genrand64(seed_num);
    //unsigned long long init[4]={0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL}, length=4;
    //init_by_array64(init, length);
    //init_genrand64(init[0]);

    //init some used variables
    upop_cross  = 1;//0exp_cross, 1bin_cross
    CR      = 1.0/pow(2,1/(CR*D));
    Np      = Np*(double)D;
    n_pop   = 5;
    elite   = 0;
    max_eval= 10000*D;
    tot_eval= 0;
    mu_v    = (double*)calloc(D, sizeof(double));
    std_v   = (double*)calloc(D, sizeof(double));
    trial   = (double*)calloc(D, sizeof(double));
    fpop    = (double*)calloc(n_pop, sizeof(double));
    pop     = (double*)calloc((D*n_pop), sizeof(double));
    ubudget = (int)floor(upop_budget*(double)max_eval);

    //initialize PV
    for(i=0; i<D; i++){
        mu_v[i] = 0.0;
        std_v[i]= 10.0;
    }
    //initialize pop
    for(i=0; i<n_pop; i++){
        for(j=0; j<D; j++){
            pop[i*D+j]  = sampleSolution(mu_v[j], std_v[j]);
        }
        fpop[i] = feval(&pop[i*D], D);
        tot_eval++;
        //Each sampled elite solution is update if is outperformed
        if(i!=elite && fpop[i] < fpop[elite]){
            elite = i;
        }
    }


#if defined(SVPLOT) && SVPLOT == 1
    printf("%e\n",fpop[elite]);
#endif

    while (tot_eval < ubudget){
        //micro population module
        for(l=0; l<num_exe; l++){
            for(i=0; i<n_pop; i++){
                //Bin crossover with the mutant vector
                i_smp[0]    = i;
                randSample(i_smp, n_pop, 3);
                j   = genrand64_int64()%D;
                for(k=0; k<D; k++){
                    rand_aux = genrand64_real2();
                    //create mutant element if needed only
                    if(rand_aux < CR_b || k==j){
                        trial[k]    = pop[i_smp[1]*D + k] + F*(pop[i_smp[2]*D + k] - pop[i_smp[3]*D + k]);
                        //fix bounds
                        while(trial[k]>1){
                            trial[k]-=2;
                        }
                        while(trial[k]<-1){
                            trial[k]+=2;
                        }
                    }else{
                        trial[k]    = pop[i*D + k];
                    }
                }
                
                //Competition and update of the model
                ftrial  = feval(trial, D);
                tot_eval++;
                if(ftrial < fpop[i]){
                    //update model
                    updatePv(mu_v, std_v, trial, &pop[i*D], Np, D);
                    for(k=0; k<D; k++){
                        if(pop[i*D+k]  != trial[k]){
                            pop[i*D+k]  = trial[k];
                        }
                    }
                    fpop[i] = ftrial;
                    //update elite position
                    if(i != elite && fpop[i]<fpop[elite]){
                        elite = i;
                    }
                }else{
                    //update model
                    updatePv(mu_v, std_v, &pop[i*D], trial, Np, D);
                }
#if defined(SVPLOT) && SVPLOT == 1
                if(tot_eval%plt==0){
                    printf("%e\n",fpop[elite]);
                }
#endif
                //break loop if max budget is reached
                if(tot_eval >= ubudget){
                    break;
                }
            }
            //break loop if max budget is reached
            if(tot_eval >= ubudget){
                break;
            }
        }

        //restart the worst elements in pop
        if(tot_eval < ubudget - 2*n_pop){
            for(i=0; i<n_pop; i++){
                if(i != elite){
                   for(j=0; j<D; j++){
                        pop[i*D+j]  = sampleSolution(mu_v[j], std_v[j]);
                    }
                    fpop[i] = feval(&pop[i*D], D);
                    tot_eval++;
                    //Update the elite if it is outperformed
                    if(fpop[i] < fpop[elite]){
                        elite = i;
                    }
#if defined(SVPLOT) && SVPLOT == 1
                    if(tot_eval%plt==0){
                        printf("%e\n",fpop[elite]);
                    }
#endif
                }
            }
        }
    }
    //printf("tot eval: %d\n", tot_eval);
    //compact scheme on the best element
    while (tot_eval < max_eval){
        j   = genrand64_int64()%D;
        k   = j;
        //Exponential crossover with mutant vector, created from the PV
        do{
            //mutant element
            xr  = sampleSolution(mu_v[k], std_v[k]);
            xs  = sampleSolution(mu_v[k], std_v[k]);
            xt  = sampleSolution(mu_v[k], std_v[k]);
            trial[k]    = xr + F*(xs - xt);
            //fix bounds
            while(trial[k]>1){
                trial[k]-=2;
            }
            while(trial[k]<-1){
                trial[k]+=2;
            }
            k++;
            if(k>=D){
                k=0;
            }
        }while(genrand64_real2() < CR && k != j);
        m   = k;
        while(k != j){
            trial[k]    = pop[elite*D + k];
            k++;
            if(k>=D){
                k=0;
            }
        }
        //Competition
        ftrial  = feval(trial, D);
        tot_eval++;
        if(ftrial < fpop[elite]){
            //update model
            updatePv(mu_v, std_v, trial, &pop[elite*D], Np, D);
            //copy the modified solution
            while(j != m){
                pop[elite*D+j]  = trial[j];
                j++;
                if(j>=D){
                    j=0;
                }
            }
            fpop[elite] = ftrial;
        }else{
            //update the model
            updatePv(mu_v, std_v, &pop[elite*D], trial, Np, D);
        }
        
#if defined(SVPLOT) && SVPLOT == 1
        if(tot_eval%plt==0){
            printf("%e\n",fpop[elite]);
        }
#endif

    }
#if defined(SVPLOT) && SVPLOT == 1
#else
    printf("%e\n",fpop[elite]);
#endif

    free(pop);
    free(fpop);
    free(mu_v);
    free(std_v);
    free(trial);
    pop     = NULL;
    mu_v    = NULL;
    std_v   = NULL;
    trial   = NULL;
    fpop    = NULL;
    return 0;
}
