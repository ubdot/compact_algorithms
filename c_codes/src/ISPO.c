#include "ISPO.h"
double ISPO(double *trial, double A, double P, double B, double S, double eps, uint64_t h, uint64_t dim, uint64_t budget, double (*feval)(double *, int)){
    double *trial1, ftrial, fbackup;
    double L, trial_backup, v;
    uint64_t tot_eval   = 0, i, j, k;

    //init solution
    for(i=0; i<dim; i++){
        trial[i]=2.0f*doubleRand()-1.0;
    }
    ftrial  = feval(trial, dim);
    tot_eval++;
    //init backup information
    fbackup = ftrial;

    while(tot_eval<budget){
        for(i=0; i< dim; i++){
            L   = 0.0f;
            trial_backup = trial[i];
            for(j=0; j<h; j++){
                //Calculate perturbation
                v           = (A/(pow((double)(j+1),P)))*(doubleRand()-0.5)+B*powf(L, (double)(j+1));
                trial[i]    +=v;
                //Repair bounds
                while (trial[i]>1.0)
                {
                    trial[i]    -= 2.0;
                }
                
                while (trial[i]<-1.0)
                {
                    trial[i]    += 2.0;
                }
                ftrial  = feval(trial, dim);
                tot_eval++;
                if(tot_eval>=budget){
                    break;
                }
            }
            if(tot_eval>=budget){
                break;
            }
            if(ftrial <= fbackup){//update backup information in case the best point is updated
                L           = v;        
                fbackup     = ftrial; 
                trial_backup= trial[i];
            }
            else{//Otherwise update the perturbation step size and back to the original point
                if (L != 0){
                    L   = L/S;
                }
                
                if (L < eps){
                    L   = 0;
                }
                trial[i]  = trial_backup;
            }
        }
    }
    return fbackup;
}