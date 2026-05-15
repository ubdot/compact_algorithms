#include "compact_lib.h"
#define SQRT_2_value 1.4142135623730950488016887242096980

double sampleSolution(double mu, double sigma){
    double x;
    if(sigma <= 1.0e-8){
        x = mu;
    }else{
        double sqrt_2_sigma, erf_util, erf_mu_plus, smp, C;
        sqrt_2_sigma  = SQRT_2_value * sigma;
        erf_mu_plus = erf((mu + 1.0)/sqrt_2_sigma);
        erf_util    = erf((mu - 1.0)/sqrt_2_sigma) - erf_mu_plus;
        smp = genrand64_real2();
        if (smp == 0 || smp == 1)
        {
            smp = genrand64_real2();
        }
        
        C   = -erf_mu_plus/(erf_util);
        C   = (smp - C)*erf_util;
        x   = mu - sqrt_2_sigma*ierf(C);
    }
    return x;
}

double denorm(double x, double low_lim, double up_lim){
    double xr;
    xr  = up_lim - low_lim;
    xr  = (x + 1.0)*xr/2.0 + low_lim;
    return xr;
}

void updatePv(double *muV, double *stdV, double *winner, double *looser, double Np, int len){
    int i;
    double win, los, mu, muvi;
    //printf("%lf \n", looser[9]);
    for(i = 0; i<len; i++){
        if(winner[i] != looser[i]){
            win = winner[i];
            los = looser[i];
            mu  = muV[i];
            if((mu > win) && (mu - win > 2.0-mu + win)){
                win = win + 2.0;
            }else if((mu < win) && (win - mu > 2.0-win + mu)){
                win = win - 2.0;
            }
            if((mu > los) && (mu - los > 2.0-mu + los)){
                los = los + 2.0;
            }else if((mu < los) && (los - mu > 2.0 - los + mu)){
                los = los - 2.0;
            }
            muvi = mu + (1.0/Np)*(win-los);
            muV[i] = muvi;
            if(muV[i] < -1.0){
                muV[i] = muV[i]+2.0;
            }
            else if (muV[i] > 1.0){
                muV[i] = muV[i]-2.0;
            }
            
            stdV[i] = fabs(stdV[i]*stdV[i] + mu*mu - muvi*muvi + (1/Np)*(win*win - los*los));
            stdV[i] = sqrt(stdV[i]);

            if(stdV[i] > 10){
                stdV[i] = 10;
            }
        }
    }
    
}

void updatePvB(double *muV, double *stdV, double *winner, double *looser, double Np, int len){
    int i;
    double win, los, mu, muvi;
    for(i = 0; i<len; i++){
        if(winner[i] != looser[i]){
            win = winner[i];
            los = looser[i];
            mu  = muV[i];
            muvi = mu + (1.0/Np)*(win-los);
            muV[i] = muvi;
            if(muV[i] < -1.0){
                muV[i] = -1.0;
            }
            else if (muV[i] > 1.0){
                muV[i] = 1.0;
            }
            stdV[i] = sqrt(fabs(stdV[i]*stdV[i] + mu*mu - muvi*muvi + (1/Np)*(win*win - los*los)));
        }
    }
}

void randSample(int *arr_smp, int n_pop, int n_smp){
    unsigned long long i, j, k, m, smp;

    for(i=1; i<n_smp+1; i++){
        smp     = genrand64_int64()%(n_pop-i);
        m       = 0;
        for(j=0; j< n_pop; j++){
            for(k=0; k<i; k++){
                if(arr_smp[k]==j){
                    break;
                }
            }
            if(k==i){
                if(m==smp)
                    break;
                else
                    m++;
            }
        }
        arr_smp[i]    = j;
    }
}
