
#include <stdint.h>
#include "compact.h"

void sampleSolution(double mu, double sigma, double *solution, int row, int col, int cols){
    int i, j;
    if(sigma <= 1.0e-8){//if the standard deviation is a small value the samplign process is avoid and the mean is set to a value
        if(row < 0){
            for(i=0; i<3; i++){
                if(i >= -row) j=i+1;
                else j=i;
                solution[(j)*cols + col] = mu;
            }
        }else{
            for(i=0; i<4; i++){
                if(i == row) continue;
                solution[i*cols + col] = mu;
            }
        }
    }else{
        double sqrt_2_sigma, inv_sqrt_2_sigma, erf_mu_plus, smp, C, erf_mu_neg;
        sqrt_2_sigma        = SQRT_2_value * sigma;//calculate values to sample a number in the error function 
        inv_sqrt_2_sigma    = 1.0/sqrt_2_sigma;
        erf_mu_plus   = myerf01B((mu + 1.0)*inv_sqrt_2_sigma); //lower range
        erf_mu_neg    = myerf01B((mu - 1.0)*inv_sqrt_2_sigma) - erf_mu_plus; //[lower ranger -  upper range]

        if(row < 0){
            if(col >= 0){//used for compact approach
                for(i=0; i<3; i++){ //A total of three values are sampled, used for DE approach
                    if(i >= -row) j=i+1; //the row is avoid so the sampled values are stored in the other memory spaces
                    else j  = i;
                    smp     = doubleRand();
                    C       = smp*erf_mu_neg + erf_mu_plus; //Precaluculated values of the error function are used each time, to avoid unnecesary calls to the error function
                    solution[(j)*cols + col] = mu - sqrt_2_sigma*myerfinv03B(C); //Each value are set in range with the inverse error function
                    if(solution[(j)*cols + col] < -1.0) { //as approximations are used this is to avoid the number falling outside the bounds
                        solution[(j)*cols + col] =-1.0;
                    }
                    else if(solution[(j)*cols + col] > 1.0){
                        solution[(j)*cols + col] =1.0;
                    }
                }
            }
            else{//used to init solution in compact approach
                for(i=0; i<cols; i++){ //this function is used to init an entire vector, when the mean and standard deviation for all values are the same
                    smp         = doubleRand();
                    C           = smp*erf_mu_neg + erf_mu_plus;
                    solution[i] = mu - sqrt_2_sigma*myerfinv03B(C);
                }
            }
        }else{
            for(i=0; i<5; i++){ //used in upop, sample 4 values corresponding the numbers in the population to be restarted
                if(i == row) continue;
                smp = doubleRand();
                C   = smp*erf_mu_neg + erf_mu_plus; //function considers the same mean and standard deviation for all values.
                solution[i*cols + col] = mu - sqrt_2_sigma*myerfinv03B(C); //each value is set in range using the inverse error function
            }
        }
    }
}


double singleSampleSolution(double mu, double sigma){
    double sqrt_2_sigma, inv_sqrt_2_sigma, erf_mu_plus, smp, C, erf_mu_neg, sampledValue;
    sqrt_2_sigma        = SQRT_2_value * sigma; //Calculate range values to somaple a number in the error function
    inv_sqrt_2_sigma    = 1.0/sqrt_2_sigma;
    erf_mu_plus   = myerf01B((mu + 1.0)*inv_sqrt_2_sigma); //lower range
    erf_mu_neg    = myerf01B((mu - 1.0)*inv_sqrt_2_sigma) - erf_mu_plus; //[lower ranger -  upper range]
    
    smp             = doubleRand();
    C               = smp*erf_mu_neg + erf_mu_plus; //Escalate, the number in the range calculated
    sampledValue    = mu - sqrt_2_sigma*myerfinv03B(C); //obtain the sampled value in the [-1, 1] range, using the inverse error function
    return sampledValue;
}

void updatePv(double *muV, double *stdV, double *winner, double *looser, double Np, int len){
    int i;
    double win, los, mu, muvi;
    for(i = 0; i<len; i++){
        if(winner[i] != looser[i]){
            win = winner[i];
            los = looser[i];
            mu  = muV[i];
            if((mu > win) && (mu - win > 1.0)){ //correct winner variable according to the center (mean)
                win = win + 2.0;
            }else if((mu < win) && (win - mu > 1.0)){
                win = win - 2.0;
            }
            if((mu > los) && (mu - los > 1.0)){//correct loser variable according to the center (mean)
                los = los + 2.0;
            }else if((mu < los) && (los - mu > 1.0)){
                los = los - 2.0;
            }
            muvi = mu + (Np)*(win-los); //update mean
            muV[i] = muvi;
            if(muV[i] < -1.0){ //reapir boundaries of the mean 
                muV[i] = muV[i]+2.0;
            }
            else if (muV[i] > 1.0){
                muV[i] = muV[i]-2.0;
            }
            
            stdV[i] = fabs(stdV[i]*stdV[i] + mu*mu - muvi*muvi + (Np)*(win*win - los*los)); //update standard deviaion value
            stdV[i] = sqrtf(stdV[i]);
            if(stdV[i] > 10.0){
                stdV[i] = 10.0;
            }
        }
    }
    
}




double doubleRand(void){
    uint64_t int_rand;
    double frand;
    int_rand = next();
    frand = ((double)int_rand)/((double)MAX_RAND);
    return frand;
}

void randSample(size_t *arr_smp, size_t n_pop, size_t n_smp){
    size_t i, j, k, m, smp;

    for(i=1; i<n_smp+1; i++){
        smp     = next()%(n_pop-i);
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

void updatePvB(double *muV, double *stdV, double *winner, double *looser, double Np, int len){
    int i;
    double win, los, mu, muvi;
    for(i = 0; i<len; i++){
        if(winner[i] != looser[i]){
            win = winner[i];
            los = looser[i];
            mu  = muV[i];
            muvi = mu + (Np)*(win-los);
            muV[i] = muvi;
            if(muV[i] < -1.0){
                muV[i] = -1.0;
            }
            else if (muV[i] > 1.0){
                muV[i] = 1.0;
            }
            
            stdV[i] = fabs(stdV[i]*stdV[i] + mu*mu - muvi*muvi + (Np)*(win*win - los*los));
            stdV[i] = sqrt(stdV[i]);
            if(stdV[i] > 10.0){
                stdV[i] = 10.0;
            }
        }
    }
}
    
    