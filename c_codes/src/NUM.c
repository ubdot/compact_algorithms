#include "NUM.h"
//Non-uniform mutation applied to a single component of the solution, the component is selected randomly
//based on a random number number determines the input parameters for the calculation of delta
//the inner function next() return a 64 bits uniform random number
//xoff_p:   the array containing the solution to be set under the mutation process
//MaxIt:    maximun number of function evaluations
//it:       actual number of function evaluations
//B:        parameter for the mutation operator, determining the degree of dependency on the iteration counter (non-uniformity)
//dim:      problem dimension
uint64_t SNUM(double *xoff_p, uint64_t MaxIt, uint64_t it, double B, uint64_t dim){
    uint64_t    i = next()%dim; //calculate the component to apply the NUM operator
    if(next() < 0x7FFFFFFFFFFFFFFF){//Determines the input for the function delta
        xoff_p[i]   +=  delNum(MaxIt, it, 1.0 - xoff_p[i], B);
    }else{
        xoff_p[i]   +=  delNum(MaxIt, it, 1.0 + xoff_p[i], B);
    }

    //If the mutated number is outside boundaries is repaired accordingly
    while(xoff_p[i] > 1.0){
        xoff_p[i] -= 2.0;
    }

    while(xoff_p[i] < -1.0){
        xoff_p[i] += 2.0;
    }

    return i;
}


//Non uniform mutation applied to all componentes in the solution (Multiple non-uniform mutation),
//based on a random number number determines the input parameters for the calculation of delta
//the inner function next() return a 64 bits uniform random number
//xoff_p:   the array containing the solution to be set under the mutation process
//MaxIt:    maximun number of function evaluations
//it:       actual number of function evaluations
//B:        parameter for the mutation operator, determining the degree of dependency on the iteration counter (non-uniformity)
//dim:      problem dimension
void MNUM(double *xoff_p, uint64_t MaxIt, uint64_t it, double B, uint64_t dim){
    uint64_t    i;
    for(i=0; i<dim; i++){
        if(next() < 0x7FFFFFFFFFFFFFFF){//Determines the input for the function delta
            xoff_p[i]   +=  delNum(MaxIt, it, 1.0 - xoff_p[i], B);
        }else{
            xoff_p[i]   +=  delNum(MaxIt, it, 1.0 + xoff_p[i], B);
        }

        //If the mutated number is outside boundaries is repaired accordingly
        while(xoff_p[i] > 1.0){
            xoff_p[i] -= 2.0;
        }

        while(xoff_p[i] < -1.0){
            xoff_p[i] += 2.0;
        }
    }
}

//calculates delta for non-uniform mutation based on he number of function evaluations, the actual value of the component in the solution and a parameter B
//MaxIt:    maximun number of function evaluations
//it:       actual number of function evaluations
//y:        value of the component previous to the mutation operator
//B:        parameter for the mutation operator, determining the degree of dependency on the iteration counter (non-uniformity)
double delNum(uint64_t MaxIt, uint64_t it, double y, double B){
    double delt;
    delt = pow((1.0-((double)it/(double)MaxIt)), B);
    delt = y*(1.0-pow(doubleRand(),delt));
    return delt;
}
