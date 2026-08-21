#include <math.h>
#include <malloc.h>
#include <stdint.h>
#include "NUM.h"


/**
 * @brief Algorithm implementation based on the paper "A single-solution–compact hybrid algorithm for continuous optimization"
 * 
 * @param elite array where it will be stored the best solution at the end of the algorithm
 * @param budget maximun number of function evaluations
 * @param Np value of 1/NP where NP is the value of the virtual population for the compact approach
 * @param B parameter determining the degree of dependency on the iteration counter (non-uniformity)
 * @param dim problem dimension
 * @param feval objective function with parameters (*solution, dimension)
 * @return double  feval(*elite, dimension)
 */
double cSM(double *elite, uint64_t budget, double Np, double B, uint64_t dim, double (*feval)(double *, int));

/*Khalfi, S., Iacca, G., & Draa, A. (2023). A single-solution–compact hybrid algorithm for continuous optimization. Memetic computing, 15(2), 155-204.*/