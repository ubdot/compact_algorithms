#include <malloc.h>
#include <stdint.h>
#include <math.h>
#include "compact.h"
#define PI_2 6.28318530717959f
/**
 * @brief Algorithm implementation based on the paper "A compact compound sinusoidal differential evolution algorithm for solving optimisation problems in memory-constrained environments"
 * 
 * @param elite array where it will be stored the best solution at the end of the algorithm
 * @param budget maximun number of function evaluations
 * @param np value of 1/NP where NP is the value of the virtual population for the compact approach
 * @param nWaves integer used to calculate the crossover and scaling factor according the paper
 * @param freq value used to evaluate the sin function
 * @param dim problem dimension
 * @param feval objective function with parameters (*solution, dimension)
 * @return double feval(*elite, dimension)
 */
double CScDE(double *elite, int budget, double np, uint64_t nWaves, double freq, uint64_t dim, double (*feval)(double *, int));

/*Khalfi, S., Draa, A., & Iacca, G. (2021). A compact compound sinusoidal differential evolution algorithm for solving optimisation problems in memory-constrained environments. Expert Systems with Applications, 186, 115705.*/