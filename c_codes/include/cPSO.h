#include <malloc.h>
#include <stdint.h>
#include <math.h>
#include "compact.h"

/**
 * @brief Algorithm implementation based on the paper "Compact Particle Swarm Optimization"
 * 
 * @param solution array where it will be stored the best solution at the end of the algorithm
 * @param c array containing the constants for PSO C[0] multiplier for the actual speed, c[1] multiplier for the local best atraction and c[2] multiplier for the local best atraction
 * @param budget maximun number of function evaluations
 * @param np value of 1/NP where NP is the value of the virtual population for the compact approach
 * @param dim problem dimension
 * @param feval objective function with parameters (*solution, dimension)
 * @return double feval(solution, dim)
 */
double cPSO(double *solution, double *c, int budget, double np, int dim, double (*feval)(double *, int));
/*Neri, F., Mininno, E., & Iacca, G. (2013). Compact particle swarm optimization. Information Sciences, 239, 96-121.*/