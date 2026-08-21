
#include <malloc.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "compact.h"

/**
 * @brief Algorithm implementation based on the paper "Compact Differential Evolution Light High performance despite limited memory requirement and modest computational overhead"
 * DE/rand/1/exp
 * 
 * @param solution array where it will be stored the best solution at the end of the algorithm
 * @param f scaling factor for the Differential Evolution 
 * @param cr crossover probablity (double)
 * @param budget maximun number of function evaluations
 * @param np value of 1/NP where NP is the value of the virtual population for the compact approach
 * @param dim problem dimension
 * @param feval objective function with parameters (*solution, dimension)
 * @return double feval(solution, dim)
 */
double cDE_light(double *solution, double f, double cr, int budget, double np, int dim, double (*feval)(double *, int));


/*Iacca, G., Caraffini, F. & Neri, F. Compact Differential Evolution Light: High Performance Despite Limited Memory Requirement and Modest 
Computational Overhead. J. Comput. Sci. Technol. 27, 1056–1076 (2012). https://doi.org/10.1007/s11390-012-1284-2*/