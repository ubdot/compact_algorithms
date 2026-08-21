
#include <malloc.h>
#include <stdint.h>
#include <math.h>
#include "compact.h"
 
/**
 * @brief Algorithm implementation based on the paper "Compact Differential Evolution"
 * DE/rand/1/exp
 * 
 * @param solution array where it will be stored the best solution at the end of the algorithm
 * @param f scaling factor for the Differential Evolution 
 * @param cr crossover probablity, in order to avoid unnecesary convertions the proportioned number must be in 64 bits considering the maximun possible number
 *          0xffffffffffffffff, so to calculate the actual value cr=0xffffffffffffffff * CR, where CR is a floating number in he range [0,1]
 * @param budget maximun number of function evaluations
 * @param np value of 1/NP where NP is the value of the virtual population for the compact approach
 * @param dim problem dimension
 * @param feval objective function with parameters (*solution, dimension
 * @return double feval(*elite, dimension)
 */
double cDE(double *solution, double f, uint64_t cr, int budget, double np, int dim, double (*feval)(double *, int));

/*E. Mininno, F. Neri, F. Cupertino and D. Naso, "Compact Differential Evolution," in IEEE Transactions on Evolutionary Computation, vol. 15, no. 1, pp. 32-54, Feb. 2011, doi: 10.1109/TEVC.2010.2058120.
keywords: {Algorithm design and analysis;Optimization;Evolutionary computation;Training;Robots;Memory management;Hardware;Adaptive systems;compact genetic algorithms;differential evolution (DE);estimation distribution algorithms},*/