
#include <malloc.h>
#include <stdint.h>
#include <stdio.h>
#include "compact.h"

#define NPOP 5
/**
 * @brief Algorithm implementation based on the paper "A restart support mechanism based on compact algorithms for micro-population optimization in resource-constrained devices"
 * 
 * @param solution array where it will be stored the best solution at the end of the algorithm
 * @param f scaling factor for the Differential Evolution 
 * @param cr_bin crossover probablity, in order to avoid unnecesary convertions the proportioned number must be in 64 bits considering the maximun possible number
 * @param num_exe number of generations of micro-population algorithm before restart the population
 * @param Dubudget number in the range [0,1], it determines the percentage of function evaluations employed in the micropopulation module and the compact module
 * @param np value of 1/NP where NP is the value of the virtual population for the compact approach
 * @param dim problem dimension
 * @param maxFE maximun number of function evaluations
 * @param feval objective function with parameters (*solution, dimension)
 * @return double feval(*elite, dimension)
 */
double ucde(double *solution, double f, uint64_t cr_bin, int num_exe, double Dubudget, double np, uint64_t dim, int maxFE, double (*feval)(double *, int));