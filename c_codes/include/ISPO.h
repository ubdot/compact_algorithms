#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include "compact.h"

//Algorithm implementation based on the paper "Simplified Intelligence Single Particle Optimization Based Neural Network for Digit Recognition"

/**
 * @brief 
 * 
 * @param trial array where it will be stored the best solution at the end of the algorithm
 * @param A ISPO acceleration
 * @param P ISPO acceleration power factor 
 * @param B ISPO learning coefficient
 * @param S shrinking factor
 * @param eps precision value
 * @param h number of perturbations
 * @param dim problem dimension
 * @param budget maximun number of function evaluations
 * @param feval objective function with parameters (*solution, dimension)
 * @return double feval(*elite, dimension)
 */
double ISPO(double *trial, double A, double P, double B, double S, double eps, uint64_t h, uint64_t dim, uint64_t budget, double (*feval)(double *, int));

/*Zhou, J., Ji, Z., & Shen, L. (2008, October). Simplified intelligence single particle optimization based neural network for digit recognition. 
In 2008 Chinese Conference on Pattern Recognition (pp. 1-5). IEEE.*/