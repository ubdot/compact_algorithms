#include <malloc.h>
#include <stdint.h>
#include "compact.h"
#include <stdio.h>
#define NPOP 5
#define pop_size 6

/**
 * @brief Algorithm implementation based on the paper "Empirical analysis of a micro-evolutionary algorithm for numerical optimization"
 * replacement and population size has been modified internally
 * 
 * @param solution array where it will be stored the best solution at the end of the algorithm
 * @param f scaling factor for the Differential Evolution 
 * @param cr crossover probablity, in order to avoid unnecesary convertions the proportioned number must be in 64 bits considering the maximun possible number
 *           0xffffffffffffffff, so to calculate the actual value cr=0xffffffffffffffff * CR, where CR is a floating number in he range [0,1]
 * @param budget maximun number of function evaluations
 * @param num_exe number of generations of micro-population algorithm before restart the population
 * @param dim problem dimension
 * @param feval objective function with parameters (*solution, dimension)
 * @return double feval(*solution, dimension)
 */
double ude_mod(double *solution, double f, uint64_t cr, uint64_t budget, uint64_t num_exe, int dim, double (*feval)(double *, int));
/*Viveros Jiménez, F., Mezura Montes, E., & Gelbukh, A. (2012). Empirical analysis of a micro-evolutionary algorithm for numerical 
optimization. Int. J. Phys. Sci, 7(8), 1235-1258.*/

/**
 * @brief Algorithm implementation based on the paper "Empirical Evaluation of an Elitist Replacement Strategy for Differential Evolution with Micro-Populationss"
 * 
 * @param solution array where it will be stored the best solution at the end of the algorithm
 * @param f scaling factor for the Differential Evolution 
 * @param cr crossover probablity, in order to avoid unnecesary convertions the proportioned number must be in 64 bits considering the maximun possible number
 *           0xffffffffffffffff, so to calculate the actual value cr=0xffffffffffffffff * CR, where CR is a floating number in he range [0,1]
 * @param budget maximun number of function evaluations
 * @param GR number of generations of micro-population algorithm before restart the population
 * @param dim problem dimension
 * @param feval objective function with parameters (*solution, dimension)
 * @return double feval(*solution, dimension)
 */
double udeerm(double *solution, double f, double cr, int budget, int GR, int dim, double (*feval)(double *, int));
/*Luna-Ortiz, I., Rodríguez-Molina, A., Villarreal-Cervantes, M. G., Aldape-Pérez, M., Rojas-López, A. G., & 
Paredes-Ballesteros, J. A. (2025). Empirical Evaluation of an Elitist Replacement Strategy for Differential 
Evolution with Micro-Populations. Biomimetics, 10(10), 685.*/


/**
 * @brief Algorithm implementation based on the paper "Empirical analysis of a micro-evolutionary algorithm for numerical optimization"
 * 
 * @param solution array where it will be stored the best solution at the end of the algorithm
 * @param f scaling factor for the Differential Evolution 
 * @param cr crossover probablity, in order to avoid unnecesary convertions the proportioned number must be in 64 bits considering the maximun possible number
 *           0xffffffffffffffff, so to calculate the actual value cr=0xffffffffffffffff * CR, where CR is a floating number in he range [0,1]
 * @param budget maximun number of function evaluations
 * @param GR number of generations of micro-population algorithm before restart the population
 * @param dim problem dimension
 * @param feval objective function with parameters (*solution, dimension)
 * @return double feval(*solution, dimension)
 */
double ude_fix(double *solution, double f, double cr, int budget, int GR, int dim, double (*feval)(double *, int));
/*Viveros Jiménez, F., Mezura Montes, E., & Gelbukh, A. (2012). Empirical analysis of a micro-evolutionary algorithm for numerical 
optimization. Int. J. Phys. Sci, 7(8), 1235-1258.*/