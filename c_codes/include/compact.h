#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "erf.h"
#include "xoshiro256plus.h"

#define SQRT_2_value 1.4142135623730950488016887242096980
#define MAX_RAND 0xffffffffffffffff
/**
 * @brief Function that samples solutions following the compact scheme, it has multiple configurations
 * if the rows values is a negative number the row in the function skips the respective row (-row), util when multiple solutions are ampled using the same pv values
 * if the column is negative the function fill the complete array, is util when init solutions
 * the function was created to perform multiple samples
 * 
 * @param mu mean of the PV
 * @param sigma standard deviation of the PV
 * @param solution array where it will be stored the solution [n * dim], where n is the number of solutions, and dim is the problem dimension
 * @param row number of rows in the array
 * @param col column that will be modified in the array
 * @param cols number of coulumns in the array containing the solutions 
 */
void        sampleSolution(double mu, double sigma, double *solution, int row, int col, int cols);
/**
 * @brief used to sample a variable from the PV, for multiple samples use sampleSolution
 * 
 * @param mu mean of the PV
 * @param sigma standard deviation of the PV
 * @return double 
 */
double      singleSampleSolution(double mu, double sigma);
/**
 * @brief Uses random uniform function (declared in header), that gives a number unifromly distributed of 64 bits to calculate a float random on the range [0,1]
 * 
 * @return double 
 */
double      doubleRand(void);
/**
 * @brief this function updates the PV considering toroidal distribution of solutions, with the mean value of each distribution being the center of the space
 * 
 * @param muV array containing the means of the PV
 * @param stdV array containing the standard deviations of the PV
 * @param winner array containg the winner solution of the competition
 * @param looser array containing the loser solution of the competition
 * @param Np result of the 1/Np, where Np is the virtual population of the compact scheme
 * @param len problem dimension, used to access the multiple arrays
 */
void        updatePv(double *muV, double *stdV, double *winner, double *looser, double Np, int len);
/**
 * @brief * @brief This function updates the PV using the compact equations; the mean values ​​are forced to the maximum or minimum value after the update if they fall outside the bounds.
 * 
 * @param muV array containing the means of the PV
 * @param stdV array containing the standard deviations of the PV
 * @param winner array containg the winner solution of the competition
 * @param looser array containing the loser solution of the competition
 * @param Np result of the 1/Np, where Np is the virtual population of the compact scheme
 * @param len problem dimension, used to access the multiple arrays
 */
void        updatePvB(double *muV, double *stdV, double *winner, double *looser, double Np, int len);
/**
 * @brief this function is used to sample positions in a population, used with Differential evolution to sample individuals to create the mutant vector
 * to avoid sampling a member of the population the position of the individual must be stored in the "arr_smp" array, in the [0] position
 * 
 * @param arr_smp array that where the integer sampled numbers will be stored
 * @param n_pop population size
 * @param n_smp number of integers to be sampled
 */
void        randSample(size_t *arr_smp, size_t n_pop, size_t n_smp);