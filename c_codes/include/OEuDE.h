#include <malloc.h>
#include <stdint.h>
#include "compact.h"
#include <stdint.h>
#define NPOP_OE 6
#define MAXINTA 0xfffffff0
#define MAXINTB 0x0fffffff
#define MASKPOS 0x0000000f

uint64_t replace_ind(uint64_t ind, uint64_t pos, uint64_t value);
uint64_t insert_ind(uint64_t ind, uint64_t pos, uint64_t value);
uint64_t move_ind(uint64_t ind, uint64_t pos0, uint64_t pos1);
uint64_t get_ind(uint64_t ind, uint64_t pos);
/**
 * @brief Algorithm implementation based on the paper "Opposition-based ensemble micro-differential evolution"
 * 
 * @param elite array where it will be stored the best solution at the end of the algorithm
 * @param budget maximun number of function evaluations
 * @param dim problem dimension
 * @param feval objective function with parameters (*solution, dimension)
 * @param CR crossover probablity, in order to avoid unnecesary convertions the proportioned number must be in 64 bits considering the maximun possible number
 *           0xffffffffffffffff, so to calculate the actual value cr=0xffffffffffffffff * CR, where CR is a floating number in he range [0,1]
 * @return double objective function with parameters (*solution, dimension)
 */
double oeude(double *elite, uint64_t budget, int dim, double (*feval)(double *, int), uint64_t CR);


/*Salehinejad, H., Rahnamayan, S., & Tizhoosh, H. R. (2017, November). Opposition-based ensemble 
micro-differential evolution. In 2017 IEEE Symposium Series on Computational Intelligence (SSCI) 
(pp. 1-8). IEEE.*/