#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "compact.h"

/*Non-uniform mutation, based on the article "Evolutionary programming based on non-uniform mutation"
for more information about the SNUM, and MNUM operators, article "A single-solution–compact hybrid 
algorithm for continuous optimization" can be consulted*/
uint64_t SNUM(double *xoff_p, uint64_t MaxIt, uint64_t it, double B, uint64_t dim);
void MNUM(double *xoff_p, uint64_t MaxIt, uint64_t it, double B, uint64_t dim);
double delNum(uint64_t MaxIt, uint64_t it, double y, double B);
