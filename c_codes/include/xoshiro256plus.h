
#ifndef XOSHIRO256PLUS_H
#define XOSHIRO256PLUS_H

#include <stdint.h>
//Function created or copied from the resource available in https://prng.di.unimi.it/
/*Blackman, D., & Vigna, S. (2021). Scrambled linear pseudorandom number generators. ACM Transactions on Mathematical Software (TOMS), 47(4), 1-32. */
uint64_t    next(void) ;
uint64_t    create_seed();
void        set_seed(uint64_t seed);
uint64_t 	nextInput(uint64_t *seed);
static uint64_t s[4];
static uint64_t x;

static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

#endif