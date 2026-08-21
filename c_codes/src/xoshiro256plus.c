#include "xoshiro256plus.h"

uint64_t next(void) {
	const uint64_t result = s[0] + s[3];

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}

uint64_t nextInput(uint64_t *seed) {
	const uint64_t result = seed[0] + seed[3];

	const uint64_t t = seed[1] << 17;

	seed[2] ^= seed[0];
	seed[3] ^= seed[1];
	seed[1] ^= seed[2];
	seed[0] ^= seed[3];

	seed[2] ^= t;

	seed[3] = rotl(seed[3], 45);

	return result;
}

void set_seed(uint64_t seed){
	x=seed;
	
	s[2] = create_seed();
	s[3] = create_seed();
	s[1] = create_seed();
	s[0] = create_seed();
	
}

uint64_t create_seed() {
	uint64_t z = (x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}
	