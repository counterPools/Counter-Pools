#pragma once

#ifndef CounterPoolsDefs
#define CounterPoolsDefs

#include <x86intrin.h>

#define STARS_AND_BARS_64_4 47905
#define STARS_AND_BARS_64_4_POOL_FAILED 47906

#define FT_SIZE 13
#define FT_SIZE_WEIGHTS 15

#define GRAD_ELEMENT_SIZE 8

#define MAX(a,b) a>b?a:b

#define CONCAT_KEYS(a,b,c) (((a) << 32) | ((b) << 16) | (c))

/* Auxilary functions - start */

inline uint64_t nChoosek(uint64_t n, uint64_t k)
{
	if (k > n) return 0;
	if (k * 2 > n) k = n - k;
	if (k == 0) return 1;

	uint64_t result = n;
	for (uint64_t i = 2; i <= k; ++i) {
		result *= (n - i + 1);
		result /= i;
	}
	return result;
}

inline uint64_t starsAndBars(uint64_t n, uint64_t k)
{
	return nChoosek(n + k - 1, k - 1);
}

using namespace std;

#endif // !CounterPoolsDefs
