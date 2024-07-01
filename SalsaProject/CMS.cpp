#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "CMS.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

CountMinBaseline::CountMinBaseline()
{
}

CountMinBaseline::~CountMinBaseline()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms[i];
	}

#ifdef USE_BOBHASH
	delete[] bobhash;
#endif // USE_BOBHASH

	delete[] baseline_cms;
}

void CountMinBaseline::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	int index = width_mask;
	log_width = 0;
	while (index >>= 1) ++log_width;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

#ifdef USE_ONE_XXHASH
	if (height > 1)
	{
		assert(width <= ((uint64_t)1 << (64 / height)) && "One XXHASH is not enough!");
	}	
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	height_plus_1_over_2 = (height + 1) / 2;
	if (height_plus_1_over_2 > 1)
	{
		assert(width < ((uint64_t)1 << (64 / height_plus_1_over_2)) && "Two XXHASHes are not enough!");
	}
#endif // USE_TWO_XXHASHES

#ifdef USE_128_BIT_XXHASH
	height_plus_1_over_2 = (height + 1) / 2;
	if (height_plus_1_over_2 > 1)
	{
		assert(width < ((uint64_t)1 << (64 / height_plus_1_over_2)) && "128 hash bits are not enough!");
	}
#endif // USE_128_BIT_XXHASH

	baseline_cms = new uint32_t*[height];

#ifdef USE_BOBHASH
	bobhash = new BOBHash[height];
#endif // USE_BOBHASH

	for (int i = 0; i < height; ++i)
	{
		baseline_cms[i] = new uint32_t[width]();

#ifdef USE_BOBHASH
		bobhash[i].initialize(seed*(7 + i) + i + 100);
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
		seeds[i] = seed * (7 + i) + i + 100;
#endif // if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES

	}

#ifdef USE_128_BIT_XXHASH
	xxhash_seed = seed * 7 + 100;
#endif // USE_128_BIT_XXHASH

}

void CountMinBaseline::increment(const char * str)
{

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
#endif // USE_TWO_XXHASHES

#ifdef USE_128_BIT_XXHASH
	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, xxhash_seed);

	int i;
	for (i = 0; i < height_plus_1_over_2; ++i) {
		uint index = hashes.low64 & width_mask;
		hashes.low64 >>= log_width;
		++baseline_cms[i][index];
	}
	for (; i < height; ++i) {
		uint index = hashes.high64 & width_mask;
		hashes.high64 >>= log_width;
		++baseline_cms[i][index];
	}
#endif // USE_128_BIT_XXHASH

#ifdef USE_TWO_XXHASHES

	int i;
	for (i = 0; i < height_plus_1_over_2; ++i) {
		uint index = hashes0 & width_mask;
		hashes0 >>= log_width;
		++baseline_cms[i][index];
	}
	for (; i < height; ++i) {
		uint index = hashes1 & width_mask;
		hashes1 >>= log_width;
		++baseline_cms[i][index];
	}
#endif
#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_BOBHASH
	for (int i = 0; i < height; ++i) {

#ifdef USE_BOBHASH
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[i])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

		++baseline_cms[i][index];
	}
#endif

}

uint64_t CountMinBaseline::query(const char * str)
{

#ifdef USE_BOBHASH
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
	uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[0])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint index = hashes & width_mask;
	hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
	uint index = hashes0 & width_mask;
	hashes0 >>= log_width;
#endif // USE_TWO_XXHASHES

#ifdef USE_128_BIT_XXHASH
	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, xxhash_seed);
	uint index = hashes.low64 & width_mask;
	hashes.low64 >>= log_width;
#endif // USE_128_BIT_XXHASH

	uint64_t min = baseline_cms[0][index];

#ifdef USE_TWO_XXHASHES

	int i;
	for (i = 1; i < height_plus_1_over_2; ++i) {
		uint index = hashes0 & width_mask;
		hashes0 >>= log_width;
		uint64_t temp = baseline_cms[i][index];
		if (min > temp)
		{
			min = temp;
		}
	}
	for (; i < height; ++i) {
		uint index = hashes1 & width_mask;
		hashes1 >>= log_width;
		uint64_t temp = baseline_cms[i][index];
		if (min > temp)
		{
			min = temp;
		}
	}

#endif

#ifdef USE_128_BIT_XXHASH
	int i;
	for (i = 1; i < height_plus_1_over_2; ++i) {
		uint index = hashes.low64 & width_mask;
		hashes.low64 >>= log_width;
		uint64_t temp = baseline_cms[i][index];
		if (min > temp)
		{
			min = temp;
		}
	}
	for (; i < height; ++i) {
		uint index = hashes.high64 & width_mask;
		hashes.high64 >>= log_width;
		uint64_t temp = baseline_cms[i][index];
		if (min > temp)
		{
			min = temp;
		}
	}
#endif // USE_128_BIT_XXHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_BOBHASH

	for (int i = 1; i < height; ++i) {

#ifdef USE_BOBHASH
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[i])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

		uint64_t temp = baseline_cms[i][index];
		if (min > temp)
		{
			min = temp;
		}
	}

#endif
	return min;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

WeightedCountMinBaseline::WeightedCountMinBaseline()
{
}

WeightedCountMinBaseline::~WeightedCountMinBaseline()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms[i];
	}
	delete[] bobhash;
	delete[] baseline_cms;
}

void WeightedCountMinBaseline::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cms = new uint64_t*[height];
	bobhash = new BOBHash[height];

	for (int i = 0; i < height; ++i)
	{
		baseline_cms[i] = new uint64_t[width]();
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}
}

void WeightedCountMinBaseline::add(const char * str, int c)
{
	for (int i = 0; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		baseline_cms[i][index] += c;
	}
}

uint64_t WeightedCountMinBaseline::query(const char * str)
{
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
	uint64_t min = baseline_cms[0][index];
	for (int i = 1; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = baseline_cms[i][index];
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

ConservativeUpdateBaseline::ConservativeUpdateBaseline()
{
}

ConservativeUpdateBaseline::~ConservativeUpdateBaseline()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cu[i];
	}
	delete[] bobhash;
	delete[] baseline_cu;

	delete[] counters;
}

void ConservativeUpdateBaseline::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cu = new uint32_t*[height];
	bobhash = new BOBHash[height];

	for (int i = 0; i < height; ++i)
	{
		baseline_cu[i] = new uint32_t[width]();
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}

	counters = new uint32_t*[height];
}

void ConservativeUpdateBaseline::increment(const char * str)
{
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
	counters[0] = baseline_cu[0] + index;
	uint64_t min = *counters[0];
	
	for (int i = 1; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		counters[i] = baseline_cu[i] + index;
		uint64_t temp = *counters[i];
		if (min > temp)
		{
			min = temp;
		}
	}

	for (int i = 0; i < height; ++i) {
		if (*counters[i] == min)
		{
			++(*counters[i]);
		}	
	}
}

uint64_t ConservativeUpdateBaseline::query(const char * str)
{
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
	uint64_t min = baseline_cu[0][index];
	for (int i = 1; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = baseline_cu[i][index];
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

ConservativeUpdateBaselineSanity::ConservativeUpdateBaselineSanity()
{
}

ConservativeUpdateBaselineSanity::~ConservativeUpdateBaselineSanity()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cu[i];
	}
	delete[] bobhash;
	delete[] baseline_cu;

	delete[] counters;
}

void ConservativeUpdateBaselineSanity::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = (width - 1) << 5;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cu = new uint32_t*[height];
	bobhash = new BOBHash[height];

	for (int i = 0; i < height; ++i)
	{
		baseline_cu[i] = new uint32_t[width]();
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}

	counters = new uint32_t*[height];
}

void ConservativeUpdateBaselineSanity::increment(const char * str)
{ 
	uint index = ((bobhash[0].run(str, FT_SIZE)) & width_mask) >> 5;
	counters[0] = baseline_cu[0] + index;
	uint64_t min = *counters[0];

	for (int i = 1; i < height; ++i) {
		uint index = ((bobhash[i].run(str, FT_SIZE)) & width_mask) >> 5;
		counters[i] = baseline_cu[i] + index;
		uint64_t temp = *counters[i];
		if (min > temp)
		{
			min = temp;
		}
	}

	for (int i = 0; i < height; ++i) {
		if (*counters[i] == min)
		{
			++(*counters[i]);
		}
	}
}

uint64_t ConservativeUpdateBaselineSanity::query(const char * str)
{
	uint index = ((bobhash[0].run(str, FT_SIZE)) & width_mask) >> 5;
	uint64_t min = baseline_cu[0][index];
	for (int i = 1; i < height; ++i) {
		uint index = ((bobhash[i].run(str, FT_SIZE)) & width_mask) >> 5;
		uint64_t temp = baseline_cu[i][index];
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
