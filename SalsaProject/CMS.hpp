#pragma once

#ifndef COUNT_MIN_SKETCH
#define COUNT_MIN_SKETCH

#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>

#include "BobHash.hpp"
#include "../xxhash.hpp"
#include "Defs.hpp"

using namespace std;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class CountMinBaseline {

	int width;
	int height;

	int width_mask;
	int log_width;

#ifdef USE_BOBHASH
	BOBHash *bobhash;
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
	int seeds[12];
#endif // if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES	

#ifdef USE_TWO_XXHASHES
	int height_plus_1_over_2;
#endif // USE_TWO_XXHASHES	

#ifdef USE_128_BIT_XXHASH
	int height_plus_1_over_2;
	int xxhash_seed;
#endif // USE_128_BIT_XXHASH

	uint32_t **baseline_cms;

public:

	CountMinBaseline();
	~CountMinBaseline();

	void initialize(int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

class WeightedCountMinBaseline {

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	uint64_t **baseline_cms;

public:

	WeightedCountMinBaseline();
	~WeightedCountMinBaseline();

	void initialize(int width, int height, int seed);
	void add(const char * str, int c);
	uint64_t query(const char * str);

};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class ConservativeUpdateBaselineSanity {

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	uint32_t **baseline_cu;

	uint32_t** counters;

public:

	ConservativeUpdateBaselineSanity();
	~ConservativeUpdateBaselineSanity();

	void initialize(int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

class ConservativeUpdateBaseline {

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	uint32_t **baseline_cu;

	uint32_t** counters;

public:

	ConservativeUpdateBaseline();
	~ConservativeUpdateBaseline();

	void initialize(int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

#endif // !COUNT_MIN_SKETCH
