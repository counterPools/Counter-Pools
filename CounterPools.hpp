#pragma once

#ifndef COUNT_CounterPools_64_4_0_1
#define COUNT_CounterPools_64_4_0_1

#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <unordered_map>

#include "xxhash.hpp"
#include "CounterPoolsDefs.h"

/* Auxilary functions - end */

class CounterPools_64_4_0_1
{

	int flag;
	int m_seed;
	int m_height_plus_1_over_2;

	int m_width;
	int m_height;

	int m_width_mask;
	int m_log_width;

	uint64_t** m_pools;
	uint16_t** m_encodings;

	static bool m_initLookupTable;
	static uint32_t m_lookup[STARS_AND_BARS_64_4];

	void initLookupTableFunc();

	uint32_t decode_sizes(uint16_t encoding, uint64_t n, uint64_t k);
	uint32_t decode_offsets(uint16_t encoding, uint64_t n, uint64_t k);

	uint16_t encode(uint32_t sizes);
	uint16_t encode_without_auxilary_map(uint32_t sizes);

	std::unordered_map< uint64_t, uint16_t> m_auxilary_map;

	inline uint64_t read_counter(uint64_t &pool, int counterBitSize, int counterBitOffset);

	virtual bool increment_h(uint64_t &pool, int counterBitSize, int counterBitOffset);
	virtual void increment_index(int row, uint32_t index);
	uint64_t deleteCounter(uint64_t &pool, uint8_t* offsets, int index);
	uint64_t chooseDeleteCounter(uint64_t pool, uint8_t* offsets, int index);
	void reInsert(char * key, uint64_t value, int raw);

	char*** keys;

public:

	CounterPools_64_4_0_1(int width, int height, int seed);
	~CounterPools_64_4_0_1();

	void increment(const char * str);
	void incrementCus(const char * str);
	void increment2(const char * str);
	virtual uint64_t query(const char * str);

};

#endif

