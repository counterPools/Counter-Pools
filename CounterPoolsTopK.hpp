#ifndef CounterPoolsTopKSim
#define CounterPoolsTopKSim

#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <unordered_map>

#include "CounterPoolsDefs.h"
#include "RobinMap/robin_map.h"

class CounterPoolsTopK
{

	int m_top_k_count;
	int m_top_k_min;
	int m_k;

	tsl::robin_map<uint64_t, tuple<string, uint32_t, uint32_t>> m_top_k_map;

public:

	CounterPoolsTopK(int k);
	~CounterPoolsTopK();

	bool increment_if_exists(uint64_t fp);
	tuple<string, uint32_t, uint32_t> add(uint64_t fp, uint32_t estimated_size, const char * str);
	uint32_t query(uint64_t fp);
	uint32_t get_min();
};

#endif

