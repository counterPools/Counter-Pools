#ifndef CounterPoolsSim
#define CounterPoolsSim

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
#include "RobinMap/robin_map.h"
#include "CounterPoolsTopK.hpp"

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

class CounterPools {

	int m_width;
	int m_height;

	int m_width_mask;
	int m_log_width;

	int m_xxhash_seed;

	uint32_t **m_counters;
	uint32_t **m_counters_fails;
	uint32_t *m_pool_fails;
	uint32_t *m_pool_fails2;
	bool **m_pool_failures;

	int m_height_plus_1_over_2;

	bool mf_check_pool_failure(uint32_t* pool_pointer);

	int m_pool_bit_size;
	int m_counters_per_pool;
	int m_initial_counter_size;
	int m_counter_bit_increase;

	int m_top_k_count;
	int m_top_k_min;
	int m_k;
	int failCount;
	int m_hash_size;

	int crushesNum;



	//############### cuckoo ################
	xxh::hash128_t **m_keys;
	//#######################################




	tsl::robin_map<uint64_t, tuple<string,uint32_t,uint32_t>> 
	
	
	
	
	
	
	m_top_k_map;

	void mf_add_to_sketch(const char * str, int weight);

public:

	CounterPools(int width, int height, int seed, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k);
	~CounterPools();

	void increment(const char * str);
	uint32_t query(const char * str);
	void increment1(const char * str);
	uint32_t query1(const char * str);
	void increment2(const char * str);
	uint32_t query2(const char * str);
	void increment3(const char * str);
	void increment_cus(const char * str);

	uint32_t query3(const char * str);
	float percent_of_crushed_pools();
	void countersPerPrint();


	//############### cuckoo ################
	void increment_cuckoo(const char * str);
	void increment_cuckoo_v(xxh::hash128_t hashes, uint32_t val);
	//#######################################


	void write_to_file(string fn);
	void read_from_file(string fn);

	uint32_t read_counter(int i, int j);
	bool read_pool_failure(int i, int j);

};

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

class CounterPoolsExternalTopK {

	int m_width;
	int m_height;

	int m_width_mask;
	int m_log_width;

	int m_xxhash_seed;

	uint32_t **m_counters;
	bool **m_pool_failures;

	int m_height_plus_1_over_2;

	bool mf_check_pool_failure(uint32_t* pool_pointer);

	int m_pool_bit_size;
	int m_counters_per_pool;
	int m_initial_counter_size;
	int m_counter_bit_increase;


	bool m_augmented;
	CounterPoolsTopK topK;

	void mf_add_to_sketch(const char * str, int weight);

public:

	CounterPoolsExternalTopK(int width, int height, int seed, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k);
	~CounterPoolsExternalTopK();

	void increment(const char * str);
	uint32_t query(const char * str);

};

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

#endif

