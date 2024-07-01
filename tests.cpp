#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string>
#include <unordered_map>

// #include "SalsaProject/CMS.hpp"

 #include "CounterPools.hpp"
#include "CounterPoolsSim.hpp"
#include "tests.hpp"
#include "Cuckoo/src/CF.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

// void test_cms_error_on_arrival_sanity_64_4_0_1(int N, int width, int height, int seed, const char* data)
// {
// 	CounterPools_64_4_0_1 cp_cms(width, height, seed);

// 	CounterPools sim_cp_cms(width, height, seed, 64, 4, 0, 1, 0);

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{

// 		if ((int64_t)cp_cms.query(data + i) != (int64_t)sim_cp_cms.query(data + i))
// 		{
// 			cout << "CP: " << (int64_t)cp_cms.query(data + i) << endl;
// 			cout << "CP_SIM: " << (int64_t)sim_cp_cms.query(data + i) << endl;
// 			cout << "i: " << i << endl;
// 			//system("pause");
// 		}

// 		cp_cms.increment(data + i);
// 		sim_cp_cms.increment(data + i);

// 		cp_cms.query(data + i);
// 		sim_cp_cms.query(data + i);
// 	}

// }

// void test_cms_error_on_arrival_64_4_0_1(int N, int width, int height, int seed, const char* data)
// {
// 	CounterPools_64_4_0_1 cp_cms(width, height, seed);

// 	unordered_map<uint64_t, uint64_t> true_sizes;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{
// 		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}

// 		long double point_error = (int64_t)cp_cms.query(data + i) - (int64_t)true_sizes[ft_key];

// 		L1e += abs(point_error);
// 		L2e += (point_error * point_error);
// 		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

// 		true_sizes[ft_key] += 1;

// 		cp_cms.increment2(data + i);
// 		cp_cms.query(data + i);
// 	}

// 	L1e /= N;
// 	L2e /= N;
// 	L2e = sqrt(L2e);

// 	ofstream results_file;
// 	string fn = "test_cms_error_on_arrival_64_4_0_1_seed_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
// }

// void test_cms_speed_64_4_0_1(int N, int width, int height, int seed, const char* data)
// {

// 	CounterPools_64_4_0_1 cp_cms(width, height, seed);

// 	int64_t stop_loop = N * FT_SIZE;

// 	auto start = chrono::steady_clock::now();
// 	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
// 	{
// 		cp_cms.increment(data + i);
// 	}
// 	auto end = chrono::steady_clock::now();

// 	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
// 	cout << "test_cms_speed_64_4_0_1: Elapsed time in milliseconds : "
// 		<< time / 1000
// 		<< " ms" << endl;

// 	ofstream results_file;
// 	string fn = "test_cms_speed_64_4_0_1_seed_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
// }

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// void test_cms_error_on_arrival_pools(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
// {
// 	CounterPools cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);

// 	unordered_map<uint64_t, uint64_t> true_sizes;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{
// 		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);


// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}

// 		long double point_error = (int64_t)cp_cms.query(data + i) - (int64_t)true_sizes[ft_key];

// 		/*
// 		if (abs(point_error) > 500)
// 		{
// 			cout << true_sizes[ft_key] << endl;
// 			cout << cp_cms.query(data + i) << endl;

// 			//system("pause");

// 			cp_cms.query(data + i);
// 		}
// 		*/

// 		L1e += abs(point_error);
// 		L2e += (point_error * point_error);
// 		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

// 		true_sizes[ft_key] += 1;

// 		cp_cms.increment(data + i);
// 		cp_cms.query(data + i);
// 	}

// 	L1e /= N;
// 	L2e /= N;
// 	L2e = sqrt(L2e);

// 	ofstream results_file;
// 	string fn = "CP_sim_test_cms_error_on_arrival_pools_seed_";
// 	fn.append(to_string(seed));
// 	fn.append("_");
// 	fn.append(to_string(pool_bit_size));
// 	fn.append("_");
// 	fn.append(to_string(counters_per_pool));
// 	fn.append("_");
// 	fn.append(to_string(initial_counter_size));
// 	fn.append("_");
// 	fn.append(to_string(counter_bit_increase));
// 	fn.append("_k_");
// 	fn.append(to_string(k));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
// }


void test_cms_error_on_arrival_pools_HH(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
{
	 CounterPools cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);


	int l=0;
	int prev = 0;


	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	unordered_map<double, vector<uint64_t>> threshold_to_hh_ft_keys;
	
	// x axis
	vector<double> thresholds = {0.0001, 0.000178, 0.000316, 0.000562, 0.001, 0.00178, 0.00316, 0.00562, 0.01};

	for (int i = 0; i < thresholds.size(); ++i)
	{
		threshold_to_hh_ft_keys[thresholds[i]] = vector<uint64_t>();
	}

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0,  NRMSE = 0;

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}
		ft_key_i_values[ft_key] = i;
		true_sizes[ft_key] += 1;

		// if(xxh::xxhash3<64>(data + i, FT_SIZE, 40) == 7861589472628683843){
		// 	l=i;
		// }

		cp_cms.increment(data + i);

		// if(l!=0){
			
		// 	if(prev != cp_cms.query(data + l)){
		// 		prev = cp_cms.query(data + l);
		// 		cout << ft_key << " " << prev << endl;
		// 	}

			
		// }
		

	}



	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = cp_cms.query(data + ft_key_i_values[it->first]) - it->second;

		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;

		
		L_max = L_max > error ? L_max : error;

		for (int i = 0; i < thresholds.size(); ++i)
		{
			if (it->second >= thresholds[i] * N)
			{
				threshold_to_hh_ft_keys[thresholds[i]].push_back(it->first);
			}
			else 
			{
				break;
			}
		}
	}

	vector<long double> HHre;
	for (int i = 0; i < thresholds.size(); ++i)
	{
		long double current_hh_rl = 0;
		for (int j = 0; j < threshold_to_hh_ft_keys[thresholds[i]].size(); ++j)
		{
			
			uint64_t current_ftkey = threshold_to_hh_ft_keys[thresholds[i]][j];
			long double current_query = cp_cms.query(data + ft_key_i_values[current_ftkey]);
			// if(i==8){
			// 	cout << "qury: " << current_query << "    size: " << true_sizes[current_ftkey] << endl;
			// }
			current_hh_rl += (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey];

			// if((long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] > 1){
			// 		cp_cms.query(data + ft_key_i_values[current_ftkey]);
			// 		cout << "ERROR: " <<  (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] <<	"  cuur: "	<< current_query << " key: " << current_ftkey <<endl;
					 		
			// }

		}

		current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		HHre.push_back(current_hh_rl);
		
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);



	NRMSE = L2e/N; // On araival
	RL1e /= N; // HH

	ofstream results_file;
	string fn = "cms_seed_";
	fn.append(to_string(seed));
	fn.append("_");
	fn.append(to_string(pool_bit_size));
	fn.append("_");
	fn.append(to_string(counters_per_pool));
	fn.append("_");
	fn.append(to_string(initial_counter_size));
	fn.append("_");
	fn.append(to_string(counter_bit_increase));
	fn.append("_k_");
	fn.append(to_string(k));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	//results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << endl;
	results_file << "OnAraival\t"<< NRMSE << "\tCrushesPer\t" << cp_cms.percent_of_crushed_pools() << endl;

	ofstream results_file_hh;
	string fn_hh = "cms_hh_seed_";
	fn_hh.append(to_string(seed));
	fn_hh.append("_");
	fn_hh.append(to_string(pool_bit_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counters_per_pool));
	fn_hh.append("_");
	fn_hh.append(to_string(initial_counter_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counter_bit_increase));
	fn_hh.append("_k_");
	fn_hh.append(to_string(k));
	fn_hh.append(".txt");
	results_file_hh.open(fn_hh, ofstream::out | ofstream::app);

	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height<<endl;
	for (int i = 0; i < HHre.size(); ++i)
	{
		results_file_hh << "\tThreshold\t" << thresholds[i] << "\tRelError\t" << HHre[i] <<endl;
	}
	results_file_hh << endl;
}




void test_cms_error_on_arrival_pools_cus_HH(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
{
	 CounterPools cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);


	int l=0;
	int prev = 0;


	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	unordered_map<double, vector<uint64_t>> threshold_to_hh_ft_keys;
	
	// x axis
	vector<double> thresholds = {0.0001, 0.000178, 0.000316, 0.000562, 0.001, 0.00178, 0.00316, 0.00562, 0.01};

	for (int i = 0; i < thresholds.size(); ++i)
	{
		threshold_to_hh_ft_keys[thresholds[i]] = vector<uint64_t>();
	}

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0,  NRMSE = 0;

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}
		ft_key_i_values[ft_key] = i;
		true_sizes[ft_key] += 1;

		// if(xxh::xxhash3<64>(data + i, FT_SIZE, 40) == 7861589472628683843){
		// 	l=i;
		// }

		cp_cms.increment_cus(data + i);

		// if(l!=0){
			
		// 	if(prev != cp_cms.query(data + l)){
		// 		prev = cp_cms.query(data + l);
		// 		cout << ft_key << " " << prev << endl;
		// 	}

			
		// }
		

	}



	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = cp_cms.query(data + ft_key_i_values[it->first]) - it->second;
		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;

		
		L_max = L_max > error ? L_max : error;

		for (int i = 0; i < thresholds.size(); ++i)
		{
			if (it->second >= thresholds[i] * N)
			{
				threshold_to_hh_ft_keys[thresholds[i]].push_back(it->first);
			}
			else 
			{
				break;
			}
		}
	}
/*
	vector<long double> HHre;
	for (int i = 0; i < thresholds.size(); ++i)
	{
		long double current_hh_rl = 0;
		for (int j = 0; j < threshold_to_hh_ft_keys[thresholds[i]].size(); ++j)
		{
			
			uint64_t current_ftkey = threshold_to_hh_ft_keys[thresholds[i]][j];
			long double current_query = cp_cms.query(data + ft_key_i_values[current_ftkey]);
			// if(i==8){
			// 	cout << "qury: " << current_query << "    size: " << true_sizes[current_ftkey] << endl;
			// }
			current_hh_rl += (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey];

			cout << current_ftkey << " " << (int64_t)cp_cms.query(data + i) << " " << (int64_t)true_sizes[current_ftkey] <<endl;

			// if((long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] > 1){
			// 		cp_cms.query(data + ft_key_i_values[current_ftkey]);
			// 		cout << "ERROR: " <<  (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] <<	"  cuur: "	<< current_query << " key: " << current_ftkey <<endl;
					 		
			// }

		}

		current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		HHre.push_back(current_hh_rl);
	
	}
*/
	
	L1e /= N;
	cout<< "sum " << L2e <<endl;
	L2e /= N;
	cout<< "Nsum " << L2e <<endl;
	L2e = sqrt(L2e);
	cout<< "sqrtNsum " << L2e <<endl;


	NRMSE = L2e/N; // On araival
	cout<< "NRMSE " << NRMSE <<endl;
	RL1e /= N; // HH


	ofstream results_file;
	string fn = "cms_seed_";
	fn.append(to_string(seed));
	fn.append("_");
	fn.append(to_string(pool_bit_size));
	fn.append("_");
	fn.append(to_string(counters_per_pool));
	fn.append("_");
	fn.append(to_string(initial_counter_size));
	fn.append("_");
	fn.append(to_string(counter_bit_increase));
	fn.append("_k_");
	fn.append(to_string(k));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	//results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << endl;
	results_file << "OnAraival\t"<< NRMSE << "\tCrushesPer\t" << cp_cms.percent_of_crushed_pools() << endl;

/*
	ofstream results_file_hh;
	string fn_hh = "cms_hh_seed_";
	fn_hh.append(to_string(seed));
	fn_hh.append("_");
	fn_hh.append(to_string(pool_bit_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counters_per_pool));
	fn_hh.append("_");
	fn_hh.append(to_string(initial_counter_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counter_bit_increase));
	fn_hh.append("_k_");
	fn_hh.append(to_string(k));
	fn_hh.append(".txt");
	results_file_hh.open(fn_hh, ofstream::out | ofstream::app);

	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height<<endl;
	for (int i = 0; i < HHre.size(); ++i)
	{
		results_file_hh << "\tThreshold\t" << thresholds[i] << "\tRelError\t" << HHre[i] <<endl;
	}
	results_file_hh << endl;
	*/
}




	


/*

void test_cms_error_on_arrival_pools_HH_cuckoo(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
{
	//CounterPools cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);
	Cuckoo_waterLevel_HH_no_FP_SIMD_256<uint32_t, uint32_t> cp_cms(1, 2, 3, 0.6,1,100000000);
	//CF<uint32_t> cp_cms(width,4,height,1);
	

	int l=0;
	int prev = 0;


	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	unordered_map<double, vector<uint64_t>> threshold_to_hh_ft_keys;
	
	// x axis
	vector<double> thresholds = {0.0001, 0.000178, 0.000316, 0.000562, 0.001, 0.00178, 0.00316, 0.00562, 0.01};
	thresholds = {0.01};

	for (int i = 0; i < thresholds.size(); ++i)
	{
		threshold_to_hh_ft_keys[thresholds[i]] = vector<uint64_t>();
	}

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0,  NRMSE = 0;

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}
		ft_key_i_values[ft_key] = i;
		true_sizes[ft_key] += 1;

		// if(xxh::xxhash3<64>(data + i, FT_SIZE, 40) == 7861589472628683843){
		// 	l=i;
		// }
		uint32_t dd = (uint8_t)atoi(data + i);
		cp_cms.insert(dd,1);

		// if(l!=0){
			
		// 	if(prev != cp_cms.query(data + l)){
		// 		prev = cp_cms.query(data + l);
		// 		cout << ft_key << " " << prev << endl;
		// 	}

			
		// }
		

	}

	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		uint32_t d = (uint8_t)atoi(data + ft_key_i_values[it->first]);
		long double error = cp_cms.query(d) - it->second;

		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;

		
		L_max = L_max > error ? L_max : error;

		for (int i = 0; i < thresholds.size(); ++i)
		{
			if (it->second >= thresholds[i] * N)
			{
				threshold_to_hh_ft_keys[thresholds[i]].push_back(it->first);
			}
			else 
			{
				break;
			}
		}
	}

	vector<long double> HHre;
	for (int i = 0; i < thresholds.size(); ++i)
	{
		long double current_hh_rl = 0;
		for (int j = 0; j < threshold_to_hh_ft_keys[thresholds[i]].size(); ++j)
		{
			
			uint64_t current_ftkey = threshold_to_hh_ft_keys[thresholds[i]][j];
			uint32_t d = (uint8_t)atoi(data + ft_key_i_values[current_ftkey]);
			long double current_query = cp_cms.query(d);
			// if(i==8){
			// 	cout << "qury: " << current_query << "    size: " << true_sizes[current_ftkey] << endl;
			// }
			current_hh_rl += (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey];

			// if((long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] > 1){
			// 		cp_cms.query(data + ft_key_i_values[current_ftkey]);
			// 		cout << "ERROR: " <<  (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] <<	"  cuur: "	<< current_query << " key: " << current_ftkey <<endl;
					 		
			// }

		}

		current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		HHre.push_back(current_hh_rl);
		
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);



	NRMSE = L2e/N; // On araival
	RL1e /= N; // HH

	ofstream results_file;
	string fn = "cms_seed_";
	fn.append(to_string(seed));
	fn.append("_");
	fn.append(to_string(pool_bit_size));
	fn.append("_");
	fn.append(to_string(counters_per_pool));
	fn.append("_");
	fn.append(to_string(initial_counter_size));
	fn.append("_");
	fn.append(to_string(counter_bit_increase));
	fn.append("_k_");
	fn.append(to_string(k));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	//results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << endl;
	results_file << "OnAraival\t"<< NRMSE << endl;

	ofstream results_file_hh;
	string fn_hh = "cms_hh_seed_";
	fn_hh.append(to_string(seed));
	fn_hh.append("_");
	fn_hh.append(to_string(pool_bit_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counters_per_pool));
	fn_hh.append("_");
	fn_hh.append(to_string(initial_counter_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counter_bit_increase));
	fn_hh.append("_k_");
	fn_hh.append(to_string(k));
	fn_hh.append(".txt");
	results_file_hh.open(fn_hh, ofstream::out | ofstream::app);

	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height<<endl;
	for (int i = 0; i < HHre.size(); ++i)
	{
		results_file_hh << "\tThreshold\t" << thresholds[i] << "\tRelError\t" << HHre[i]<<endl;
	}
	results_file_hh << endl;
}







*/
void test_cms_error_on_arrival_pools_HH_NoFail1(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
{
	 CounterPools cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);


	int l=0;
	int prev = 0;


	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	unordered_map<double, vector<uint64_t>> threshold_to_hh_ft_keys;
	
	// x axis
	vector<double> thresholds = {0.0001, 0.000178, 0.000316, 0.000562, 0.001, 0.00178, 0.00316, 0.00562, 0.01};
	//thresholds = {0.0001};

	for (int i = 0; i < thresholds.size(); ++i)
	{
		threshold_to_hh_ft_keys[thresholds[i]] = vector<uint64_t>();
	}

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0,  NRMSE = 0;

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}
		ft_key_i_values[ft_key] = i;
		true_sizes[ft_key] += 1;

		// if(xxh::xxhash3<64>(data + i, FT_SIZE, 40) == 7861589472628683843){
		// 	l=i;
		// }

		cp_cms.increment1(data + i);
		// cout << "got: "<<cp_cms.query1(data + i) << endl;
		// cout << "true: "<< true_sizes[ft_key]<<endl<<endl;

		// if(l!=0){
			
			// if(prev != cp_cms.query1(data + l)){
			// 	prev = cp_cms.query1(data + l);
			// 	cout << ft_key << " " << prev << endl;
			// }

			
		// }
		

	}

	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = cp_cms.query1(data + ft_key_i_values[it->first]) - it->second;

		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;

		
		L_max = L_max > error ? L_max : error;

		for (int i = 0; i < thresholds.size(); ++i)
		{
			if (it->second >= thresholds[i] * N)
			{
				threshold_to_hh_ft_keys[thresholds[i]].push_back(it->first);
			}
			else 
			{
				break;
			}
		}
	}

	vector<long double> HHre;
	for (int i = 0; i < thresholds.size(); ++i)
	{
		long double current_hh_rl = 0;
		for (int j = 0; j < threshold_to_hh_ft_keys[thresholds[i]].size(); ++j)
		{
			
			uint64_t current_ftkey = threshold_to_hh_ft_keys[thresholds[i]][j];
			long double current_query = cp_cms.query1(data + ft_key_i_values[current_ftkey]);
			// if(i==8){
			// 	cout << "qury: " << current_query << "    size: " << true_sizes[current_ftkey] << endl;
			// }
			current_hh_rl += (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey];

			// if((long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] > 1){
			// 		cp_cms.query(data + ft_key_i_values[current_ftkey]);
			// 		cout << "ERROR: " <<  (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] <<	"  cuur: "	<< current_query << " key: " << current_ftkey <<endl;
					 		
			// }

		}

		current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		HHre.push_back(current_hh_rl);
		
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);



	NRMSE = L2e/N; // On araival
	RL1e /= N; // HH

	ofstream results_file;
	string fn = "cms_seed_";
	fn.append(to_string(seed));
	fn.append("_");
	fn.append(to_string(pool_bit_size));
	fn.append("_");
	fn.append(to_string(counters_per_pool));
	fn.append("_");
	fn.append(to_string(initial_counter_size));
	fn.append("_");
	fn.append(to_string(counter_bit_increase));
	fn.append("_k_NoFail1_");
	fn.append(to_string(k));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	//results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << endl;
	results_file << "OnAraival\t"<< NRMSE << endl;

	ofstream results_file_hh;
	string fn_hh = "cms_hh_seed_";
	fn_hh.append(to_string(seed));
	fn_hh.append("_");
	fn_hh.append(to_string(pool_bit_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counters_per_pool));
	fn_hh.append("_");
	fn_hh.append(to_string(initial_counter_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counter_bit_increase));
	fn_hh.append("_k__NoFail1_");
	fn_hh.append(to_string(k));
	fn_hh.append(".txt");
	results_file_hh.open(fn_hh, ofstream::out | ofstream::app);

	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height<<endl;
	for (int i = 0; i < HHre.size(); ++i)
	{
		results_file_hh << "\tThreshold\t" << thresholds[i] << "\tRelError\t" << HHre[i]<<endl;
	}
	results_file_hh << endl;
}





void test_cms_error_on_arrival_pools_HH_NoFail2(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
{
	 CounterPools cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);


	int l=0;
	int prev = 0;


	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	unordered_map<double, vector<uint64_t>> threshold_to_hh_ft_keys;
	
	// x axis
	vector<double> thresholds = {0.0001, 0.000178, 0.000316, 0.000562, 0.001, 0.00178, 0.00316, 0.00562, 0.01};

	for (int i = 0; i < thresholds.size(); ++i)
	{
		threshold_to_hh_ft_keys[thresholds[i]] = vector<uint64_t>();
	}

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0,  NRMSE = 0;

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}
		ft_key_i_values[ft_key] = i;
		true_sizes[ft_key] += 1;

		// if(xxh::xxhash3<64>(data + i, FT_SIZE, 40) == 7861589472628683843){
		// 	l=i;
		// }

		cp_cms.increment2(data + i);
		// cout << "got: "<<cp_cms.query1(data + i) << endl;
		// cout << "true: "<< true_sizes[ft_key]<<endl<<endl;

		// if(l!=0){
			
			// if(prev != cp_cms.query1(data + l)){
			// 	prev = cp_cms.query1(data + l);
			// 	cout << ft_key << " " << prev << endl;
			// }

			
		// }
		

	}

	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = cp_cms.query2(data + ft_key_i_values[it->first]) - it->second;

		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;

		
		L_max = L_max > error ? L_max : error;

		for (int i = 0; i < thresholds.size(); ++i)
		{
			if (it->second >= thresholds[i] * N)
			{
				threshold_to_hh_ft_keys[thresholds[i]].push_back(it->first);
			}
			else 
			{
				break;
			}
		}
	}

	vector<long double> HHre;
	for (int i = 0; i < thresholds.size(); ++i)
	{
		long double current_hh_rl = 0;
		for (int j = 0; j < threshold_to_hh_ft_keys[thresholds[i]].size(); ++j)
		{
			
			uint64_t current_ftkey = threshold_to_hh_ft_keys[thresholds[i]][j];
			long double current_query = cp_cms.query2(data + ft_key_i_values[current_ftkey]);
			// if(i==8){
			// 	cout << "qury: " << current_query << "    size: " << true_sizes[current_ftkey] << endl;
			// }
			current_hh_rl += (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey];

			// if((long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] > 1){
			// 		cp_cms.query(data + ft_key_i_values[current_ftkey]);
			// 		cout << "ERROR: " <<  (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] <<	"  cuur: "	<< current_query << " key: " << current_ftkey <<endl;
					 		
			// }

		}

		current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		HHre.push_back(current_hh_rl);
		
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);



	NRMSE = L2e/N; // On araival
	RL1e /= N; // HH


	
	ofstream results_file;
	string fn = "cms_seed_";
	fn.append(to_string(seed));
	fn.append("_");
	fn.append(to_string(pool_bit_size));
	fn.append("_");
	fn.append(to_string(counters_per_pool));
	fn.append("_");
	fn.append(to_string(initial_counter_size));
	fn.append("_");
	fn.append(to_string(counter_bit_increase));
	fn.append("_k_NoFail2_");
	fn.append(to_string(10));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	//results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << endl;
	results_file << "OnAraival\t"<< NRMSE << endl;



	ofstream results_file_hh;
	string fn_hh = "cms_hh_seed_";
	fn_hh.append(to_string(seed));
	fn_hh.append("_");
	fn_hh.append(to_string(pool_bit_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counters_per_pool));
	fn_hh.append("_");
	fn_hh.append(to_string(initial_counter_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counter_bit_increase));
	fn_hh.append("_k__NoFail2_poolCounter_");
	fn_hh.append(to_string(10));
	fn_hh.append(".txt");
	results_file_hh.open(fn_hh, ofstream::out | ofstream::app);

	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height<<endl;
	for (int i = 0; i < HHre.size(); ++i)
	{
		results_file_hh << "\tThreshold\t" << thresholds[i] << "\tRelError\t" << HHre[i]<<endl;
	}
	results_file_hh << endl;
}




/*



// void test_cms_error_on_arrival_pools_etk(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
// {
// 	CounterPoolsExternalTopK cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);

// 	unordered_map<uint64_t, uint64_t> true_sizes;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{
// 		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}

// 		long double point_error = (int64_t)cp_cms.query(data + i) - (int64_t)true_sizes[ft_key];

// 		if (abs(point_error) > 500)
// 		{
// 			//cout << true_sizes[ft_key] << endl;
// 			//cout << cp_cms.query(data + i) << endl;

// 			////system("pause");

// 			cp_cms.query(data + i);
// 		}

// 		L1e += abs(point_error);
// 		L2e += (point_error * point_error);
// 		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

// 		true_sizes[ft_key] += 1;

// 		cp_cms.increment(data + i);
// 		cp_cms.query(data + i);
// 	}

// 	L1e /= N;
// 	L2e /= N;
// 	L2e = sqrt(L2e);

// 	ofstream results_file; 
// 	string fn = "sim_test_cms_error_on_arrival_pools_etk_seed_";
// 	fn.append(to_string(seed));
// 	fn.append("_");
// 	fn.append(to_string(pool_bit_size));
// 	fn.append("_");
// 	fn.append(to_string(counters_per_pool));
// 	fn.append("_");
// 	fn.append(to_string(initial_counter_size));
// 	fn.append("_");
// 	fn.append(to_string(counter_bit_increase));
// 	fn.append("_hh_");
// 	fn.append(to_string(k));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
// }

// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////

// void sim_test_cms_error_on_arrival_sanity(int N, int width, int height, int seed, const char* data)
// {

// 	CounterPools sim_cp_cms(width, height, seed, 64, 4, 0, 1, 0);
// 	CounterPools sim_cp_cms_2(width, height, seed, 64, 4, 7, 4, 0);

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{

// 		if ((int64_t)sim_cp_cms.query(data + i) > (int64_t)sim_cp_cms_2.query(data + i))
// 		{
// 			cout << "CP_SIM: " << (int64_t)sim_cp_cms.query(data + i) << endl;
// 			cout << "CP_SIM_2: " << (int64_t)sim_cp_cms_2.query(data + i) << endl;
// 			cout << "i: " << i << endl;
// 			//system("pause");
// 		}

// 		sim_cp_cms.increment(data + i);
// 		sim_cp_cms_2.increment(data + i);

// 		sim_cp_cms.query(data + i);
// 		sim_cp_cms_2.query(data + i);

// 	}

// }

// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////

// void test_cms_error_on_arrival_pools_sanity(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
// {
// 	CounterPools cp_cms_0(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, 0);
// 	CounterPools cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);

// 	unordered_map<uint64_t, uint64_t> true_sizes;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{
// 		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}

// 		long double point_error_0 = (int64_t)cp_cms_0.query(data + i) - (int64_t)true_sizes[ft_key];
// 		long double point_error = (int64_t)cp_cms.query(data + i) - (int64_t)true_sizes[ft_key];

		
// 		if (abs(point_error) > abs(point_error_0))
// 		{
// 			cout << true_sizes[ft_key] << endl;
// 			cout << ft_key << endl;
// 			cout << cp_cms_0.query(data + i) << endl;
// 			cout << cp_cms.query(data + i) << endl;
// 			cout << i/13 << endl;

// 			//system("pause");
// 		}

// 		/*
// 		if (ft_key == 4391600616352669936)
// 		{
// 			cout << i << endl;
// 			cout << cp_cms_0.query(data + i) << endl;
// 			cout << cp_cms.query(data + i) << endl;
// 			//system("pause");
// 		}
// 		*/

// 		/*
// 		int bug = 5798;
// 		if ((cp_cms_0.query(data + bug) >= 2) && (cp_cms.query(data + bug) <= 1))
// 		{
// 			cout << i << endl;
// 			cout << cp_cms_0.query(data + bug) << endl;
// 			cout << cp_cms.query(data + bug) << endl;
// 			//system("pause");
// 		}
// 		*/
		

// 		L1e += abs(point_error);
// 		L2e += (point_error * point_error);
// 		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

// 		true_sizes[ft_key] += 1;

// 		cp_cms.increment(data + i);
// 		cp_cms_0.increment(data + i);
// 		cp_cms.query(data + i);
// 		cp_cms_0.query(data + i);
// 	}

// }

// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////

// void test_cms_error_on_arrival_top_k_sanity(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
// {
// 	CounterPools cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);
// 	CounterPoolsExternalTopK cp_cms_etk(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);

// 	unordered_map<uint64_t, uint64_t> true_sizes;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{
// 		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}

// 		long double point_error_etk = (int64_t)cp_cms_etk.query(data + i) - (int64_t)true_sizes[ft_key];
// 		long double point_error = (int64_t)cp_cms.query(data + i) - (int64_t)true_sizes[ft_key];

// 		/*
// 		if ((point_error_etk != point_error) && (abs(point_error_etk - point_error) > 10))
// 		{
// 			cout << true_sizes[ft_key] << endl;
// 			cout << ft_key << endl;
// 			cout << cp_cms_etk.query(data + i) << endl;
// 			cout << cp_cms.query(data + i) << endl;
// 			cout << i / 13 << endl;

// 			//system("pause");
// 		}
// 		*/

// 		/*
// 		if (ft_key == 15964518976569952023)
// 		{
// 			cout << i << endl;
// 			cout << cp_cms_etk.query(data + i) << endl;
// 			cout << cp_cms.query(data + i) << endl;
// 			cout << endl;
// 			////system("pause");
// 		}
// 		*/
		

// 		/*
// 		int bug = 2116231;
// 		if (cp_cms_etk.query(data + bug)!=cp_cms.query(data + bug))
// 		{
// 			cout << i << endl;
// 			cout << cp_cms_etk.query(data + bug) << endl;
// 			cout << cp_cms.query(data + bug) << endl;
// 			//system("pause");
// 		}
// 		*/


// 		L1e += abs(point_error);
// 		L2e += (point_error * point_error);
// 		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

// 		true_sizes[ft_key] += 1;

// 		cp_cms.increment(data + i);
// 		cp_cms_etk.increment(data + i);
// 		cp_cms.query(data + i);
// 		cp_cms_etk.query(data + i);
// 	}

// }

// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////

// void test_cms_rw_files_pools(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
// {
// 	CounterPools cp_cms_write(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);
// 	CounterPools cp_cms_read(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{
// 		cp_cms_write.increment(data + i);
// 	}

// 	string fn = "test_cms_rw_files_pools_seed_";
// 	fn.append(to_string(seed));
// 	fn.append("_");
// 	fn.append(to_string(pool_bit_size));
// 	fn.append("_");
// 	fn.append(to_string(counters_per_pool));
// 	fn.append("_");
// 	fn.append(to_string(initial_counter_size));
// 	fn.append("_");
// 	fn.append(to_string(counter_bit_increase));
// 	fn.append("_hh_");
// 	fn.append(to_string(k));
// 	fn.append(".txt");

// 	cp_cms_write.write_to_file(fn);
// 	cp_cms_read.read_from_file(fn);

// 	for (int i = 0; i < height; ++i)
// 	{
// 		for (int j = 0; j < width; ++j)
// 		{
// 			if (cp_cms_write.read_pool_failure(i, j) != cp_cms_write.read_pool_failure(i, j))
// 			{
// 				cout << "Pool Failure Error" << endl;
// 				//system("pause");
// 			}
// 			if (!cp_cms_write.read_pool_failure(i, j))
// 			{
// 				if (cp_cms_write.read_counter(i, j) != cp_cms_write.read_counter(i, j))
// 				{
// 					cout << "Counter Value Error" << endl;
// 					//system("pause");
// 				}
// 			}

// 		}
// 	}
// }



void test_cms_error_on_arrival_pools_HH_NoFail3(int N, int width, int height, int seed, const char* data, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k)
{
	 CounterPools cp_cms(width, height, seed, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);


	int l=0;
	int prev = 0;


	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	unordered_map<double, vector<uint64_t>> threshold_to_hh_ft_keys;
	
	// x axis
	vector<double> thresholds = {0.0001, 0.000178, 0.000316, 0.000562, 0.001, 0.00178, 0.00316, 0.00562, 0.01};

	for (int i = 0; i < thresholds.size(); ++i)
	{
		threshold_to_hh_ft_keys[thresholds[i]] = vector<uint64_t>();
	}

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0,  NRMSE = 0;

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = xxh::xxhash3<64>(data + i, FT_SIZE, seed);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}
		ft_key_i_values[ft_key] = i;
		true_sizes[ft_key] += 1;



		cp_cms.increment3(data + i);

		

	}

	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = cp_cms.query3(data + ft_key_i_values[it->first]) - it->second;

		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;

		
		L_max = L_max > error ? L_max : error;

		for (int i = 0; i < thresholds.size(); ++i)
		{
			if (it->second >= thresholds[i] * N)
			{
				threshold_to_hh_ft_keys[thresholds[i]].push_back(it->first);
			}
			else 
			{
				break;
			}
		}
	}

	vector<long double> HHre;
	for (int i = 0; i < thresholds.size(); ++i)
	{
		long double current_hh_rl = 0;
		for (int j = 0; j < threshold_to_hh_ft_keys[thresholds[i]].size(); ++j)
		{
			
			uint64_t current_ftkey = threshold_to_hh_ft_keys[thresholds[i]][j];
			long double current_query = cp_cms.query3(data + ft_key_i_values[current_ftkey]);
			// if(i==8){
			// 	cout << "qury: " << current_query << "    size: " << true_sizes[current_ftkey] << endl;
			// }
			current_hh_rl += (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey];

			// if((long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] > 1){
			// 		cp_cms.query(data + ft_key_i_values[current_ftkey]);
			// 		cout << "ERROR: " <<  (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey] <<	"  cuur: "	<< current_query << " key: " << current_ftkey <<endl;
					 		
			// }

		}

		current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		HHre.push_back(current_hh_rl);
		
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);



	NRMSE = L2e/N; // On araival
	RL1e /= N; // HH


	
	ofstream results_file;
	string fn = "cms_seed_";
	fn.append(to_string(seed));
	fn.append("_");
	fn.append(to_string(pool_bit_size));
	fn.append("_");
	fn.append(to_string(counters_per_pool));
	fn.append("_");
	fn.append(to_string(initial_counter_size));
	fn.append("_");
	fn.append(to_string(counter_bit_increase));
	fn.append("_k_NoFail3");
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	//results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << endl;
	results_file << "OnAraival\t"<< NRMSE << endl;



	ofstream results_file_hh;
	string fn_hh = "cms_hh_seed_";
	fn_hh.append(to_string(seed));
	fn_hh.append("_");
	fn_hh.append(to_string(pool_bit_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counters_per_pool));
	fn_hh.append("_");
	fn_hh.append(to_string(initial_counter_size));
	fn_hh.append("_");
	fn_hh.append(to_string(counter_bit_increase));
	fn_hh.append("_k__NoFail3_");
	fn_hh.append(".txt");
	results_file_hh.open(fn_hh, ofstream::out | ofstream::app);

	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height<<endl;
	for (int i = 0; i < HHre.size(); ++i)
	{
		results_file_hh << "\tThreshold\t" << thresholds[i] << "\tRelError\t" << HHre[i]<<endl;
	}
	results_file_hh << endl;

	
}
