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

#include "../Defs.hpp"

#include "PCMSketch.h"
#include "PCUSketch.h"
#include "PCSketch.h"


#include "PyramidTests.hpp"
#include "../SalsaCMSBaseline.hpp"

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_pyramid_cms_error_on_arrival(int N, int width, int height, int seed, const char* data)
{
	int word_num = (width * height) >> 2;
	int word_size = 64;

	PCMSketch cms_baseline(word_num, height, word_size);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		long double point_error = (int64_t)cms_baseline.Query(data + i) - (int64_t)true_sizes[ft_key];

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		cms_baseline.Insert(data + i);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string fn = "test_pyramid_cms_error_on_arrival_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}






void test_pyramid_cms_error_on_arrival_HH(int N, int width, int height, int seed, const char* data)
{
	int word_num = (width * height) >> 2;
	int word_size = 64;
	PCMSketch cms_baseline(word_num, height, word_size);

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

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0, NRMSE = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	int count =0;
	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}
		ft_key_i_values[ft_key] = i;
		true_sizes[ft_key] += 1;
		cms_baseline.Insert(data + i);
	}

	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = ((long double)cms_baseline.Query(data + ft_key_i_values[it->first])) - it->second;
			long double current_query = cms_baseline.Query(data + ft_key_i_values[it->first]);
			if((long double)abs(current_query - true_sizes[it->first]) / (long double)true_sizes[it->first] > 1 &&  (long double)abs(current_query - true_sizes[it->first]) > 10000000){
				cms_baseline.Query(data + ft_key_i_values[it->first]);
				// cout << "ERROR: " <<  (long double)abs(current_query - true_sizes[it->first]) / (long double)true_sizes[it->first] <<	"  cuur: "	<< current_query << " key: " << it->first <<endl;		
			}

		L1e += abs(error);
		L2e += error * error;
		if(it->second != 0){
			RL1e += abs(error) / it->second;
		}
		long double prev = L_max;
		L_max = L_max > error ? L_max : error;
		// if(prev < 100000000 && L_max > 100000000){
		// 	cout << "true: " <<  true_sizes[it->first] <<	"  query: "	<< current_query << "  error:" << error  << " key: " << it->first <<endl;		
		// }


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

		// cout << "threshold: " <<  thresholds[i] << "   HHNUM: " << threshold_to_hh_ft_keys[thresholds[i]].size() <<endl;

		long double current_hh_rl = 0;
		for (int j = 0; j < threshold_to_hh_ft_keys[thresholds[i]].size(); ++j)
		{
			uint64_t current_ftkey = threshold_to_hh_ft_keys[thresholds[i]][j];
			long double current_query = cms_baseline.Query(data + ft_key_i_values[current_ftkey]);



			current_hh_rl += (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey];
			//  if(i==8){
			//  	cout << "qury: " << current_query << "    size: " << true_sizes[current_ftkey] << endl;
			//  }
		}
		current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		HHre.push_back(current_hh_rl);
		
		// if(threshold_to_hh_ft_keys[thresholds[i]].size() == 0){
		// 	HHre.push_back(0);
		// }
		// else{
		// 	current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		// 	HHre.push_back(current_hh_rl);
		// }
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);
	NRMSE = L2e/N; // On araival
	RL1e /= N; // HH

	ofstream results_file;
	string fn = "pyramid_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	//results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << endl;
	results_file << "OnAraival\t"<< NRMSE << endl;

	
	ofstream results_file_hh;
	string fn_hh = "pyramid_hh_seed_";
	fn_hh.append(to_string(seed));
	fn_hh.append(".txt");
	results_file_hh.open(fn_hh, ofstream::out | ofstream::app);

	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << endl;
	for (int i = 0; i < HHre.size(); ++i)
	{
		results_file_hh << "\tThreshold\t" << thresholds[i] << "\tRelError\t" << HHre[i] << endl;
	}
	results_file_hh << endl;
}


void test_pyramid_cms_speed(int N, int width, int height, int seed, const char* data)
{

	int word_num = (width * height) >> 2;
	int word_size = 64;

	PCMSketch cms_baseline(word_num, height, word_size);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.Insert(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_pyramid_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_pyramid_cms_speed_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////

// void test_pyramid_cus_error_on_arrival(int N, int width, int height, int seed, const char* data)
// {
// 	int word_num = (width * height) >> 2;
// 	int word_size = 64;

// 	PCUSketch cms_baseline(word_num, height, word_size);

// 	unordered_map<uint64_t, uint64_t> true_sizes;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	BOBHash ft_to_bobkey_1(seed);
// 	BOBHash ft_to_bobkey_2(seed + 17);

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{

// 		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}

// 		long double point_error = (int64_t)cms_baseline.Query(data + i) - (int64_t)true_sizes[ft_key];

// 		L1e += abs(point_error);
// 		L2e += (point_error * point_error);
// 		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

// 		true_sizes[ft_key] += 1;
// 		cms_baseline.Insert(data + i);
// 	}

// 	L1e /= N;
// 	L2e /= N;
// 	L2e = sqrt(L2e);

// 	ofstream results_file;
// 	string fn = "test_pyramid_cus_error_on_arrival_seed_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
// }

// void test_pyramid_cus_speed(int N, int width, int height, int seed, const char* data)
// {

// 	int word_num = width >> 1;
// 	int word_size = 64;

// 	PCUSketch cms_baseline(word_num, height, word_size);

// 	int64_t stop_loop = N * FT_SIZE;

// 	auto start = chrono::steady_clock::now();
// 	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
// 	{
// 		cms_baseline.Insert(data + i);
// 	}
// 	auto end = chrono::steady_clock::now();

// 	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
// 	cout << "test_pyramid_cus_speed: Elapsed time in milliseconds : "
// 		<< time / 1000
// 		<< " ms" << endl;

// 	ofstream results_file;
// 	string fn = "test_pyramid_cus_speed_seed_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
// }

// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////

// void test_pyramid_cs_error_on_arrival(int N, int width, int height, int seed, const char* data)
// {
// 	int word_num = (width * height) >> 2;
// 	int word_size = 64;

// 	PCSketch cms_baseline(word_num, height, word_size);

// 	unordered_map<uint64_t, uint64_t> true_sizes;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0;

// 	BOBHash ft_to_bobkey_1(seed);
// 	BOBHash ft_to_bobkey_2(seed + 17);

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{

// 		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}

// 		long double point_error = (int64_t)cms_baseline.Query(data + i) - (int64_t)true_sizes[ft_key];

// 		L1e += abs(point_error);
// 		L2e += (point_error * point_error);
// 		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

// 		true_sizes[ft_key] += 1;
// 		cms_baseline.Insert(data + i);
// 	}

// 	L1e /= N;
// 	L2e /= N;
// 	L2e = sqrt(L2e);

// 	ofstream results_file;
// 	string fn = "test_pyramid_cs_error_on_arrival_seed_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
// }

// void test_pyramid_cs_speed(int N, int width, int height, int seed, const char* data)
// {

// 	int word_num = width >> 1;
// 	int word_size = 64;

// 	PCSketch cms_baseline(word_num, height, word_size);

// 	int64_t stop_loop = N * FT_SIZE;

// 	auto start = chrono::steady_clock::now();
// 	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
// 	{
// 		cms_baseline.Insert(data + i);
// 	}
// 	auto end = chrono::steady_clock::now();

// 	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
// 	cout << "test_pyramid_cs_speed: Elapsed time in milliseconds : "
// 		<< time / 1000
// 		<< " ms" << endl;

// 	ofstream results_file;
// 	string fn = "test_pyramid_cs_speed_seed_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
// }

// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////

// void test_pyramid_cms_final_error(int N, int width, int height, int seed, const char* data)
// {

// 	int word_num = (width * height) >> 2;
// 	int word_size = 64;

// 	PCMSketch cms_baseline(word_num, height, word_size, seed);

// 	unordered_map<uint64_t, uint64_t> true_sizes;
// 	unordered_map<uint64_t, uint64_t> ft_key_i_values;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0;

// 	BOBHash ft_to_bobkey_1(seed);
// 	BOBHash ft_to_bobkey_2(seed + 17);

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{

// 		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}
// 		ft_key_i_values[ft_key] = i;
// 		true_sizes[ft_key] += 1;
// 		cms_baseline.Insert(data + i);
// 	}

// 	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
// 	{
// 		long double error = cms_baseline.Query(data + ft_key_i_values[it->first]) - it->second;
// 		L1e += abs(error);
// 		L2e += error * error;
// 		RL1e += abs(error) / it->second;
// 		L_max = L_max > error ? L_max : error;
// 	}

// 	L1e /= true_sizes.size();
// 	L2e /= true_sizes.size();
// 	L2e = sqrt(L2e);
// 	RL1e /= true_sizes.size();

// 	ofstream results_file;
// 	string fn = "test_pyramid_cms_final_error_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;

// }

// void test_salsa_cms_final_error(int N, int width, int height, int seed, const char* data)
// {

// 	MaximumSalsaCMSBaseline cms_baseline;
// 	cms_baseline.initialize(width, height, seed);

// 	unordered_map<uint64_t, uint64_t> true_sizes;
// 	unordered_map<uint64_t, uint64_t> ft_key_i_values;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0;

// 	BOBHash ft_to_bobkey_1(seed);
// 	BOBHash ft_to_bobkey_2(seed + 17);

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{

// 		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}
// 		ft_key_i_values[ft_key] = i;
// 		true_sizes[ft_key] += 1;
// 		cms_baseline.increment(data + i);
// 	}

// 	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
// 	{
// 		long double error = cms_baseline.query(data + ft_key_i_values[it->first]) - it->second;
// 		L1e += abs(error);
// 		L2e += error * error;
// 		RL1e += abs(error) / it->second;
// 		L_max = L_max > error ? L_max : error;
// 	}

// 	L1e /= true_sizes.size();
// 	L2e /= true_sizes.size();
// 	L2e = sqrt(L2e);
// 	RL1e /= true_sizes.size();

// 	ofstream results_file;
// 	string fn = "test_salsa_cms_final_error_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;

// }

// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////

// void test_pyramid_cus_final_error(int N, int width, int height, int seed, const char* data)
// {

// 	int word_num = (width * height) >> 2;
// 	int word_size = 64;

// 	PCUSketch cms_baseline(word_num, height, word_size);

// 	unordered_map<uint64_t, uint64_t> true_sizes;
// 	unordered_map<uint64_t, uint64_t> ft_key_i_values;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0;

// 	BOBHash ft_to_bobkey_1(seed);
// 	BOBHash ft_to_bobkey_2(seed + 17);

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{

// 		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}
// 		ft_key_i_values[ft_key] = i;
// 		true_sizes[ft_key] += 1;
// 		cms_baseline.Insert(data + i);
// 	}

// 	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
// 	{
// 		long double error = cms_baseline.Query(data + ft_key_i_values[it->first]) - it->second;
// 		L1e += abs(error);
// 		L2e += error * error;
// 		RL1e += abs(error) / it->second;
// 		L_max = L_max > error ? L_max : error;
// 	}

// 	L1e /= true_sizes.size();
// 	L2e /= true_sizes.size();
// 	L2e = sqrt(L2e);
// 	RL1e /= true_sizes.size();

// 	ofstream results_file;
// 	string fn = "test_pyramid_cus_final_error_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;

// }

// void test_salsa_cus_final_error(int N, int width, int height, int seed, const char* data)
// {

// 	MaximumSalsaCUSBaseline cms_baseline;
// 	cms_baseline.initialize(width, height, seed);

// 	unordered_map<uint64_t, uint64_t> true_sizes;
// 	unordered_map<uint64_t, uint64_t> ft_key_i_values;

// 	int64_t stop_loop = N * FT_SIZE;

// 	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0;

// 	BOBHash ft_to_bobkey_1(seed);
// 	BOBHash ft_to_bobkey_2(seed + 17);

// 	for (int64_t i = 0; i < stop_loop; i += 13)
// 	{

// 		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

// 		if (true_sizes.find(ft_key) == true_sizes.end())
// 		{
// 			true_sizes[ft_key] = 0;
// 		}
// 		ft_key_i_values[ft_key] = i;
// 		true_sizes[ft_key] += 1;
// 		cms_baseline.increment(data + i);
// 	}

// 	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
// 	{
// 		long double error = cms_baseline.query(data + ft_key_i_values[it->first]) - it->second;
// 		L1e += abs(error);
// 		L2e += error * error;
// 		RL1e += abs(error) / it->second;
// 		L_max = L_max > error ? L_max : error;
// 	}

// 	L1e /= true_sizes.size();
// 	L2e /= true_sizes.size();
// 	L2e = sqrt(L2e);
// 	RL1e /= true_sizes.size();

// 	ofstream results_file;
// 	string fn = "test_salsa_cus_final_error_";
// 	fn.append(to_string(seed));
// 	fn.append(".txt");
// 	results_file.open(fn, ofstream::out | ofstream::app);
// 	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;

// }

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////