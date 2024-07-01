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

#include "Defs.hpp"

#include "PyramidSketch/PCMSketch.h"
#include "PyramidSketch/PyramidTests.hpp"

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

void test_pyramid_cms_speed(int N, int width, int height, int seed, const char* data)
{

	int word_num = width >> 1;
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

