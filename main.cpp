#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <string>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include <filesystem>
#include <unistd.h>

// #include <windows.h>

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// From SALSA

#include "SalsaProject/Defs.hpp"
#include "SalsaProject/CMSTests.hpp"
#include "SalsaProject/CountSketchTests.hpp"
#include "SalsaProject/PyramidSketch/PyramidTests.hpp"
#include "SalsaProject/BobHash.hpp"


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#include "tests.hpp"

#include "GradCountSketch.h"

using namespace std;
namespace fs = std::filesystem;
using namespace fs;

int main(int argc, char* argv[])
{

	if (argc < 8) {
		cout << "Usage Error:\n";
		cout << "argv[1]: int N\n";
		cout << "argv[2]: int Seed\n";
		cout << "argv[3]: int Alg\n";
		cout << "argv[4]: bool weighted_experiment\n";
		cout << "argv[5]: String Trace\n";
		cout << "argv[6]: int RowCounterNum\n";
		cout << "argv[7]: int RowsNum\n";
		//system("pause");
		return 0;
	}

	// Read arguments
	int N = atoi(argv[1]);
	int Seed = atoi(argv[2]);
	int Alg = atoi(argv[3]);
	bool weighted_experiment = (atoi(argv[4]) == 0) ? false : true;
	string Trace = string(argv[5]);
	int RowCounterNum = atoi(argv[6]);
	int RowsNum = atoi(argv[7]);

	// Pools configuration
	int pool_bit_size = 64;
	int counters_per_pool = 4;
	int initial_counter_size = 0;
	int counter_bit_increase = 1;

	int k = 0;

	if (argc > 8)
	{
		pool_bit_size = atoi(argv[8]);
		counters_per_pool = atoi(argv[9]);
		initial_counter_size = atoi(argv[10]);
		counter_bit_increase = atoi(argv[11]);
		cout << "Using pools configlll:" << pool_bit_size << "_" << counters_per_pool << "_" << initial_counter_size << "_" << counter_bit_increase << endl;

		k = atoi(argv[12]);
		cout << "Using HH table of: " << k << endl;
	}

	// Input sanity checks
	assert((weighted_experiment == false) && "no weighted simulations");

	// Print input arguments
	cout << "Received input arguments:" << endl;
	cout << "argv[1]: int N = " << N << endl;
	cout << "argv[2]: int Seed = " << Seed << endl;
	cout << "argv[3]: int Alg = " << Alg << endl;
	cout << "argv[4]: bool weighted_experiment = " << weighted_experiment << endl;
	cout << "argv[5]: String Trace = \"" << Trace << "\"" << endl;
	cout << "argv[6]: Counters in row = \"" << RowCounterNum << "\"" << endl;
	cout << "argv[7]: Number of rows = \"" << RowsNum << "\"" << endl;
	cout << endl;

	// ML or networking?
	bool network_traces = true;

	string path = "./";
	path.append(Trace);
	

	// if (Trace.rfind("youtube", 0) == 0)
	// {
	// 	fstream fileStream;
	// 	path.append("youtube/");
	// 	path.append(Trace);
	// }
	// else if (Trace.rfind("zipf", 0) == 0)
	// {
	// 	fstream fileStream;
	// 	path.append("zipf/");
	// 	path.append(Trace);
	// 	fileStream.open(path);
	// 	if (fileStream.fail()) {
	// 		// file could not be opened, return error exit status
	// 		return 1;
	// 	}
	// 	else
	// 	{
	// 		fileStream.close();
	// 	}
	// }
	// else if (Trace.rfind("grad", 0) == 0)
	// {
	// 	network_traces = false;
	// 	fstream fileStream;
	// 	path.append("gradients/");
	// 	path.append(Trace);
	// 	fileStream.open(path);
	// 	if (fileStream.fail()) {
	// 		// file could not be opened, return error exit status
	// 		return 1;
	// 	}
	// 	else
	// 	{
	// 		fileStream.close();
	// 	}
	// }
	// else
	// {
	// 	path.append(Trace);
	// }

	// if (weighted_experiment)
	// {
	// 	path.append("_weighted");
	// }

	// Asserts
	assert((Seed < 1229) && "Seed too large for BobHash");

	// Seed 
	srand(Seed);

	printf("bb:path ");
	path.append(".bin");
	printf("%s\n",path.c_str());

	// Read file?
	ifstream f(path, ios::binary);

	char* data = NULL;
	uint16_t* data_weights = NULL;

	if (network_traces)
	{
		int64_t read_so_far = 0;
		int64_t remaining = N;

		data = new char[FT_SIZE * N]();

		if (weighted_experiment)
		{
			assert((N <= 98000000) && "N too large for used traces");

			data_weights = new uint16_t[N]();
			char* temp_data_buffer = new char[100000 * FT_SIZE_WEIGHTS]();

			while (remaining > 0)
			{
				int64_t to_read = remaining > 100000 ? 100000 : remaining;
				f.read(temp_data_buffer, to_read * FT_SIZE_WEIGHTS);

				for (int i = 0; i < to_read; ++i)
				{
					memcpy(data + (i + read_so_far)* FT_SIZE, temp_data_buffer + i * FT_SIZE_WEIGHTS, FT_SIZE);

					uint8_t* p1 = (uint8_t*)(temp_data_buffer + FT_SIZE_WEIGHTS * i + FT_SIZE);
					uint8_t* p2 = (uint8_t*)(temp_data_buffer + FT_SIZE_WEIGHTS * i + FT_SIZE + 1);
					data_weights[read_so_far + i] = (*p1) * 256 + *p2;
				}

				remaining -= to_read;
				read_so_far += to_read;
			}

			delete[] temp_data_buffer;
		}

		else
		{
			while (remaining > 0)
			{
				int64_t to_read = remaining > 100000 ? 100000 : remaining;
				f.read(data + read_so_far * FT_SIZE, to_read * FT_SIZE);
				remaining -= to_read;
				read_so_far += to_read;
			}
		}

		assert((read_so_far == N) && "Read trace file error");
		assert((remaining == 0) && "Read trace file error");
	}
	else
	{
		int64_t read_so_far = 0;
		int64_t remaining = N;

		data = new char[GRAD_ELEMENT_SIZE * N]();

		while (remaining > 0)
		{
			int64_t to_read = remaining > 100000 ? 100000 : remaining;
			f.read(data + read_so_far * GRAD_ELEMENT_SIZE, to_read * GRAD_ELEMENT_SIZE);
			remaining -= to_read;
			read_so_far += to_read;
		}

		assert((read_so_far == N) && "Read trace file error");
		assert((remaining == 0) && "Read trace file error");
	}
	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	if (is_directory(Trace.c_str()) || create_directory(Trace.c_str()))
	{
		// Enter
		if (chdir(Trace.c_str())!=-1)
		{
			cout << "Directory change ok" << endl;
		}	
		else
		{
			cout << "Directory change failed" << endl;
			return 0;
		}			
	}
	else
	{
		// Failed to create directory
		cout << "Failed to create trace directory" << endl;
		return 0;
	}

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	//test_cms_error_on_arrival_sanity_64_4_0_1(N, RowCounterNum << 1, RowsNum, Seed, data);
	
	//test_cms_error_on_arrival_pools(N, RowCounterNum << 1, RowsNum, Seed, data, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, 0);
	//test_cms_error_on_arrival_pools(N, RowCounterNum << 1, RowsNum, Seed, data, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);
	//test_cms_error_on_arrival_pools_etk(N, RowCounterNum << 1, RowsNum, Seed, data, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);

	// we have a bug!!! TODO: run with sanity with univ
	//test_cms_error_on_arrival_pools_sanity(N, RowCounterNum << 1, RowsNum, Seed, data, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);
	//test_cms_error_on_arrival_top_k_sanity(N, RowCounterNum << 1, RowsNum, Seed, data, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);

	//test_cms_rw_files_pools(N, RowCounterNum << 1, RowsNum, Seed, data, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);



	// GradCountSketchFixedPoint<double, int64_t> cs_int(RowCounterNum, RowsNum, Seed, N);
	// cs_int.process_or_read_from_file(path.append("_sketched"), path);
	// cs_int.compareTopK(1000, (double*)data);


	//GradCountSketchFixedPointFL<double, int64_t> test(RowCounterNum, RowsNum, Seed, N, 10, 1, "prefix");

	/*
	cout << "%%%%%%%%%%%%%%%%%%%%%%" << endl;
	cout << "%%%%%%% cs_int %%%%%%%" << endl;
	cout << "%%%%%%%%%%%%%%%%%%%%%%" << endl;

	GradCountSketchFixedPoint<double, int64_t> cs_int(RowCounterNum, RowsNum, Seed);
	cs_int.process((double*)data, N);
	cs_int.compareTopK(1000, N, (double*)data);
	cs_int.write_to_file("cs_int_test");

	cout << "%%%%%%%%%%%%%%%%%%%%%%" << endl;

	GradCountSketchFixedPoint<double, int64_t> cs_int_2(RowCounterNum, RowsNum, Seed);
	cs_int_2.read_from_file("cs_int_test");
	cs_int_2.compareTopK(1000, N, (double*)data);
	*/

	/*
	cout << "%%%%%%%%%%%%%%%%%%%%%%" << endl;
	cout << "%%%%%%%% cs %%%%%%%%%%" << endl;
	cout << "%%%%%%%%%%%%%%%%%%%%%%" << endl;

	GradCountSketch<double> cs(RowCounterNum, RowsNum, Seed);
	cs.process((double*)data, N);
	cs.compareTopK(1000, N, (double*)data);
	*/

	//system("pause");


	if (Alg == 0)
	{
		/* error */

		// CMS baseline (32-bit counters)
		//test_cms_error_on_arrival_pools_cus_HH(N, RowCounterNum, RowsNum, Seed, data, pool_bit_size,  counters_per_pool,  initial_counter_size,  counter_bit_increase,k);
		//test_cms_error_on_arrival_64_4_0_1(N, RowCounterNum, RowsNum, Seed, data);
		//test_cms_error_on_arrival_pools_HH(N, RowCounterNum, RowsNum, Seed, data, pool_bit_size,  counters_per_pool,  initial_counter_size,  counter_bit_increase,k);
		//test_cms_error_on_arrival_pools_HH_NoFail1(N, RowCounterNum, RowsNum, Seed, data, pool_bit_size,  counters_per_pool,  initial_counter_size,  counter_bit_increase,k);
		test_cms_error_on_arrival_pools_HH_NoFail2(N, RowCounterNum, RowsNum, Seed, data, pool_bit_size,  counters_per_pool,  initial_counter_size,  counter_bit_increase,k);
		//test_cms_error_on_arrival_pools_HH_NoFail3(N, RowCounterNum, RowsNum, Seed, data, pool_bit_size,  counters_per_pool,  initial_counter_size,  counter_bit_increase,k);

		// CP
		// test_cms_error_on_arrival_64_4_0_1(N, RowCounterNum << 1, RowsNum, Seed, data);

		// Salsa max (start with 8-bit counters)
		// test_maximum_salsa_baseline_cms_error_on_arrival(N, RowCounterNum << 2, RowsNum, Seed, data);
	}
	// else if (Alg == 1)
	// {
	// 	/* speed */

	
	// 	// CP
	 	//test_cms_speed_64_4_0_1(N, RowCounterNum, RowsNum, Seed, data);

	// 	// Salsa max (start with 8-bit counters)
	 	//test_maximum_salsa_baseline_cms_speed(N, RowCounterNum, RowsNum, Seed, data);



		// test_aee_cms_speed(N, RowCounterNum, RowsNum, Seed, data);
		// // 	// CMS baseline (32-bit counters)
		// test_baseline_cms_speed(N, RowCounterNum, RowsNum, Seed, data);
		// // 	// CMS salsa baseline
		// test_salsa_baseline_cms_speed(N, RowCounterNum, RowsNum, Seed, data);



	// }
	// else if (Alg == 2)
	// {
	// 	// CP-SIM
	// 	test_cms_error_on_arrival_pools(N, RowCounterNum << 1, RowsNum, Seed, data, pool_bit_size, counters_per_pool, initial_counter_size, counter_bit_increase, k);
	// }
	////system("pause");

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	delete[] data;
	delete[] data_weights;

	return 0;

}





