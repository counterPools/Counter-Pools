#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>

#include "xxhash.hpp"
#include "CounterPoolsDefs.h"
#include "SalsaProject/topK.hpp"

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

template<class Ctype>
class GradCountSketch
{

private:

	int m_width;
	int m_height;

	int m_seed;

	int m_width_and_sign_mask;

	Ctype **m_counters;
	Ctype *m_counter_values;

	Ctype MedianCounterValue()
	{
		sort(m_counter_values, m_counter_values + m_height);
		return m_counter_values[m_height >> 1];
	}

public:

	GradCountSketch(int width, int height, int seed);
	~GradCountSketch();

	void process(const Ctype * sequence, int len);
	Ctype query(int element);

	vector<pair<int, Ctype>> queryTopK(int k, int len);

	// for debug
	inline Ctype read_row_value(int index, int row)
	{
		return m_counters[row][index];
	}
	inline Ctype read_row_value_element(int element, int row)
	{
		uint64_t hash = xxh::xxhash3<64>(&element, sizeof(element), (m_seed << 5) + row);
		int index = (hash & m_width_and_sign_mask) >> 1;
		return m_counters[row][index] * (1 - 2 * (hash & 0b1));
	}

	// testing
	void compareTopK(int k, int len, const Ctype * sequence);

};

///////////////////////////////////////////////////////////////////////////////////

template<class Ctype>
GradCountSketch<Ctype>::GradCountSketch(int width, int height, int seed)
{
	m_width = width;
	m_height = height;

	m_width_and_sign_mask = (width << 1) - 1;

	assert((m_height % 2) == 1 && "We assume that height is odd");
	assert((m_width & (m_width - 1)) == 0 && "We assume that width is a power of 2");
	assert(m_width >= 4 && "We assume that (width % 4 == 0)");
	assert(m_width < (1 << 30) && "We assume that you are not a database lunatic");

	m_counters = new Ctype*[height];

	m_counter_values = new Ctype[height];

	for (int i = 0; i < height; ++i)
	{
		m_counters[i] = new Ctype[width]();
	}
}

template<class Ctype>
GradCountSketch<Ctype>::~GradCountSketch()
{
	for (int i = 0; i < m_height; ++i)
	{
		delete[] m_counters[i];
	}
	delete[] m_counters;

	delete[] m_counter_values;
}

template<class Ctype>
void GradCountSketch<Ctype>::process(const Ctype * sequence, int len)
{
	for (int element = 0; element < len; ++element)
	{
		for (int row = 0; row < m_height; ++row) 
		{
			uint64_t hash = xxh::xxhash3<64>(&element, sizeof(element), (m_seed << 5) + row);

			int index = (hash & m_width_and_sign_mask) >> 1;
			int sign = 1 - 2 * (hash & 0b1);

			m_counters[row][index] += sign * sequence[element];

			/*
			if (read_row_value_element(8042, row) != 0)
			{
				cout << row << endl;
				cout << index << endl;
				cout << (hash & m_width_and_sign_mask) << endl;
				int nelement = 8042;
				cout << (xxh::xxhash3<64>(&nelement, sizeof(nelement), (m_seed << 5) + row) & m_width_and_sign_mask) << endl;
				system("pause");
			}
			*/

		}
	}
}

template<class Ctype>
Ctype GradCountSketch<Ctype>::query(int element)
{
	for (int row = 0; row < m_height; ++row)
	{

		uint64_t hash = xxh::xxhash3<64>(&element, sizeof(element), (m_seed << 5) + row);

		int index = (hash & m_width_and_sign_mask) >> 1;
		int sign = 1 - 2 * (hash & 0b1);

		m_counter_values[row] = sign * m_counters[row][index];

	}
	return MedianCounterValue();
}

template<class Ctype>
inline vector<pair<int, Ctype>> GradCountSketch<Ctype>::queryTopK(int k, int len)
{
	orderedMapTopK<int, Ctype> topK(k);
	for (int element = 0; element < len; ++element)
	{
		Ctype value = query(element);
		topK.update(element, abs(value));
	}
	return topK.items();
}

template<class Ctype>
inline void GradCountSketch<Ctype>::compareTopK(int k, int len, const Ctype * sequence)
{
	long double l2s = 0;
	orderedMapTopK<int, Ctype> topK(k);
	for (int element = 0; element < len; ++element)
	{
		l2s += sequence[element] * sequence[element];
		Ctype value = sequence[element];
		topK.update(element, abs(value));
	}
	vector<pair<int, Ctype>> true_topK = topK.items();
	vector<pair<int, Ctype>> sketch_topK = queryTopK(k, len);

	auto result = *std::min_element(true_topK.cbegin(), true_topK.cend(), [](const auto& lhs, const auto& rhs) {
		return lhs.second < rhs.second;
	});

	double topk_threshold = result.second - 0.0001*sqrtl(l2s);
	double hh_threshold = 0.01*sqrtl(l2s);
	cout << "topk_threshold: " << topk_threshold << endl;
	cout << "result.second: " << result.second << endl;
	cout << "sqrtl(l2s): " << sqrtl(l2s) << endl;

	int hits = 0;
	long double sse = 0;
	long double sre = 0;
	long double hh_sre = 0;
	int hh_hits = 0;
	for (int i = 0; i < k; ++i)
	{
		long double error = abs(sequence[sketch_topK[i].first]) - abs(sketch_topK[i].second);
		//cout << sequence[sketch_topK[i].first] << "\t" << sketch_topK[i].second << endl;
		if (abs(sequence[sketch_topK[i].first]) >= topk_threshold)
		{
			++hits;
			if (abs(sequence[true_topK[i].first]) >= hh_threshold)
			{
				++hh_hits;
				hh_sre += abs(error / sequence[sketch_topK[i].first]);
			}
		}
		sse += error * error;
		sre += abs(error / sequence[sketch_topK[i].first]);
	}

	cout << "hits: " << hits << endl;
	cout << "mse: " << sse / k << endl;
	cout << "rmse: " << sqrtl(sse / k) << endl;
	cout << "are: " << sre / k << endl;
	cout << "hh are: " << hh_sre / hh_hits << endl;

	/*
	int min_estimated_hh_index = abs(sketch_topK.begin()->first);
	for (int i = 1; i < k; ++i)
	{
		if (abs(sequence[sketch_topK[i].first]) < abs(sequence[min_estimated_hh_index]))
		{
			min_estimated_hh_index = sketch_topK[i].first;
		}
	}
	cout << result.second << "\t" << abs(sequence[min_estimated_hh_index]) << endl;
	cout << sqrtl(l2s) << endl;
	*/

}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

template<class CoordinateType, class CounterType>
class GradCountSketchFixedPoint
{

private:

	int m_width;
	int m_height;

	int m_seed;

	int m_width_and_sign_mask;

	bool m_processed;

	CounterType **m_counters;
	CounterType *m_counter_values;

	CoordinateType MedianCounterValue()
	{
		sort(m_counter_values, m_counter_values + m_height);
		CounterType median = m_counter_values[m_height >> 1];
		CoordinateType result = ((((CoordinateType)median) - 0.5) / ((uint64_t)1 << m_counter_type_bits))*m_range;
		return result;
	}

	CoordinateType m_min;
	CoordinateType m_max;
	CoordinateType m_range;
	int m_counter_type_bits;

	int m_gradient_dimension;

public:

	GradCountSketchFixedPoint(int width, int height, int seed, int gradient_dimension);
	~GradCountSketchFixedPoint();

	void process(const CoordinateType * sequence);
	CoordinateType query(int element);

	vector<pair<int, CoordinateType>> queryTopK(int k);

	// for debug
	inline CounterType read_row_value(int index, int row)
	{
		return m_counters[row][index];
	}
	inline CounterType read_row_value_element(int element, int row)
	{
		uint64_t hash = xxh::xxhash3<64>(&element, sizeof(element), (m_seed << 5) + row);
		int index = (hash & m_width_and_sign_mask) >> 1;
		return m_counters[row][index] * (1 - 2 * (hash & 0b1));
	}

	// testing
	void compareTopK(int k, const CoordinateType * sequence);

	void write_sketch_to_file(string fn);
	void read_sketch_from_file(string fn);

	void read_gradient_from_file(string fn, CoordinateType* buff);

	void process_or_read_from_file(string sketch_file_name, string gradient_file_name);

	void reset();

};

///////////////////////////////////////////////////////////////////////////////////

template<class CoordinateType, class CounterType>
GradCountSketchFixedPoint<CoordinateType, CounterType>::GradCountSketchFixedPoint(int width, int height, int seed, int gradient_dimension)
{
	m_width = width;
	m_height = height;

	m_seed = seed;

	m_gradient_dimension = gradient_dimension;

	m_processed = false;

	m_width_and_sign_mask = (width << 1) - 1;

	assert((m_height % 2) == 1 && "We assume that height is odd");
	assert((m_width & (m_width - 1)) == 0 && "We assume that width is a power of 2");
	assert(m_width >= 4 && "We assume that (width % 4 == 0)");
	assert(m_width < (1 << 30) && "We assume that you are not a database lunatic");

	m_counters = new CounterType*[height];

	m_counter_values = new CounterType[height];

	for (int i = 0; i < height; ++i)
	{
		m_counters[i] = new CounterType[width]();
	}

	m_counter_type_bits = 48;
}

template<class CoordinateType, class CounterType>
GradCountSketchFixedPoint<CoordinateType, CounterType>::~GradCountSketchFixedPoint()
{
	for (int i = 0; i < m_height; ++i)
	{
		delete[] m_counters[i];
	}
	delete[] m_counters;

	delete[] m_counter_values;
}

template<class CoordinateType, class CounterType>
void GradCountSketchFixedPoint<CoordinateType, CounterType>::process(const CoordinateType * sequence)
{

	if (m_processed)
	{
		cout << "m_processed is true" << endl;
		exit(1);
	}
	m_processed = true;

	m_min = sequence[0];
	m_max = sequence[0];
	for (int element = 1; element < m_gradient_dimension; ++element)
	{
		if (m_min > sequence[element])
		{
			m_min = sequence[element];
		}
		else if (m_max < sequence[element])
		{
			m_max = sequence[element];
		}
	}

	m_range = abs(m_max) > abs(m_min) ? abs(m_max) : abs(m_min);
	for (int element = 0; element < m_gradient_dimension; ++element)
	{
		// random rounding?
		CounterType weight = (CounterType)(((sequence[element] / m_range)) * ((uint64_t)1 << m_counter_type_bits) + 0.5);

		for (int row = 0; row < m_height; ++row)
		{
			uint64_t hash = xxh::xxhash3<64>(&element, sizeof(element), (m_seed << 5) + row);

			int index = (hash & m_width_and_sign_mask) >> 1;
			int sign = 1 - 2 * (hash & 0b1);

			m_counters[row][index] += sign * weight;
		}
	}
}

template<class CoordinateType, class CounterType>
CoordinateType GradCountSketchFixedPoint<CoordinateType, CounterType>::query(int element)
{
	for (int row = 0; row < m_height; ++row)
	{

		uint64_t hash = xxh::xxhash3<64>(&element, sizeof(element), (m_seed << 5) + row);

		int index = (hash & m_width_and_sign_mask) >> 1;
		int sign = 1 - 2 * (hash & 0b1);

		m_counter_values[row] = sign * m_counters[row][index];

	}
	return MedianCounterValue();
}

template<class CoordinateType, class CounterType>
inline vector<pair<int, CoordinateType>> GradCountSketchFixedPoint<CoordinateType, CounterType>::queryTopK(int k)
{
	if (!m_processed)
	{
		cout << "m_processed is false" << endl;
		exit(1);
	}
	orderedMapTopK<int, CoordinateType> topK(k);
	for (int element = 0; element < m_gradient_dimension; ++element)
	{
		CoordinateType value = query(element);
		topK.update(element, abs(value));
	}
	return topK.items();
}

template<class CoordinateType, class CounterType>
inline void GradCountSketchFixedPoint<CoordinateType, CounterType>::compareTopK(int k, const CoordinateType * sequence)
{
	long double l2s = 0;
	orderedMapTopK<int, CoordinateType> topK(k);
	for (int element = 0; element < m_gradient_dimension; ++element)
	{
		l2s += sequence[element] * sequence[element];
		CoordinateType value = sequence[element];
		topK.update(element, abs(value));
	}
	vector<pair<int, CoordinateType>> true_topK = topK.items();
	vector<pair<int, CoordinateType>> sketch_topK = queryTopK(k);

	auto result = *std::min_element(true_topK.cbegin(), true_topK.cend(), [](const auto& lhs, const auto& rhs) {
		return lhs.second < rhs.second;
	});

	double topk_threshold = result.second - 0.0001*sqrtl(l2s);
	double hh_threshold = 0.01*sqrtl(l2s);
	cout << "topk_threshold: " << topk_threshold << endl;
	cout << "result.second: " << result.second << endl;
	cout << "sqrtl(l2s): " << sqrtl(l2s) << endl;

	int hits = 0;
	long double sse = 0;
	long double sre = 0;
	long double hh_sre = 0;
	int hh_hits = 0;
	for (int i = 0; i < k; ++i)
	{
		long double error = abs(sequence[sketch_topK[i].first]) - abs(sketch_topK[i].second);

		if (abs(sequence[sketch_topK[i].first]) >= topk_threshold)
		{
			++hits;
			if (abs(sequence[true_topK[i].first]) >= hh_threshold)
			{
				++hh_hits;
				hh_sre += abs(error / sequence[sketch_topK[i].first]);
			}
		}
		sse += error * error;
		sre += abs(error / sequence[sketch_topK[i].first]);
	}

	cout << "hits: " << hits << endl;
	cout << "mse: " << sse / k << endl;
	cout << "rmse: " << sqrtl(sse / k) << endl;
	cout << "are: " << sre / k << endl;
	cout << "hh are: " << hh_sre / hh_hits << endl;

}

///////////////////////////////////////////////////////////////////////////////////

template<class CoordinateType, class CounterType>
void GradCountSketchFixedPoint<CoordinateType, CounterType>::write_sketch_to_file(string fn)
{
	ofstream results_file;
	results_file.open(fn, ofstream::out | ios::binary);
	for (int i = 0; i < m_height; ++i)
	{
		results_file.write((char*)(m_counters[i]), m_width * sizeof(CounterType));
	}
	results_file.write((char*)(&m_min), sizeof(m_min));
	results_file.write((char*)(&m_max), sizeof(m_max));
	results_file.write((char*)(&m_range), sizeof(m_range));
	results_file.close();
}

template<class CoordinateType, class CounterType>
void GradCountSketchFixedPoint<CoordinateType, CounterType>::read_sketch_from_file(string fn)
{
	ifstream results_file;
	results_file.open(fn, ifstream::in | ios::binary);
	for (int i = 0; i < m_height; ++i)
	{
		results_file.read((char*)(m_counters[i]), m_width * sizeof(CounterType));
	}
	results_file.read((char*)(&m_min), sizeof(m_min));
	results_file.read((char*)(&m_max), sizeof(m_max));
	results_file.read((char*)(&m_range), sizeof(m_range));
	if (!results_file)
	{
		cout << "failed to read file " << fn << endl;
		exit(1);
	}
	results_file.close();
	m_processed = true;
}

template<class CoordinateType, class CounterType>
inline void GradCountSketchFixedPoint<CoordinateType, CounterType>::read_gradient_from_file(string fn, CoordinateType* buff)
{
	ifstream gradient_file;
	gradient_file.open(fn, ifstream::in | ios::binary);
	gradient_file.read((char*)(buff), sizeof(CoordinateType)*m_gradient_dimension);
	if (!gradient_file)
	{
		cout << "failed to read file " << fn << endl;
		exit(1);
	}
	gradient_file.close();
}

template<class CoordinateType, class CounterType>
void GradCountSketchFixedPoint<CoordinateType, CounterType>::process_or_read_from_file(string sketch_file_name, string gradient_file_name)
{
	ifstream sketch_file;
	sketch_file.open(sketch_file_name, ifstream::in | ios::binary);
	if (!sketch_file.is_open())
	{
		CoordinateType* buff = new CoordinateType[m_gradient_dimension];
		read_gradient_from_file(gradient_file_name, buff);
		process(buff);
		cout << "call: process(buff)" << endl;
		delete[] buff;
		write_sketch_to_file(sketch_file_name);
	}
	else
	{
		sketch_file.close();
	}
	read_sketch_from_file(sketch_file_name);
}

template<class CoordinateType, class CounterType>
inline void GradCountSketchFixedPoint<CoordinateType, CounterType>::reset()
{
	m_processed = false;
	for (int row = 0; row < m_height; ++row)
	{
		memset(m_counters[row], 0, m_width * sizeof(CounterType));
	}

	m_min = 0;
	m_max = 0;
	m_range = 0;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

template<class CoordinateType, class CounterType>
class GradCountSketchFixedPointFL
{

private:

	int m_width;
	int m_height;
	int m_seed;

	int m_num_workers;
	int m_num_iterations;

	string m_file_prefix;

	GradCountSketchFixedPoint<CoordinateType, CounterType>** workers;


public:

	GradCountSketchFixedPointFL(int width, int height, int seed, int gradient_dimension, int num_workers, int num_iterations, string file_prefix);
	~GradCountSketchFixedPointFL();

};

template<class CoordinateType, class CounterType>
inline GradCountSketchFixedPointFL<CoordinateType, CounterType>::GradCountSketchFixedPointFL(int width, int height, int seed, int gradient_dimension, int num_workers, int num_iterations, string file_prefix)
{
	m_width = width;
	m_height = height;
	m_seed = seed;

	m_num_workers = num_workers;
	m_num_iterations = num_iterations;

	m_file_prefix = file_prefix;

	workers = new GradCountSketchFixedPoint<CoordinateType, CounterType>*[m_num_workers];
	for (int i = 0; i < m_num_workers; ++i)
	{
		workers[i] = new GradCountSketchFixedPoint<CoordinateType, CounterType>(width, height, seed, gradient_dimension);
	}
}

template<class CoordinateType, class CounterType>
inline GradCountSketchFixedPointFL<CoordinateType, CounterType>::~GradCountSketchFixedPointFL()
{
	for (int i = 0; i < m_num_workers; ++i)
	{
		delete workers[i];
	}
	delete[] workers;
}

// finished process_or_read_from_file... now we want to run it for all workers...