#include <sstream> 
#include <iostream> 
#include<stdio.h> 
#include <math.h> 

#include "CounterPoolsSim.hpp"

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

CounterPools::CounterPools(int width, int height, int seed, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k) :
	m_top_k_map(4*k)
{
	m_width = width;
	m_height = height;

	m_width_mask = width - 1;

	int index = m_width_mask;
	m_log_width = 0;
	while (index >>= 1) ++m_log_width;

	assert(width > 0 && "We assume too much!");
	//assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	//assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	m_height_plus_1_over_2 = (height + 1) / 2;
	if (m_height_plus_1_over_2 > 1)
	{
		assert(width < ((uint64_t)1 << (64 / m_height_plus_1_over_2)) && "128 hash bits are not enough!");
	}

	m_counters = new uint32_t*[height];
	m_counters_fails = new uint32_t*[height];
	m_pool_failures = new bool*[height];

	// fail array size:
	m_hash_size = 0.1*width*height;

	m_pool_fails = new uint32_t[m_hash_size];
	m_pool_fails2 = new uint32_t[height*width*2];

	crushesNum =0;
	for (int i = 0; i < height; ++i)
	{
		m_counters[i] = new uint32_t[width]();
		m_counters_fails[i] = new uint32_t[width]();
		m_pool_failures[i] = new bool[width / counters_per_pool + 1]();
	}

	m_xxhash_seed = seed * 7 + 100;

	m_pool_bit_size = pool_bit_size;
	m_counters_per_pool = counters_per_pool;
	m_initial_counter_size = initial_counter_size;
	m_counter_bit_increase = counter_bit_increase;

	m_k = k;
	m_top_k_count = 0;
	m_top_k_min = 0;
	failCount=0;


	//############### cuckoo ################



	//#######################################
}

CounterPools::~CounterPools()
{
	for (int i = 0; i < m_height; ++i)
	{
		delete[] m_counters[i];
		delete[] m_counters_fails[i];
		delete[] m_pool_failures[i];
	}
	delete[] m_counters;
	delete[] m_counters_fails;
	delete[] m_pool_failures;
	delete[] m_pool_fails;
	delete[] m_pool_fails2;
}




void CounterPools::increment(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	
	int i;

	/*if (m_k > 0)
	{
		uint64_t fp = hashes.low64 ^ hashes.high64;

		auto it = m_top_k_map.find(fp);
		if (it != m_top_k_map.end()) {
			tuple<string, uint32_t, uint32_t> &value = it.value();
			++get<2>(value);
			
			if (m_top_k_min == get<1>(value) + get<2>(value) - 1)
			{
				m_top_k_min = INT32_MAX;
				for (auto it = m_top_k_map.begin(); it != m_top_k_map.end(); ++it)
				{
					auto it_value = it.value();
					if (get<1>(it_value) + get<2>(it_value) < m_top_k_min)
					{
						m_top_k_min = get<1>(it_value) + get<2>(it_value);
					}
				}
			}
			return;
		}

		uint32_t estimated_size = query(str);
		if (estimated_size >= m_top_k_min)
		{

			tuple<string, uint32_t, uint32_t> new_value(string(str, FT_SIZE), estimated_size, 1);
			if (m_top_k_count == m_k)
			{

				auto value = m_top_k_map.begin().value();
				for (auto it = m_top_k_map.begin(); it != m_top_k_map.end(); ++it) 
				{
					value = it.value();
					if (get<1>(value) + get<2>(value) == m_top_k_min)
					{
						m_top_k_map.erase(it);
						break;
					}
				}
				m_top_k_map[fp] = new_value;
				mf_add_to_sketch(get<0>(value).c_str(), get<2>(value));
			}
			else
			{
				m_top_k_map[fp] = new_value;
				++m_top_k_count;
			}

			m_top_k_min = INT32_MAX;
			for (auto it = m_top_k_map.begin(); it != m_top_k_map.end(); ++it)
			{
				auto it_value = it.value();
				if (get<1>(it_value) + get<2>(it_value) < m_top_k_min)
				{
					m_top_k_min = get<1>(it_value) + get<2>(it_value);
				}
			}

			return;
		}
	}*/


	for (i = 0; i < m_height_plus_1_over_2; ++i) {

		//int index1 = hashes.low64 & m_width_mask;
		


		int index = hashes.low64 % m_width;
		//  if((index != index1) || ((hashes.low64 >>m_log_width) != (hashes.low64 / (m_width/2)))){
		//  	cout <<"low:" <<hashes.low64 << "  width:"<<m_width <<   "  width_mask:"<<m_width_mask << "    m_log_width:"<< m_log_width  <<endl;
		//  	cout <<"index:"<<index<<"    index1: "<<index1 <<endl;
		// 	cout <<"hash:"<<(hashes.low64 / (m_width/2))<<"    hash1: "<< (hashes.low64 >>= m_log_width) <<endl;
		//  	exit(0);
		//  }


		hashes.low64 /= (m_width/2);
		//hashes.low64 >>= m_log_width;
		

		int pool_index = index / m_counters_per_pool;

	

		if (!m_pool_failures[i][pool_index])
		{
			++m_counters[i][index];
			if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
			{
				m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
				crushesNum++;
			}
		}		
	}
	for (; i < m_height; ++i) {

		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;

		if (!m_pool_failures[i][pool_index])
		{
			++m_counters[i][index];
			if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
			{
				m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
				crushesNum++;
			}
		}
	}
}






void CounterPools::increment1(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	
	int i;

	for (i = 0; i < m_height_plus_1_over_2; ++i) {

		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;

	


		++m_counters[i][index];
		if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
		{
			//m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
				--m_counters[i][index];
				m_counters_fails[i][index] = -1;
				++failCount;
			}
		}
			
	}
	for (; i < m_height; ++i) {

		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;


		++m_counters[i][index];
		if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
		{
			//m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
				--m_counters[i][index];
				m_counters_fails[i][index] = -1;
				++failCount;
			}
		}
		
	}
}


/* FAIL COUNTERS

with failing
1 counter
4 counters
array= [width*hight*2] 

*/


// with a separate no fail counters array 
void CounterPools::increment2(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	
	int i;
	for (i = 0; i < m_height_plus_1_over_2; ++i) {

		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;

	


		++m_counters[i][index];
		if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
		{
			if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
				--m_counters[i][index];
				m_counters_fails[i][index] = -1;
				m_pool_fails[(i+pool_index)%m_hash_size] +=1;

			}
		}
			
	}
	for (; i < m_height; ++i) {

		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;


		++m_counters[i][index];
		if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
		{
			//m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
				--m_counters[i][index];
				m_counters_fails[i][index] = -1;
				m_pool_fails[(i+pool_index)%m_hash_size] +=1;
			}
		}
		
	}
}




// failed arrays turn into two counters
void CounterPools::increment3(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	
	int i;
	for (i = 0; i < m_height_plus_1_over_2; ++i) {

		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;




		++m_counters[i][index];
		if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
		{
			if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
				--m_counters[i][index];
				m_counters_fails[i][index] = -1;
				// if the index is odd use the lower high*width counter, if its even use the higher counters
				int high_low_choose = ((index%2)*m_height*m_width);
				m_pool_fails2[(i+m_height*pool_index)+high_low_choose] +=1;

			}
		}
			
	}
	for (; i < m_height; ++i) {

		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;


		++m_counters[i][index];
		if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
		{
			//m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
				--m_counters[i][index];
				m_counters_fails[i][index] = -1;
				// if the index is odd use the lower high*width counter, if its even use the higher counters
				int high_low_choose = ((index%2)*m_height*m_width);
				m_pool_fails2[(i+m_height*pool_index)+high_low_choose] +=1;
			}
		}
		
	}
}


void CounterPools::increment_cus(const char * str)
{

	

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	xxh::hash128_t hashes1 = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);

	
	int i;
	int j;

	
	
	uint32_t min = INT32_MAX; // if fails...

	for (i = 0; i < m_height_plus_1_over_2; ++i) {
		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);



		uint32_t temp = m_counters[i][index];
		if (min > temp)
		{
			min = temp;
		}
		
	}
	for (; i < m_height; ++i) {

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		uint32_t temp = m_counters[i][index];
		if (min > temp)
		{
			min = temp;
		}
		
	}




	hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	hashes1 = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);







	for (i = 0; i < m_height_plus_1_over_2; ++i) {

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;

		

		if(m_counters[i][index] == min){

			++m_counters[i][index];
			if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
			{
				if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
					--m_counters[i][index];
					m_counters_fails[i][index] = -1;
					// if the index is odd use the lower high*width counter, if its even use the higher counters
					int high_low_choose = ((index%2)*m_height*m_width);
					m_pool_fails2[(i+m_height*pool_index)+high_low_choose] +=1;

				}
			}
		}
			
	}
	for (; i < m_height; ++i) {
		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;

		if(m_counters[i][index] == min){
			++m_counters[i][index];
			if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
			{
				//m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
				if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
					--m_counters[i][index];
					m_counters_fails[i][index] = -1;
					// if the index is odd use the lower high*width counter, if its even use the higher counters
					int high_low_choose = ((index%2)*m_height*m_width);
					m_pool_fails2[(i+m_height*pool_index)+high_low_choose] +=1;
				}
			}
		}
		
	}
}


void CounterPools::increment_cuckoo(const char * str)
{
	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	increment_cuckoo_v(hashes, 1);
}


void CounterPools::increment_cuckoo_v(xxh::hash128_t hash, uint32_t val)
{
	int i;
	xxh::hash128_t hashes = hash;

	for (i = 0; i < m_height_plus_1_over_2; ++i) {

		uint64_t currHash = hashes.low64;

		hashes.low64 /= m_width;
		int index = currHash % m_width;
		int index2 = (currHash >> 32) % m_width;
		//hashes.low64 >>= m_log_width;
		//int index = currHash & m_width_mask;
		//int index2 = (currHash >> 32) & m_width_mask;

	

		int keyExist=0;
		if(m_keys[i][index] == hash){
			m_counters[i][index]+=val;

			if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
			{
				int pool_index = index / m_counters_per_pool;
				if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
					xxh::hash128_t tmpK = m_keys[i][index2];
					uint64_t tmpV = m_counters[i][index2];
					m_keys[i][index2] = m_keys[i][index];
					m_counters[i][index2] = m_counters[i][index];
					m_keys[i][index].low64 =0;
					m_keys[i][index].high64 =0;
					increment_cuckoo_v(tmpK, tmpV);
				}
			}
		}
		else if(m_keys[i][index2] == hash || m_keys[i][index2].low64 == 0 && m_keys[i][index2].high64==0){
			m_counters[i][index2]+=val;
			if ((m_counters[i][index2] & (m_counters[i][index2] - 1)) == 0)
			{
				int pool_index = index2 / m_counters_per_pool;
				if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
					xxh::hash128_t tmpK = m_keys[i][index];
					uint64_t tmpV = m_counters[i][index];
					m_keys[i][index] = m_keys[i][index2];
					m_counters[i][index] = m_counters[i][index2];
					m_keys[i][index2].low64 =0;
					m_keys[i][index2].high64 =0;
					increment_cuckoo_v(tmpK, tmpV);
				}
			}
		}
		else if(m_keys[i][index].low64 == 0 && m_keys[i][index].high64==0){
			m_counters[i][index]+=val;

			if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
			{
				int pool_index = index / m_counters_per_pool;
				if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
					xxh::hash128_t tmpK = m_keys[i][index2];
					uint64_t tmpV = m_counters[i][index2];
					m_keys[i][index2] = m_keys[i][index];
					m_counters[i][index2] = m_counters[i][index];
					m_keys[i][index].low64 =0;
					m_keys[i][index].high64 =0;
					increment_cuckoo_v(tmpK, tmpV);
				}
			}
		}
			
	}
	for (; i < m_height; ++i) {

		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;


		++m_counters[i][index];
		if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
		{
			//m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			if(mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool)){
				--m_counters[i][index];
				m_counters_fails[i][index] = -1;
				++failCount;
			}
		}
		
	}
}

void CounterPools::countersPerPrint(){

ofstream results_file;
string fn = "countersPer.txt";
results_file.open(fn, ofstream::out | ofstream::app);
results_file<<"width: " <<m_width<<endl;
int sum[33] = {0};

for(int i=0; i<m_height;i++){
	for(int j=0; j<m_width;j++){
		if(int(m_counters[i][j] !=0)){
			sum [ (int)log2(int(m_counters[i][j]))]++;
		}
		else{
			sum[0]++;
		}
	}
}

for(int k=0;k<=32;k++){
	results_file<< "k=" <<k << "->  "<<((float) sum[k])/((float)(m_height*m_width))<<endl;

}
	results_file<<" "<<endl;
	results_file<<" "<<endl;
}







uint32_t CounterPools::query(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);

	if (m_k > 0)
	{
		uint64_t fp = hashes.low64 ^ hashes.high64;

		auto it = m_top_k_map.find(fp);
		if (it != m_top_k_map.end()) {
			auto value = it.value();
			return get<1>(value) + get<2>(value);
		}
	}

	uint32_t min = INT32_MAX / 2; // if fails...
	int i;
	for (i = 0; i < m_height_plus_1_over_2; ++i) {
		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);


		if (!m_pool_failures[i][index / m_counters_per_pool])
		{
			uint32_t temp = m_counters[i][index];
			if (min > temp)
			{
				min = temp;
			}
		}
	}
	for (; i < m_height; ++i) {
		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);
		if (!m_pool_failures[i][index / m_counters_per_pool])
		{
			uint32_t temp = m_counters[i][index];
			if (min > temp)
			{
				min = temp;
			}
		}
	}
	return min;
}



uint32_t CounterPools::query1(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);

	uint32_t min = INT32_MAX / 2; // if fails...
	int i;
	for (i = 0; i < m_height_plus_1_over_2; ++i) {
		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);
		

		int pool_index = index / m_counters_per_pool;


		uint32_t temp = m_counters[i][index];


		if(m_counters_fails[i][index] == -1){
			temp += failCount;
		}

		if (min > temp)
		{
			min = temp;
		}
	}
	for (; i < m_height; ++i) {
		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;


		uint32_t temp = m_counters[i][index];


		if(m_counters_fails[i][index] == -1){
			temp += failCount;
		}

		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}





uint32_t CounterPools::query2(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);

	uint32_t min = INT32_MAX / 2; // if fails...
	int i;
	for (i = 0; i < m_height_plus_1_over_2; ++i) {
		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;


		uint32_t temp = m_counters[i][index];


		if(m_counters_fails[i][index] == -1){
			temp += m_pool_fails[(i+pool_index)%m_hash_size];
		}

		if (min > temp)
		{
			min = temp;
		}
	}
	for (; i < m_height; ++i) {
		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;


		uint32_t temp = m_counters[i][index];


		if(m_counters_fails[i][index] == -1){
			temp += m_pool_fails[(i+pool_index)%m_hash_size];
		}

		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}





uint32_t CounterPools::query3(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);

	uint32_t min = INT32_MAX / 2; // if fails...
	int i;
	for (i = 0; i < m_height_plus_1_over_2; ++i) {
		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;


		uint32_t temp = m_counters[i][index];


		if(m_counters_fails[i][index] == -1){
			// if the index is odd use the lower high*width counter, if its even use the higher counters
			int high_low_choose = ((index%2)*m_height*m_width);
			temp += m_pool_fails2[(i+m_height*pool_index)+high_low_choose];
		}

		if (min > temp)
		{
			min = temp;
		}
	}
	for (; i < m_height; ++i) {
		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;


		uint32_t temp = m_counters[i][index];


		if(m_counters_fails[i][index] == -1){
			// if the index is odd use the lower high*width counter, if its even use the higher counters
			int high_low_choose = ((index%2)*m_height*m_width);
			temp += m_pool_fails2[(i+m_height*pool_index)+high_low_choose];
		}

		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}

float CounterPools::percent_of_crushed_pools(){
	
	return crushesNum/(m_height*(m_width / m_counters_per_pool + 1));
}

void CounterPools::write_to_file(string fn)
{
	ofstream results_file;
	results_file.open(fn, ofstream::out);
	for (int i = 0; i < m_height; ++i) 
	{
		for (int j = 0; j < m_width; ++j)
		{
			int pool_index = j / m_counters_per_pool;
			if (!m_pool_failures[i][pool_index])
			{
				results_file << m_counters[i][j];
			}
			else
			{
				results_file << INT32_MAX;
			}
			if (j < m_width - 1)
			{
				results_file << ",";
			}
		}
		if (i < m_height - 1)
		{
			results_file << endl;
		}		
	}
}

void CounterPools::read_from_file(string fn)
{
	std::ifstream data(fn);
	std::string line;

	int i = 0;
	while (getline(data, line))
	{

		stringstream lineStream(line);
		string cell;

		uint32_t value;

		int j = 0;
		while (getline(lineStream, cell, ','))
		{
			int pool_index = j / m_counters_per_pool;

			stringstream scell(cell);
			scell >> value;
			
			if (value != INT32_MAX)
			{
				m_counters[i][j] = value;
				m_pool_failures[i][pool_index] = false;
			}
			else
			{
				m_counters[i][j] = INT32_MAX / 2;
				m_pool_failures[i][pool_index] = true;
			}
			++j;
		}
		++i;
	}
}

bool CounterPools::mf_check_pool_failure(uint32_t* pool_pointer)
{
	int remaining_bits = m_pool_bit_size - m_counters_per_pool * m_initial_counter_size;
	for (int offset = 0; offset < m_counters_per_pool; ++offset)
	{
		int required_bits = ceil(log2(pool_pointer[offset] + 1));
		int additional_required_bits = required_bits - m_initial_counter_size;
		if (additional_required_bits > 0)
		{
			int allocated_additional_bits = ceil(((double)additional_required_bits) / m_counter_bit_increase)*m_counter_bit_increase;
			remaining_bits -= allocated_additional_bits;
		}	
	}
	return remaining_bits < 0;
}

void CounterPools::mf_add_to_sketch(const char * str, int weight)
{
	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	int i;

	for (i = 0; i < m_height_plus_1_over_2; ++i) {

		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);


		int pool_index = index / m_counters_per_pool;
		if (!m_pool_failures[i][pool_index])
		{
			uint32_t leading_zeros_1 = _lzcnt_u32(m_counters[i][index]);
			m_counters[i][index] += weight;
			uint32_t leading_zeros_2 = _lzcnt_u32(m_counters[i][index]);
			if (leading_zeros_1 != leading_zeros_2)
			{
				m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			}
		}
	}
	for (; i < m_height; ++i) {

		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;
		if (!m_pool_failures[i][pool_index])
		{
			uint32_t leading_zeros_1 = _lzcnt_u32(m_counters[i][index]);
			m_counters[i][index] += weight;
			uint32_t leading_zeros_2 = _lzcnt_u32(m_counters[i][index]);
			if (leading_zeros_1 != leading_zeros_2)
			{
				m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			}
		}
	}
}

uint32_t CounterPools::read_counter(int i, int j)
{
	return m_counters[i][j];
}

bool CounterPools::read_pool_failure(int i, int j)
{
	return m_pool_failures[i][j / m_counters_per_pool];
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

CounterPoolsExternalTopK::CounterPoolsExternalTopK(int width, int height, int seed, int pool_bit_size, int counters_per_pool, int initial_counter_size, int counter_bit_increase, int k) :
	topK(k)
{
	m_width = width;
	m_height = height;

	m_width_mask = width - 1;

	int index = m_width_mask;
	m_log_width = 0;
	while (index >>= 1) ++m_log_width;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	m_height_plus_1_over_2 = (height + 1) / 2;
	if (m_height_plus_1_over_2 > 1)
	{
		assert(width < ((uint64_t)1 << (64 / m_height_plus_1_over_2)) && "128 hash bits are not enough!");
	}

	m_counters = new uint32_t*[height];
	m_pool_failures = new bool*[height];
	for (int i = 0; i < height; ++i)
	{
		m_counters[i] = new uint32_t[width]();
		m_pool_failures[i] = new bool[width / counters_per_pool + 1]();
	}

	m_xxhash_seed = seed * 7 + 100;

	m_pool_bit_size = pool_bit_size;
	m_counters_per_pool = counters_per_pool;
	m_initial_counter_size = initial_counter_size;
	m_counter_bit_increase = counter_bit_increase;

	m_augmented = k > 0;
}

CounterPoolsExternalTopK::~CounterPoolsExternalTopK()
{
	for (int i = 0; i < m_height; ++i)
	{
		delete[] m_counters[i];
		delete[] m_pool_failures[i];
	}
	delete[] m_counters;
	delete[] m_pool_failures;
}

void CounterPoolsExternalTopK::increment(const char * str)
{
	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	int i;

	if (m_augmented)
	{
		uint64_t fp = hashes.low64 ^ hashes.high64;

		if (topK.increment_if_exists(fp))
		{
			return;
		}

		uint32_t estimated_size = query(str);
		if (estimated_size >= topK.get_min())
		{
			auto value = topK.add(fp, estimated_size, str);
			if (get<2>(value))
			{
				mf_add_to_sketch(get<0>(value).c_str(), get<2>(value));
			}
			return;
		}			
	}

	for (i = 0; i < m_height_plus_1_over_2; ++i) {

		int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;


		int index1 = hashes.low64 % m_width;
		if((index != index1) || ((hashes.low64 >>m_log_width) != (hashes.low64 / m_width))){
			printf("@@@@@@@@@@@@@@@@@");
		}
		hashes.low64 /= m_width;


		int pool_index = index / m_counters_per_pool;
		if (!m_pool_failures[i][pool_index])
		{
			++m_counters[i][index];
			if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
			{
				m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			}
		}
	}
	for (; i < m_height; ++i) {

		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;
		if (!m_pool_failures[i][pool_index])
		{
			++m_counters[i][index];
			if ((m_counters[i][index] & (m_counters[i][index] - 1)) == 0)
			{
				m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			}
		}
	}
}

uint32_t CounterPoolsExternalTopK::query(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);

	if (m_augmented)
	{
		uint64_t fp = hashes.low64 ^ hashes.high64;
		int top_k_val = topK.query(fp);
		if (top_k_val) {
			return top_k_val;
		}
	}

	uint32_t min = INT32_MAX / 2; // if fails...
	int i;
	for (i = 0; i < m_height_plus_1_over_2; ++i) {
		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);

		if (!m_pool_failures[i][index / m_counters_per_pool])
		{
			uint32_t temp = m_counters[i][index];
			if (min > temp)
			{
				min = temp;
			}
		}
	}
	for (; i < m_height; ++i) {
		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;


		int index = hashes.high64 % m_width;
		hashes.high64 /= m_width;

		if (!m_pool_failures[i][index / m_counters_per_pool])
		{
			uint32_t temp = m_counters[i][index];
			if (min > temp)
			{
				min = temp;
			}
		}
	}
	return min;
}

void CounterPoolsExternalTopK::mf_add_to_sketch(const char * str, int weight)
{
	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_xxhash_seed);
	int i;

	for (i = 0; i < m_height_plus_1_over_2; ++i) {

		//int index = hashes.low64 & m_width_mask;
		//hashes.low64 >>= m_log_width;

		int index = hashes.low64 % m_width;
		hashes.low64 /= (m_width/2);


		int pool_index = index / m_counters_per_pool;
		if (!m_pool_failures[i][pool_index])
		{
			uint32_t leading_zeros_1 = _lzcnt_u32(m_counters[i][index]);
			m_counters[i][index] += weight;
			uint32_t leading_zeros_2 = _lzcnt_u32(m_counters[i][index]);
			if (leading_zeros_1 != leading_zeros_2)
			{
				m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			}
		}
	}
	for (; i < m_height; ++i) {

		//int index = hashes.high64 & m_width_mask;
		//hashes.high64 >>= m_log_width;

		int index = hashes.high64 % m_width;
		hashes.high64 /= (m_width/2);

		int pool_index = index / m_counters_per_pool;
		if (!m_pool_failures[i][pool_index])
		{
			uint32_t leading_zeros_1 = _lzcnt_u32(m_counters[i][index]);
			m_counters[i][index] += weight;
			uint32_t leading_zeros_2 = _lzcnt_u32(m_counters[i][index]);
			if (leading_zeros_1 != leading_zeros_2)
			{
				m_pool_failures[i][pool_index] = mf_check_pool_failure(m_counters[i] + pool_index * m_counters_per_pool);
			}
		}
	}
}

bool CounterPoolsExternalTopK::mf_check_pool_failure(uint32_t* pool_pointer)
{
	int remaining_bits = m_pool_bit_size - m_counters_per_pool * m_initial_counter_size;
	for (int offset = 0; offset < m_counters_per_pool; ++offset)
	{
		int required_bits = ceil(log2(pool_pointer[offset] + 1));
		int additional_required_bits = required_bits - m_initial_counter_size;
		if (additional_required_bits > 0)
		{
			int allocated_additional_bits = ceil(((double)additional_required_bits) / m_counter_bit_increase)*m_counter_bit_increase;
			remaining_bits -= allocated_additional_bits;
		}
	}
	return remaining_bits < 0;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////