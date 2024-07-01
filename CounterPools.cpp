#include <algorithm>
#include "CounterPools.hpp"

bool CounterPools_64_4_0_1::m_initLookupTable = false;
uint32_t CounterPools_64_4_0_1::m_lookup[STARS_AND_BARS_64_4] = { 0 };


// second transition is to choose the smallest counter thats larger than us


CounterPools_64_4_0_1::CounterPools_64_4_0_1(int width, int height, int seed)
{
	initLookupTableFunc();

	/*
	for (uint16_t encoding = 0; encoding < STARS_AND_BARS_64_4; ++encoding)
	{
		if ((encoding != encode(m_lookup[encoding])) || (encoding != encode_without_auxilary_map(m_lookup[encoding])))
		{
			std::cout << encoding << "\t" << encode(m_lookup[encoding]) << std::endl;
		}
	}
	*/
	flag=0;
	m_width = width;
	m_height = height;

	m_width_mask = width - 1;

	int index = m_width_mask;
	m_log_width = 0;
	while (index >>= 1) ++m_log_width;

	assert(width >= 4 && "We assume at least 4 counters!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	m_height_plus_1_over_2 = (height + 1) / 2;
	m_seed = seed * (7 + 0) + 0 + 100;
	if (height > 1)
	{
		assert(width < ((uint64_t)1 << (64 / m_height_plus_1_over_2)) && "128 hash bits are not enough!");
	}

	m_pools = new uint64_t*[height];
	m_encodings = new uint16_t*[height];

	for (int i = 0; i < height; ++i)
	{
		m_pools[i] = new uint64_t[width >> 2]();
		m_encodings[i] = new uint16_t[width >> 2]();

		// changing starting number of bits to 14
		/*
		for (int j = 0; j < width >> 2; ++j)
		{
			m_encodings[i][j] = 15832;
		}
		*/
	}
}

CounterPools_64_4_0_1::~CounterPools_64_4_0_1()
{
	for (int i = 0; i < m_height; ++i)
	{
		delete[] m_pools[i];
		delete[] m_encodings[i];
	}
	delete[] m_pools;
	delete[] m_encodings;
}

void CounterPools_64_4_0_1::initLookupTableFunc()
{
	if (m_initLookupTable)
	{
		return;
	}
	m_initLookupTable = true;

	/*
	X = sum(seq) + 2
	Y = len(seq)
	SABsumDict = dict([((n, x, y), sum([nSAB(n - _, y) for _ in range(x)])) for n in range(X) for x in range(X) for y in range(Y)])
	*/

	uint64_t X = 66; 
	uint64_t Y = 4;

	for (uint64_t n = 0; n < X; ++n)
	{
		for (uint64_t x = 0; x < X; ++x)
		{
			for (uint64_t y = 1; y < Y; ++y)
			{
				uint64_t res = 0;
				for (uint64_t u = 0; u < x; ++u)
				{
					res += starsAndBars(n - u, y);
				}
				uint64_t key = ((n << 32) | (x << 16) | y);
				m_auxilary_map[key] = res;
			}
		}
	}

	for (uint16_t encoding = 0; encoding < STARS_AND_BARS_64_4; ++encoding)
	{
		m_lookup[encoding] = decode_offsets(encoding, 64, 4);
	}

	// comment out for testing at constructor
	m_auxilary_map.clear();
}

uint32_t CounterPools_64_4_0_1::decode_sizes(uint16_t encoding, uint64_t n, uint64_t k)
{
	/*
	def SABdecode(conf,n,k):
		if k == 1:
			return [n]
		f = 0
		while SABsumDict[(n,f+1,k-1)] <= conf:
			f += 1
		return [f] + SABdecode(conf-SABsumDict[(n,f,k-1)],n-f,k-1)
	*/

	if (k == 1)
	{
		return n;
	}

	uint64_t f = 0;
	while (m_auxilary_map[CONCAT_KEYS(n, f + 1, k - 1)] <= encoding)
	{
		++f;
	}

	return (f << (8 * (k - 1))) | decode_sizes(encoding - m_auxilary_map[CONCAT_KEYS(n, f, k - 1)], n - f, k - 1);
}

uint32_t CounterPools_64_4_0_1::decode_offsets(uint16_t encoding, uint64_t n, uint64_t k)
{
	uint32_t sizes = decode_sizes(encoding, n, k);
	uint8_t* sizes_array = (uint8_t*)&sizes;
	uint32_t offsets = 0;
	uint32_t offset = sizes_array[k - 1];

	for (int i = 1; i < k - 1; ++i)
	{
		offsets |= offset;
		offsets <<= 8;
		offset += sizes_array[k - 1 - i];
	}
	offsets |= offset;

	uint8_t* offsets_array = (uint8_t*)&offsets;
	std::reverse(offsets_array, offsets_array + k);
	return offsets;
}

uint16_t CounterPools_64_4_0_1::encode(uint32_t sizes)
{
	/*
	def SABencode(sizes):
		k = len(sizes)
		n = sum(sizes)
		f = sizes[0]
		if k == 1:
			return 0
		prevCombs = SABsumDict.get((n,f,k-1),0) #0 if f == 0 else SABsumDict[(n,f,k-1)]
		return prevCombs+SABencode(sizes[1:])
	*/

	uint64_t n = 64;
	uint64_t k = 4;
	uint16_t encoding = 0;
	for (int i = 0; i < k - 1; ++i)
	{
		uint64_t size = ((sizes & (0xFF << ((k - i - 1) << 3))) >> ((k - i - 1) << 3));
		encoding += m_auxilary_map[CONCAT_KEYS(n, size, k - i - 1)];
		n -= size;
	}
	return encoding;
}

uint16_t CounterPools_64_4_0_1::encode_without_auxilary_map(uint32_t sizes)
{
	uint64_t n = 64;
	uint64_t k = 4;

	uint16_t encoding = 0;
	for (int i = 0; i < k - 1; ++i)
	{
		uint64_t size = ((sizes & (0xFF << ((k - i - 1) << 3))) >> ((k - i - 1) << 3));

		uint64_t res = 0;
		for (uint64_t u = 0; u < size; ++u)
		{
			res += starsAndBars(n - u, k - i - 1);
		}
		encoding += res;

		n -= size;
	}

	return encoding;
}

void CounterPools_64_4_0_1::increment(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_seed);

	for (int row = 0; row < m_height; ++row) {

		uint32_t index;
		if (row < m_height_plus_1_over_2)
		{
			index = hashes.low64 & m_width_mask;
			hashes.low64 >>= m_log_width;
		}
		else
		{
			index = hashes.high64 & m_width_mask;
			hashes.high64 >>= m_log_width;
		}
		
		uint32_t pool_index = index >> 2;

		uint16_t &pool_encoding = m_encodings[row][pool_index];
		uint64_t &pool = m_pools[row][pool_index];

		uint8_t* offsets = (uint8_t*)(m_lookup + pool_encoding);
		uint64_t k = 4;

		int poolCounterIndex = index & 0b11;
		int counterBitOffset = offsets[poolCounterIndex];

		
		int counterBitSize;
		uint64_t mask;
		uint64_t value = pool >> counterBitOffset;

		if (poolCounterIndex < 3)
		{
			counterBitSize = offsets[poolCounterIndex + 1] - counterBitOffset;
			mask = (((uint64_t)1 << (counterBitSize)) - 1);
			value &= mask;

			if (value < mask)
			{
				pool += ((uint64_t)1 << counterBitOffset);
				continue;
			}
		}
		else
		{
			counterBitSize = 64 - counterBitOffset;

			if (value + 1 < ((uint64_t)1 << counterBitSize))
			{
				pool += ((uint64_t)1 << counterBitOffset);
				continue;
			}
			else if (counterBitSize == 64)
			{
				++pool;
				continue;
			}
		}
				
		if (pool_encoding != STARS_AND_BARS_64_4_POOL_FAILED)
		{
			uint8_t sizes_array[4] = { 0 };
			sizes_array[0] = offsets[1];
			for (int i = 1; i < k - 1; ++i)
			{
				sizes_array[i] = offsets[i + 1] - offsets[i];
			}
			sizes_array[k - 1] = 64 - offsets[k - 1];

			if (pool & (((uint64_t)1) << 63))
			{
				pool_encoding = STARS_AND_BARS_64_4_POOL_FAILED;
				continue;
			}

			++sizes_array[poolCounterIndex];
			--sizes_array[k - 1];

			uint32_t new_sizes = 0;
			for (int i = 0; i < k - 1; ++i)
			{
				new_sizes |= sizes_array[i];
				new_sizes <<= 8;
			}
			new_sizes |= sizes_array[k - 1];
			pool_encoding = encode_without_auxilary_map(new_sizes);

			uint64_t pool_lsb = pool & (((uint64_t)1 << (counterBitOffset + counterBitSize)) - 1);
			uint64_t pool_msb = pool & (~(((uint64_t)1 << (counterBitOffset + counterBitSize)) - 1));
			pool = pool_lsb | (pool_msb << 1);

			++counterBitSize;
			if (value + 1 < ((uint64_t)1 << counterBitSize))
			{
				pool += ((uint64_t)1 << counterBitOffset);
			}
			else
			{
				++pool;
			}
		}
		
	}
}


















void CounterPools_64_4_0_1::incrementCus(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_seed);
	uint64_t minVal = std::numeric_limits<uint64_t>::max();

	uint32_t minValIndexes[m_height] = {0};
	int minValrows[m_height] = {0};
	int arrayIndex=-1;



	for (int row = 0; row < m_height; ++row) {
		uint32_t index;
		if (row < m_height_plus_1_over_2)
		{
			index = hashes.low64 & m_width_mask;
			hashes.low64 >>= m_log_width;
		}
		else
		{
			index = hashes.high64 & m_width_mask;
			hashes.high64 >>= m_log_width;
		}
		
		uint32_t pool_index = index >> 2;

		uint16_t &pool_encoding = m_encodings[row][pool_index];
		uint64_t &pool = m_pools[row][pool_index];

		uint8_t* offsets = (uint8_t*)(m_lookup + pool_encoding);

		int poolCounterIndex = index & 0b11;
		int counterBitOffset = offsets[poolCounterIndex];

		

		
		uint64_t value = pool >> counterBitOffset;
		uint64_t mask;
		int counterBitSize;

		if (poolCounterIndex < 3)
		{
			counterBitSize = offsets[poolCounterIndex + 1] - counterBitOffset;
			mask = (((uint64_t)1 << (counterBitSize)) - 1);
			value &= mask;
		}


		if(value <= minVal){
			minVal = value;
			arrayIndex=0;
			minValIndexes[arrayIndex] = index;
			minValrows[arrayIndex] = row;
		}
		else if(value == minVal){
			arrayIndex++;
			minValIndexes[arrayIndex] = index;
			minValrows[arrayIndex] = row;

		}
	}








	for (int i = 0; i <= arrayIndex; ++i) {

		uint32_t index = minValIndexes[i];
		int row = minValrows[i];

		uint32_t pool_index = index >> 2;		
		uint16_t &pool_encoding = m_encodings[row][pool_index];
		uint64_t &pool = m_pools[row][pool_index];
		uint8_t* offsets = (uint8_t*)(m_lookup + pool_encoding);
		uint64_t k = 4;
		int poolCounterIndex = index & 0b11;
		int counterBitOffset = offsets[poolCounterIndex];

		

		
		uint64_t value = pool >> counterBitOffset;


		int counterBitSize;
		uint64_t mask;
		if (poolCounterIndex < 3)
			{
				counterBitSize = offsets[poolCounterIndex + 1] - counterBitOffset;
				mask = (((uint64_t)1 << (counterBitSize)) - 1);
				value &= mask;

				if (value < mask)
				{
					pool += ((uint64_t)1 << counterBitOffset);
					continue;
				}
			}
			else
			{
				counterBitSize = 64 - counterBitOffset;

				if (value + 1 < ((uint64_t)1 << counterBitSize))
				{
					pool += ((uint64_t)1 << counterBitOffset);
					continue;
				}
				else if (counterBitSize == 64)
				{
					++pool;
					continue;
				}
			}

				
		if (pool_encoding != STARS_AND_BARS_64_4_POOL_FAILED)
		{
			uint8_t sizes_array[4] = { 0 };
			sizes_array[0] = offsets[1];
			for (int i = 1; i < k - 1; ++i)
			{
				sizes_array[i] = offsets[i + 1] - offsets[i];
			}
			sizes_array[k - 1] = 64 - offsets[k - 1];

			if (pool & (((uint64_t)1) << 63))
			{
				pool_encoding = STARS_AND_BARS_64_4_POOL_FAILED;
				continue;
			}

			++sizes_array[poolCounterIndex];
			--sizes_array[k - 1];

			uint32_t new_sizes = 0;
			for (int i = 0; i < k - 1; ++i)
			{
				new_sizes |= sizes_array[i];
				new_sizes <<= 8;
			}
			new_sizes |= sizes_array[k - 1];
			pool_encoding = encode_without_auxilary_map(new_sizes);

			uint64_t pool_lsb = pool & (((uint64_t)1 << (counterBitOffset + counterBitSize)) - 1);
			uint64_t pool_msb = pool & (~(((uint64_t)1 << (counterBitOffset + counterBitSize)) - 1));
			pool = pool_lsb | (pool_msb << 1);

			++counterBitSize;
			if (value + 1 < ((uint64_t)1 << counterBitSize))
			{
				pool += ((uint64_t)1 << counterBitOffset);
			}
			else
			{
				++pool;
			}
		}
	}
		
}




uint64_t CounterPools_64_4_0_1::chooseDeleteCounter(uint64_t pool, uint8_t* offsets, int index ){
	if(index<3){
		return deleteCounter(pool,offsets,index+1);
	}
	else{
		return deleteCounter(pool,offsets,0);
	}
}



uint64_t CounterPools_64_4_0_1::deleteCounter(uint64_t &pool, uint8_t* offsets, int index ){
	if(index <3){
		int bitsNum = offsets[index];
		int fromBit=1;

		// extract "bitsNum"  bits from "fromBit" position
		uint64_t first = (((1 << bitsNum) - 1) & (pool >> fromBit-1));

		
		int fromBit1 = (offsets[index+1]) + 1;
		int bitsNum1 = 64 - fromBit1;
		

		// extract "bitsNum1"  bits from "fromBit1" position
		uint64_t last = (((1 << bitsNum1) - 1) & (pool >> (fromBit1 - 1)));


		uint64_t diff = offsets[index+1] - offsets[index];


		index++;
		while(index<4){
			offsets[index] = offsets[index] - diff;
			index++;
		}
		pool = first|(last << bitsNum);
	}
	else{
		int p = 1;
		int k = offsets[index];
		pool = (((1 << k) - 1) & (pool >> p-1));

	}
	return 1;
}



// void CounterPools_64_4_0_1::increment2(const char * str)
// {

// 	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_seed);

// 	for (int row = 0; row < m_height; ++row) {

// 		uint32_t index;
// 		if (row < m_height_plus_1_over_2)
// 		{
// 			index = hashes.low64 & m_width_mask;
// 			hashes.low64 >>= m_log_width;
// 		}
// 		else
// 		{
// 			index = hashes.high64 & m_width_mask;
// 			hashes.high64 >>= m_log_width;
// 		}
		
// 		uint32_t pool_index = index >> 2;

// 		uint16_t &pool_encoding = m_encodings[row][pool_index];
// 		uint64_t &pool = m_pools[row][pool_index];
// 		char* key =  keys[row][pool_index];

// 		uint8_t* offsets = (uint8_t*)(m_lookup + pool_encoding);
// 		uint64_t k = 4;

// 		int poolCounterIndex = index & 0b11;
// 		int counterBitOffset = offsets[poolCounterIndex];
// 		int counterBitSize = offsets[poolCounterIndex+1] - offsets[poolCounterIndex];

// 		if ( counterBitSize!=0 && strcmp(key, str)!=0 ){

// 			xxh::hash128_t hashes2 = xxh::xxhash3<128>(str, FT_SIZE, m_seed+1);

// 			uint32_t index2;
// 			if (row < m_height_plus_1_over_2)
// 			{
// 				index2 = hashes2.low64 & m_width_mask;
// 				hashes2.low64 >>= m_log_width;
// 			}
// 			else
// 			{
// 				index2 = hashes2.high64 & m_width_mask;
// 				hashes2.high64 >>= m_log_width;
// 			}
			
// 			uint32_t pool_index2 = index2 >> 2;

// 			uint16_t &pool_encoding2 = m_encodings[row][pool_index2];
// 			uint64_t &pool2 = m_pools[row][pool_index2];
// 			char* key2 =  keys[row][pool_index2];

// 			uint8_t* offsets2 = (uint8_t*)(m_lookup + pool_encoding2);
// 			uint64_t k2 = 4;

// 			int poolCounterIndex2 = index2 & 0b11;
// 			int counterBitOffset2 = offsets2[poolCounterIndex2];
// 			int counterBitSize2 = offsets2[poolCounterIndex2+1] - offsets[poolCounterIndex];

// 			if ( counterBitSize!=0 && strcmp(key, str)!=0 ){
// 				pool2 = deleteCounter(pool2, offsets2, index2 );
// 				uint64_t value = pool2 >> counterBitOffset2;
// 				if (poolCounterIndex < 3)
// 				{
// 					uint64_t mask = (((uint64_t)1 << (counterBitSize2)) - 1);
// 					value &= mask;
// 				}
// 				reInsert(key2, value, row);
// 			}
// 		}

			
// 		int counterBitSize1;
// 		uint64_t mask;
// 		uint64_t value = pool >> counterBitOffset;

// 		if (poolCounterIndex < 3)
// 		{
// 			counterBitSize1 = offsets[poolCounterIndex + 1] - counterBitOffset;
// 			mask = (((uint64_t)1 << (counterBitSize1)) - 1);
// 			value &= mask;
// 		if (value < mask)
// 			{	
// 				pool += ((uint64_t)1 << counterBitOffset);
// 				continue;
// 			}
// 		}
// 		else
// 		{
// 			counterBitSize1 = 64 - counterBitOffset;

// 			if (value + 1 < ((uint64_t)1 << counterBitSize1))
// 			{
// 				pool += ((uint64_t)1 << counterBitOffset);
// 				continue;
// 			}
// 			else if (counterBitSize1 == 64)
// 			{
// 				++pool;
// 				continue;
// 			}
// 		}
				

// 		uint8_t sizes_array[4] = { 0 };
// 		sizes_array[0] = offsets[1];
// 		for (int i = 1; i < k - 1; ++i)
// 		{
// 			sizes_array[i] = offsets[i + 1] - offsets[i];
// 		}
// 		sizes_array[k - 1] = 64 - offsets[k - 1];

// 		if (pool & (((uint64_t)1) << 63))
// 		{
// 			//relocateCounter(row, index);
// 			increment2(str);
// 			continue;
// 		}

// 		++sizes_array[poolCounterIndex];
// 		--sizes_array[k - 1];

// 		uint32_t new_sizes = 0;
// 		for (int i = 0; i < k - 1; ++i)
// 		{
// 			new_sizes |= sizes_array[i];
// 			new_sizes <<= 8;
// 		}
// 		new_sizes |= sizes_array[k - 1];
// 		pool_encoding = encode_without_auxilary_map(new_sizes);

// 		uint64_t pool_lsb = pool & (((uint64_t)1 << (counterBitOffset + counterBitSize1)) - 1);
// 		uint64_t pool_msb = pool & (~(((uint64_t)1 << (counterBitOffset + counterBitSize1)) - 1));
// 		pool = pool_lsb | (pool_msb << 1);

// 		++counterBitSize1;
// 		if (value + 1 < ((uint64_t)1 << counterBitSize1))
// 		{
// 			pool += ((uint64_t)1 << counterBitOffset);
// 		}
// 		else
// 		{
// 			++pool;
// 		}
	
		
// 	}
// }

/*
void CounterPools_64_4_0_1::relocateCounter(int row, uint32_t index){
	int counterNum = ChooseCounterToRelocate(row, index);

}


int CounterPools_64_4_0_1::ChooseCounterToRelocate(int row, uint32_t index){
	int poolCounterIndex = index & 0b11;
	if(poolCounterIndex < 3){
		return poolCounterIndex+1;
	}
	else{
		return 0;
	}
}
*/

void CounterPools_64_4_0_1::increment_index(int row, uint32_t index)
{
	uint32_t pool_index = index >> 2;
	uint16_t &pool_encoding = m_encodings[row][pool_index];

	if (pool_encoding != STARS_AND_BARS_64_4_POOL_FAILED)
	{
		uint64_t &pool = m_pools[row][pool_index];
		//uint32_t sizes = m_lookup[pool_encoding];
		//uint8_t* sizes = (uint8_t*)(m_lookup + pool_encoding);
		uint8_t* offsets = (uint8_t*)(m_lookup + pool_encoding);
		uint64_t k = 4;

		int poolCounterIndex = index & 0b11;
		//int counterBitOffset = 0;
		int counterBitOffset = offsets[poolCounterIndex];
		//for (int i = 3; i > 3 - poolCounterIndex; --i)
		//{
			//uint64_t size = ((sizes & (0xFF << ((k - i - 1) << 3))) >> ((k - i - 1) << 3));
		//	counterBitOffset += sizes[i];
		//}
		//int counterBitSize = ((sizes & (0xFF << ((k - poolCounterIndex - 1) << 3))) >> ((k - poolCounterIndex - 1) << 3));
		//int counterBitSize = sizes[3 - poolCounterIndex];
		int counterBitSize = poolCounterIndex < 3 ? offsets[poolCounterIndex + 1] - counterBitOffset : 64 - counterBitOffset;

		if (increment_h(pool, counterBitSize, counterBitOffset))
		{
			uint8_t sizes_array[4] = { 0 };
			sizes_array[0] = offsets[1];
			for (int i = 1; i < k - 1; ++i)
			{
				//sizes_array[i] = ((sizes & (0xFF << ((k - i - 1) << 3))) >> ((k - i - 1) << 3));
				sizes_array[i] = offsets[i + 1] - offsets[i];
			}
			sizes_array[k - 1] = 64 - offsets[k - 1];

			// for signed counters - the sign bit is the LSB.
			if (pool & (((uint64_t)1) << 63))
			{
				pool_encoding = STARS_AND_BARS_64_4_POOL_FAILED;
			}

			++sizes_array[poolCounterIndex];
			--sizes_array[k - 1];

			uint32_t new_sizes = 0;
			for (int i = 0; i < k - 1; ++i)
			{
				new_sizes |= sizes_array[i];
				new_sizes <<= 8;
			}
			new_sizes |= sizes_array[k - 1];
			pool_encoding = encode_without_auxilary_map(new_sizes);

			uint64_t pool_lsb = pool & (((uint64_t)1 << (counterBitOffset + counterBitSize)) - 1);
			uint64_t pool_msb = pool & (~(((uint64_t)1 << (counterBitOffset + counterBitSize)) - 1));
			pool = pool_lsb | (pool_msb << 1);

			increment_h(pool, counterBitSize + 1, counterBitOffset);
		}
	}
}

inline uint64_t CounterPools_64_4_0_1::read_counter(uint64_t & pool, int counterBitSize, int counterBitOffset)
{
	if (counterBitSize + counterBitOffset == 64)
	{
		return pool >> counterBitOffset;
	}
	return (pool & (((uint64_t)1 << (counterBitOffset + counterBitSize)) - 1)) >> counterBitOffset;
}

bool CounterPools_64_4_0_1::increment_h(uint64_t & pool, int counterBitSize, int counterBitOffset)
{
	
	uint64_t mask = (counterBitSize + counterBitOffset != 64) ? (((uint64_t)1 << (counterBitOffset + counterBitSize)) - 1) : -1;
	uint64_t value = (pool & mask) >> counterBitOffset;

	if (value + 1 < ((uint64_t)1 << counterBitSize))
	{
		pool += ((uint64_t)1 << counterBitOffset);
		return false;
	}
	else if (counterBitSize == 64)
	{
		++pool;
		return false;
	}
	return true;
}

uint64_t CounterPools_64_4_0_1::query(const char * str)
{

	xxh::hash128_t hashes = xxh::xxhash3<128>(str, FT_SIZE, m_seed);
	uint64_t min = ~(uint64_t)0;
	for (int row = 0; row < m_height; ++row)
	{

		uint32_t index;
		if (row < m_height_plus_1_over_2)
		{
			index = hashes.low64 & m_width_mask;
			hashes.low64 >>= m_log_width;
		}
		else
		{
			index = hashes.high64 & m_width_mask;
			hashes.high64 >>= m_log_width;
		}

		uint32_t pool_index = index >> 2;
		uint16_t &pool_encoding = m_encodings[row][pool_index];
		
		if (pool_encoding == STARS_AND_BARS_64_4_POOL_FAILED)
		{
			continue;
		}

		uint64_t &pool = m_pools[row][pool_index];
		uint8_t* offsets = (uint8_t*)(m_lookup + pool_encoding);

		uint64_t k = 4;

		int poolCounterIndex = index & 0b11;
		int counterBitOffset = offsets[poolCounterIndex];
		int counterBitSize = poolCounterIndex < 3 ? offsets[poolCounterIndex + 1] - counterBitOffset : 64 - counterBitOffset;
		uint64_t row_counter = read_counter(pool, counterBitSize, counterBitOffset);

		min = (min <= row_counter) ? min : row_counter;
	}

	return min;
}

