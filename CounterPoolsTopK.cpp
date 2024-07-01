#include "CounterPoolsTopK.hpp"

CounterPoolsTopK::CounterPoolsTopK(int k) :
	m_top_k_map(4 * k)
{
	m_top_k_count = 0;
	m_top_k_min = 0;
	m_k = k;
}

CounterPoolsTopK::~CounterPoolsTopK()
{
}

bool CounterPoolsTopK::increment_if_exists(uint64_t fp)
{
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
		return true;
	}
	return false;
}

tuple<string, uint32_t, uint32_t> CounterPoolsTopK::add(uint64_t fp, uint32_t estimated_size, const char * str)
{
	tuple<string, uint32_t, uint32_t> new_value(string(str, FT_SIZE), estimated_size, 1);
	if ((m_top_k_count == m_k) && (m_top_k_count))
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

		m_top_k_min = INT32_MAX;
		for (auto it = m_top_k_map.begin(); it != m_top_k_map.end(); ++it)
		{
			auto it_value = it.value();
			if (get<1>(it_value) + get<2>(it_value) < m_top_k_min)
			{
				m_top_k_min = get<1>(it_value) + get<2>(it_value);
			}
		}

		return value;
	}
	else
	{
		m_top_k_map[fp] = new_value;
		++m_top_k_count;
	}
	return tuple<string, uint32_t, uint32_t>("", 0, 0);
}

uint32_t CounterPoolsTopK::query(uint64_t fp)
{
	auto it = m_top_k_map.find(fp);
	if (it != m_top_k_map.end()) {
		auto value = it.value();
		return get<1>(value) + get<2>(value);
	}
	return 0;
}

uint32_t CounterPoolsTopK::get_min()
{
	return m_top_k_min;
}


