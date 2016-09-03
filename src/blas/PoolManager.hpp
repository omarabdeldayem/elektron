#ifndef POOL_MANAGER_HPP_
#define POOL_MANAGER_HPP_

#define POOL_LEN_DEFAULT 512

#include <array>
#include <memory>

namespace nlib
{

template <typename T>
class PoolManager
{
public:
	PoolManager(bool ow) : overwrite_(ow) { };

private:

#if (POOL_LEN && POOL_LEN > POOL_LEN_DEFAULT) 
	// Requested pool size likely too large for stack - force into data or .bss
	static std::array<T, POOL_LEN> pool_;
#elif (POOL_LEN && POOL_LEN < POOL_LEN_DEFAULT)
	std::array<T, POOL_LEN> pool_; 
#else
	std::array<T, POOL_LEN_DEFAULT> pool_;
#endif
	
	std::unique_ptr<T> pool_head_;

	bool is_full_ = false;
	bool overwrite_ = false;

};

}
#endif
