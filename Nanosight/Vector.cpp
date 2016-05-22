#include "Vector.hpp"

template <typename T> 
nanos::Vector<T>::Vector(std::vector<T> elems)
{
}

template <typename T> 
nanos::Vector<T>::Vector(T* T_array, int len)
{
}

template <typename T> 
inline uint_fast16_t nanos::Vector<T>::get_dim()
{
	return dim;
}
