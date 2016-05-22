#ifndef VECTOR_H
#define VECTOR_H
#endif

#include <stdio.h>
#include <vector>

namespace nanos
{
	template <typename T>
	class Vector
	{
	public:
		Vector(std::vector<T> elems);
		Vector(T* T_array, int len);

		inline std::vector<T>* get_rform();
		inline std::vector<T>* get_cform();
		inline uint_fast16_t get_dim();

	private:
		uint_fast16_t dim;
		std::vector<T> elements;
	};
}