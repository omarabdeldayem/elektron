#ifndef MATRIX_H
#define	MATRIX_H
#endif 

#include "Vector.hpp"

namespace nanos
{
	template <typename T>
	class Matrix
	{
	public:
		Matrix(T** T_matrix, int r, int c);

		inline uint_fast16_t rdim();
		inline uint_fast16_t cdim();

		Matrix transpose();
		Matrix ctranspose();
		Matrix LUD();
		Matrix QRD();
		Matrix EigenD();
		Matrix CholeskyD();

	private:
		uint_fast16_t r_dim, c_dim;
	};

}
