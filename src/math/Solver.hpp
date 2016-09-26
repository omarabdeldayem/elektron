#ifndef ELEKTRON_SOLVER_HPP_
#define ELEKTRON_SOLVER_HPP_

#include <cstddef>

#include "Matrix.hpp"

namespace elektron
{

// forward-back substitution process for LU decomp
template <typename T, std::size_t ROWS, std::size_t COLS>
void lu_solve(const Matrix<T, ROWS, COLS>& L, const Matrix<T, ROWS, COLS>& U, const Matrix<T, ROWS, 1> b, Matrix<T, ROWS, 1> x)
{
	Matrix<T, ROWS, 1> y_tmp;

	for (int i = 0; i < L.r(); i++)
	{
		y_tmp(i, 0) = b(i, 0);

		for (int j = 0; j < i; j++)
		{
			y_tmp(i, 0) -= L(i, j) * y_tmp(j, 0);
		}

		y_tmp(i, 0) /= L(i, i);
	}

	for (int i = L.r() - 1; i >= 0; i--)
	{
		x(i, 0) = y_tmp(i, 0);

		for (int j = i + 1; j < L.r(); j++)
		{
			x(i, 0) -= U(i, j) * x(j, 0);
		}

		x(i, 0) /= U(i, i);
	}
}

} // End of namespace elektron

#endif
