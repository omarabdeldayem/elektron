#ifndef ELEKTRON_DERIVATIVES_H_
#define ELEKTRON_DERIVATIVES_H_

#include "Matrix.hpp"

#include <cstddef>

namespace elektron
{

template <typename T_>
inline double deriv_newton(T_ y1, T_ y2, double delta)
{
	return static_cast<double>((y2 - y1)/delta);
}

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<double, R_, C_> jacob(const Matrix<T_, R_, C_>& f1, const Matrix<T_, R_, C_>& f2, const RVec<T_, R_>& deltas)
{
	Matrix<T_, R_, C_> J = f2 - f1;

	for (int i = 0; i < f2.rdim(); i++)
	{
		for (int j = 0; j < f2.cdim(); j++)	
		{
			J(i, j) = J(i, j) / deltas(0, j); 
		}
	}

	return J;
}

} // End of namespace elektron
#endif
