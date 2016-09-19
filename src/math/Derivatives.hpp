#ifndef ELEKTRON_DERIVATIVES_H_
#define ELEKTRON_DERIVATIVES_H_

#include "Matrix.hpp"

#include <cstddef>

namespace elektron
{

template <typename T_>
inline double deriv_newton(T_ y1, T_ y2, double delta)
{
	return static_cast<double>(y2 - y1)/delta;
}

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, R_, C_> jacob(const Matrix<T_, R_, C_>& M)
{
	Matrix<T_, R_, C_> J = M;
}

} // End of namespace elektron
#endif
