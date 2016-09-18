#ifndef INTEGRAL_IMAGE_H_
#define INTEGRAL_IMAGE_H_

#include "../math/Matrix.hpp"

#include <cstddef>

namespace elektron 
{
namespace cv
{

// Computes the integral image (summed-area table) of matrix m
template <typename T, std::size_t ROWS, std::size_t COLS>
void compute_iimg(const Matrix<T, ROWS, COLS>& m, Matrix<T, ROWS, COLS>& iimg_m)
{
	// Required initial values	
	iimg_m(0, 0) = m(0, 0);
	iimg_m(0, 0) = m(0, 0) + m(0, 1);
	iimg_m(0, 0) = m(0, 0) + m(1, 0);

	for (int i = 1; i < m.rdim(); i++)
	{
		for (int j = 1; j < m.cdim(); j++)
		{
			iimg_m(i, j) = m(i, j) + iimg_m(i-1, j) + iimg_m(i, j-1) -
				iimg_m(i-1, j-1);
		}
	}
}

} // end of namespace cv
} // end of namespace nlib

#endif
