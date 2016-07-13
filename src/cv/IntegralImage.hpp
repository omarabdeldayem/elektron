#ifndef INTEGRAL_IMAGE_H_
#define INTEGRAL_IMAGE_H_

#include "../blas/Matrix.hpp"

namespace nlib 
{

template <typename T>
void compute_iimg(const Matrix<T>& mat, const Matrix<T>& iim_mat);

}

#endif
