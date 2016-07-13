#ifndef LINREGRESSION_H_
#define LINREGRESSION_H_

#include "Regression.hpp"

namespace nlib
{

template <typename T>
class LinRegression : Regression 
{ 
public:
	LinRegression(Matrix<T> X, Matrix<T> y)
		: X(X)
		, Y(X)
		, is_fit(false)
	{ }
			 
private:
}

}

#endif
