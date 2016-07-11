#include <vector.h>

#include "../blas/Matrix.hpp"

namespace nlib 
{
	template <typename T>
	class Regression 
	{
	public:
		virtual void fit() = 0;
		virtual void predict(T x) = 0;
		
		inline void set_weights(Matrix<T> new_W) { weights = new_W; };
		inline void set_X(Matrix<T> new_X) { X = new_X; is_fit = false; };
		inline void set_Y(Matrix<T> new_Y) { Y = new_y; is_fit = false; };

	protected:
		// Regression model of the form
		// E(Y | X) = f(X, weights)
		Matrix<T> weights;
		Matrix<T> X;
		Matrix<T> Y;

		bool is_fit;

		virtual void grad(Matrix<T> X, Matrix<T> Y, Matrix<T> w) = 0;
	}
}
