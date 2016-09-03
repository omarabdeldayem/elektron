#ifndef MATRIX_H_
#define	MATRIX_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <numeric>
#include <cmath>
#include <cstddef>

namespace nlib
{

const double QR_THRESHOLD = 0.0000000001;

template<typename T, std::size_t ROWS, std::size_t COLS>
class Matrix
{
public:
	// CONSTRUCTORS
	Matrix();
	
	// OPERATOR OVERLOADS
	template<std::size_t MROWS, std::size_t MCOLS>
	Matrix<T, ROWS, MCOLS> operator*(Matrix<T, MROWS, MCOLS> M);
	
	Matrix<T, ROWS, COLS> operator*(T scalar);
	Matrix<T, ROWS, COLS> operator/(T scalar);
	Matrix<T, ROWS, COLS> operator+(Matrix<T, ROWS, COLS> M);
	Matrix<T, ROWS, COLS> operator-(Matrix<T, ROWS, COLS> M);
	T& operator()(int r, int c) { return mat_[r * cols_ + c]; };

	template<std::size_t MROWS, std::size_t MCOLS>
	void sub(Matrix<T, MROWS, MCOLS> M, int r_i, int r_f, int c_i, int c_f);
	
	// PRIVATE MEMBER ACCESS METHODS
	inline int rdim() { return rows_; };
	inline int cdim() { return cols_; };

	// MATRIX OPERATIONS
	T trace();
	double norm();
	Matrix<T, COLS, ROWS> tpose();

	// MATRIX DECOMPOSITIONS
	void lud(Matrix<T, ROWS, COLS>& L, Matrix<T, ROWS, COLS>& U);
	void qrd(Matrix<T, ROWS, ROWS>& Q, Matrix<T, ROWS, COLS>& R);		
//	void svd(Matrix<T>& u, Matrix<T>& sigma, Matrix<T>& v);

	// UTILITIES
	void zeros();
	void ones();
	void eye();
	void print();

private:
	int rows_;
	int cols_;

	// Default: store matrix in array, row-major form 
	std::array<T, ROWS * COLS> mat_;
};

template <typename T, std::size_t ROWS, std::size_t COLS>
Matrix<T, ROWS, COLS>::Matrix()
{
	rows_ = static_cast<int>(ROWS);
	cols_ = static_cast<int>(COLS);

	this->zeros();
}

template <typename T, std::size_t ROWS, std::size_t COLS>
template <std::size_t MROWS, std::size_t MCOLS>
Matrix<T, ROWS, MCOLS> Matrix<T, ROWS, COLS>::operator*(Matrix<T, MROWS, MCOLS> M)
{
	Matrix<T, ROWS, MCOLS> res;

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < res.cdim(); j++)
		{
			for (int k = 0; k < cols_; k++)
			{
				res(i, j) += mat_[i * cols_ + k] * M(k, j);
			}
		}
	}

	return res;
}

template <typename T, std::size_t ROWS, std::size_t COLS>
Matrix<T, ROWS, COLS> Matrix<T, ROWS, COLS>::operator*(T scalar)
{
	Matrix<T, ROWS, COLS> res = Matrix<T, ROWS, COLS>();

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * cols_ + j] * scalar;
		}
	}
	
	return res;
}

template <typename T, std::size_t ROWS, std::size_t COLS>
Matrix<T, ROWS, COLS> Matrix<T, ROWS, COLS>::operator/(T scalar)
{
	Matrix<T, ROWS, COLS> res = Matrix<T, ROWS, COLS>();

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * cols_ + j] / scalar;
		}
	}

	return res;
}

template <typename T, std::size_t ROWS, std::size_t COLS>
Matrix<T, ROWS, COLS> Matrix<T, ROWS, COLS>::operator+(Matrix<T, ROWS, COLS> M)
{
	if (rows_ != M.rdim() && cols_ != M.cdim()) {
		// TODO: throw error	
	}

	Matrix<T, ROWS, COLS> res = Matrix<T, ROWS, COLS>();
	
	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * cols_ + j] + M(i, j);
		}
	}

	return res;
}

template <typename T, std::size_t ROWS, std::size_t COLS>
Matrix<T, ROWS, COLS> Matrix<T, ROWS, COLS>::operator-(Matrix<T, ROWS, COLS> M)
{
	if (rows_ != M.rdim() && cols_ != M.cdim()) {
		// TODO: throw error	
	}

	Matrix<T, ROWS, COLS> res = Matrix<T, ROWS, COLS>();
	
	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * cols_ + j] - M(i, j);
		}
	}

	return res;

}

// Return submatrix from row i to row f, col i to col f exclusive
// 
template <typename T, std::size_t ROWS, std::size_t COLS>
template <std::size_t MROWS, std::size_t MCOLS>
void Matrix<T, ROWS, COLS>::sub(Matrix<T, MROWS, MCOLS> M, int r_i, int r_f, int c_i, int c_f)
{
	if (static_cast<int>(MROWS) != (r_f - r_i) || static_cast<int>(MCOLS) != (c_f - c_i))
	{
		return;
	}

	for (int i = r_i; i < r_f; i++)
	{
		for (int j = c_i; j < c_f; j++)
		{
			M(i, j) = (*this)(i, j);
		}
	}
}

// Computes L_(2, 1) matrix norm 
template <typename T, std::size_t ROWS, std::size_t COLS>
double Matrix<T, ROWS, COLS>::norm()
{
	double norm = 0;
	
	for (int i = 0; i < rows_; i++)
	{
		double col_sum = 0;
		for (int j = 0; j < cols_; j++)
		{
			col_sum += pow(mat_[i * cols_ + j], 2.0);
		}
		norm += sqrt(col_sum);
	}
	return norm;
}
		
template <typename T, std::size_t ROWS, std::size_t COLS>
Matrix<T, COLS, ROWS> Matrix<T, ROWS, COLS>::tpose()
{
	Matrix<T, COLS, ROWS> mat_T = Matrix<T, COLS, ROWS>();
	for (int i = 0; i < mat_T.rdim(); i++)
	{
		for (int j = 0; j < mat_T.cdim(); j++)
		{
			mat_T(i, j) = (*this)(j, i);
		}
	}

	return mat_T;
}

template <typename T, std::size_t ROWS, std::size_t COLS>
T Matrix<T, ROWS, COLS>::trace()
{
	if (rows_ == cols_)
	{
		T sum = 0;
		for (auto it_ = mat_.begin(); it_ < it_.end(); it_+=cols_+1)
		{
			sum += *it_;
		}
		return sum;
	}
	return NULL;
}

// QR Decomposition using householder reflections
template <typename T, std::size_t ROWS, std::size_t COLS>
void Matrix<T, ROWS, COLS>::qrd(Matrix<T, ROWS, ROWS>& Q, Matrix<T, ROWS, COLS>& R)
{
	// For numerical stability
	double epsilon = 0.0;
	double alpha = 0.0;
	
	Q.eye();
	R = *this;

	Matrix<T, ROWS, 1> u; 
	Matrix<T, ROWS, 1> v;
	Matrix<T, ROWS, ROWS> P;
   	Matrix<T, ROWS, ROWS> I;

	I.eye();
		
	for (int j = 0; j < cols_; j++) 
	{
		u.zeros();
		v.zeros();
		
		epsilon = 0.0;

		for (int i = j; i < rows_; i++)
		{
			u(i, 0) = R(i, j);
			epsilon += u(i, 0) * u(i, 0);   	
		}
		
		epsilon = sqrt(epsilon);
		alpha = copysign(epsilon, -u(j, 0)); // If you replace epsilon here with u.norm(), the Q matrix result is symmetric... figure out why
		epsilon = 0.0;	
		
		for (int i = j; i < rows_; i++)
		{
			v(i, 0) = i == j ? u(i, 0) + alpha : u(i, 0);
			epsilon += v(i, 0) * v(i, 0);
		}

		epsilon = sqrt(epsilon);

		if (epsilon > QR_THRESHOLD)
		{
			for (int i = j; i < rows_; i++) v(i, 0) /= epsilon;

			P = I - (v * v.tpose()) * 2.0;
			R = P * R;
			Q = Q * P;
		}	
	}
}

// LU Decomposition via Dolittle algorithm
template <typename T, std::size_t ROWS, std::size_t COLS>
void Matrix<T, ROWS, COLS>::lud(Matrix<T, ROWS, COLS>& L, Matrix<T, ROWS, COLS>& U)
{
	for (int i = 0; i < rows_; i++) 
	{
		for (int j = 0; j < rows_; j++)
		{
			if (j < i)
			{
				L(j, i) = 0;
			}
			else 
			{
				L(j, i) = (*this)(j, i);
				for (int k = 0; k < i; k++)
				{
					L(j, i) = L(j, i) - L(j, k) * U(k, i);
				}
			}
		}
		for (int j = 0; j < rows_; j++)
		{
			if (j < i)
			{
				U(i, j) = 0;
			}
			else if (j == i)
			{
				U(i, j) = 1;
			}
			else
			{
				U(i, j) = (*this)(i, j) / L(i, i);
				for (int k = 0; k < i; k++)
				{
					U(i, j) = U(i, j) - (L(i, k) * U(k , j)) / L(i, i);
				}
			}
		}
	}
}

//template <typename T, std::size_t ROWS, std::size_t COLS>
//void Matrix<T>::svd(Matrix<T>& u, Matrix<T>& sigma, Matrix<T>& v)
//{
	// TODO: Implement
//}

template <typename T, std::size_t ROWS, std::size_t COLS>
void Matrix<T, ROWS, COLS>::zeros()
{
	for (auto it_ = mat_.begin(); it_ < mat_.end(); it_++)
	{
		*it_ = 0;
	}
}

template <typename T, std::size_t ROWS, std::size_t COLS>
void Matrix<T, ROWS, COLS>::ones()
{
	for (auto it_ = mat_.begin(); it_ < mat_.end(); it_++)
	{
		*it_ = 1;
	}
}

template <typename T, std::size_t ROWS, std::size_t COLS>
void Matrix<T, ROWS, COLS>::eye()
{
	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			mat_[i * cols_ + j] =  (i == j) ? 1 : 0;
		}
	}
}

// ------ DEBUG ------ //
template <typename T, std::size_t ROWS, std::size_t COLS>
void Matrix<T, ROWS, COLS>::print()
{
	int i = 0;

	for (T val : mat_)
	{
		i += 1;	
		std::cout << std::setw(16) << val << " ";
		
		if (i % cols_ == 0) {
			std::cout << "\n";
		}
	}
	std::cout << std::endl;
}


}

#endif 
