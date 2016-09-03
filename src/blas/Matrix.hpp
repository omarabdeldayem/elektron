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
//	Matrix(int r, int c);
//	Matrix(T def_val, int r, int c)
//		: rows_(r)
//		, cols_(c)
//	{ mat_.resize(rows_ * cols_, def_val); }
	
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
//	void lud(Matrix<T>& l, Matrix<T>& u, Matrix<T>& p);
//	void svd(Matrix<T>& u, Matrix<T>& sigma, Matrix<T>& v);
	void qrd(Matrix<T, ROWS, ROWS>& Q, Matrix<T, ROWS, COLS>& R);	
	
	// UTILITIES
	void zeros();
	void ones();
	void eye();
	void print();

private:
	int rows_;
	int cols_;

	// Default: store matrix in vector, row-major form 
//	std::vector<T> mat_;
	std::array<T, ROWS * COLS> mat_;
//	Matrix<T> pivot();
//	Matrix<T> bidiag();
};

// Creates r x c identity mat_rix
//template <typename T, std::size_t ROWS, std::size_t COLS>
//Matrix<T>::Matrix(int r, int c) : rows_(r), cols_(c)
//{
//	mat_.resize(r * c, 0);
//	(*this).eye();
//}
//


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

/**
template <typename T>
Matrix<T> Matrix<T>::pivot()
{
	if (rows_ == cols_)
	{
		Matrix<T> id = Matrix<T>(rows_, cols_);
		for (int i = 0; i < rows_; i++)
		{
			T maxv = (*this)(i, i);
			int row = i;

			for (int j = i; j < rows_; j++)
			{
				if ((*this)(j, i) > maxv)
				{
					maxv = (*this)(j, i);
					row = j;
				}
			}

			if (i != row)
			{
				std::vector<T> i_tmp(mat_.begin() + i * rows_, mat_.begin() + (i * rows_) + cols_);
				std::vector<T> row_tmp(mat_.begin() + row * rows_, mat_.begin() + (row * rows_) + cols_);

				std::swap_ranges(mat_.begin() + i * rows_, mat_.begin() + (i * rows_) + cols_,
					   	row_tmp.begin());
				std::swap_ranges(mat_.begin() + row * rows_, mat_.begin() + (row * rows_) + cols_, 
						i_tmp.begin());
			}
		}
		return id;
	}
}


// Golub-Kahan-Lanczos Bidiagonalization
template <typename T>
Matrix<T> Matrix<T>::bidiag() 
{
	// TODO: Implement
	float beta_0 = 0;
	return NULL;
}

template <typename T>
void Matrix<T>::lud(Matrix<T>& l, Matrix<T>& u, Matrix<T>& p)
{
	p = p.pivot();
	Matrix<T> m2 = p * (*this);

	for (int j = 0; j < rows_; j++)
	{
		l(j, j) = 1.0;
		for (int i = 0; i < j + 1; i++)
		{
			double s = 0.0;
			for (int k = 0; k < i; k++)
			{
				s += u(k, j) * l(i, k);
			}
			u(i, j) = m2(i, j) - s;
		}
		for (int i = j; i < rows_; i++)
		{
			double s = 0.0;
			for (int k = 0; k < j; k++)
			{
				s += u(k, j) * l(i, k);
			}
			l(i, j) = (m2(i, j) - s) / u(j, j);
		}
	}
}

//template <typename T, std::size_t ROWS, std::size_t COLS>
//void Matrix<T>::svd(Matrix<T>& u, Matrix<T>& sigma, Matrix<T>& v)
//{
	// TODO: Implement
//}
**/
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
