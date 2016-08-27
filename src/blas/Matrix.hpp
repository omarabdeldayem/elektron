#ifndef MATRIX_H_
#define	MATRIX_H_

#include <stdexcept>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>

namespace nlib
{

template<typename T>
class Matrix
{
public:
	// CONSTRUCTORS
	Matrix(int r, int c);

	Matrix(T def_val, int r, int c)
		: rows_(r)
		, cols_(c)
	{ mat_.resize(rows_ * cols_, def_val); }
	
	// OPERATOR OVERLOADS
	Matrix<T> operator*(Matrix<T> m);
	Matrix<T> operator*(T scalar);
	Matrix<T> operator+(Matrix<T> m);
	T& operator()(int r, int c) { return mat_.at(r * rows_ + c); };
	
	// PRIVATE MEMBER ACCESS METHODS
	inline int rdim() { return rows_; };
	inline int cdim() { return cols_; };

	// MATRIX OPERATIONS
	T tr();
	T norm();
	Matrix<T> tpose();
	void lud(Matrix<T>& l, Matrix<T>& u, Matrix<T>& p);
	void svd(Matrix<T>& u, Matrix<T>& sigma, Matrix<T>& v);	
	//Matrix QRD();
	//Matrix eigenD();
	//Matrix choleskyD();
	void print();

private:
	int rows_;
	int cols_;

	// Default store mat_rix in row-major form
	std::vector<T> mat_;
	typename std::vector<T>::iterator iter_;

	Matrix<T> pivot();
	Matrix<T> bidiag();
};

// Creates r x c identity mat_rix
template <typename T>
Matrix<T>::Matrix(int r, int c) : rows_(r), cols_(c)
{
	mat_.resize(r * c, 0);
	for (iter_ = mat_.begin(); iter_ < mat_.end(); iter_+=cols_+1)
	{
		*iter_ = 1;
	}
}

template <typename T>
Matrix<T> Matrix<T>::operator*(Matrix<T> m)
{
	Matrix<T> res = Matrix<T>(0, rows_, m.cdim());

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < m.cdim(); j++)
		{
			for (int k = 0; k < cols_; k++)
			{
				res(i, j) += mat_[i * rows_ + k] * m(i, j);
			}
		}
	}

	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(T scalar)
{
	Matrix<T> res = Matrix<T>(rows_, cols_);

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * rows_ + j] * scalar;
		}
	}
	
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(Matrix<T> m)
{
	if (rows_ != m.rdim() && cols_ != m.cdim()) {
		// TODO: throw error	
	}

	Matrix<T> res = Matrix<T>(rows_, cols_);
	
	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * rows_ + j] + m(i, j);
		}
	}

	return res;
}

// Computes L_(2, 1) matrix norm 
template <typename T>
T Matrix<T>::norm()
{
	T norm = 0;
	
	for (int i = 0; i < rows_; i++)
	{
		T col_sum = 0;
		for (int j = 0; j < cols_; j++)
		{
			col_sum += mat_(i, j);
		}
		norm += sqrt(col_sum);
	}

	return norm;
}
		
template <typename T>
Matrix<T> Matrix<T>::tpose()
{
	Matrix<T> mat_T = Matrix<T>(NULL, cols_, rows_);

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			mat_T(i, j) = mat_(j, i);
		}
	}

	return mat_T;
}

template <typename T>
T Matrix<T>::tr()
{
	if (rows_ == cols_)
	{
		T sum = 0;
		for (iter_ = mat_.begin(); iter_ < iter_.end(); iter_+=cols_+1)
		{
			sum += *iter_;
		}
		return sum;
	}
	return NULL;
}

template <typename T>
Matrix<T> Matrix<T>::pivot()
{
	if (rows_ == cols_)
	{
		Matrix<T> id = Matrix<T>(rows_, cols_);
		for (int i = 0; i < rows_; i++)
		{
			T maxv = mat_(i, i);
			int row = i;

			for (int j = i; j < rows_; j++)
			{
				if (mat_(j, i) > maxv)
				{
					maxv = mat_(j, i);
					row = j;
				}
			}

			if (i != row)
			{
				std::vector<T> i_tmp(iter_.begin() + i * rows_, iter_.begin() + (i * rows_) + cols_);
				std::vector<T> row_tmp(iter_.begin() + row * rows_, iter_.begin() + (row * rows_) + cols_);

				std::swap_ranges(iter_.begin() + i * rows_, iter_.begin() + (i * rows_) + cols_,
					   	row_tmp.begin());
				std::swap_ranges(iter_.begin() + row * rows_, iter_.begin() + (row * rows_) + cols_, 
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
	Matrix<T> m2 = p * mat_;

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

template <typename T>
void Matrix<T>::svd(Matrix<T>& u, Matrix<T>& sigma, Matrix<T>& v)
{
	// TODO: Implement
}

// ------ DEBUG ------ //
template <typename T>
void Matrix<T>::print()
{
	int i = 0;

	for (T val : mat_)
	{
		i += 1;	
		std::cout << val << " ";
		
		if (i % rows_ == 0) {
			std::cout << "\n";
		}
	}
	std::cout.flush();
}


}

#endif 
