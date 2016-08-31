#ifndef MATRIX_H_
#define	MATRIX_H_
#define QR_THRESHOLD 0.0000000001

#include <stdexcept>
#include <iostream>
#include <iomanip>
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
	Matrix()  { };
	Matrix(int r, int c);
	Matrix(T def_val, int r, int c)
		: rows_(r)
		, cols_(c)
	{ mat_.resize(rows_ * cols_, def_val); }
	
	// OPERATOR OVERLOADS
	Matrix<T> operator*(Matrix<T> m);
	Matrix<T> operator*(T scalar);
	Matrix<T> operator/(T scalar);
	Matrix<T> operator+(Matrix<T> m);
	Matrix<T> operator-(Matrix<T> m);
	T& operator()(int r, int c) { return mat_.at(r * cols_ + c); };
	Matrix<T> operator()(int r_i, int r_f, int c_i, int c_f);
	
	// PRIVATE MEMBER ACCESS METHODS
	inline int rdim() { return rows_; };
	inline int cdim() { return cols_; };

	// MATRIX OPERATIONS
	T tr();
	double norm();
	Matrix<T> tpose();

	// MATRIX DECOMPOSITIONS
	void lud(Matrix<T>& l, Matrix<T>& u, Matrix<T>& p);
	void svd(Matrix<T>& u, Matrix<T>& sigma, Matrix<T>& v);
	void qrd(Matrix<T>& Q, Matrix<T>& R);	
	
	// UTILITIES
	void print();

private:
	int rows_;
	int cols_;

	// Default: store matrix in row-major form
	std::vector<T> mat_;

	Matrix<T> pivot();
	Matrix<T> bidiag();
};

// Creates r x c identity mat_rix
template <typename T>
Matrix<T>::Matrix(int r, int c) : rows_(r), cols_(c)
{
	mat_.resize(r * c, 0);
	for (auto iter_ = mat_.begin(); iter_ < mat_.end(); iter_+=cols_+1)
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
		for (int j = 0; j < res.cdim(); j++)
		{
			for (int k = 0; k < cols_; k++)
			{
				res(i, j) += mat_[i * cols_ + k] * m(k, j);
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
Matrix<T> Matrix<T>::operator/(T scalar)
{
	Matrix<T> res = Matrix<T>(rows_, cols_);

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * rows_ + j] / scalar;
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

template <typename T>
Matrix<T> Matrix<T>::operator-(Matrix<T> m)
{
	if (rows_ != m.rdim() && cols_ != m.cdim()) {
		// TODO: throw error	
	}

	Matrix<T> res = Matrix<T>(rows_, cols_);
	
	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * rows_ + j] - m(i, j);
		}
	}

	return res;

}

// Return submatrix from row i to row f, col i to col f exclusive
template <typename T>
Matrix<T> Matrix<T>::operator()(int r_i, int r_f, int c_i, int c_f)
{
	Matrix<T> sub = Matrix<T>(r_f - r_i, c_f - c_i);

	for (int i = r_i; i < r_f; i++)
	{
		for (int j = c_i; j < c_f; j++)
		{
			sub(i, j) = (*this)(i, j);
		}
	}

	return sub;
}

// Computes L_(2, 1) matrix norm 
template <typename T>
double Matrix<T>::norm()
{
	double norm = 0;
	
	for (int i = 0; i < rows_; i++)
	{
		double col_sum = 0;
		for (int j = 0; j < cols_; j++)
		{
			col_sum += pow(mat_[i * rows_ + j], 2.0);
		}
		norm += sqrt(col_sum);
	}
	return norm;
}
		
template <typename T>
Matrix<T> Matrix<T>::tpose()
{
	Matrix<T> mat_T = Matrix<T>(0, cols_, rows_);
	for (int i = 0; i < mat_T.rdim(); i++)
	{
		for (int j = 0; j < mat_T.cdim(); j++)
		{
			mat_T(i, j) = (*this)(j, i);
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
		for (auto iter_ = mat_.begin(); iter_ < iter_.end(); iter_+=cols_+1)
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

// QR Decomposition using householder reflections
template <typename T>
void Matrix<T>::qrd(Matrix<T>& Q, Matrix<T>& R)
{
	// For numerical stability
	double epsilon, alpha;
	
	Q = Matrix<T>(rows_, rows_);
	R = *this;

	Matrix<T> u, v;
	Matrix<T> P(rows_, rows_), I(rows_, rows_);
	
	for (int j = 0; j < cols_; j++) 
	{
		u = Matrix<T>(0, rows_, 1);
		v = Matrix<T>(0, rows_, 1);
		
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

		if (epsilon > 0.0000000001)
		{
			for (int i = j; i < rows_; i++) v(i, 0) /= epsilon;

			P = I - (v * v.tpose()) * 2.0;
			R = P * R;
			Q = Q * P;
		}	
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
		std::cout << std::setw(16) << val << " ";
		
		if (i % cols_ == 0) {
			std::cout << "\n";
		}
	}
	std::cout << std::endl;
}


}

#endif 
