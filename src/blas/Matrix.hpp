#ifndef MATRIX_H_
#define	MATRIX_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <vector>

namespace nlib
{

const int SVD_NUM_ITERATIONS = 50;
const double QR_EPSILON = 0.0000000001;
const double SVD_EPSILON = 1e-17;


enum MatInit { z, o, i };

template<typename T_, std::size_t R_, std::size_t C_>
class Matrix
{
public:
	// CONSTRUCTORS
	Matrix();
	Matrix(MatInit m_init);
	Matrix(const std::array<T_, R_ * C_>& arr);
	Matrix(const T_* arr);
//	Matrix(const Matrix<T_, R_, C_>& M); Leave until fully implemented

	// OPERATOR OVERLOADS
	template<std::size_t MR_, std::size_t MC_>
	Matrix<T_, R_, MC_> operator*(Matrix<T_, MR_, MC_> M);
	
	Matrix<T_, R_, C_> operator*(T_ scalar);
	Matrix<T_, R_, C_> operator/(T_ scalar);
	Matrix<T_, R_, C_> operator+(Matrix<T_, R_, C_> M);
	Matrix<T_, R_, C_> operator-(Matrix<T_, R_, C_> M);
	T_& operator()(int r, int c);

	template<std::size_t MR_, std::size_t MC_>
	void sub(Matrix<T_, MR_, MC_> M, int r_i, int r_f, int c_i, int c_f);
	
	// PRIVATE MEMBER ACCESS METHODS
	inline int rdim() { return rows_; };
	inline int cdim() { return cols_; };

	// MATRIX OPERATIONS
	T_ trace();
	double norm();
	Matrix<T_, C_, R_> tpose();

	// MATRIX DECOMPOSITIONS
	void lud(Matrix<T_, R_, C_>& L, Matrix<T_, R_, C_>& U);
	void qrd(Matrix<T_, R_, R_>& Q, Matrix<T_, R_, C_>& R);		
	void svd(Matrix<T_, R_, C_>& U, Matrix<T_, C_, C_>& S, Matrix<T_, C_, C_>& V_T);

	// UTILITIES
	void zeros();
	void ones();
	void eye();
	void print();

private:
	int rows_;
	int cols_;

#ifdef NLIB_USE_HEAP 
	std::vector<T_> mat_;
#else
	std::array<T_, R_ * C_> mat_;
#endif

	// SVD UTILITIES
	double pythagorean(double a, double b);
	void svd_reord(Matrix<T_, R_, C_>& U, Matrix<T_, C_, C_>& S, Matrix<T_, C_, C_>& V_T);

};

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, R_, C_>::Matrix()
{
	rows_ = static_cast<int>(R_);
	cols_ = static_cast<int>(C_);

#ifdef NLIB_USE_HEAP
	mat_.resize(rows_ * cols_);
#endif

	this->zeros();
}

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, R_, C_>::Matrix(MatInit m_init)
{
	rows_ = static_cast<int>(R_);
	cols_ = static_cast<int>(C_);

#ifdef NLIB_USE_HEAP
	mat_.resize(rows_ * cols_);
#endif

	switch(m_init)
	{
		case z : this->zeros(); break;
		case o : this->ones(); 	break;
		case i : this->eye();	break;
	}
}

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, R_, C_>::Matrix(const std::array<T_, R_ * C_>& arr)
{
	rows_ = static_cast<int>(R_);
	cols_ = static_cast<int>(C_);

	mat_ = arr;
}

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, R_, C_>::Matrix(const T_* arr)
{
	rows_ = static_cast<int>(R_);
	cols_ = static_cast<int>(C_);

	std::copy(std::begin(arr), std::end(arr), std::begin(mat_));
}

// Don't uncomment this until fully implemented, otherwise you'll segfault
//template <typename T_, std::size_t R_, std::size_t C_>
//Matrix<T_, R_, C_>::Matrix(const Matrix<T_, R_, C_>& M)
//{
	// TODO
//}

template <typename T_, std::size_t R_, std::size_t C_>
template <std::size_t MR_, std::size_t MC_>
Matrix<T_, R_, MC_> Matrix<T_, R_, C_>::operator*(Matrix<T_, MR_, MC_> M)
{
	Matrix<T_, R_, MC_> res;

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

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, R_, C_> Matrix<T_, R_, C_>::operator*(T_ scalar)
{
	Matrix<T_, R_, C_> res = Matrix<T_, R_, C_>();
	
	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * cols_ + j] * scalar;
		}
	}
	
	return res;
}

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, R_, C_> Matrix<T_, R_, C_>::operator/(T_ scalar)
{
	Matrix<T_, R_, C_> res = Matrix<T_, R_, C_>();

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * cols_ + j] / scalar;
		}
	}

	return res;
}

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, R_, C_> Matrix<T_, R_, C_>::operator+(Matrix<T_, R_, C_> M)
{
	if (rows_ != M.rdim() && cols_ != M.cdim()) {
		// T_ODO: throw error	
	}

	Matrix<T_, R_, C_> res = Matrix<T_, R_, C_>();
	
	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * cols_ + j] + M(i, j);
		}
	}

	return res;
}

template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, R_, C_> Matrix<T_, R_, C_>::operator-(Matrix<T_, R_, C_> M)
{
	if (rows_ != M.rdim() && cols_ != M.cdim()) {
		// T_ODO: throw error	
	}

	Matrix<T_, R_, C_> res = Matrix<T_, R_, C_>();
	
	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			res(i, j) = mat_[i * cols_ + j] - M(i, j);
		}
	}

	return res;

}

template <typename T_, std::size_t R_, std::size_t C_>
T_& Matrix<T_, R_, C_>::operator()(int r, int c)
{
	assert((r>=0) && (c>=0) && (r<rows_) && (c<cols_));
	return mat_[r * cols_ + c];
}

// Return submatrix from row i to row f, col i to col f exclusive
template <typename T_, std::size_t R_, std::size_t C_>
template <std::size_t MR_, std::size_t MC_>
void Matrix<T_, R_, C_>::sub(Matrix<T_, MR_, MC_> M, int r_i, int r_f, int c_i, int c_f)
{
	if (static_cast<int>(MR_) != (r_f - r_i) || static_cast<int>(MC_) != (c_f - c_i))
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
template <typename T_, std::size_t R_, std::size_t C_>
double Matrix<T_, R_, C_>::norm()
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
		
template <typename T_, std::size_t R_, std::size_t C_>
Matrix<T_, C_, R_> Matrix<T_, R_, C_>::tpose()
{
	Matrix<T_, C_, R_> mat_T_ = Matrix<T_, C_, R_>();
	for (int i = 0; i < mat_T_.rdim(); i++)
	{
		for (int j = 0; j < mat_T_.cdim(); j++)
		{
			mat_T_(i, j) = (*this)(j, i);
		}
	}

	return mat_T_;
}

template <typename T_, std::size_t R_, std::size_t C_>
T_ Matrix<T_, R_, C_>::trace()
{
	if (rows_ == cols_)
	{
		T_ sum = 0;
		for (auto it_ = mat_.begin(); it_ < it_.end(); it_+=cols_+1)
		{
			sum += *it_;
		}
		return sum;
	}
	return NULL;
}

// QR Decomposition using householder reflections
template <typename T_, std::size_t R_, std::size_t C_>
void Matrix<T_, R_, C_>::qrd(Matrix<T_, R_, R_>& Q, Matrix<T_, R_, C_>& R)
{
	// For numerical stability
	double epsilon = 0.0;
	double alpha = 0.0;
	
	Q.eye();
	R = *this;

	Matrix<T_, R_, 1> u; 
	Matrix<T_, R_, 1> v;
	Matrix<T_, R_, R_> P;
   	Matrix<T_, R_, R_> I;

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

		if (epsilon > QR_EPSILON)
		{
			for (int i = j; i < rows_; i++) v(i, 0) /= epsilon;

			P = I - (v * v.tpose()) * 2.0;
			R = P * R;
			Q = Q * P;
		}	
	}
}

// LU Decomposition via Dolittle algorithm
template <typename T_, std::size_t R_, std::size_t C_>
void Matrix<T_, R_, C_>::lud(Matrix<T_, R_, C_>& L, Matrix<T_, R_, C_>& U)
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

template <typename T_, std::size_t R_, std::size_t C_>
void Matrix<T_, R_, C_>::svd(Matrix<T_, R_, C_>& U, Matrix<T_, C_, C_>& S, Matrix<T_, C_, C_>& V_T)
{
	int i = 0;
	int j = 0;
	int jj = 0;
	int k = 0;
	int its = 0;
	int flag = 0;
	int l = 0;	
	int nm = 0;

	double scale = 0.0;
	double anorm = 0.0;

	double c = 0.0;
	double f = 0.0;
	double g = 0.0;
	double h = 0.0;
	double s = 0.0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	Matrix<T_, 1, C_> rv1;
	
	U = *this;
	
	// Householder reduction into bidiagonal form
	for (i = 0; i < cols_; i++)
	{
		l = i + 2;

		rv1(0, i) = scale * g;
		
		// R_eset g, scale, s
		scale = 0.0;
		g = 0.0;
		s = 0.0;

		if (i < rows_)
		{
			for (k = i; k < rows_; k++)
			{
				scale += std::fabs(U(k, i));
			}

			if (scale)
			{
				for (k = i; k < rows_; k++)
				{
					U(k, i) = U(k, i) / scale;
					s += U(k, i) * U(k, i);
				}

				f = U(i, i);
				g = -std::copysign(sqrt(s), f);
				h = (f * g) - s;
				U(i, i) = f - g;

				for (j = l-1; j < cols_; j++)
				{
					s = 0.0;
					for (k = i; k < rows_; k++)
					{
						s += U(k, i) * U(k, j); 
					}
					
					f = s / h;
					
					for (k = i; k < rows_; k++)
					{
						U(k, j) += f * U(k, i);
					}
				}

				for (k = i; k < rows_; k++)
				{
					U(k, i) = U(k, i) * scale;
				}
			}
		}
		S(i, i) = g * scale;
	
		scale = 0.0;
		g = 0.0;
		s = 0.0;
		
		if (i + 1 <= rows_ && i + 1 != cols_)
		{
			for (k = l - 1; k < cols_; k++)
			{
				scale += std::fabs(U(i, k));
			}

			if (scale)
			{
				for (k = l - 1; k < cols_; k++)
				{
					U(i, k) = U(i, k) / scale;
					s += U(i, k) * U(i, k);
				}
			
				f = U(i, l - 1);
				g = -std::copysign(sqrt(s), f);
				h = (f * g) - s;
				U(i, l - 1) = f - g;

				for (k = l - 1; k < cols_; k++)
				{
					rv1(0, k) = U(i, k) / h;
				}

				for (j = l - 1; j < rows_; j++)
				{
					s = 0.0;
					for (k = l-1; k < cols_; k++)
					{
						s += U(j, k) * U(i, k);
					}
					for (k = l-1; k < cols_; k++)
					{
						U(j, k) += s * rv1(0, k);
					}
				}

				for (int k = l - 1; k < cols_; k++)
				{
					U(i, k) = U(i, k) * scale;
				}
			}
		}
		anorm = std::max(anorm, (std::fabs(S(i, i)) + std::fabs(rv1(0, i))));
	}
	// -------------------------------------------------------------- //
	
	// Accumulate left-hand / right-hand transformations	
	for (i = cols_ - 1; i >= 0; i--)
	{
		if (i < cols_ -1)
		{
			if (g)
			{
				for (j = l; j < cols_; j++)
				{
					V_T(j, i) = (U(i, j) / U(i, l)) / g;
				}
				for (j = l; j < cols_; j++)
				{
					s = 0.0;
					for (k = l; k < cols_; k++)
					{
						s += U(i, k) * V_T(k, j);
					}
					for (k = l; k < cols_; k++)
					{
						V_T(k, j) += s * V_T(k, i);
					}
				}
			}
			for (j = l; j < cols_; j++)
			{
				V_T(i, j) = 0.0;
				V_T(j, i) = 0.0;
			}
		}
		V_T(i, i) = 1.0;
		g = rv1(0, i);
		l = i;
	}
	
	for (i = std::min(rows_, cols_) - 1; i >= 0; i--)
	{
		l = i + 1;
		g = S(i, i);

		for (j = l; j < cols_; j++)
		{
			U(i, j) = 0.0;	
		}
		
		if (g)
		{
			g = 1.0 / g;
			
			for (j = l; j < cols_; j++)
			{
				s = 0.0;
				for (s = 0.0, k = l; k < rows_; k++)
				{
					s += U(k, i) * U(k, j);
				}
					
				f = (s / U(i, i)) * g;
					
				for (k = i; k < rows_; k++)
				{
					U(k, j) += f * U(k, i);
				}
			}
			
			for (j = i; j < rows_; j++)
			{
				U(j, i) = U(j, i) * g;
			}
		}
		else
		{
			for (j = i; j < rows_; j++)
			{
				U(j, i) = 0.0;
			}
		}	
		++U(i, i);
	}
	// ----------------------------------------------------------------//
	
	// Diagonalize
	for (k = cols_ - 1; k >= 0; k--)
	{
		for (its = 0; its <= SVD_NUM_ITERATIONS; its++)
		{
			flag = 1;

			for (l = k; l >= 0; l--)
			{
				nm = l - 1;

				if (l == 0 || std::fabs(rv1(0, l)) <= SVD_EPSILON * anorm)
				{
					flag = 0;
					break;
				}

				if (std::fabs(S(nm, nm)) <= SVD_EPSILON * anorm) 
				{ 
					break;
			   	}
			}
		
		
			if (flag)
			{
				c = 0.0;
				s = 1.0;

				for (i = l; i < k + 1; i++)
				{
					f = s * rv1(0, i);
					rv1(0, i) = c * rv1(0, i);
			
					if (std::fabs(f) <= SVD_EPSILON * anorm) { break; }
						
					g = S(i, i);
					h = pythagorean(f, g);
					S(i, i) = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;

					for (j = 0; j < rows_; j++)
					{
						y = U(j, nm);
						z = U(j, i);

						U(j, nm) = (y * c) + (z * s);
						U(j, i) = (z * c) - (y * s);
					}
				}
			}
			z = S(k, k);

			if (l == k)
			{
				// Make negative singular value non-negative
				if (z < 0.0)
				{
					S(k , k) = -z;
					for (j = 0; j < cols_; j++)
					{
						V_T(j, k) = (-V_T(j, k));
					}
				}
				break;
			}

			if (its == SVD_NUM_ITERATIONS)
			{
				// NO CONVERGENCE
				std::cout << "Early return, no convergence\n";
				return;
			}
		
			x = S(l, l);
			nm = k - 1;
			y = S(nm, nm);
			g = rv1(0, nm);
			h = rv1(0, k);
			f = (((y - z) * (y + z)) + ((g - h) * (g + h))) / (2.0 * h * y);
			g = pythagorean(f, 1.0);
			f = (((x - z) * (x + z)) + (h * ((y / (f + std::copysign(g, f))) - h))) / x;

			// QR_ transformation
			c = 1.0;
			s = 1.0;

			for (j = l; j <= nm; j++)
			{
				i = j + 1;
				g = rv1(0, i);
				y = S(i, i);
				h = s * g;
				g = c * g;
				z = pythagorean(f, h);
				rv1(0, j) = z;
				c = f / z;
				s = h / z;
				f = (x * c) + (g * s);
				g = (g * c) - (x * s);
				h = y * s;
				y = y * c;
				
				for (jj = 0; jj < cols_; jj++)
				{
					x = V_T(jj, j);
					z = V_T(jj, i);
					V_T(jj, j) = (x * c) + (z * s);
					V_T(jj, i) = (z * c) - (x * s);
				}
				
				z = pythagorean(f, h);
				S(j , j) = z;

				if (z)
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}

				f = c * g + s * y;
				x = c * y - s * g;
	
				for (jj = 0; jj < rows_; jj++)
				{
					y = U(jj, j);
					z = U(jj, i);
					U(jj, j) = (y * c) + (z * s);
					U(jj, i) = (z * c) - (y * s);
				}
			}
			rv1(0, l) = 0.0;
			rv1(0, k) = f;
			S(k, k) = x;
		}
	}

	svd_reord(U, S, V_T);
}

template <typename T_, std::size_t R_, std::size_t C_>
double Matrix<T_, R_, C_>::pythagorean(double a, double b)
{
	double absa = fabs(a);
	double absb = fabs(b);
	double ct = 0.0;
	double result = 0.0;

	if (absa > absb)
	{
		ct = absb / absa;
		result = absa * sqrt(1.0 + (ct * ct));
	}
	else if (absb > 0.0)
	{
		ct = absa / absb;
		result = absb * sqrt(1.0 + (ct * ct));
	}

	return result;
}

// Sorts S s.t. s_0 > s_1 > ... > s_k
// Also sorts corresponding cols of U and V by dec mag
// Maximizes number of positive elements
template <typename T_, std::size_t R_, std::size_t C_>
void Matrix<T_, R_, C_>::svd_reord(Matrix<T_, R_, C_>& U, Matrix<T_, C_, C_>& S, Matrix<T_, C_, C_>& V_T)
{
	int inc = 1;
	
	Matrix<T_, R_, 1> su;
	Matrix<T_, 1, C_> sv;

	do { inc *= 3; inc++; } while (inc <= cols_);

	do 
	{
		inc /= 3;

		for (int i = inc; i < cols_; i++)
		{
			double sw = S(i, i);
			
			for (int k = 0; k < rows_; k++)
			{
				su(k, 0) = U(k, i);
			}
			for (int k = 0; k < cols_; k++)
			{
				sv(0, k) = V_T(k, i);
			}
			
			int j = i;

			while(S(j-inc, j-inc) < sw)
			{
				S(j, j) = S(j-inc, j-inc);

				for (int k = 0; k < rows_; k++)
				{
					U(k, j) = U(k, j-inc);
				}

				for (int k = 0; k < cols_; k++)
				{
					V_T(k, j) = V_T(k, j-inc);
				}

				j -= inc;
				
				if (j < inc) break;
			}
			
			S(j, j) = sw;

			for (int k = 0; k < rows_; k++)
			{
				U(k, j) = su(k, 0);
			}
			for (int k = 0; k < cols_; k++)
			{
				V_T(k, j) = sv(0, k);
			}
		}
	} while (inc > 1);
	for (int k = 0; k < cols_; k++)
	{
		int s = 0;

		for (int i = 0; i < rows_; i++)
		{
			if (U(i, k) < 0.0) { s++; }
		}
		for (int j = 0; j < cols_; j++)
		{
			if (V_T(j, k) < 0.0) { s++; }
		}

		if (s > (rows_ + cols_) / 2)
		{
			for (int i = 0; i < rows_; i++)
			{
				U(i, k) = - U(i, k);
			}
			for (int j = 0; j < cols_; j++)
			{
				V_T(j, k) = -V_T(j, k);
			}
		}
	}
}

template <typename T_, std::size_t R_, std::size_t C_>
void Matrix<T_, R_, C_>::zeros()
{
	for (auto it_ = mat_.begin(); it_ < mat_.end(); it_++)
	{
		*it_ = 0;
	}
}

template <typename T_, std::size_t R_, std::size_t C_>
void Matrix<T_, R_, C_>::ones()
{
	for (auto it_ = mat_.begin(); it_ < mat_.end(); it_++)
	{
		*it_ = 1;
	}
}

template <typename T_, std::size_t R_, std::size_t C_>
void Matrix<T_, R_, C_>::eye()
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
template <typename T_, std::size_t R_, std::size_t C_>
void Matrix<T_, R_, C_>::print()
{
	int i = 0;

	for (T_ val : mat_)
	{
		i += 1;	
		std::cout << std::setw(16) << val << " ";
		
		if (i % cols_ == 0) {
			std::cout << "\n";
		}
	}
	std::cout << std::endl;
}


} // end of namespace nlib

#endif 
