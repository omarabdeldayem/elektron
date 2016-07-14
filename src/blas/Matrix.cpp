#include "Matrix.hpp"

namespace nlib
{

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
std::vector<T> Matrix<T>::operator[] (int i)
{
	if (i < 0 || i > mat_.size())
	{
		throw std::out_of_range("Index out of bounds.");
	}
	typename std::vector<T>::const_iterator begin = mat_.begin() + (i * cols_);
	typename std::vector<T>::const_iterator end = begin + cols_;
	std::vector<T> new_vec(begin, end);	
	return new_vec;
}

template <typename T>
Matrix<T> Matrix<T>::operator* (Matrix<T> a)
{
	Matrix<T> res = Matrix<T>(0, rows_, a.cdim());

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < a.cdim(); j++)
		{
			for (int k = 0; k < cols_; k++)
			{
				res[i][j] += mat_[i * rows_ + cols_] * a[k][j];
			}
		}
	}

	return res;
}

template <typename T>
void Matrix<T>::operator= (std::vector<std::vector<T>>& mat_)
{
	mat_ = m;
	rows_ = m.size();
	cols_ = m[0].size();
}

template <typename T>
Matrix<T> Matrix<T>::tpose()
{
	Matrix<T> mat_T = Matrix<T>(NULL, cols_, rows_);

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			mat_T[i][j] = m[j][i];
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
			T maxv = m[i][i];
			int row = i;

			for (int j = i; j < rows_; j++)
			{
				if (m[j][i] > maxv)
				{
					maxv = m[j][i];
					row = j;
				}
			}

			if (i != row)
			{
				std::vector<T> tmp = id[i];
				id[i] = id[row];
				id[row] = tmp;
			}
		}
		return id;
	}
	
}

template <typename T>
void Matrix<T>::luD(Matrix<T>& l, Matrix<T>& u, Matrix<T>& p)
{
	p = p.pivot();
	Matrix<T> m2 = p * m;

	for (int j = 0; j < rows_; j++)
	{
		l[j][j] = 1.0;
		for (int i = 0; i < j + 1; i++)
		{
			double s = 0.0;
			for (int k = 0; k < i; k++)
			{
				s += u[k][j] * l[i][k];
			}
			u[i][j] = m2[i][j] - s;
		}
		for (int i = j; i < rows_; i++)
		{
			double s = 0.0;
			for (int k = 0; k < j; k++)
			{
				s += u[k][j] * l[i][k];
			}
			l[i][j] = (m2[i][j] - s) / u[j][j];
		}
	}
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
