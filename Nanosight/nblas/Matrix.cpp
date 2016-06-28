#include "Matrix.hpp"

namespace nanos
{
	template class Matrix<uint_fast8_t>;
	template class Matrix<int>;
	template class Matrix<float>;
	template class Matrix<long>;

	template <typename T>
	Matrix<T>::Matrix()
	{

	}

	template <typename T>
	Matrix<T>::Matrix(std::vector<std::vector<T>>& T_matrix)
	{
		m = T_matrix;
		r_dim = T_matrix.size();
		c_dim = T_matrix[0].size();
	}

	template <typename T>
	Matrix<T>::Matrix(T def_val, int r, int c)
	{
		for (int i = 0; i < r; i++)
		{
			std::vector<T> t;
			m.push_back(t);

			for (int j = 0; j < c; j++)
			{
				m[i].push_back(def_val);
			}
		}

		r_dim = r;
		c_dim = c;
	}

	template <typename T>
	std::vector<T>& Matrix<T>::operator[] (int i)
	{
		if (i < 0 || i > m.size())
		{
			throw std::out_of_range("Index out of bounds.");
		}
		return m[i];
	}

	template <typename T>
	Matrix<T> Matrix<T>::operator* (Matrix<T> a)
	{
		Matrix<T> res = Matrix<T>(0, r_dim, a.cdim());

		for (int i = 0; i < r_dim; i++)
		{
			for (int j = 0; j < a.cdim(); j++)
			{
				for (int k = 0; k < c_dim; k++)
				{
					res[i][j] += m[i][k] * a[k][j];
				}
			}
		}

		return res;
	}

	template <typename T>
	Matrix<T> Matrix<T>::tpose()
	{
		Matrix<T> m_T = Matrix<T>(NULL, c_dim, r_dim);

		for (int i = 0; i < r_dim; i++)
		{
			for (int j = 0; j < c_dim; j++)
			{
				m_T[i][j] = m_T[j][i];
			}
		}

		return m_T;
	}

	template <typename T>
	T Matrix<T>::tr()
	{
		if (is_sqr)
		{
			T sum = 0;
			for (int i = 0; i < r_dim; i++)
			{
				sum += m[i][i];
			}
			return sum;
		}
		return NULL;
	}

	template <typename T>
	void Matrix<T>::LUD(Matrix<T>& l, Matrix<T>& u)
	{
		for (int j = 0; j < r_dim; j++)
		{
			for (int i = 0; i < r_dim; i++)
			{
				if (i <= j)
				{
					u[i][j] = m[i][j];
					for (int k = 0; k < i - 1; k++)
					{
						u[i][j] -= l[i][k] * u[k][j];
					}
					if (i == j)
					{
						l[i][j] = 1;
					}
					else
					{
						l[i][j] = 0;
					}
				}
				else
				{
					l[i][j] = m[i][j];
					for (int k = 0; k <= j - 1; k++)
					{
						l[i][j] -= l[i][k] * u[k][j];
					}
					l[i][j] /=	u[j][j];
					u[i][j] = 0;
				}
			}
		}
	}

	template <typename T>
	void Matrix<T>::print()
	{
		for (std::vector<T> i : m)
		{
			for (T j : i)
			{
				std::cout << j << " ";
			}
			std::cout << "\n";
		}
		std::cout.flush();

		std::system("pause");	// DEBUG
	}
}