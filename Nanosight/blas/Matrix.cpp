#include "Matrix.hpp"

namespace nanos
{
	template class Matrix<uint_fast8_t>;
	template class Matrix<int>;
	template class Matrix<float>;
	template class Matrix<double>;
	template class Matrix<long>;

	template <typename T>
	Matrix<T>::Matrix()
	{
		r_dim = 0;
		c_dim = 0;
	}

	// Creates r x c identity matrix
	template <typename T>
	Matrix<T>::Matrix(int r, int c)
	{
		for (int i = 0; i < r; i++)
		{
			std::vector<T> t;
			m.push_back(t);

			for (int j = 0; j < c; j++)
			{
				if (i == j)
				{
					m[i].push_back(1);
				}
				else
				{
					m[i].push_back(0);
				}
			}
		}

		r_dim = r;
		c_dim = c;
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

			for (auto& mat: m)
			{
				mat.push_back(def_val);
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
				m_T[i][j] = m[j][i];
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
	Matrix<T> Matrix<T>::pivot()
	{
		if (is_sqr)
		{
			Matrix<T> id = Matrix<T>(r_dim, c_dim);
			for (int i = 0; i < r_dim; i++)
			{
				T maxv = m[i][i];
				int row = i;

				for (int j = i; j < r_dim; j++)
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
	void Matrix<T>::luD(Matrix<T>& l, Matrix<T>& u, Matrix<T>&p)
	{
		p = p.pivot();
		Matrix<T> m2 = p * m;

		for (int j = 0; j < r_dim; j++)
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
			for (int i = j; i < r_dim; i++)
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