#include "Matrix.hpp"

namespace nanos
{
	template class Matrix<int>;
	template class Matrix<float>;
	template class Matrix<long>;

	template <typename T>
	Matrix<T>::Matrix(std::vector<std::vector<T>> T_matrix)
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
	Matrix<T> Matrix<T>::transpose()
	{
		Matrix<T> m_T = Matrix<T>(NULL, c_dim, r_dim);

		for (int i = 0; i < r_dim; i++)
		{
			for (int j = 0; j < c_dim; j++)
			{
				m_T.set_elem(m[i][j], j, i);
			}
		}

		return m_T;
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