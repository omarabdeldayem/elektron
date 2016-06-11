#include "Matrix.hpp"

template class nanos::Matrix<int>;
template class nanos::Matrix<float>;
template class nanos::Matrix<long>;

template<typename T>
nanos::Matrix<T>::Matrix()
{

}

template <typename T> 
nanos::Matrix<T>::Matrix(std::vector<std::vector<T>> T_matrix)
{
	m = T_matrix;
	r_dim = T_matrix.size();
	c_dim = T_matrix[0].size();
}

template <typename T>
nanos::Matrix<T>::Matrix(T def_val, int r, int c)
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
nanos::Matrix<T>::Matrix(int r, int c)
{ 
}

template <typename T>
inline int nanos::Matrix<T>::rdim()
{
	return r_dim;
}

template <typename T>
inline int nanos::Matrix<T>::cdim()
{
	return c_dim;
}

template <typename T>
inline T nanos::Matrix<T>::elem(int r, int c)
{
	return m[r][c];
}

template <typename T>
void nanos::Matrix<T>::print()
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