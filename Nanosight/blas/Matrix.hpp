#ifndef MATRIX_HPP
#define	MATRIX_HPP

#include <stdexcept>
#include <iostream>
#include <vector>
#include <numeric>
#include <math.h>

namespace nanos
{
	template<typename T>
	class Matrix
	{
	public:
		// CONSTRUCTORS
		Matrix();
		Matrix(int r, int c);
		Matrix(std::vector<std::vector<T>>& T_matrix);
		Matrix(T def_val, int r, int c);
		
		// OPERATOR OVERLOADS
		std::vector<T>& operator[](int i);
		Matrix<T> operator*(Matrix<T> a);
		void operator=(std::vector<std::vector<T>>& mat);

		// PRIVATE MEMBER ACCESS METHODS
		inline int rdim() { return r_dim; };
		inline int cdim() { return c_dim; };
		inline std::vector<std::vector<T>> data() { return m; };

		// MATRIX OPERATIONS
		T tr();
		Matrix<T> tpose();
		void luD(Matrix<T>& l, Matrix<T>& u, Matrix<T>& p);
		
		//Matrix QRD();
		//Matrix eigenD();
		//Matrix choleskyD();
		//int det();
		void print();

	private:
		int r_dim;
		int c_dim;
		bool is_sqr = r_dim == c_dim ? true : false;

		std::vector<std::vector<T>> m;		

		Matrix<T> pivot();
	};

}

#endif 