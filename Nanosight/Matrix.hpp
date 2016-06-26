#ifndef MATRIX_HPP
#define	MATRIX_HPP

#include <stdexcept>
#include <iostream>
#include <vector>
#include <math.h>

namespace nanos
{
	template<typename T>
	class Matrix
	{
	public:
		Matrix(std::vector<std::vector<T>>& T_matrix);
		Matrix(T def_val, int r, int c);
		
		std::vector<T>& operator[](int i);

		inline int rdim() { return r_dim; };
		inline int cdim() { return c_dim; };

		T tr();
		
		Matrix tpose();
		
		//Matrix LUD();
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
	};

}
#endif 