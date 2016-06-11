#ifndef MATRIX_H
#define	MATRIX_H
#endif 

#include <iostream>
#include <vector>

namespace nanos
{
	template<typename T>
	class Matrix
	{
	public:
		Matrix(std::vector<std::vector<T>> T_matrix);
		Matrix(T def_val, int r, int c);
		
		inline int rdim() { return r_dim; };
		inline int cdim() { return c_dim; };
		inline T elem(int r, int c) { return m[r][c]; };
		inline void set_elem(T val, int r, int c) { m[r][c] = val; };

		Matrix transpose();
		//Matrix LUD();
		//Matrix QRD();
		//Matrix eigenD();
		//Matrix choleskyD();

		void print();

	private:
		std::vector<std::vector<T>> m;

		int r_dim;
		int c_dim;
	};

}
