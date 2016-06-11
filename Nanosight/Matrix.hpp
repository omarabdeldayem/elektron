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
		Matrix();
		Matrix(std::vector<std::vector<T>> T_matrix);
		Matrix(T def_val, int r, int c);
		Matrix(int r, int c);
		


		inline int rdim();
		inline int cdim();
		inline T elem(int r, int c);

		//Matrix transpose();
		//Matrix ctranspose();
		//Matrix LUD();
		//Matrix QRD();
		//Matrix EigenD();
		//Matrix CholeskyD();

		void print();

	private:
		std::vector<std::vector<T>> m;

		int r_dim;
		int c_dim;
	};

}
