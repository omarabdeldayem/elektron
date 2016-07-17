#ifndef MATRIX_H_
#define	MATRIX_H_

#include <stdexcept>
#include <iostream>
#include <vector>
#include <numeric>
#include <math.h>

namespace nlib
{

template<typename T>
class Matrix
{
public:
	// CONSTRUCTORS
	Matrix(int r, int c);

	Matrix(std::vector<std::vector<T>>& T_matrix)
		: mat_(T_matrix)
		, rows_(T_matrix.size())
		, cols_(T_matrix[0].size())
	{ }
	
	Matrix(T def_val, int r, int c)
		: rows_(r)
		, cols_(c)
	{ mat_.resize(rows_ * cols_, def_val); }
	
	// OPERATOR OVERLOADS
	Matrix<T> operator*(Matrix<T> a);
	T& operator()(int r, int c) { return mat_.at(r * rows_ + c); };
	
	// PRIVATE MEMBER ACCESS METHODS
	inline int rdim() { return rows_; };
	inline int cdim() { return cols_; };

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
	int rows_;
	int cols_;

	// Default store mat_rix in row-major form
	std::vector<T> mat_;
	typename std::vector<T>::iterator iter_;

	Matrix<T> pivot();
};

}

#include "Matrix.cpp"
#endif 
