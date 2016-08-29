#include "blas/Matrix.hpp"
#include <iostream>

int main()
{
	// TEST 
	nlib::Matrix<double> m = nlib::Matrix<double>(8.0, 3, 3);
	nlib::Matrix<double> m_l = nlib::Matrix<double>(0.0, 3, 3),
						m_u = nlib::Matrix<double>(0.0, 3, 3),
						m_p = nlib::Matrix<double>(3, 3);

	m(0, 0) = 3;
	m(2, 1) = 5;

	m.lud(m_l, m_u, m_p);
	std::cout << "LUDecomposition Start: \n" << "----------\n";
	std::cout << "Matrix M: \n";
	m.print();
	std::cout << "Matrix L: \n";
	m_l.print();
	std::cout << "Matrix U: \n";
	m_u.print();
	std::cout << "Matrix P: \n";
	m_p.print();
	std::cout << "Matrix L * U: \n";
	(m_l * m_u).print();
	
	std::cout << "Scalar Operations: \n" << "----------\n";
	nlib::Matrix<int> x = nlib::Matrix<int>(3, 3);
	x.print();
	std::cout << "Multiply by 5\n";
	x = x * 5;
	x.print();
	std::cout << "Divide by 5\n";
	x = x / 5;
	x.print();

	std::cout << "Matrix Operations: \n" << "----------\n";
	nlib::Matrix<int> y = nlib::Matrix<int>(9, 4, 6);
	y(0, 2) = 3;
	y(2, 3) = 5;
	y.print();
	std::cout << "Tranposed: \n";
	y = y.tpose();
	y.print();	
//	nlib::Matrix<double> m_T = m.tpose();
//	m_T.print();

	return 0;
}
