#include "blas/Matrix.hpp"
#include "blas/Solver.hpp"

#include <iostream>


int main()
{
	// Test
	std::cout << "Scalar Operations: \n" << "----------\n";
	nlib::Matrix<int, 4, 4> x(nlib::o);
	x.print();
	std::cout << "Multiply by 5\n";
	x = x * 5;
	x.print();
	std::cout << "L_(2, 1) norm: " << x.norm() << "\n" << "\n"; 
	std::cout << "Divide by 5\n";
	x = x / 5;
	x.print();


	std::cout << "Matrix Operations: \n" << "----------\n";
	nlib::Matrix<int, 6, 4> y;
	nlib::Matrix<int, 4, 6> yt;
	y.ones();
	y(0, 2) = 3;
	y(2, 3) = 5;
	y.print();
	std::cout << "Tranposed: \n";
	yt = y.tpose();
	yt.print();
	std::cout << "Multiplied: \n";
	nlib::Matrix<int, 6, 6> prod = y * yt;
	prod.print();
	std::cout << "Sub-matrix: \n";
	nlib::Matrix<int, 1, 4> y_sub;
	y.sub(y_sub, 0, 1, 0, 4);
	y_sub.print();	
	
	std::cout << "LUDecomposition Start: \n" << "----------\n";
	nlib::Matrix<double, 3, 3> m;
	nlib::Matrix<double, 3, 3> m_l;
	nlib::Matrix<double, 3, 3> m_u;
	m.ones();
//	m_l.ones();
//	m_u.ones();
	m = m * 8;
	
	m(0, 0) = 3;
	m(2, 1) = 5;
	m.lud(m_l, m_u);
	std::cout << "Matrix M: \n";
	m.print();
	std::cout << "Matrix L: \n";
	m_l.print();
	std::cout << "Matrix U: \n";
	m_u.print();
	std::cout << "Matrix L * U: \n";
	(m_l * m_u).print();

	std::cout << "QRDecomposition Start: \n" << "----------\n";
	nlib::Matrix<double, 3, 3> A;
	nlib::Matrix<double, 3, 3> Q;
   	nlib::Matrix<double, 3, 3> R;
	A(0, 0) = 12; 
	A(0, 1) = -51;
	A(0, 2) = 4;
	A(1, 0) = 6;
	A(1, 1) = 167;
	A(1, 2) = -68;
	A(2, 0) = -4;
	A(2, 1) = 24;
	A(2, 2) = -41;
	std::cout << "Matrix A: \n";
	A.print();
	A.qrd(Q, R);
	std::cout << "Matrix Q: \n";
	Q.print();
	std::cout << "Matrix R: \n";
	R.print();
	std::cout << "Matrix Q * R: \n";	
	(Q * R).print(); 


	std::cout << "SVDecomposition Start: \n" << "----------\n";
	nlib::Matrix<double, 3, 2> a_svd;
	nlib::Matrix<double, 3, 2> a_u;
	nlib::Matrix<double, 2, 2> a_s;
	nlib::Matrix<double, 2, 2> a_vt;
	
	a_svd.ones();
	a_svd = a_svd * 2;
	a_svd(0, 0) = 1;
	a_svd(2, 1) = 1;

	std::cout << "Matrix A: \n";
	a_svd.print();
	a_svd.svd(a_u, a_s, a_vt);
	
	std::cout << "Matrix U: \n";
	a_u.print();
	std::cout << "Matrix S: \n";
	a_s.print();
	std::cout << "Matrix V_T: \n";
	a_vt.print();
	std::cout << "U * S * V_T: \n";
	(a_u * a_s * a_vt).print();
	return 0;
}
