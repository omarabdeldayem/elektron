#include "blas/Matrix.hpp"

int main()
{
	// TEST 
	nlib::Matrix<double> m = nlib::Matrix<double>(8.0, 3, 3);
	nlib::Matrix<double> m_l = nlib::Matrix<double>(0.0, 3, 3),
						m_u = nlib::Matrix<double>(0.0, 3, 3),
						m_p = nlib::Matrix<double>(3, 3);

	m[0][0] = 1;
	m[0][1] = 2;
	m[0][2] = 4;
	m[1][0] = 3;
	m[1][1] = 8;
	m[1][2] = 14;
	m[2][0] = 2;
	m[2][1] = 6;
	m[2][2] = 13;

//	m.luD(m_l, m_u, m_p);
	m.print();
	m_l.print();
	m_u.print();
	m_p.print();

//	nlib::Matrix<double> m_T = m.tpose();
//	m_T.print();

	return 0;
}
