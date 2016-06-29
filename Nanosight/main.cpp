#include "nblas/Matrix.hpp"

int main()
{
	// TEST 
	nanos::Matrix<double> m_int = nanos::Matrix<double>(8, 3, 3);
	nanos::Matrix<double> m_int_l = nanos::Matrix<double>(0, 3, 3),
					   m_int_u = nanos::Matrix<double>(0, 3, 3);

	m_int.LUD(m_int_l, m_int_u);
	m_int.print();
	m_int_l.print();
	m_int_u.print();

	nanos::Matrix<double> m_int_T = m_int.tpose();
	m_int_T.print();

	return 0;
}