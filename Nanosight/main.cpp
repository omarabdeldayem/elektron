#include "nblas/Matrix.hpp"

int main()
{
	// TEST 
	nanos::Matrix<int> m_int = nanos::Matrix<int>(8, 3, 3);
	nanos::Matrix<int> m_int_l = nanos::Matrix<int>(0, 3, 3), 
					   m_int_u = nanos::Matrix<int>(0, 3, 3);

	nanos::Matrix<int> mult_res = m_int * m_int_l;
	mult_res.print();

	m_int.LUD(m_int_l, m_int_u);
	m_int.print();
	m_int_l.print();
	m_int_u.print();

	nanos::Matrix<int> m_int_T = m_int.tpose();
	m_int_T.print();
	
	m_int[2][0] = 4;
	m_int.print();

	return 0;
}