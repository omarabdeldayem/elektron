#include "Matrix.hpp"

int main()
{
	// TEST 
	nanos::Matrix<int> m_int = nanos::Matrix<int>(8, 3, 4);
	m_int.print();
	nanos::Matrix<int> m_int_T = m_int.tpose();
	m_int_T.print();
	
	m_int[2][0] = 4;
	m_int.print();

	return 0;
}