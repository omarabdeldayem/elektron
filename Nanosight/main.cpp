#include "Matrix.hpp"

int main()
{
	// TEST 
	nanos::Matrix<int> m_int = nanos::Matrix<int>(8, 3, 4);
	m_int.print();
	nanos::Matrix<int> m_int_T = m_int.transpose();
	m_int_T.print();

	return 0;
}