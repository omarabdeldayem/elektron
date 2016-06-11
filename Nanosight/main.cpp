#include "Matrix.hpp"

int main()
{
	// TEST 
	nanos::Matrix<int> m_int = nanos::Matrix<int>(8, 3, 4);
	m_int.print();

	return 0;
}