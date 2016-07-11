#include "IntegralImage.hpp"

namespace nlib 
{
	template <typename T>
	void compute_iimg(const Matrix<T>& m, const Matrix<T>& iimg_m)
	{
		// TODO: Add size checking: (m.r_dim == iimg_m.r_dim) &&
		//							(m.c_dim == iimg_m.c_dim)

		// Required initial values	
		iimg_m[0][0] = m[0][0];
		iimg_m[0][1] = m[0][0] + m[0][1];
		iimg_m[1][0] = m[0][1] + m[1][0];

		for (int i = 1; i < m.rdim(); i++)
		{
			for (int j = 1; j < m.cdim(); j++)
			{
				iimg_m[i][j] = m[i][j] + iimg_m[i-1][j] + iimg_m[i][j-1] -
					iimg_m[i-1][j-1];
			}
		}
	}
}
