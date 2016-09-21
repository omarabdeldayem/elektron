#ifndef ELEKTRON_IMAGE_H_
#define ELEKTRON_IMAGE_H_


#include "../math/Matrix.hpp"

#include <cstddef>

namespace elektron
{

enum ImType = {RGB, HSV, HSL, HSI, BW};

tempalte <typename T_, std::size_t R_, std::size_t C_>
class Image : public Matrix<T_, R_, C_>
{
public:
	Image();
	Image(const std::array<T_, R_*S_>& ch1, const std::array<T_, R_*S_>& ch2, const std::array<T_, R_*S_>& ch3, ImType t);
	Image(const Matrix<T_, R_, S_>& ch1, const Matrix<T_, R_, S_>& ch2, const Matrix<T_, R_, S_>& ch3, ImType t);
	
	void to_RGB();
	void to_HSV();
	void to_HSL();
	void to_HIS();
	void to_BW();
	
	ImType type();

	inline int r() { return rows_ };
	inline int c() { return cols_ };

private:	
	Matrix<T_, R_, S_> ch1_;
	Matrix<T_, R_, S_> ch2_;
	Matrix<T_, R_, S_> ch3_;

	ImType im_t_;

	int rows_;
	int cols_;

}

template <typename T_, std::size_t R_, std::size_t C_>
Image<T_, R_, C_>::Image()
{
	ch1_ = Matrix<T_, R_, C_>(z);
	ch2_ = Matrix<T_, R_, C_>(z);
	ch3_ = Matrix<T_, R_, C_>(z);

	im_t_ = RGB;
}

tempalte <typename T_, std::size_t R_, std::size_t C_>
Image<T_, R_, C_>::Image(const std::array<T_, R_*S_>& ch1, const std::array<T_, R_*C_>& ch2, const std::array<T_, R_*C_>& ch3, ImType t)
	: im_t_(t)
{

}

template <typename T_, std::size_t R_, std::size_t C_>
Image<T_, R_, C_>::Image(const Matrix<T_, R_, S_>& ch1, const Matrix<T_, R_, S_>& ch2, const Matrix<T_, R_, S_>& ch3, ImType t)
	: im_t_(t),
	  ch1_(ch1),
	  ch2_(ch2),
	  ch3_(ch3)
{

}

}

#endif
