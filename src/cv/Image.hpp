#ifndef ELEKTRON_IMAGE_H_
#define ELEKTRON_IMAGE_H_


#include "../math/Matrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace elektron
{

const double HSV_H_THRESHOLD = 0.0001;
const double RGB2GRAY_RCOEFF = 0.299;
const double RGB2GRAY_GCOEFF = 0.587;
const double RGB2GRAY_BCOEFF = 0.114;

enum ImType {RGB, HSV, HSL, HSI, GRAY};

template <typename T_, std::size_t R_, std::size_t C_>
class Image
{
public:
	// CONSTRUCTORS
	Image();
	Image(const std::array<T_, R_*C_>& ch1, const std::array<T_, R_*C_>& ch2, const std::array<T_, R_*C_>& ch3, ImType t);
	Image(const Matrix<T_, R_, C_>& ch1, const Matrix<T_, R_, C_>& ch2, const Matrix<T_, R_, C_>& ch3, ImType t);

	// CHANNEL CONVERSION INTERFACE
	void to_RGB();
	void to_HSV();
	void to_HSL();
	void to_HSI();
	void to_GRAY();

	// IMAGE UTILITIES
	void iimg(Image<T_, R_, C_>& iimg_m);

	// PRIVATE MEMBER ACCESS METHODS
	Matrix<T_, R_, C_>& ch1() { return ch1_; };
	Matrix<T_, R_, C_>& ch2() { return ch2_; };
	Matrix<T_, R_, C_>& ch3() { return ch3_; };

	ImType type() { return im_t_; };
	inline int r() { return rows_; };
	inline int c() { return cols_; };

private:
	//  sizeof(T_)-bit IMAGE CHANNELS
	Matrix<T_, R_, C_> ch1_;
	Matrix<T_, R_, C_> ch2_;
	Matrix<T_, R_, C_> ch3_;

	// IMAGE TYPE
	ImType im_t_;

	int rows_;
	int cols_;

	// CHANNEL CONVERSIONS
	void hsl2rgb();
	void hsi2rgb();
	void hsv2rgb();
	void rgb2hsv();
	void rgb2hsi();
	void rgb2hsl();
	void rgb2gray();

	void sumtable(const Matrix<T_, R_, C_>& m, Matrix<T_, R_, C_>& sum);
};

template <typename T_, std::size_t R_, std::size_t C_>
Image<T_, R_, C_>::Image()
{
	ch1_ = Matrix<T_, R_, C_>(z);
	ch2_ = Matrix<T_, R_, C_>(z);
	ch3_ = Matrix<T_, R_, C_>(z);

	im_t_ = RGB;
}

template <typename T_, std::size_t R_, std::size_t C_>
Image<T_, R_, C_>::Image(const std::array<T_, R_*C_>& ch1, const std::array<T_, R_*C_>& ch2, const std::array<T_, R_*C_>& ch3, ImType t)
	: im_t_(t)
{
	ch1_ = Matrix<T_, R_, C_>(ch1);
	ch2_ = Matrix<T_, R_, C_>(ch2);
	ch3_ = Matrix<T_, R_, C_>(ch3);
}

template <typename T_, std::size_t R_, std::size_t C_>
Image<T_, R_, C_>::Image(const Matrix<T_, R_, C_>& ch1, const Matrix<T_, R_, C_>& ch2, const Matrix<T_, R_, C_>& ch3, ImType t)
	: im_t_(t),
	  ch1_(ch1),
	  ch2_(ch2),
	  ch3_(ch3)
{ }

template <typename T_, std::size_t R_, std::size_t C_>
void Image<T_, R_, C_>::to_HSV()
{
	switch(im_t_)
	{
		case RGB:
			rgb2hsv();
			break;
		case HSL:
			hsl2rgb();
			rgb2hsv();
			break;
		case HSI:
			hsi2rgb();
			rgb2hsv();
			break;
		case GRAY:
			break;
	}
}

// Computes the integral image (summed-area table) of image
template <typename T_, std::size_t R_, std::size_t C_>
void Image<T_, R_, C_>::iimg(Image<T_, R_, C_>& iimg_m)
{
	if (im_t_ == GRAY)
	{
		sumtable(ch1_, iimg_m.ch1());
	}
	else
	{
		sumtable(ch1_, iimg_m.ch1());
		sumtable(ch2_, iimg_m.ch2());
		sumtable(ch3_, iimg_m.ch3());
	}
}

// NTS: This doesn't belong here - move
template <typename T_, std::size_t R_, std::size_t C_>
void Image<T_, R_, C_>::sumtable(const Matrix<T_, R_, C_>& m, Matrix<T_, R_, C_>& sum)
{
	// Required initial values
	sum(0, 0) = m(0, 0);
	sum(0, 0) = m(0, 0) + m(0, 1);
	sum(0, 0) = m(0, 0) + m(1, 0);
	for (int i = 1; i < m.r(); i++)
	{
		for (int j = 1; j < m.c(); j++)
		{
			sum(i, j) = m(i, j) + sum(i-1, j) + sum(i, j-1) -
				sum(i-1, j-1);
		}
	}
}


template <typename T_, std::size_t R_, std::size_t C_>
void Image<T_, R_, C_>::rgb2hsv()
{
	double M;
	double m;
	double C;

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			// Max channel value
			M = std::max(std::max(ch1_(i, j), ch2_(i, j)), ch3_(i, j));
			// Min channel value
			m = std::min(std::min(ch1_(i, j), ch2_(i, j)), ch3_(i, j));
			// Delta
			C = M - m;

			if (C < HSV_H_THRESHOLD || M < 0.0)
			{
				ch1_(i, j) = static_cast<T_>(0);
				ch2_(i, j) = static_cast<T_>(0);
			}
			else
			{
				if (ch1_(i, j) >= M)
				{
					ch1_(i, j) = 60 * ((ch2_(i, j) - ch3_(i, j)) / C);
				}
				else if (ch2_(i, j) >= M)
				{
					ch1_(i, j) = 60 * (2.0 + (ch3_(i, j) - ch1_(i, j)) / C);
				}
				else
				{
					ch1_(i, j) = 60 * (4.0 + (ch1_(i, j) - ch2_(i, j)) / C);
				}

				if (ch1_(i, j) < 0.0) ch1_(i, j) += 360.0;

				ch2_(i, j) = static_cast<T_>(C / M);
			}

			// Set V channel to Max
			ch3_(i, j) = M;
		}
	}
}

template <typename T_, std::size_t R_, std::size_t C_>
void Image<T_, R_, C_>::hsv2rgb()
{
	double C;
	double H;
	double X;
	double M;

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			C = ch2_(i, j) * ch3_(i, j);
			H = std::fmod(ch1_(i, j) / 60.0, 6);
			X = C * (1 - std::fabs(std::fmod(H, 2) - 1));
			M = ch3_(i, j) - C;

			if (0.0 <= H && H < 1.0)
			{
				ch1_(i, j) = C + M;
				ch2_(i, j) = X + M;
				ch3_(i, j) = 0;
			}
			else if (1.0 <= H && H < 2.0)
			{
				ch1_(i, j) = X + M;
				ch2_(i, j) = C + M;
				ch3_(i, j) = 0;
			}
			else if (2.0 <= H && H < 3.0)
			{
				ch1_(i, j) = 0;
				ch2_(i, j) = C + M;
				ch3_(i, j) = X + M;
			}
			else if (3.0 <= H && H < 4.0)
			{
				ch1_(i, j) = 0;
				ch2_(i, j) = X + M;
				ch3_(i, j) = C + M;
			}
			else if (4.0 <= H && H < 5.0)
			{
				ch1_(i, j) = X + M;
				ch2_(i, j) = 0;
				ch3_(i, j) = C + M;
			}
			else if (5.0 <= H && H < 6.0)
			{
				ch1_(i, j) = C + M;
				ch2_(i, j) = 0;
				ch3_(i, j) = X + M;
			}
			else
			{
				ch1_(i, j) = M;
				ch2_(i, j) = M;
				ch3_(i, j) = M;
			}
		}
	}
}

template <typename T_, std::size_t R_, std::size_t C_>
void Image<T_, R_, C_>::rgb2hsl()
{
	double M;
	double m;
	double C;

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			// Max channel value
			M = std::max(std::max(ch1_(i, j), ch2_(i, j)), ch3_(i, j));
			// Min channel value
			m = std::min(std::min(ch1_(i, j), ch2_(i, j)), ch3_(i, j));
			// Delta
			C = M - m;

			if (C < HSV_H_THRESHOLD || M < 0.0)
			{
				ch1_(i, j) = static_cast<T_>(0);
				ch2_(i, j) = static_cast<T_>(0);
			}
			else
			{
				if (ch1_(i, j) >= M)
				{
					ch1_(i, j) = 60 * ((ch2_(i, j) - ch3_(i, j)) / C);
				}
				else if (ch2_(i, j) >= M)
				{
					ch1_(i, j) = 60 * (2.0 + (ch3_(i, j) - ch1_(i, j)) / C);
				}
				else
				{
					ch1_(i, j) = 60 * (4.0 + (ch1_(i, j) - ch2_(i, j)) / C);
				}

				if (ch1_(i, j) < 0.0) ch1_(i, j) += 360.0;

				ch2_(i, j) = static_cast<T_>(C / M);
			}

			// Set V channel to Max
			ch3_(i, j) = static_cast<T_>((0.5 * (M + m)));
		}
	}
}

template <typename T_, std::size_t R_, std::size_t C_>
void Image<T_, R_, C_>::hsl2rgb()
{
	double C;
	double H;
	double X;
	double M;

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			C = (1 - std::fabs(2.0  * ch3_(i, j) - 1)) * ch2_(i, j);
			H = std::fmod(ch1_(i, j) / 60.0, 6);
			X = C * (1 - std::fabs(std::fmod(H, 2) - 1));
			M = ch3_(i, j) - C;

			if (0.0 <= H && H < 1.0)
			{
				ch1_(i, j) = C + M;
				ch2_(i, j) = X + M;
				ch3_(i, j) = 0;
			}
			else if (1.0 <= H && H < 2.0)
			{
				ch1_(i, j) = X + M;
				ch2_(i, j) = C + M;
				ch3_(i, j) = 0;
			}
			else if (2.0 <= H && H < 3.0)
			{
				ch1_(i, j) = 0;
				ch2_(i, j) = C + M;
				ch3_(i, j) = X + M;
			}
			else if (3.0 <= H && H < 4.0)
			{
				ch1_(i, j) = 0;
				ch2_(i, j) = X + M;
				ch3_(i, j) = C + M;
			}
			else if (4.0 <= H && H < 5.0)
			{
				ch1_(i, j) = X + M;
				ch2_(i, j) = 0;
				ch3_(i, j) = C + M;
			}
			else if (5.0 <= H && H < 6.0)
			{
				ch1_(i, j) = C + M;
				ch2_(i, j) = 0;
				ch3_(i, j) = X + M;
			}
			else
			{
				ch1_(i, j) = M;
				ch2_(i, j) = M;
				ch3_(i, j) = M;
			}
		}
	}
}


template <typename T_, std::size_t R_, std::size_t C_>
void Image<T_, R_, C_>::rgb2hsi()
{
	double M;
	double m;
	double C;

	T_ r;
	T_ g;
	T_ b;

	for (int i = 0; i < rows_; i++)
	{
		for (int j = 0; j < cols_; j++)
		{
			r = ch1_(i, j);
			g = ch2_(i, j);
			b = ch3_(i, j);

			// Max channel value
			M = std::max(std::max(r, g), b);
			// Min channel value
			m = std::min(std::min(r, g), b);
			// Delta
			C = M - m;

			if (C < HSV_H_THRESHOLD || M < 0.0)
			{
				ch1_(i, j) = static_cast<T_>(0);
				ch2_(i, j) = static_cast<T_>(0);
			}
			else
			{
				if (ch1_(i, j) >= M)
				{
					ch1_(i, j) = 60 * ((g - b) / C);
				}
				else if (ch2_(i, j) >= M)
				{
					ch1_(i, j) = 60 * (2.0 + (b - r) / C);
				}
				else
				{
					ch1_(i, j) = 60 * (4.0 + (r - b) / C);
				}

				if (ch1_(i, j) < 0.0) ch1_(i, j) += 360.0;

				ch2_(i, j) = static_cast<T_>(C / M);
			}

			ch3_(i, j) = static_cast<T_>((1/3) * (r + g + b));
		}
	}
}

// Fills channel 1 with the grayscale values
// Leaves channels 2 and 3 unchanged
template <typename T_, std::size_t R_, std::size_t C_>
void Image<T_, R_, C_>::rgb2gray()
{
	ch1_ = (RGB2GRAY_RCOEFF * ch1_) + (RGB2GRAY_GCOEFF * ch2_) + (RGB2GRAY_BCOEFF * ch3_);
}

} // End of namespace elektron

#endif
