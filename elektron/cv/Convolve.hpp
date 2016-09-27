#ifndef ELEKTRON_CONVOLVE_H_
#define ELEKTRON_CONVOLVE_H_

#include "../math/Matrix.hpp"

#include <cstddef>

namespace elektron
{

template <std::size_t R_, std::size_t C_>
struct Kernel
{
    Matrix<int, R_, C_> m;
    double prescale;
};


// Overwrites input image
// No border processing for now
template <typename T_, std::size_t R_, std::size_t C_>
void conv(Image<T_, R_, C_>& im, const Kernel<3, 3>& kern)
{
    for (int i = 1; i < im.r() - 1; i++)
    {
        for (int j = 1; j < im.c() - 1; i++)
        {
            T_ ch1_accum = 0;
            T_ ch2_accum = 0;
            T_ ch3_accum = 0;

            for (int k = 0; k < kern.m.r(); k++)
            {
                for (int l = 0; l < kern.m.c(); l++)
                {
                    ch1_accum += kern.prescale * kern.m(k, l) * im.ch1()(i - k - 1, j - l - 1);
                    ch2_accum += kern.prescale * kern.m(k, l) * im.ch2()(i - k - 1, j - l - 1);
                    ch3_accum += kern.prescale * kern.m(k, l) * im.ch3()(i - k - 1, j - l - 1);
                }
            }

            im.ch1()(i, j) = ch1_accum;
            im.ch2()(i, j) = ch2_accum;
            im.ch3()(i, j) = ch3_accum;
        }
    }
}

// Convolved image is put into out
// No border processing for now
template <typename T_, std::size_t R_, std::size_t C_>
void conv(const Image<T_, R_, C_>& im, const Kernel<3, 3>& kern, Image<T_, R_, C_>& out)
{
    for (int i = 1; i < im.r() - 1; i++)
    {
        for (int j = 1; j < im.c() - 1; i++)
        {
            T_ ch1_accum = 0;
            T_ ch2_accum = 0;
            T_ ch3_accum = 0;

            for (int k = 0; k < kern.m.r(); k++)
            {
                for (int l = 0; l < kern.m.c(); l++)
                {
                    ch1_accum += kern.prescale * kern.m(k, l) * im.ch1()(i - k - 1, j - l - 1);
                    ch2_accum += kern.prescale * kern.m(k, l) * im.ch2()(i - k - 1, j - l - 1);
                    ch3_accum += kern.prescale * kern.m(k, l) * im.ch3()(i - k - 1, j - l - 1);
                }
            }

            out.ch1()(i, j) = ch1_accum;
            out.ch2()(i, j) = ch2_accum;
            out.ch3()(i, j) = ch3_accum;
        }
    }
}

}

#endif
