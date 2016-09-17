#ifndef ELEKTRON_KALMAN_H
#define ELEKTRON_KALMAN_H

#include <cstddef>

namespace elektron
{
namespace filters
{

template <std::size_t S_, std::size_t M_>
class Kalman 
{
public:
	Kalman(double dt) : delta_t(dt) { };

	void predict(CVec<double, S_>& x, const Matrix<double, S_, S_>& F, const Matrix<double, S_, S_>& B, const Matrix<double, S_, S_>& P, const Matrix<double, M_, M_>& Q);
	void update();
private:
	// Kalman gain
	Matrix<double, S_, M_> K;

	// Elapsed time
	double delta_t;
};

// Matrix F - state transition model
// Matrix B - control model
// Matrix P - error covariance matrix
// Matrix Q - process noise covariance
template <std::size_t S_, std::size_t M_>
void Kalman<S_, M_>::predict(CVec<double, S_>& x, const Matrix<double, S_, S_>& F, const Matrix<double, S_, S_>& B, const Matrix<double, S_, S_>& P, const Matrix<double, M_, M_>& Q)
{
	x = (F * x) + B;
	P = (F * P * F.tpose()) + Q;
}

template <std::size_t S_, std::size_t M_>
void Kalman<S_, M_>::update()
{
}

}	
}

#endif
