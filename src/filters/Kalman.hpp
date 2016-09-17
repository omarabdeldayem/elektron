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
	Kalman();

	void predict(CVec<double, S_>& x, const CVec<double, S_>& u, const Matrix<double, S_, S_>& F, const Matrix<double, S_, S_>& B, const Matrix<double, S_, S_>& P, const Matrix<double, M_, M_>& Q);
	void update(CVec<double, S_>& x, const CVec<double, M_>& z, const CVec<double, M_>& v, const Matrix<double, S_, S_>& H, const Matrix<double, S_, S_>& P, const Matrix<double, M_, M_>& R);

private:
	// Kalman gain
	Matrix<double, S_, M_> K;
	Matrix<double, M_, M_> I;
};

template <std::size_t S_, std::size_t M_>
Kalman<S_, M_>::Kalman() 
{
	I = Matrix<double, M_, M_>(i);
	K = Matrix<double, S_, M_>();
}

// Matrix F - state transition model
// Matrix B - control model
// Matrix P - error estimate covariance
// Matrix Q - process noise covariance
// CVec x - state vector
// CVec u - control vector
template <std::size_t S_, std::size_t M_>
void Kalman<S_, M_>::predict(CVec<double, S_>& x, const CVec<double, S_>& u, const Matrix<double, S_, S_>& F, const Matrix<double, S_, S_>& B, const Matrix<double, S_, S_>& P, const Matrix<double, M_, M_>& Q)
{
	// Predict new state based on state transition model and applied control
	x = (F * x) + (B * u);
	// Error in state estimamte 
	P = (F * P * F.tpose()) + Q;
}

// Matrix H - observation model
// Matrix P - error estimate covariance
// Matrix R - observation noise covariance
// CVec x - state vector
// CVec z - observation vector
// CVec v - observation noise
template <std::size_t S_, std::size_t M_>
void Kalman<S_, M_>::update(CVec<double, S_>& x, const CVec<double, M_>& z, const CVec<double, M_>& v, const Matrix<double, S_, S_>& H, const Matrix<double, S_, S_>& P, const Matrix<double, M_, M_>& R)
{
	K = P * H.tpose() * ((H * P * H.tpose()) + R).inverse(); 
	x = x + K * (z - H * x);
	P = (I - K * H) * P; 
}

} // End of namespace filters	
} // Eng of namespace elektron

#endif
