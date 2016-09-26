#ifndef ELEKTRON_EKALMAN_H
#define ELEKTRON_EKALMAN_H

#include "../math/Matrix.hpp"
#include "../math/Derivatives.hpp"

#include <cstddef>

namespace elektron
{

template <std::size_t S_, std::size_t M_>
class EKalman 
{
public:
	EKalman();

	void predict(CVec<double, S_>& x, 
				const CVec<double, S_>& u,
			   	const Matrix<double, S_, S_>& F, 
				const Matrix<double, S_, S_>& B,
			   	const Matrix<double, S_, S_>& P,
			   	const Matrix<double, M_, M_>& Q);

	void update(CVec<double, S_>& x,
		   		const CVec<double, M_>& z, 
				const CVec<double, M_>& v,
			   	const Matrix<double, S_, S_>& H,
			   	const Matrix<double, S_, S_>& P,
			   	const Matrix<double, M_, M_>& R);
	
	void update_jacobians(const Matrix<double, S_, S_>& F1,
		   		const Matrix<double, S_, S_>& F2,
			   	const Matrix<double, S_, S_>& H1,
			   	const Matrix<double, S_, S_>& H2,
			   	const RVec<double, S_>& deltas);

private:
	// Kalman gain
	Matrix<double, S_, M_> K;
	Matrix<double, M_, M_> I;
	Matrix<double, S_, S_> F_jacobian;
	Matrix<double, S_, S_> H_jacobian;
	double delta_t;
};

template <std::size_t S_, std::size_t M_>
EKalman<S_, M_>::EKalman() 
{
	I = Matrix<double, M_, M_>(i);
	K = Matrix<double, S_, M_>(r);
}

template <std::size_t  S_, std::size_t M_>
void EKalman<S_, M_>::update_jacobians(const Matrix<double, S_, S_>& F1, const Matrix<double, S_, S_>& F2, const Matrix<double, S_, S_>& H1, const Matrix<double, S_, S_>& H2, const RVec<double, S_>& deltas)
{
	F_jacobian = update_jacobians(F1, F2, deltas);
	H_jacobian = update_jacobians(H1, H2, deltas);
}	

// Matrix F - state transition model
// Matrix B - control model
// Matrix P - error estimate covariance
// Matrix Q - process noise covariance
// CVec x - state vector
// CVec u - control vector
template <std::size_t S_, std::size_t M_>
void EKalman<S_, M_>::predict(CVec<double, S_>& x, const CVec<double, S_>& u, const Matrix<double, S_, S_>& F, const Matrix<double, S_, S_>& B, const Matrix<double, S_, S_>& P, const Matrix<double, M_, M_>& Q)
{
	// Predict new state based on state transition model and applied control
	x = (F * x) + (B * u);
	// Error in state estimamte 
	P = (F_jacobian * P * F_jacobian.tpose()) + Q;
}

// Matrix H - observation model
// Matrix P - error estimate covariance
// Matrix R - observation noise covariance
// CVec x - state vector
// CVec z - observation vector
// CVec v - observation noise
template <std::size_t S_, std::size_t M_>
void EKalman<S_, M_>::update(CVec<double, S_>& x, const CVec<double, M_>& z, const CVec<double, M_>& v, const Matrix<double, S_, S_>& H, const Matrix<double, S_, S_>& P, const Matrix<double, M_, M_>& R)
{
	K = P * H_jacobian.tpose() * ((H * P * H_jacobian.tpose()) + R).inverse(); 
	x = x + K * (z - H_jacobian * x);
	P = (I - K * H_jacobian) * P; 
}

} // Eng of namespace elektron

#endif
