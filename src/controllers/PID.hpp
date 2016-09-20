#ifndef ELEKTRON_PCTRL_H_
#define ELKETRON_PCTRL_H_

#include "../math/CVec.hpp"

#include <cstddef>

namespace elektron
{

template <typename T_, std::size_t R_>
class PID
{
public:
	PID(const CVec<T_, R_>& s);
	PID(const CVec<T_, R_>& s, const CVec<T_, R_>& k_p, const CVec<T_, R_>& k_i, const CVec<T_, R_>& k_d);

	void ctrl(CVec<T_, R_>& U, const CVec<T_, R_>& x) const;

private:
	CVec<T_, R_> set_point_;
	CVec<T_, R_> K_p_;
	CVec<T_, R_> K_i_;
	CVec<T_, R_> K_d_;
}


template <typename T_, std::size_t R_>
PID<T_, R_>::PID(const CVec<T_, R_>& s) : set_point_(s)
{
	K_p_ = CVec<T_, R_>(i);
	K_i_ = CVec<T_, R_>(i);
	K_d_ = CVec<T_, R_>(i);
}

template <typename T_, std::size_t R_>
PID<T_, R_>::PID(const CVec<T_, R_>& s, const CVec<T_, R_>& k_p, const CVec<T_, R_>& k_i, const CVec<T_, R_>& k_d)
	: set_point_(s), 
	  K_p_(k_p),
	  K_i_(k_i),
	  K_d_(k_d)
{ }


template <typename T_, std::size_t R_>
void PID<T_, R_>::ctrl(CVec<T_, R_>& U, const CVec<T_, R_>& x1, const CVec<T_, R_>& x2, const double dt) const
{
	CVec<T_, R_> err = set_point_ - x;
	CVec<T_, R_> err_i = 0;
	CVec<T_, R_> err_d = 0;
	
	U = (K_p_.tpose() * err) + (K_i_.tpose() * err_i) + (K_d_.tpose() * err_d);
}

} // End of namespace elektron
#endif
