#ifndef ELEKTRON_INTEGRALS_H_
#define ELEKTRON_INTEGRALS_H_

namespace elektron
{

template <typename T_>
inline double integral_rect(T_ f, double delta)
{
	return static_cast<double>(f * delta);
}

template <typename T_>
inline double integral_trapez(T_ f1, T_ f2, double delta)
{
	return static_cast<double>(delta * (f2-f1)/2);
}

} // End of namespace elektron
#endif
