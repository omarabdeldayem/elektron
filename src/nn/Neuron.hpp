#ifndef ELEKTRON_NODE_H_
#define ELEKTRON_NODE_H_

#include "../math/Matrix.hpp"

#include <cstddef>

namespace elektron
{

template <std::size_t I_>
struct Neuron
{
	Neuron();
	Neuron(const CVec<double, I_>& wts);

	CVec<double, I_> wts; // Weights
}

template <std::size_t I_>
Neuron<I_>::Neuron()
{ 
	wts = CVec<double, I_>(r);
}

template <std::size_t I_>
Neuron<I_>::Neuron(CVec<double, I_>) :
	wts(wts)
{
}

} // End of namespace elektron

#endif

