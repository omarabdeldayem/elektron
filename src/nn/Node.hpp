#ifndef ELEKTRON_NODE_H_
#define ELEKTRON_NODE_H_

#include "../blas/Matrix.hpp"

#include <cstddef>

namespace elektron
{

template <typename T_, std::size_t INPUTS_>
struct Neuron
{
	Neuron();
	Neuron(Matrix<T, INPUTS_, 1> wts);

	Marix<T_, INPUTS_, 1> wts; // Weights
}

template <typename T_, std::size_t INPUTS_>
Neuron<T_, INPUTS_> Neuron<T_, INPUTS_>::Neuron<T_, INPUTS_>()
{ 
	// DONOTHING FOR NOW
}

template <typename T_, std::size_t INPUTS_>
Neuron<T_, INPUTS_> Neuron<T_, INPUTS_>::Neuron<T_, INPUTS_>(Matrix<T_, INPUTS_, 1>) :
	wts(wts)
{
   // DONOTHING FOR NOW	
}

} // End of namespace elektron

#endif

