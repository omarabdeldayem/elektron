#ifndef ELEKTRON_LAYER_H_
#define ELEKTRON_LAYER_H_

#include "Neuron.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace elektron
{

template <std::size_t H_>
class Layer
{
public:

#ifdef ELEKTRON_USE_HEAP
#else
	void add_neuron(Neuron* n);
#endif
	
private:

#ifdef ELEKTRON_USE_HEAP
	std::array<Neuron*, H_> neurons_;
#else
	std::vector<Neuron*, H_> neurons_;
#endif	
	
}	

}
#endif
