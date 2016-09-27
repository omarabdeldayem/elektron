#ifndef ELEKTRON_PARTICLE_H_
#define ELEKTRON_PARTICLE_H_

#include <random>

namespace elektron
{

struct Sample
{
	double wt;
}

template <std::size_t S_> 
class Particle
{
public:
	Particle();
	

	sample();
	select();
	belief();

private:
	int num_samples_;

#ifdef ELEKTRON_USE_HEAP
	std::vector<Sample> samples_;
#else
	std::array<Sample, S_> samples_;
#endif

}

template <std::size_t S_>
Particle::Particle()
{
	for (auto it_ = samples_.begin(); it_ < samples_.end(); ++it_)
	{
		// Sample particles
	}
}

}

#endif
