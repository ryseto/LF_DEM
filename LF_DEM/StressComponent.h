#ifndef __LF_DEM__StressComponent__
#define __LF_DEM__StressComponent__
#include <vector>
#include <string>
#include "Sym2Tensor.h"
#include "RateDependence.h"

enum class StressType : unsigned {
	velocity,
	brownian,
	velocitygrad,
	xf
};

struct StressComponent
{
	StressType type;
	RateDependence rate_dependence;
	std::string group;
	std::vector<Sym2Tensor> particle_stress;

	StressComponent(){};
	StressComponent(StressType _type,
					std::size_t size,
					RateDependence _rate_dependence,
					const std::string &_group):
	type(_type),
	rate_dependence(_rate_dependence), group(_group)
	{
		particle_stress.resize(size);
		reset();
	}

	Sym2Tensor getTotalStress() const
	{
		Sym2Tensor total_stress;
		for (const auto &s: particle_stress) {
			total_stress += s;
		}
		return total_stress;
	}

	void reset()
	{
		for (auto &s: particle_stress) {
			s.reset();
		}
	}
};

#endif/* defined(__LF_DEM__StressComponent__) */
