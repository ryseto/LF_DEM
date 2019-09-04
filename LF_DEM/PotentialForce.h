#ifndef __LF_DEM__PotentialForce__
#define __LF_DEM__PotentialForce__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "PairwiseInteraction.h"

namespace Interactions
{

class PotentialForce {
protected:
	PairwiseInteraction* interaction;
	vec3d force_vector; // normal contact force
	double force_norm;
	
public:

	PotentialForce(PairwiseInteraction* interaction_) : interaction(interaction_) {};
	virtual void calcForce() = 0;
	void addUpForce(std::vector<vec3d> &force) const;
	inline double getForceNorm() const
	{
		return force_norm;
	}
	vec3d getForceVector() const
	{
		return force_vector;
	}
	void addUpStressXF(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1);
	virtual double calcEnergy() const = 0;
};

} // namespace Interactions

#endif /* defined(__LF_DEM__PotentialForce__) */
