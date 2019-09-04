#include "PairwiseInteraction.h"
#include "PairwiseConfig.h"

namespace Interactions
{

PairwiseInteraction::PairwiseInteraction(unsigned i, unsigned j, 
										 double a_i, double a_j, 
										 vec3d sep)
										 // const Geometry::PairwiseConfig &pconf)
										 // const struct PairVelocityGrad &vgrad)
{
	if (j > i) {
		p0 = i, p1 = j;
		a0 = a_i, a1 = a_j;
	} else {
		p0 = j, p1 = i;
		a0 = a_j, a1 = a_i;
	}
	setSeparation(sep);
	// setVelocities(pconf);
	// setVelocityGrad(vgrad);
}

void PairwiseInteraction::setSeparation(const vec3d &sep)
{
	rvec = sep;
	r = rvec.norm();
	nvec = rvec/r;
	reduced_gap = 2*r/(a0+a1)-2;
}

// void PairwiseInteraction::setSeparation(const Geometry::PairwiseConfig &pconf)
// {
// 	rvec = pconf.getSeparation(p0, p1);
// 	r = rvec.norm();
// 	nvec = rvec/r;
// 	reduced_gap = 2*r/(a0+a1)-2;
// }

// void PairwiseInteraction::setVelocities(const Geometry::PairwiseConfig &pconf)
// {
// 	pconf.getVelocities(p0, p1, vel);
// 	// pconf.getBackgroundVelocities(p0, p1, velinf);
// }

// void PairwiseInteraction::setVelocityGrad(const struct PairVelocityGrad &vgrad)
// {
// 	velgrad = vgrad;
// }

// double PairwiseInteraction::separation2gap(double sep)
// {
// }

/* observation */
double PairwiseInteraction::getNormalVelocity(const struct PairVelocity &vel) const
{
	return dot(vel.U[1] - vel.U[0], nvec);
}

}