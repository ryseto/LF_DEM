#include "PairwiseInteraction.h"
#include "PairwiseConfig.h"

namespace Interactions
{

PairwiseInteraction::PairwiseInteraction(const PairId &data, vec3d sep)
{
	if (data.p1 > data.p0) {
		p0 = data.p0, p1 = data.p1;
		a0 = data.a0, a1 = data.a1;
	} else {
		p0 = data.p1, p1 = data.p0;
		a0 = data.a1, a1 = data.a0;
	}
	setSeparation(sep);
}

void PairwiseInteraction::setSeparation(const vec3d &sep)
{
	rvec = sep;
	r = rvec.norm();
	nvec = rvec/r;
	reduced_gap = 2*r/(a0+a1)-2;
}

/* observation */
double PairwiseInteraction::getNormalVelocity(const struct PairVelocity &vel) const
{
	return dot(vel.U[1] - vel.U[0], nvec);
}

}