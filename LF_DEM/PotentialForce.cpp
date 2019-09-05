#include "PotentialForce.h"

namespace Interactions
{

void PotentialForce::addUpForce(std::vector<vec3d> &force) const
{
	force[interaction->p0] += force_vector;
	force[interaction->p1] -= force_vector;
}

void PotentialForce::addUpStressXF(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1)
{
	/**
	 \brief Compute the XF stress associated with the repulsive force.

	 */
	/* force_vector is force acting on particle 0
	 * rvec is from particle 0 to particle 1
	 */
	auto sc = 0.5*outer_sym(interaction->rvec, force_vector);
	/* NOTE:
		As the repulsive force is not a contact force, there is an ambiguity defining the stress per particle. Here we make the choice of attributing 1/2 of the interaction stress to each particle.
	*/
	stress_p0 += sc;
	stress_p1 += sc;
}

} // namespace Interactions