//
//
//  Created by Ryohei Seto and Romain Mari on 02/11/15.
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <vector>
#include "VanDerWaals.h"

namespace Interactions
{

vanDerWaalsForce::vanDerWaalsForce(PairwiseInteraction* interaction_, struct vanDerWaalsForceParams params) :
PotentialForce(interaction_),
p(params)
{
	geometric_factor = interaction->a0*interaction->a1/(interaction->a0+interaction->a1);
	calcForce();
}


void vanDerWaalsForce::calcReducedForceNorm()
{
	/**
		\brief Compute the repulsive force in its own units.
	*/
	double gap = interaction->getGap();
	if (gap > 0) {
		/* van del Waals attraction */
		double hh = gap+p.singularity_cutoff;
		reduced_force_norm = -p.coefficient/hh/hh;
		double f_tmp = 1-p.coefficient/p.singularity_cutoff/p.singularity_cutoff;
		reduced_force_norm /= f_tmp;
		reduced_force_norm *= geometric_factor;
	}
}

void vanDerWaalsForce::calcScaledForce()
{
	/**
		\brief Computes the force from a previously computed reduced force.
	*/
	force_norm = p.amplitude*reduced_force_norm;
	/* nvec is from particle 0 to particle 1.
	 * force_vector is force acting on particle 0
	 */
	force_vector = -force_norm*interaction->nvec;
}

void vanDerWaalsForce::calcForce()
{
	/**
		\brief Compute the repulsive force.

		The force is normal and has an amplitude \f$ f_{R} = f_{R}^0
		\exp(-h/\lambda) \f$ if \f$h>0\f$ and \f$ f_{R} = f_{R}^0 \f$
		if \f$h<0\f$, where \f$h\f$ is the interparticle gap.
	*/
	calcReducedForceNorm();
	calcScaledForce();
}

double vanDerWaalsForce::calcEnergy() const
{
	double energy;
	// ????
	throw std::runtime_error("calcEnergy not yet implemented for VanDerWaals");
	return energy;
}

} // namespace Interactions