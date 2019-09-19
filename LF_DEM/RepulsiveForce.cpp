//
//
//  Created by Ryohei Seto and Romain Mari on 02/11/15.
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#include <vector>
#include "RepulsiveForce.h"

namespace Interactions
{

RepulsiveForce::RepulsiveForce(PairwiseInteraction* interaction_, struct RepulsiveForceParams params) :
PotentialForce(interaction_),
p(params)
{
	geometric_factor = interaction->a0*interaction->a1/(interaction->a0+interaction->a1);
	calcForce();
}

void RepulsiveForce::calcReducedForceNorm()
{
	/**
		\brief Compute the repulsive force in its own units.

		The force is normal and has an amplitude \f$ f_{R} = f_{R}^0
		\exp(-h/\lambda) \f$ if \f$h>0\f$ and \f$ f_{R} = f_{R}^0 \f$
		if \f$h<0\f$, where \f$h\f$ is the interparticle gap.

		This method returns an amplitude \f$ \hat{f}_{R} = \exp(-h/\lambda)\f$.
	*/
	double gap = interaction->getGap();
	if (gap > 0) {
		/* separating */
		if (p.smoothing == 0) {
			reduced_force_norm = exp(-gap/p.screening_length);
		} else {
			reduced_force_norm = exp(-gap/p.screening_length)*0.5*(1+tanh(-(gap-p.max_gap+3*p.smoothing)/p.smoothing));
		}
		reduced_force_norm *= geometric_factor;
	} else {
		/* contacting */
		reduced_force_norm = geometric_factor;
	}
}

void RepulsiveForce::calcScaledForce()
{
	/**
		\brief Computes the force from a previously computed reduced force.
	*/
	force_norm = p.repulsion*reduced_force_norm;
	/* nvec is from particle 0 to particle 1.
	 * force_vector is force acting on particle 0
	 */
	force_vector = -force_norm*interaction->nvec;
}

void RepulsiveForce::calcForce()
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

double RepulsiveForce::calcEnergy() const
{
	double energy;
	double gap = interaction->getGap();
	if (gap > 0) {
		/* separating */
		if (p.smoothing == 0) {
			energy = geometric_factor*p.screening_length*exp(-gap/p.screening_length);
		} else {
			// This energy is not exact one to generate the cut-offed repulsive force.
			energy = geometric_factor*p.screening_length*exp(-gap/p.screening_length)*0.5*(1+tanh(-(gap-p.max_gap)/p.smoothing));
		}
		//		reduced_force_vector = -force_norm*interaction->nvec;
	} else {
		/* contacting */
		energy = geometric_factor*p.screening_length*(1-gap);
	}
	return energy;
}

} // namespace Interactions