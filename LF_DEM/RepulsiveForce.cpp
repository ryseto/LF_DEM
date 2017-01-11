//
//
//  Created by Ryohei Seto and Romain Mari on 02/11/15.
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "RepulsiveForce.h"
#include "Interaction.h"
#include "System.h"

void RepulsiveForce::init(System* sys_, Interaction* interaction_)
{
	sys = sys_;
	interaction = interaction_;
	force_norm = 0;
}

void RepulsiveForce::activate()
{
	std::tie(p0, p1) = interaction->get_par_num();
	/*
	 * The size dependence of repulsive force:
	 * a0*a1/(a1+a2)
	 */
	double a0 = sys->radius[p0];
	double a1 = sys->radius[p1];
	geometric_factor = a0*a1/(a0+a1);
	screening_length = sys->p.repulsive_length;
	max_length = sys->p.repulsive_max_length;
	force_vector.reset();
	force_norm = 0;
	reduced_force_norm = 0;
	stresslet_XF.reset();
	cutoff_roundlength = 1e-3;
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
	double gap = interaction->get_gap();
	if (gap > 0) {
		/* separating */
		if (max_length == -1) {
			reduced_force_norm = geometric_factor*exp(-gap/screening_length);
		} else {
			reduced_force_norm = geometric_factor*exp(-gap/screening_length)*0.5*(1+tanh(-(gap-max_length)/cutoff_roundlength));
		}
	} else {
		/* contacting */
		reduced_force_norm = geometric_factor;
	}
}

void RepulsiveForce::calcScaledForce()
{
	/**
		\brief Computes the force in the System class from a previously computed reduced force.
	*/
	force_norm = sys->p.repulsion*reduced_force_norm;
	/* nvec is from particle 0 to particle 1.
	 * force_vector is force acting on particle 0
	 */
	force_vector = -force_norm*interaction->nvec;
}

void RepulsiveForce::calcForce()
{
	/**
		\brief Compute the repulsive force in the System class units.

		The force is normal and has an amplitude \f$ f_{R} = f_{R}^0
		\exp(-h/\lambda) \f$ if \f$h>0\f$ and \f$ f_{R} = f_{R}^0 \f$
		if \f$h<0\f$, where \f$h\f$ is the interparticle gap.
	*/

	calcReducedForceNorm();
	calcScaledForce();
}

void RepulsiveForce::addUpForce(std::vector<vec3d> &force) const
{
	force[p0] += force_vector;
	force[p1] -= force_vector;
}

void RepulsiveForce::calcStressXF()
{
	/**
	 \brief Compute the XF stress associated with the repulsive force.

	 \b NOTE: this method does not recompute the reduced force, this force must be first computed by calcReducedForce().
	 This method however converts the force in the System units from the reduced force.
	 */
	calcScaledForce();
	/* force_vector is force acting on particle 0
	 * rvec is from particle 0 to particle 1
	 */
	stresslet_XF = outer_sym(interaction->rvec, force_vector);
}

void RepulsiveForce::addUpStressXF(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1)
{
	calcStressXF();
	/* NOTE:
		As the repulsive force is not a contact force, there is an ambiguity defining the stress per particle. Here we make the choice of attributing 1/2 of the interaction stress to each particle.
	*/
	Sym2Tensor sc = 0.5*stresslet_XF;
	stress_p0 += sc;
	stress_p1 += sc;
}

double RepulsiveForce::calcEnergy() const
{
	double energy;
	double gap = interaction->get_gap();
	if (gap > 0) {
		/* separating */
		if (max_length == -1) {
			energy = geometric_factor*screening_length*exp(-gap/screening_length);
		} else {
			// This energy is not exact one to generate the cut-offed repulsive force.
			energy = geometric_factor*screening_length*exp(-gap/screening_length)*0.5*(1+tanh(-(gap-max_length)/cutoff_roundlength));
		}
		//		reduced_force_vector = -force_norm*interaction->nvec;
	} else {
		/* contacting */
		energy = geometric_factor*screening_length*(1-gap);
	}
	return energy;
}
