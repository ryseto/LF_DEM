//
//
//  Created by Ryohei Seto and Romain Mari on 02/11/15.
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "RepulsiveForce.h"
#include "Interaction.h"

void RepulsiveForce::init(System *sys_, Interaction *interaction_)
{
	sys = sys_;
	interaction = interaction_;
}

void RepulsiveForce::activate()
{
	interaction->get_par_num(p0, p1);
	/*
	 * The size dependence of repulsive force:
	 * a0*a1/(a1+a2)/2
	 */
	geometric_factor = interaction->a0*interaction->a1/interaction->ro;
	length = sys->get_repulsiveforce_length();
	force_vector.reset();
	force_norm = 0;
	reduced_force_norm = 0;
	stresslet_XF = 0;
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
	//	if(gap>0){
	if (interaction->contact.state == 0) { // why not testing for gap? When particles separate, there is a time step for which gap>0 and contact.state>0, is that the reason?
		// ---> I forgot why I did so:-)
		/* separating */
		reduced_force_norm = geometric_factor*exp(-gap/length);
		//		reduced_force_vector = -force_norm*interaction->nvec;
	} else {
		/* contacting */
		reduced_force_norm = geometric_factor;
		//		reduced_force_vector = -force_norm*interaction->nvec;
	}


}
	
void RepulsiveForce::calcScaledForce()
{
	/** 
		\brief Computes the force in the System class from a previously computed reduced force.
	*/

	force_norm = sys->amplitudes.repulsion*reduced_force_norm;
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

void RepulsiveForce::addUpForce()
{
	sys->repulsive_force[p0] += force_vector;
	sys->repulsive_force[p1] -= force_vector;
}

void RepulsiveForce::calcStressXF()
{
	/**
	   \brief Compute the XF stress associated with the repulsive force.

	   \b NOTE: this method does not recompute the reduced force, this force must be first computed by calcReducedForce().
	   This method however converts the force in the System units from the reduced force.
	*/
	calcScaledForce();
	stresslet_XF.set(interaction->rvec, force_vector);
}

double RepulsiveForce::calcEnergy()
{
	double energy;
	
	double gap = interaction->get_gap();
	if (interaction->contact.state == 0) {
		/* separating */
		energy = geometric_factor*length*exp(-gap/length);
		//		reduced_force_vector = -force_norm*interaction->nvec;
	} else {
		/* contacting */
		energy = geometric_factor*length*(1-gap);
	}
	return energy;
}
