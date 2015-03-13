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
	amplitude = interaction->a0*interaction->a1/interaction->ro;
	length = sys->get_repulsiveforce_length();
	force_vector.reset();
	force_norm = 0;
	stresslet_XF = 0;
}

void RepulsiveForce::calcForce()
{
	/** 
		\brief Compute the repulsive force.
		
		The force is normal and has an amplitude \f$ f_{R} = f_{R}^0
		\exp(-h/\lambda) \f$ if \f$h>0\f$ and \f$ f_{R} = f_{R}^0 \f$
		if \f$h<0\f$, where \f$h\f$ is the interparticle gap.
	*/
	double gap = interaction->get_gap();
	//	if(gap>0){
	if (interaction->contact.state == 0) { // why not testing for gap? When particles separate, there is a time step for which gap>0 and contact.state>0, is that the reason?
		// ---> I forgot why I did so:-)
		/* separating */
		force_norm = sys->amplitudes.repulsion*amplitude*exp(-gap/length);
		force_vector = -force_norm*interaction->nvec;
	} else {
		/* contacting */
		force_norm = sys->amplitudes.repulsion*amplitude;
		force_vector = -force_norm*interaction->nvec;
	}
}

void RepulsiveForce::addUpForce()
{
	sys->repulsive_force[p0] += force_vector;
	sys->repulsive_force[p1] -= force_vector;
}

void RepulsiveForce::calcStressXF()
{
	calcForce();
	stresslet_XF.set(interaction->rvec, force_vector);
}
