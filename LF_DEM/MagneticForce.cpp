//
//  MagneticForce.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 4/17/15.
//  Copyright (c) 2015 Ryohei Seto. All rights reserved.
//

#include "MagneticForce.h"
#include "Interaction.h"

void MagneticForce::init(System *sys_, Interaction *interaction_)
{
	sys = sys_;
	interaction = interaction_;
}

void MagneticForce::activate()
{
	interaction->get_par_num(p0, p1);
	/*
	 * The size dependence of repulsive force:
	 * a0*a1/(a1+a2)/2
	 */
	force_vector0.reset();
	torque0.reset();
	torque1.reset();
}

void MagneticForce::calcForceToruqe()
{
	/**
		\brief Compute the magnetic force in the System class units.
		
		The force is normal and has an amplitude \f$ f_{R} = f_{R}^0
		\exp(-h/\lambda) \f$ if \f$h>0\f$ and \f$ f_{R} = f_{R}^0 \f$
		if \f$h<0\f$, where \f$h\f$ is the interparticle gap.
	 */
	double r = interaction->get_r();
	vec3d nvec = interaction->get_nvec();
	vec3d m0 = sys->magnetic_moment[p0];
	vec3d m1 = sys->magnetic_moment[p1];
	force_vector0 = pow(r,-4)*(dot(m0, nvec)*m1
							  +dot(m1,nvec)*m0
							  +dot(m0,m1)*nvec
							  -5*dot(m0,nvec)*dot(m1,nvec)*nvec);
	vec3d b0 = pow(r,-3)*(3*dot(nvec, m0)*nvec - m0);
	vec3d b1 = pow(r,-3)*(3*dot(nvec, m1)*nvec - m1);
	torque0 = cross(m0, b1);
	torque1 = cross(m1, b0);
}

void MagneticForce::addUpForceTorque()
{
	sys->magnetic_force[p0] += force_vector0;
	sys->magnetic_force[p1] -= force_vector0;
	sys->magnetic_torque[p0] += torque0;
	sys->magnetic_torque[p1] += torque1;
}
