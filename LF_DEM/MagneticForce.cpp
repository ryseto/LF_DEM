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
	force_vector0.reset();
	torque0.reset();
	torque1.reset();
	coeffient = sys->magnetic_coeffient;
}

void MagneticForce::calcForceToruqe()
{
	/**
		\brief Compute the magnetic force in the System class units.
	 */
	double r = interaction->get_r();
	double r_cubic = r*r*r;
	const vec3d &nvec = interaction->nvec;
	const vec3d &m0 = sys->magnetic_moment[p0];
	const vec3d &m1 = sys->magnetic_moment[p1];
	/*
	 * magnetic_coeffient = (3*mu0)/(4*M_PI)
	 */
	force_vector0 = -(coeffient/(r_cubic*r))*(dot(m1, nvec)*m0+dot(m0, nvec)*m1
											  +dot(m1, m0)*nvec
											  -5*dot(m0,nvec)*dot(m1,nvec)*nvec);
	if (sys->p.magnetic_type == 1) {
		vec3d b0 = (coeffient/r_cubic)*(3*dot(nvec, m0)*nvec-m0);
		vec3d b1 = (coeffient/r_cubic)*(3*dot(nvec, m1)*nvec-m1);
		torque0 = cross(m0, b1);
		torque1 = cross(m1, b0);
	}
}

void MagneticForce::addUpForceTorque()
{
	sys->magnetic_force[p0] += force_vector0;
	sys->magnetic_force[p1] -= force_vector0;
	if (sys->p.magnetic_type == 1) {
		sys->magnetic_torque[p0] += torque0;
		sys->magnetic_torque[p1] += torque1;
	}
}

double MagneticForce::calcEnergy()
{
	double energy;
	double r = interaction->get_r();
	double r_cubic = r*r*r;
	const vec3d &m0 = sys->magnetic_moment[p0];
	const vec3d &m1 = sys->magnetic_moment[p1];
	const vec3d &nvec = interaction->nvec;
	/*
	 * magnetic_coeffient/3 = mu0/(4*M_PI)
	 */
	energy = -(coeffient/(3*r_cubic))*(3*dot(m0, nvec)*dot(m1, nvec)-dot(m0,m1));
	return energy;
}
