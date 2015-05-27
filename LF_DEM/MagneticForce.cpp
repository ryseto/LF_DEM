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
	if (sys->magnetic_moment_norm[p0] != 0
		&& sys->magnetic_moment_norm[p1] != 0) {
		active = true;
		coeffient = sys->magnetic_coeffient;
		if (sys->magnetic_rotation_active == false){
			chi0chi1 = sys->magnetic_susceptibility[p0]*sys->magnetic_susceptibility[p1];
		}
	} else {
		active = false;
	}
	force_vector0.reset();
	torque0.reset();
	torque1.reset();
}

void MagneticForce::calcForceToruqe()
{
	/**
		\brief Compute the magnetic force in the System class units.
	 */
	if (active == false) {
		return;
	}
	double r = interaction->get_r();
	double r_cubic = r*r*r;
	const vec3d &nvec = interaction->nvec;
	/*
	 * magnetic_coeffient = (3*mu0)/(4*M_PI)
	 */
	if (sys->magnetic_rotation_active == false){
		double H_dot_n = dot(sys->p.external_magnetic_field, nvec);
		force_vector0 = -(coeffient/(r_cubic*r))*chi0chi1*(2*H_dot_n*sys->p.external_magnetic_field
														   +(sys->magnetic_field_square-5*H_dot_n*H_dot_n)*nvec);
	} else {
		const vec3d &m0 = sys->magnetic_moment[p0];
		const vec3d &m1 = sys->magnetic_moment[p1];
		force_vector0 = -(coeffient/(r_cubic*r))*(dot(m1, nvec)*m0+dot(m0, nvec)*m1
												  +dot(m1, m0)*nvec
												  -5*dot(m0,nvec)*dot(m1,nvec)*nvec);
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
	if (sys->magnetic_rotation_active) {
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
