//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Interaction.h"

void
Interaction::init(System *sys_){
	sys = sys_;
	active = false;
	lubrication.init(sys);
	contact.init(sys, this);
}

void
Interaction::setResistanceCoeff(double normal_rc, double tangent_rc){
	lubrication.setResistanceCoeff(normal_rc, tangent_rc);
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1. ( i --> j)
 * pd_z : Periodic boundary condition
 */
void
Interaction::calcNormalVectorDistanceGap(){
	rvec = sys->position[p1]-sys->position[p0];
	sys->periodize_diff(rvec, zshift);
	r = rvec.norm();
	nvec = rvec/r;
	nxnx = nvec.x*nvec.x;
	nxny = nvec.x*nvec.y;
	nxnz = nvec.x*nvec.z;
	nynz = nvec.y*nvec.z;
	nyny = nvec.y*nvec.y;
	nznz = nvec.z*nvec.z;
	gap_nondim = r/ro_12-2;
	if (contact.active) {
		double overlap_12 = 0.5*(a0+a1-r);
		a0_dash = a0-overlap_12;
		a1_dash = a1-overlap_12;
	} else {
		double lub_coeff = 1/(gap_nondim+sys->lub_reduce_parameter);
		lubrication.setResistanceCoeff(lub_coeff, log(lub_coeff));
		a0_dash = a0;
		a1_dash = a1;
	}
}

/* Activate interaction between particles i and j.
 * Always j>i is satisfied.
 */
void
Interaction::activate(int i, int j){
	active = true;
	if (j > i) {
		p0 = i;
		p1 = j;
	} else {
		p0 = j;
		p1 = i;
	}
	// tell it to particles i and j
	sys->interaction_list[i].insert(this);
	sys->interaction_list[j].insert(this);
	// tell them their new partner
	sys->interaction_partners[i].insert(j);
	sys->interaction_partners[j].insert(i);
	a0 = sys->radius[p0];
	a1 = sys->radius[p1];
	set_ro(a0+a1); // ro=a0+a1
	interaction_range_scaled = ro_12*sys->get_lub_max();
	/* NOTE:
	 * lub_coeff_contact includes kn.
	 * If the scaled kn is used there,
	 * particle size dependence appears in the simulation.
	 * I don't understand this point yet.
	 * lub_coeff_contact_scaled = 4*kn_scaled*sys->contact_relaxation_time;
	 */
	/*
	 * The size dependence of colloidal force:
	 * a0*a1/(a1+a2)/2
	 * Is
	 */
	colloidalforce_amplitude = sys->get_colloidalforce_amplitude()*a0*a1/ro;
	colloidalforce_length = sys->get_colloidalforce_length();
	calcNormalVectorDistanceGap();
	// deal with contact
	contact.getInteractionData();
	if (gap_nondim <= 0) {
		contact.activate();
	} else {
		contact.deactivate();
	}
	lubrication.getInteractionData();
	lubrication.calcLubConstants();
}

void
Interaction::deactivate(){
	// r > lub_max
	contact.deactivate();
	active = false;
	sys->interaction_list[p0].erase(this);
	sys->interaction_list[p1].erase(this);
	sys->interaction_partners[p0].erase(p1);
	sys->interaction_partners[p1].erase(p0);
}

void
Interaction::updateState(bool &deactivated){
	/* update tangential displacement: we do it before updating nvec as:
     *  - it should be along the tangential vector defined in the previous time step
	 *  - (VERY IMPORTANT) it must compute the relative velocities with PBC at time t. 
	 *    This is encoded in variable zshift which takes care of Lees-Edwards in the z direction.
	 *    zshift is updated for time t+1 in calcNormalVectorDistanceGap(), 
	 *    so this function should be called after calcRelativeVelocities();
	 */
	if (sys->friction && is_contact()) {
		calcRelativeVelocities();
		contact.incrementTangentialDisplacement();
	}
	
	calcNormalVectorDistanceGap();
	
	if (contact.active) {
		if (gap_nondim > 0){
			contact.deactivate();
		}
	} else {
		if (gap_nondim <= 0) {
			contact.activate();
		} else if (r > interaction_range_scaled) {
			/* all interaction is switched off. */
			deactivate();
			deactivated = true;
			return;
		}
	}
	if (contact.active) {
		contact.calcContactInteraction();
	}
	if (sys->colloidalforce) {
		if (contact.active) {
			/* For continuity, the colloidal force is kept as constant for h < 0.
			 * This force does not affect the friction law,
			 * i.e. it is separated from Fc_normal_norm.
			 */
			f_colloidal_norm = colloidalforce_amplitude;
			f_colloidal = -f_colloidal_norm*nvec;
		} else {
			/* separating */
			f_colloidal_norm = colloidalforce_amplitude*exp(-(r-ro)/colloidalforce_length);
			f_colloidal = -f_colloidal_norm*nvec;
		}
	}
}

/* Relaxation to generate initial configuration.
 * This process should be reconsidered.
 */
void
Interaction::updateStateRelax(bool &deactivated){
	deactivated = false;
	if (active == false) {
		return;
	}
	calcNormalVectorDistanceGap();
	if (contact.active) {
		contact.calcContactInteractionRelax();
		if (gap_nondim > 0) {
			contact.deactivate();
		}
		f_colloidal_norm = colloidalforce_amplitude;
		f_colloidal = -f_colloidal_norm*nvec;
	} else {
		f_colloidal_norm = colloidalforce_amplitude*exp(-(r-ro)/colloidalforce_length);
		f_colloidal = -f_colloidal_norm*nvec;
		if (gap_nondim <= 0) {
			contact.activate();
		} else if (r > interaction_range_scaled) {
			deactivate();
			deactivated = true;
		}
	}
}

/*
 * Colloidal stabilizing force
 */
void
Interaction::addUpColloidalForce(){
	sys->colloidal_force[p0] += f_colloidal;
	sys->colloidal_force[p1] -= f_colloidal;
}

/* Relative velocity of particle 1 from particle 0.
 *
 * Use:
 *  sys->velocity and ang_velocity
 *
 */
void
Interaction::calcRelativeVelocities(){
	/* relative velocity particle 1 from particle 0.
	 */
	/*
	 * v1' = v1 - Lz = v1 - zshift*lz;
	 */
	/******************************************************
	 * if p1 is upper, zshift = -1.
	 * zshift = -1; //  p1 (z ~ lz), p0 (z ~ 0)
	 *
	 ******************************************************/
	relative_velocity = sys->velocity[p1]-sys->velocity[p0]; //true velocity, in predictor and in corrector
	relative_velocity.x += zshift*sys->vel_difference;
	relative_surface_velocity = relative_velocity-cross(a0*sys->ang_velocity[p0]+a1*sys->ang_velocity[p1], nvec);
	relative_surface_velocity -= dot(relative_surface_velocity, nvec)*nvec;
}

void
Interaction::addColloidalStress(){
	colloidal_stresslet_XF.set(rvec, f_colloidal);
}

/* observation */
double
Interaction::getContactVelocity(){
	if (contact.active == false) {
		return 0;
	}
	return relative_surface_velocity.norm();
}

/* observation */
double
Interaction::getNormalVelocity(){
	vec3d d_velocity = sys->velocity[p1]-sys->velocity[p0];
	if (zshift != 0) {
		d_velocity.x += zshift*sys->vel_difference;
	}
	return dot(d_velocity, nvec);
}
