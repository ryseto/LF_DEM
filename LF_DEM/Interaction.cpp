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
	f_repulsive_norm = 0;
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
}

void
Interaction::updateResistanceCoeff(){
	if (contact.state > 0) {
		//		double overlap_12 = 0.5*(a0+a1-r);
		//		if (overlap_12 < 0) {
		//			overlap_12 = 0;
		//		}
		////		a0_dash = a0-overlap_12;
		//		a1_dash = a1-overlap_12;
		if (!contact_state_changed_after_predictor) {
			lubrication.setResistanceCoeff(sys->lub_coeff_contact,
										   sys->log_lub_coeff_contact_tan_total);
		} else {
			/* This is to avoid discontinous change.
			 * Before the predictor, particles are apart.
			 * The displacement in the predictor makes the particles in contact.
			 * In the corrector of the same time step,
			 * the resistance coeffient is set to the maximum value of separating state.
			 * Thus, no drift force is generated.
			 */
			double lub_coeff = 1/sys->lub_reduce_parameter;
			lubrication.setResistanceCoeff(lub_coeff, log(lub_coeff));
		}
	} else {
		//		a0_dash = a0;
		//		a1_dash = a1;
		if (!contact_state_changed_after_predictor) {
			double lub_coeff = 1/(gap_nondim+sys->lub_reduce_parameter);
			lubrication.setResistanceCoeff(lub_coeff, log(lub_coeff));
		} else {
			/* This is to avoid discontinous change.
			 * Before the predictor, the particles are in contact.
			 * The displacement in the predictor makes particles apart.
			 * In the corrector for the same time step,
			 * the resistance coeffient is set to the ones used in contact state.
			 * Thus, no drift force is generated.
			 */
			lubrication.setResistanceCoeff(sys->lub_coeff_contact,
										   sys->log_lub_coeff_contact_tan_total);
		}
	}
}

/* Activate interaction between particles i and j.
 * Always j>i is satisfied.
 */
void
Interaction::activate(unsigned short i, unsigned short j){
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
	 * The size dependence of repulsive force:
	 * a0*a1/(a1+a2)/2
	 */
	if (sys->repulsiveforce) {
		repulsiveforce_amplitude = sys->get_repulsiveforce_amplitude()*a0*a1/ro;
		repulsiveforce_length = sys->get_repulsiveforce_length();
	}
	calcNormalVectorDistanceGap();
	// deal with contact
	contact.getInteractionData();
	if (gap_nondim <= 0) {
		contact.activate();
	} else {
		contact.deactivate();
	}
	updateResistanceCoeff();
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
		if (sys->rolling_friction) {
			calcRollingVelocities();
			contact.incrementRollingDisplacement();
		}
	}
	calcNormalVectorDistanceGap();
	contact_state_changed_after_predictor = false;
	if (contact.state > 0) {
		// contacting in previous step
		//		if (contact.kn_scaled*gap_nondim > sys->cohesive_force){
		if (gap_nondim > 0){
			contact.deactivate();
			if (sys->in_predictor) {
				contact_state_changed_after_predictor = true;
			}
		}
	} else {
		// not contacting in previous step
		if (gap_nondim <= 0) {
			// now contact
			contact.activate();
			if (sys->in_predictor) {
				contact_state_changed_after_predictor = true;
			}
		} else if (r > interaction_range_scaled) {
			/* all interaction is switched off. */
			deactivate();
			deactivated = true;
			return;
		}
	}
	updateResistanceCoeff();
	if (contact.state > 0) {
		contact.calcContactInteraction();
	}
	if (sys->repulsiveforce) {
		if (contact.state > 0) {
			/* For continuity, the repulsive force is kept as constant for h < 0.
			 * This force does not affect the friction law,
			 * i.e. it is separated from Fc_normal_norm.
			 */
			f_repulsive_norm = repulsiveforce_amplitude;
			f_repulsive = -f_repulsive_norm*nvec;
		} else {
			/* separating */
			f_repulsive_norm = repulsiveforce_amplitude*exp(-(r-ro)/repulsiveforce_length);
			f_repulsive = -f_repulsive_norm*nvec;
		}
	}
}

/*
 * Repulsive stabilizing force
 */
void
Interaction::addUpRepulsiveForce(){
	sys->repulsive_force[p0] += f_repulsive;
	sys->repulsive_force[p1] -= f_repulsive;
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
	 * [note]
	 * velocity[] = velocity_predictor[] in predictor
	 * and
	 * velocity[] = 0.5*(velocity_predictor[]+velocity_corrector[]) in corrector.
	 * Therefore, relative_velocity and relative_surface_velocity
	 * are true velocities.
	 * relative_surface_velocity is used to update `disp_tan' in contact object.
	 *
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
Interaction::calcRollingVelocities(){
	rolling_velocity = -cross(a0*sys->ang_velocity[p0]-a1*sys->ang_velocity[p1], nvec);
}

void
Interaction::calcRepulsiveStress(){
	repulsive_stresslet_XF.set(rvec, f_repulsive);
}

/* observation */
double
Interaction::getContactVelocity(){
	if (contact.state == 0) {
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
