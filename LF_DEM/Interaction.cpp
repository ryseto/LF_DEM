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
	rvec = sys->position[par_num[1]]-sys->position[par_num[0]];
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
		par_num[0] = i;
		par_num[1] = j;
	} else {
		par_num[0] = j;
		par_num[1] = i;
	}
	// tell it to particles i and j
	sys->interaction_list[i].insert(this);
	sys->interaction_list[j].insert(this);
	// tell them their new partner
	sys->interaction_partners[i].insert(j);
	sys->interaction_partners[j].insert(i);
	a0 = sys->radius[par_num[0]];
	a1 = sys->radius[par_num[1]];
	set_ro(a0+a1); // ro=a0+a1
	interaction_range_scaled = ro_12*sys->get_lub_max();
	/* NOTE:
	 * lub_coeff_contact includes kn.
	 * If the scaled kn is used there,
	 * particle size dependence appears in the simulation.
	 * I don't understand this point yet.
	 * lub_coeff_contact_scaled = 4*kn_scaled*sys->contact_relaxzation_time;
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
	contact.resetObservables();
	
	lubrication.getInteractionData();
	strain_lub_start = sys->get_shear_strain(); // for output
	lubrication.calcLubConstants();
}

void
Interaction::deactivate(){
	// r > lub_max
	outputSummary();
	contact.deactivate();
	active = false;
	sys->interaction_list[par_num[0]].erase(this);
	sys->interaction_list[par_num[1]].erase(this);
	sys->interaction_partners[par_num[0]].erase(par_num[1]);
	sys->interaction_partners[par_num[1]].erase(par_num[0]);
}

void
Interaction::updateState(bool &deactivated){
	/* update tangential displacement: we do it before updating nvec
	 * as it should be along the tangential vector defined in the previous time step
	 */
	calcNormalVectorDistanceGap();
	if (contact.active) {
		if (gap_nondim > 0){
			contact.deactivate();
		}
	} else {
		if (gap_nondim <= 0) {
			contact.activate();
		} else if (r > interaction_range_scaled) {
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

void
Interaction::updateStateRelax(bool &deactivated){
	deactivated = false;
	if (active == false) {
		return;
	}
	/* update tangential displacement: we do it before updating nvec
	 * as it should be along the tangential vector defined in the previous time step
	 */
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
	sys->colloidal_force[par_num[0]] += f_colloidal;
	sys->colloidal_force[par_num[1]] -= f_colloidal;
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
	/**** NOTE ********************************************
	 * In the Corrector, this relative_surface_velocity
	 * is also the correcting velocity.
	 * This correcting velocity should not involve the
	 * velocity diffrence due to crossing the z boundary.
	 * fix_interaction_status = true : in the Predictor
	 * fix_interaction_status = false : in the Corrector
	 *
	 * if p1 is upper, zshift = -1.
	 * zshift = -1; //  p1 (z ~ lz), p0 (z ~ 0)
	 *
	 ******************************************************/
	vec3d dv = sys->velocity[par_num[1]]-sys->velocity[par_num[0]];
	if (sys->in_predictor && zshift != 0) {
		dv.x += zshift*sys->vel_difference;
	}
	relative_surface_velocity = dv-cross(a0*sys->ang_velocity[par_num[0]]+a1*sys->ang_velocity[par_num[1]], nvec);
	relative_surface_velocity -= dot(relative_surface_velocity, nvec)*nvec;
}

void
Interaction::addColloidalStress(){
	colloidal_stresslet_XF.set(rvec, f_colloidal);
}

void
Interaction::outputSummary(){
	duration = sys->get_shear_strain()-strain_lub_start;
	sys->fout_int_data << strain_lub_start << ' '; // 1
	sys->fout_int_data << duration << ' '; // 2
	sys->fout_int_data << contact.get_duration() << ' '; // 3
	sys->fout_int_data << endl;
}

double
Interaction::getContactVelocity(){
	if (contact.active == false) {
		return 0;
	}
	return relative_surface_velocity.norm();
}

double
Interaction::getNormalVelocity(){
	sys->in_predictor = true;
	vec3d d_velocity = sys->velocity[par_num[1]]-sys->velocity[par_num[0]];
	if (zshift != 0) {
		d_velocity.x += zshift*sys->vel_difference;
	}
	return dot(d_velocity, nvec);
}

double
Interaction::getPotentialEnergy(){
	double energy;
	if (gap_nondim < 0) {
		energy = 0.5*sys->get_kn()*gap_nondim*gap_nondim;
		energy += -colloidalforce_amplitude*gap_nondim;
	} else {
		energy = colloidalforce_length*colloidalforce_amplitude*(exp(-(r-ro)/colloidalforce_length)-1);
	}
	return energy;
}
