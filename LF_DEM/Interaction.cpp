//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Interaction.h"

void Interaction::init(System *sys_)
{
	sys = sys_;
	active = false;
	lubrication.init(sys);
	contact.init(sys, this);
	repulsion.init(sys, this);
	magneticforce.init(sys, this);
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1. ( i --> j)
 * pd_z : Periodic boundary condition
 */
void Interaction::calcNormalVectorDistanceGap()
{
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
	reduced_gap = r/ro_12-2;
}

/* Activate interaction between particles i and j.
 * Always j>i is satisfied.
 */
void Interaction::activate(unsigned short i, unsigned short j, double range)
{
	active = true;
	if (j > i) {
		p0 = i, p1 = j;
	} else {
		p0 = j, p1 = i;
	}
	// tell it to particles i and j
	sys->interaction_list[i].insert(this);
	sys->interaction_list[j].insert(this);
	// tell them their new partner
	sys->interaction_partners[i].insert(j);
	sys->interaction_partners[j].insert(i);
	a0 = sys->radius[p0];
	a1 = sys->radius[p1];
	/* [note]
	 * The reduced (or effective) radius is defined as
	 * 1/a_reduced = 1/a0 + 1/a1
	 * This definition comes from sphere vs half-plane geometory of contact mechanics.
	 * If sphere contacts with a half-plane (a1 = infinity), a_reduced = a0.
	 * For equal sized spheres, a_reduced = 0.5*a0 = 0.5*a1
	 */
	a_reduced = a0*a1/(a0+a1);
	set_ro(a0+a1); // ro=a0+a1
	interaction_range = range;
	/* NOTE:
	 * lub_coeff_contact includes kn.
	 * If the scaled kn is used there,
	 * particle size dependence appears in the simulation.
	 * I don't understand this point yet.
	 * lub_coeff_contact_scaled = 4*kn_scaled*sys->contact_relaxation_time;
	 */
	if (sys->repulsiveforce) {
		repulsion.activate();
	}
	if (sys->magnetic) {
		magneticforce.activate();
	}
	
//	double lub_range = (sys->*System::calcInteractionRange)(p0, p1);
	

	
//	sq_lub_range = lub_range*lub_range;

	
	calcNormalVectorDistanceGap();
	// deal with contact
	contact.setInteractionData();
	if (reduced_gap <= 0) {
		contact.activate();
	} else {
		contact.deactivate();
	}
	
	
	
	
	
	contact_state_changed_after_predictor = false;
	lubrication.getInteractionData();
	lubrication.updateResistanceCoeff();
	lubrication.calcLubConstants();
}

void Interaction::deactivate()
{
	// r > interaction_range
	contact.deactivate();
	active = false;
	sys->interaction_list[p0].erase(this);
	sys->interaction_list[p1].erase(this);
	sys->interaction_partners[p0].erase(p1);
	sys->interaction_partners[p1].erase(p0);
}

void Interaction::updateState(bool &deactivated)
{
	if (is_contact()) {
		// (VERY IMPORTANT): we increment displacements BEFORE updating the normal vector not to mess up with Lees-Edwards PBC
		contact.incrementDisplacements();
	}
	calcNormalVectorDistanceGap();
	updateContactState(deactivated);

	lubrication.updateResistanceCoeff();
	if (contact.state > 0) {
		contact.calcContactInteraction();
	}
	if (sys->repulsiveforce) {
		repulsion.calcForce();
	}
	if (sys->magnetic) {
		magneticforce.calcForceToruqe();
	}
}

void Interaction::updateContactState(bool &deactivated)
{
	contact_state_changed_after_predictor = false;
	if (contact.state > 0) {
		// contacting in previous step
		bool breakup_contact_bond = false;
		if (!sys->cohesion) {
			if (reduced_gap > 0) {
				breakup_contact_bond = true;
			}
		} else {
			/*
			 * Checking cohesive bond breaking.
			 */
			if (contact.get_normal_load() < 0) {
				breakup_contact_bond = true;
			}
		}
		if (breakup_contact_bond) {
			contact.deactivate();
			if (sys->in_predictor && sys->brownian) {
				contact_state_changed_after_predictor = true;
			}
		}
	} else {
		// not contacting in previous step
		if (reduced_gap <= sys->new_contact_gap) {
			// now contact
			contact.activate();
			if (reduced_gap < -0.1){
				cerr << "new contact may have problem\n";
				cerr << "gap = " << reduced_gap << endl;
				//exit(1);
			}
			if (sys->in_predictor && sys->brownian) {
				contact_state_changed_after_predictor = true;
			}
		} else if (r > interaction_range) {
			/* all interaction is switched off. */
			deactivate();
			deactivated = true;
			return;
		}
	}
}


/* Relative velocity of particle 1 from particle 0.
 *
 * Use:
 *  sys->velocity and ang_velocity
 *
 */
void Interaction::calcRelativeVelocities()
{
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
	if (zshift != 0) {
		relative_velocity.x += zshift*sys->vel_difference;
	}
	relative_surface_velocity = relative_velocity-cross(a0*sys->ang_velocity[p0]+a1*sys->ang_velocity[p1], nvec);
	relative_surface_velocity -= dot(relative_surface_velocity, nvec)*nvec;
}

void Interaction::calcRollingVelocities()
{
	/**
	 Calculate rolling velocity
	 We follow Luding(2008).
	 The factor 2 is added.
	 cf. Book by Marshall and Li
	 equation 3.6.13 ??
	 */
	rolling_velocity = 2*a_reduced*cross(sys->ang_velocity[p1]-sys->ang_velocity[p0], nvec);
}

/* observation */
double Interaction::getContactVelocity()
{
	if (contact.state == 0) {
		return 0;
	}
	return relative_surface_velocity.norm();
}

/* observation */
double Interaction::getNormalVelocity()
{
	vec3d d_velocity = sys->velocity[p1]-sys->velocity[p0];
	if (zshift != 0) {
		d_velocity.x += zshift*sys->vel_difference;
	}
	return dot(d_velocity, nvec);
}
