//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "Interaction.h"
#include "System.h"

using namespace std;

void Interaction::init(System* sys_)
{
	sys = sys_;
	contact.init(sys, this);
	if (sys->lubrication) {
		lubrication.init(sys, this);
	}
	if (sys->repulsiveforce) {
		repulsion.init(sys, this);
	}
	if (sys->p.lubrication_model == "normal") {
	 	RFU_DBlocks_lub = &Lubrication::RFU_DBlocks_squeeze;
		RFU_ODBlock_lub = &Lubrication::RFU_ODBlock_squeeze;
	} else if (sys->p.lubrication_model == "tangential") {
	 	RFU_DBlocks_lub = &Lubrication::RFU_DBlocks_squeeze_tangential;
		RFU_ODBlock_lub = &Lubrication::RFU_ODBlock_squeeze_tangential;
	}
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1. ( i --> j)
 * pd_z : Periodic boundary condition
 */
void Interaction::calcNormalVectorDistanceGap()
{
	rvec = sys->position[p1]-sys->position[p0];
	z_offset = sys->periodize_diff(rvec);
	r = rvec.norm();
	nvec = rvec/r;
	reduced_gap = 2*r/ro-2;
}

/* Activate interaction between particles i and j.
 * Always j>i is satisfied.
 */
void Interaction::activate(unsigned int i, unsigned int j,
						   double interaction_range_)
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
	sys->interaction_partners[i].push_back(j);
	sys->interaction_partners[j].push_back(i);
	sys->updateNumberOfInteraction(p0, p1, 1);

	set_ro(sys->radius[p0]+sys->radius[p1]); // ro=a0+a1
	interaction_range = interaction_range_;
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
	calcNormalVectorDistanceGap();
	// deal with contact
	contact.setInteractionData();
	if (reduced_gap <= 0) {
		contact.activate();
	}
	contact_state_changed_after_predictor = false;
	if (sys->lubrication) {
		lubrication.setParticleData();
		lubrication.updateActivationState();
		if (lubrication.is_active()) {
			lubrication.updateResistanceCoeff();
		}
	}
}

void Interaction::deactivate()
{
	// r > interaction_range
	if (contact.is_active()) {
		contact.deactivate();
	}
	if (sys->lubrication) {
		if (lubrication.is_active()) {
			lubrication.deactivate();
		}
	}
	active = false;
	sys->interaction_list[p0].erase(this);
	sys->interaction_list[p1].erase(this);
	sys->removeNeighbors(p0,p1);
	sys->updateNumberOfInteraction(p0, p1, -1);
}

void Interaction::updateState(bool& deactivated)
{
	if (contact.is_active()) {
		// (VERY IMPORTANT): we increment displacements BEFORE updating the normal vector not to mess up with Lees-Edwards PBC
		contact.incrementDisplacements();
	}

	calcNormalVectorDistanceGap();

	if (r > interaction_range) {
		/* all interactions are switched off. */
		deactivate();
		deactivated = true;
		return;
	}

	updateContactState();
	if (contact.is_active()) {
		contact.calcContactSpringForce();
	}
	if (sys->lubrication) {
		lubrication.updateActivationState();
		if (lubrication.is_active()) {
			lubrication.updateResistanceCoeff();
		}
	}
	if (sys->repulsiveforce) {
		repulsion.calcForce();
	}
}

void Interaction::updateContactState()
{
	contact_state_changed_after_predictor = false;
	if (contact.is_active()) {
		// contacting in previous step
		bool breakup_contact_bond = false;
		if (!sys->cohesion) {
			// no cohesion: breakup based on distance
			if (reduced_gap > 0) {
				breakup_contact_bond = true;
			}
		} else {
			/*
			 * Checking cohesive bond breaking.
			 * breakup based on force
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
		if (reduced_gap <= 0) {
			// now contact
			contact.activate();
			if (sys->in_predictor && sys->brownian) {
				contact_state_changed_after_predictor = true;
			}
		}
	}
}

bool Interaction::hasPairwiseResistance()
{
	if (sys->lubrication) {
		return contact.dashpot.is_active() || lubrication.is_active();
	} else {
		return contact.dashpot.is_active();
	}
}

struct ODBlock Interaction::RFU_ODBlock()
{
	// This is a bit complex to read, we should find a better way to write a piece of code doing the same thing.
	if (!sys->lubrication && contact.dashpot.is_active()) {
		return contact.dashpot.RFU_ODBlock();
	}
	if (contact.dashpot.is_active()) {
		if (!contact_state_changed_after_predictor) {
			return contact.dashpot.RFU_ODBlock();
		} else {
			/*
			 * Brownian ONLY (contact_state_changed_after_predictor == false in other cases):
			 * we take the resistance provided by the lubrication, i.e. the resistance as it was
			 * BEFORE the first step (predictor) of the mid-point algorithm).
			 * This is to avoid a discontinuous change of the resistance between the two steps of the
			 * midpoint algo if the contact state changed during the first step.
			 * A discontinuous change has to be avoided because the Brownian force contains div(Mobility),
			 * and this is estimated as a difference in the resistance between the two steps.
			 * Allowing a discontinuity here would lead to a non-physical drift.
			 */
			 return (lubrication.*RFU_ODBlock_lub)();
		}
	}
	if (lubrication.is_active()) {
		if (!contact_state_changed_after_predictor) {
			return (lubrication.*RFU_ODBlock_lub)();
		} else {
			/*
			 * Brownian ONLY (contact_state_changed_after_predictor == false in other cases):
			 * we take the resistance provided by the contact dashpot, i.e. the resistance as it was
			 * BEFORE the first step (predictor) of the mid-point algorithm).
			 * See above for rationale.
			 */
			return contact.dashpot.RFU_ODBlock();
		}
	}
	struct ODBlock b;
	resetODBlock(b);
	return b;
}

std::pair<struct DBlock, struct DBlock> Interaction::RFU_DBlocks()
{
	// This is a bit complex to read, we should find a better way to write a piece of code doing the same thing.
	if (!sys->lubrication && contact.dashpot.is_active()) {
		return contact.dashpot.RFU_DBlocks();
	}
	if (contact.dashpot.is_active()) {
		if (!contact_state_changed_after_predictor) { // used in Brownian only see above for rationale
			return contact.dashpot.RFU_DBlocks();
		} else {
			return (lubrication.*RFU_DBlocks_lub)();
		}
	}
	if (lubrication.is_active()) {
		if (!contact_state_changed_after_predictor) { // used in Brownian only see above for rationale
			return (lubrication.*RFU_DBlocks_lub)();
		} else {
			return contact.dashpot.RFU_DBlocks();
		}
	}
	struct DBlock b;
	resetDBlock(b);
	return std::make_pair(b, b);
}

/* observation */
double Interaction::getNormalVelocity()
{
	vec3d vel_offset = z_offset*sys->get_vel_difference();
	vec3d d_velocity = sys->velocity[p1]-sys->velocity[p0] + vel_offset;
	return dot(d_velocity, nvec);
}
