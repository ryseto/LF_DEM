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

Interaction::Interaction(System* sys_,
						 unsigned int i, unsigned int j,
						 double interaction_range_):
sys(sys_),
reduced_gap(0),
r(0),
interaction_range(interaction_range_),
contact_state_changed_after_predictor(false),
rvec(0),
nvec(0),
z_offset(0)
{
	if (j > i) {
		p0 = i, p1 = j;
	} else {
		p0 = j, p1 = i;
	}
	ro = sys->radius[p0]+sys->radius[p1]; // ro=a0+a1
	init();
	// tell it to particles i and j
	sys->interaction_list[i].insert(this);
	sys->interaction_list[j].insert(this);
	activateForceMembers();
	if (sys->p.output.recording_interaction_history) {
		if (sys->get_cumulated_strain() > sys->p.output.recording_start) {
			record = true;
			birth_strain = sys->get_cumulated_strain();
			strain_history.clear();
			angle_history.clear();
			normalforce_history.clear();
			gap_history.clear();
		} else {
			record = false;
		}
	}
}

Interaction::Interaction(const Interaction &other):
sys(other.sys),
p0(other.p0),
p1(other.p1),
ro(other.ro),
reduced_gap(other.reduced_gap),
r(other.r),
interaction_range(other.interaction_range),
contact_state_changed_after_predictor(other.contact_state_changed_after_predictor),
rvec(other.rvec),
nvec(other.nvec),
record(other.record),
birth_strain(other.birth_strain),
strain_history(other.strain_history),
angle_history(other.angle_history),
normalforce_history(other.normalforce_history),
gap_history(other.gap_history),
z_offset(other.z_offset)
{
	init();
	sys->interaction_list[p0].insert(this);
	sys->interaction_list[p1].insert(this);
	activateForceMembers();
}

Interaction &Interaction::operator = (const Interaction &inter)
{
	Interaction tmp(inter);
	swap(tmp);
	init();
	sys->interaction_list[p0].insert(this);
	sys->interaction_list[p1].insert(this);
	activateForceMembers();
	return *this;
}

Interaction::~Interaction()
{
	sys->interaction_list[p0].erase(this);
	sys->interaction_list[p1].erase(this);
}

void Interaction::swap(Interaction& other)
{
	sys->interaction_list[p0].erase(this);
	sys->interaction_list[p1].erase(this);
	sys->interaction_list[other.p0].erase(&other);
	sys->interaction_list[other.p1].erase(&other);
	std::swap(sys, other.sys);
	std::swap(p0, other.p0);
	std::swap(p1, other.p1);
	std::swap(ro, other.ro);
	std::swap(reduced_gap, other.reduced_gap);
	std::swap(r, other.r);
	std::swap(interaction_range, other.interaction_range);
	std::swap(contact_state_changed_after_predictor, other.contact_state_changed_after_predictor);
	std::swap(rvec, other.rvec);
	std::swap(nvec, other.nvec);
	std::swap(z_offset, other.z_offset);
	std::swap(contact, other.contact);
	std::swap(lubrication, other.lubrication);
	std::swap(repulsion, other.repulsion);
	delayed_adhesion.swap(other.delayed_adhesion);
	std::swap(label, other.label);
	sys->interaction_list[p0].insert(this);
	sys->interaction_list[p1].insert(this);
	sys->interaction_list[other.p0].insert(&other);
	sys->interaction_list[other.p1].insert(&other);
	std::swap(record, other.record);
	if (sys->p.output.recording_interaction_history) {
		std::swap(birth_strain, other.birth_strain);
		std::swap(strain_history, other.strain_history);
		std::swap(angle_history, other.angle_history);
		std::swap(normalforce_history, other.normalforce_history);
		std::swap(gap_history, other.gap_history);
	}
}

void Interaction::init()
{
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
	if (sys->delayed_adhesion) {
		delayed_adhesion = \
		std::unique_ptr<TActAdhesion::TimeActivatedAdhesion>(new TActAdhesion::TimeActivatedAdhesion(sys->p.TA_adhesion,
																									 p0, p1,
																									 sys->radius[p0],
																									 sys->radius[p1]));
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
	if (!sys->ext_flow) {
		z_offset = sys->periodizeDiff(rvec);
	} else {
		sys->periodizeDiffExtFlow(rvec, pd_shift, p0, p1);
	}
	r = rvec.norm();
	nvec = rvec/r;
	reduced_gap = 2*r/ro-2;
}

void Interaction::activateForceMembers()
{
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

	if (sys->delayed_adhesion) {
		delayed_adhesion->update(sys->get_time(), reduced_gap, nvec);
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
	if (sys->delayed_adhesion) {
		delayed_adhesion->deactivate();
	}
	if (sys->p.output.recording_interaction_history) {
		outputHisotry();
	}
}

void Interaction::outputHisotry()
{
	if (record) {
		unsigned dk = 20;
		for (unsigned k=0; k < strain_history.size(); k += dk) {
			double ang = angle_history[k];
			if (sys->ext_flow) {
				ang += sys->p.magic_angle;
			} else {
				ang -= M_PI/4;
			}
			if (ang < 0) {
				ang += M_PI;
			} else if (ang > M_PI) {
				ang -= M_PI;
			}
			cerr << strain_history[k] << ' ';
			cerr << ang << ' ';
			cerr << normalforce_history[k] <<' ';
			cerr << gap_history[k] << endl;
		}
		cerr << endl;
	}
	record = false;
	strain_history.clear();
	angle_history.clear();
	normalforce_history.clear();
	gap_history.clear();
}

void Interaction::updateState(bool& deactivated)
{
	if (contact.is_active()) {
		// (VERY IMPORTANT): we increment displacements BEFORE updating the normal vector not to mess up with Lees-Edwards PBC
		contact.incrementDisplacements();
	}
	calcNormalVectorDistanceGap();
	if (r > interaction_range) {
		/* all forces are switched off, but NOT the interaction itself
		This has to be done by the caller, based on the value of deactivated */
		deactivate();
		deactivated = true;
		return;
	}
	updateContactState();
	contact.calcContactSpringForce();
	if (sys->lubrication) {
		lubrication.updateActivationState();
		if (lubrication.is_active()) {
			lubrication.updateResistanceCoeff();
		}
	}
	if (sys->repulsiveforce) {
		repulsion.calcForce();
	}
	if (sys->delayed_adhesion) {
		delayed_adhesion->update(sys->get_time(), reduced_gap, nvec);
	}
}

void Interaction::updateContactState()
{
	contact_state_changed_after_predictor = false;
	if (contact.is_active()) {
		// contacting in previous step
		bool breakup_contact_bond = false;
		if (!sys->adhesion) {
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
double Interaction::getNormalVelocity() const
{
	vec3d vel_offset = z_offset*sys->get_vel_difference();
	vec3d d_velocity = sys->velocity[p1]-sys->velocity[p0]+vel_offset;
	return dot(d_velocity, nvec);
}

void Interaction::recordHistory()
{
	if (record) {
		double total_normal_force = 0;
		if (contact.is_active()) {
			total_normal_force += contact.getNormalForceValue();
		}
		if (sys->lubrication) {
			if (lubrication.is_active()) {
				total_normal_force += -lubrication.force;
			}
		}
		double interaction_strain = sys->get_cumulated_strain()-birth_strain;
		strain_history.push_back(interaction_strain);
		angle_history.push_back(nvec.angle_0_pi());
		normalforce_history.push_back(total_normal_force);
		gap_history.push_back(reduced_gap);
	}
}
