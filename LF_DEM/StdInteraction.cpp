//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#include "StdInteraction.h"
#include "PairwiseResistanceVelocitySolver.h"

namespace Interactions {

StdInteraction::StdInteraction(const PairId &pairid,
							   vec3d sep,
							   double interaction_range_,
							   struct StdInteractionParams params,
							   Dynamics::PairwiseResistanceVelocitySolver *vel_solver = nullptr) :
PairwiseInteraction(pairid, sep),
interaction_range(interaction_range_),
p(std::move(params)),
solver(vel_solver)
{
	setSeparation(sep);

	if (p.contp && reduced_gap <= 0) {
		switchOnContact();
	}

	// Lub
	if (p.lubp) {
		updateLubricationState();
	}

	if (p.repp) {
		p.repp->max_gap = calcRepulsiveForceMaxGap(*(p.repp));
		updateRepulsionState();
	}
	
	if (p.actadhp) {
		updateActivatedAdhesionState(0);
	}

	if (!vel_solver && (p.lubp || has_dashpot(*(p.contp)))) { // should not happen as it is screened in InteractionManager
		throw std::runtime_error(" StdIntercation:: Error: Dashpot and lubrication require a pairwise velocity solver.");
	}
}

void StdInteraction::switchOnContact()
{
	auto dash_coeffs = calcContactDashpotResistanceCoeffs();
	contact = std::unique_ptr<Contact>(new Contact (this, *(p.contp), dash_coeffs.first, dash_coeffs.second));
	if (!lubrication) {
		solver->declareResistance(p0, p1);
	}
}

void StdInteraction::switchOffContact()
{
	contact.reset(nullptr);
	if (!lubrication) {
		solver->eraseResistance(p0, p1);
	}
}


std::pair<double, double> StdInteraction::calcContactDashpotResistanceCoeffs()
{
	double normal_coeff, tangential_coeff;
	if (p.contp->relaxation_time >= 0) {
		/* t = beta/kn
		 *  beta = t*kn
		 * normal_coeff = 4*beta = 4*kn*rtime_normal
		 */
		if (p.contp->relaxation_time == 0) {
			std::cerr << "WARNING: The relaxation time is set to zero, which set the dashpot resistance zero.\n" ;
		}
		normal_coeff = 4*p.contp->kn*p.contp->relaxation_time;
	} else {
		if (p.lubp && p.lubp->model != "none") { // take the same resistance as lubrication
			// 1/(h+c) --> 1/c
			normal_coeff = 1/p.lubp->regularization_length;
		} else {
			throw std::runtime_error(" ContactDashpot:: Error: normal relaxation time set <=0, but no lubrication.");
		}
	}

	if (p.contp->friction_model == FrictionModel::frictionless && p.lubp && p.lubp->model != "tangential") {
		tangential_coeff = 0;
	} else {
		if (p.contp->relaxation_time_tan >= 0) {
			if (p.contp->relaxation_time_tan == 0) {
				std::cerr << "WARNING: The relaxation time is set to zero, which set the dashpot resistance zero.\n" ;
			} else {
				if (p.lubp && p.lubp->model != "none") { // the contact can get unstable if the tangential resistance difference is too big between with and wihout contact
					throw std::runtime_error(" ContactDashpot:: Error: with lubrication, tangential relaxation time cannot be set positive.");
				}
			}
			tangential_coeff = 6*p.contp->kt*p.contp->relaxation_time_tan;
		} else {
			if (p.lubp && p.lubp->model != "none") {// take the same resistance as lubrication
				// 1/(h+c) --> 1/c
				if (p.lubp->regularization_length < 1) {
					tangential_coeff = log(1/p.lubp->regularization_length);
				} else {
					tangential_coeff = 0;
				}
			} else {
				throw std::runtime_error(" ContactDashpot:: Error: tangential relaxation time set <=0, but no lubrication.");
			}
		}
	}
	return std::make_pair(normal_coeff, tangential_coeff);
}

void StdInteraction::updateState(const struct PairVelocity &vel,
								vec3d sep,
	// const Geometry::PairwiseConfig &pconf, 
							 	double dt,
							 	bool& deactivated)
{
	// setVelocities(pconf);
	if (contact) {
		// (VERY IMPORTANT): we increment displacements BEFORE updating the normal vector not to mess up with Lees-Edwards PBC
		contact->incrementDisplacements(dt, vel);
	}

	setSeparation(sep);
	// setSeparation(pconf);
	if (r > interaction_range) {
		deactivated = true;
		return;
	}

	if (p.contp) {
		updateContactState();
	}
	if (p.lubp) {
		// Lub
		updateLubricationState();
	}

	if (p.repp) {
		updateRepulsionState();
	}
	if (p.actadhp) {
		updateActivatedAdhesionState(dt);
	}
}

void StdInteraction::updateLubricationState()
{
	if (lubrication) {
		// check if it is still alive
		if (reduced_gap > p.lubp->max_gap || contact) {
			lubrication.reset(nullptr);
			if (!(contact && contact->dashpot)) {
				solver->eraseResistance(p0, p1);
			}
		} else {
			lubrication->updateResistanceCoeff();
		}
	} else {
		if (reduced_gap <= p.lubp->max_gap && !contact) {
			lubrication = std::unique_ptr<Lub::Lubrication>(new Lub::Lubrication (this, *(p.lubp)));
			solver->declareResistance(p0, p1);
		}
	}
}

void StdInteraction::updateContactState()
{
	if (contact) {
		// contacting in previous step
		bool breakup_contact_bond = false;
		if (!has_adhesion(*(p.contp))) {
			// no cohesion: breakup based on distance
			if (reduced_gap > 0) {
				breakup_contact_bond = true;
			}
		} else {
			/*
			 * Checking cohesive bond breaking.
			 * breakup based on force
			 */
			if (-contact->getNormalSpringForce() > p.contp->adhesion) {
				breakup_contact_bond = true;
			}
		}
		if (breakup_contact_bond) {
			switchOffContact();
		}
	} else {
		// not contacting in previous step
		if (reduced_gap <= 0) {
			// now contact
			switchOnContact();
		}
	}
	if (contact) {
		contact->calcContactSpringForce();
	}
}

void StdInteraction::updateRepulsionState()
{
	if (repulsion) {
		// check if it is still alive
		if (r > p.repp->max_gap) {
			repulsion.reset(nullptr);
		}
	} else {
		if (r <= p.repp->max_gap) {
			repulsion = std::unique_ptr<RepulsiveForce>(new RepulsiveForce (this, *(p.repp)));
		}
	}
	if (repulsion) {
		repulsion->calcForce();
	}
}

void StdInteraction::updateActivatedAdhesionState(double dt)
{
	if (act_adhesion) {
		// check if it is still alive
		if (reduced_gap > p.actadhp->range) {
			act_adhesion.reset(nullptr);
		}
	} else {
		if (reduced_gap <= p.actadhp->range) {
			act_adhesion = std::unique_ptr<ActAdhesion::ActivatedAdhesion>(
								new ActAdhesion::ActivatedAdhesion(this, *(p.actadhp)));
		}
	}
	if (act_adhesion) {
		act_adhesion->update(dt);
	}
}


bool StdInteraction::hasPairwiseResistance()
{
	return lubrication || (contact && contact->dashpot) ;
}

struct ODBlock StdInteraction::RFU_ODBlock()
{
	// This is a bit complex to read, we should find a better way to write a piece of code doing the same thing.
	if (!lubrication && (contact && contact->dashpot)) {
		return contact->dashpot->RFU_ODBlock();
	}
	if (contact && contact->dashpot) {
		return contact->dashpot->RFU_ODBlock();
	}
	if (lubrication) {
		return lubrication->RFU_ODBlock();
	}
	struct ODBlock b;
	resetODBlock(b);
	return b;
}

std::pair<struct DBlock, struct DBlock> StdInteraction::RFU_DBlocks()
{
	// This is a bit complex to read, we should find a better way to write a piece of code doing the same thing.
	if (!lubrication && (contact && contact->dashpot)) {
		return contact->dashpot->RFU_DBlocks();
	}
	if (contact && contact->dashpot) {
		return contact->dashpot->RFU_DBlocks();
	}
	if (lubrication) {
		return lubrication->RFU_DBlocks();
	}
	struct DBlock b;
	resetDBlock(b);
	return std::make_pair(b, b);
}

void StdInteraction::saveState()
{
	if (contact) {
		contact->saveState();
	}
}

void StdInteraction::restoreState()
{
	if (contact) {
		contact->restoreState();
	}
}


// void Interaction::recordHistory()
// {
// 	if (record) {
// 		double total_normal_force = 0;
// 		if (contact.is_active()) {
// 			total_normal_force += contact.getNormalForceValue();
// 		}
// 		if (sys->lubrication) {
// 			if (lubrication.is_active()) {
// 				total_normal_force += -lubrication.force;
// 			}
// 		}
// 		double interaction_strain = sys->get_cumulated_strain()-birth_strain;
// 		strain_history.push_back(interaction_strain);
// 		angle_history.push_back(nvec.angle_0_pi());
// 		normalforce_history.push_back(total_normal_force);
// 		gap_history.push_back(reduced_gap);
// 	}
// }

// void Interaction::outputHisotry()
// {
// 	if (record) {
// 		unsigned dk = 20;
// 		for (unsigned k=0; k < strain_history.size(); k += dk) {
// 			double ang = angle_history[k];
// 			if (sys->simu_type == sys->SimulationType::extensional_flow) {
// 				ang += sys->p.magic_angle;
// 			} else {
// 				ang -= M_PI/4;
// 			}
// 			if (ang < 0) {
// 				ang += M_PI;
// 			} else if (ang > M_PI) {
// 				ang -= M_PI;
// 			}
// 			cerr << strain_history[k] << ' ';
// 			cerr << ang << ' ';
// 			cerr << normalforce_history[k] <<' ';
// 			cerr << gap_history[k] << endl;
// 		}
// 		cerr << endl;
// 	}
// 	record = false;
// 	strain_history.clear();
// 	angle_history.clear();
// 	normalforce_history.clear();
// 	gap_history.clear();
// }

} // namespace Interactions