//
//  Interaction.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class Interaction
 \brief Interaction object, master class holding any interaction between two particles
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__StdInteraction__
#define __LF_DEM__StdInteraction__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
#include "vec3d.h"
#include "Contact.h"
#include "Lubrication.h"
#include "RepulsiveForce.h"
#include "VanDerWaals.h"
#include "ActivatedAdhesion.h"
#include "ActivatedAdhesion_io.h"
#include "Sym2Tensor.h"
#include "PairwiseInteraction.h"
#include "StdInteractionParams.h"

namespace Dynamics {
	class PairwiseResistanceVelocitySolver;
}
namespace Interactions {

class StdInteraction : public PairwiseInteraction {

private:
	/*********************************
	 *        Members                *
	 *********************************/
	//===== forces and stresses ==================== //
	double interaction_range;  // max distance
	bool record;
	/*********************************
	 *       Private Methods         *
	 *********************************/
	//===== forces and stresses computations =====//
	void integrateStress();
	std::pair<double, double> calcContactDashpotResistanceCoeffs();
	void switchOnContact();
	void switchOffContact();
	void updateContactState();
	void updateRepulsionState();
	void updateLubricationState();
	void updateActivatedAdhesionState(double dt);

	//===== forces/stresses  ========================== //
	double birth_strain;
	std::vector<double> strain_history;
	std::vector<double> angle_history;
	std::vector<double> normalforce_history;
	std::vector<double> gap_history;
	
	struct StdInteractionParams p;
	friend void ActAdhesion::setupInteractions(Interactions::StdInteractionManager &interactions, 
											    const std::vector <struct ActAdhesion::State> &adhesion_states);
public:
	std::unique_ptr<Contact> contact;
	std::unique_ptr<Lub::Lubrication> lubrication;
	std::unique_ptr<RepulsiveForce> repulsion;
	std::unique_ptr<ActAdhesion::ActivatedAdhesion> act_adhesion;

	Dynamics::PairwiseResistanceVelocitySolver *solver;

	 
	/*********************************
	 *       Public Methods          *
	 *********************************/
	StdInteraction(const PairId &pairid,
				   vec3d sep,	
				   double interaction_range_,
				   struct StdInteractionParams params,
				   Dynamics::PairwiseResistanceVelocitySolver *vel_solver);
	//======= state updates  ====================//
	void updateState(const struct PairVelocity &vel,
				    vec3d sep,
					double dt,
					bool& deactivated);
	double separation_distance() const {return r;}
	// void recordHistory();
	// void outputHisotry();

	bool hasPairwiseResistance();
	struct ODBlock RFU_ODBlock();
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks();

	void saveState();
	void restoreState();
};

} // namespace Interactions

#endif /* defined(__LF_DEM__StdInteraction__) */
