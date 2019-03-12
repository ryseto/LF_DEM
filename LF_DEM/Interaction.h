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

#ifndef __LF_DEM__Interaction__
#define __LF_DEM__Interaction__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
#include "vec3d.h"
#include "Contact.h"
#include "Lubrication.h"
#include "RepulsiveForce.h"
#include "TimeActivatedAdhesion.h"
#include "Sym2Tensor.h"

class System;

class Interaction{

private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	//======= internal state =====================//
	unsigned int p0;
	unsigned int p1;
	double ro; // ro = a0+a1;
	//======= relative position/velocity data  =========//
	double reduced_gap; // gap between particles (dimensionless gap = s - 2, s = 2r/(a1+a2) )
	double r; // center-center distance
	
	//===== forces and stresses ==================== //
	double interaction_range;  // max distance
	bool record;
	/*********************************
	 *       Private Methods         *
	 *********************************/
	//===== forces and stresses computations =====//
	void calcNormalVectorDistanceGap();
	void integrateStress();
	void updateContactState();
	void activateForceMembers();
	void deactivate();
	void init();
	void swap(Interaction& other); // used by assignment operator

	struct ODBlock (Lubrication::*RFU_ODBlock_lub)() const;
	std::pair<struct DBlock, struct DBlock> (Lubrication::*RFU_DBlocks_lub)() const;
	//===== forces/stresses  ========================== //
	/* To avoid discontinous change between predictor and corrector,
	 * the change of contact state is informed in updateResiCoeff.
	 */
	bool contact_state_changed_after_predictor;
	double birth_strain;
	std::vector<double> strain_history;
	std::vector<double> angle_history;
	std::vector<double> normalforce_history;
	std::vector<double> gap_history;
	
public:
	Contact contact;
	Lubrication lubrication;
	RepulsiveForce repulsion;
	std::unique_ptr<TActAdhesion::TimeActivatedAdhesion> delayed_adhesion;
	 
	vec3d rvec; // vector center to center
	vec3d nvec; // normal vector
	int z_offset;
	vec3d pd_shift;

	/*********************************
	 *       Public Methods          *
	 *********************************/
	Interaction(System *sys_,
				unsigned int i,
				unsigned int j,
				double interaction_range_);
	// copy ctor
	Interaction (const Interaction &);
	// delete move ctor to avoid implicit implementation that would not inform sys->interaction_list.
	Interaction (Interaction &&) = delete;
	~Interaction();
	Interaction & operator = (const Interaction &inter); // assignement by copy-swap
	//======= state updates  ====================//
	/* Update the follow items:
	 * - r_vec, z_offset, _r, and nr_vec
	 * - contact_velocity_tan
	 * - disp_tan
	 * - Fc_normal and Fc_tan
	 * - check breakup of static friction
	 * - State (deactivation, contact)
	 */
	void updateState(bool& deactivated);
	double separation_distance() const {return r;}
	void recordHistory();
	void outputHisotry();

	//======= particles data  ====================//
	int partner(unsigned int i) const {return (i == p0 ? p1 : p0);}
	std::pair<unsigned int, unsigned int> get_par_num() const {return std::make_pair(p0, p1);}
	double get_reduced_gap() const {return reduced_gap;}
	double get_gap() const {return r-ro;}
	unsigned int get_p0() const {return p0;}
	unsigned int get_p1() const {return p1;}
	bool hasPairwiseResistance();
	double getNormalVelocity() const;
	struct ODBlock RFU_ODBlock();
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks();

	unsigned int label; // deprecated, only left for GenerateInitConfig, please don't use for new code
};

struct compare_interaction {
	bool operator() (Interaction *inter1, Interaction *inter2) const
	{
		auto ij1 = inter1->get_par_num();
		auto ij2 = inter2->get_par_num();

		bool equal = (ij1.first == ij2.first) && (ij1.second == ij2.second);
		if (equal) {
			return inter1 < inter2;
		} else {
			return (ij1.first + ij1.second) < (ij2.first + ij2.second);
		}
	}
};
#endif /* defined(__LF_DEM__Interaction__) */
