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
#include "vec3d.h"
#include "Contact.h"
#include "Lubrication.h"
#include "RepulsiveForce.h"
#include "StressTensor.h"

class System;

class Interaction{
	friend class Contact;
	friend class RepulsiveForce;
	friend class Lubrication;

private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	double a0; // radii
	double a1; // second raddi > a0
	double ro_12; // ro_12 = ro/2
	double a_reduced; // 1/a_reduced = 1/a0 + 1/a1
	//======= internal state =====================//
	bool active;
	unsigned int label;
	unsigned int p0;
	unsigned int p1;
	//======= relative position/velocity data  =========//
	int zshift;
	double reduced_gap; // gap between particles (dimensionless gap = s - 2, s = 2r/(a1+a2) )
	vec3d relative_velocity;
	vec3d rolling_velocity;
	//===== forces and stresses ==================== //
	double interaction_range;  // max distance
	/*********************************
	 *       Private Methods         *
	 *********************************/
	//===== forces and stresses computations =====//
	inline void set_ro(double val)
	{
		ro = val; // ro = a0 + a1
		ro_12 = ro/2;
	};
	void calcNormalVectorDistanceGap();
	void calcRelativeVelocities();
	void calcRollingVelocities();
	void integrateStress();
	void updateContactState();
	struct ODBlock (Lubrication::*RFU_ODBlock_lub)();
	std::pair<struct DBlock, struct DBlock> (Lubrication::*RFU_DBlocks_lub)();
	//===== forces/stresses  ========================== //
	/* To avoid discontinous change between predictor and corrector,
	 * the change of contact state is informed in updateResiCoeff.
	 */
	bool contact_state_changed_after_predictor;

public:
	Contact contact;
	Lubrication lubrication;
	RepulsiveForce repulsion;
	vec3d relative_surface_velocity;
	double ro; // ro = a0+a1;
	double r; // center-center distance
	vec3d rvec; // vector center to center
	vec3d nvec; // normal vector

	/*********************************
	 *       Public Methods          *
	 *********************************/
	Interaction(){};
	void init(System *sys_);
	//======= state updates  ====================//
	/* Update the follow items:
	 * - r_vec, zshift, _r, and nr_vec
	 * - contact_velocity_tan
	 * - disp_tan
	 * - Fc_normal and Fc_tan
	 * - check breakup of static friction
	 * - State (deactivation, contact)
	 */
	void updateState(bool& deactivated);
	void activate(unsigned int i, unsigned int j, double interaction_range_);
	void deactivate();

	inline bool is_active()
	{
		return active;
	}
	//======= particles data  ====================//
	inline int partner(unsigned int i)
	{
		return (i == p0 ? p1 : p0);
	}
	inline std::pair<unsigned int, unsigned int>	get_par_num()
	{
		return std::make_pair(p0, p1);
	}
	inline void set_label(unsigned int val)
	{
		label = val;
	}
	inline unsigned int get_label()
	{
		return label;
	}
	inline double get_a_reduced()
	{
		return a_reduced;
	}
	inline double get_reduced_gap()
	{
		return reduced_gap;
	}
	inline double get_gap()
	{
		return r-ro;
	}
	bool hasPairwiseResistance();
	double getNormalVelocity();
	double getRelativeVelocity()
	{
		return relative_velocity.norm();
	}
	double getContactVelocity();
	struct ODBlock RFU_ODBlock();


	std::pair<struct DBlock, struct DBlock> RFU_DBlocks();
};
#endif /* defined(__LF_DEM__Interaction__) */
