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
#include "System.h"
#include "Contact.h"
#include "Lubrication.h"
#include "RepulsiveForce.h"
#include "StressTensor.h"

using namespace std;
class System;
class Lubrication;
class Contact;
class RepulsiveForce;

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
	double ro; // ro = a0+a1;
	double ro_12; // ro_12 = ro/2
	double a_reduced;
	double c_rolling_veolocity;
	bool a0_eq_a1;
	//======= internal state =====================//
	bool active;
	unsigned int label;
	unsigned short p0;
	unsigned short p1;
	//======= relative position/velocity data  =========//
	double r; // center-center distance
	int zshift;
	double reduced_gap; // gap between particles (dimensionless gap = s - 2, s = 2r/(a1+a2) )
	vec3d rvec; // vector center to center
	vec3d nvec; // normal vector
	double nxnx;
	double nxny;
	double nxnz;
	double nynz;
	double nyny;
	double nznz;
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
	/*********************************
	 *       Public Methods          *
	 *********************************/
	Interaction(): contact(), lubrication(Lubrication(this)) {}
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
	void updateState(bool &deactivated);
	void activate(unsigned short i, unsigned short j, double range);
	void deactivate();
	inline bool is_overlap()
	{
		return r < ro;
	}
	inline bool is_contact()
	{
		return contact.state >= 1;
	}
	inline bool is_friccontact()
	{
		return contact.state >= 2;
	}
	inline bool is_active()
	{
		return active;
	}
	//======= particles data  ====================//
	inline int partner(unsigned int i)
	{
		return (i == p0 ? p1 : p0);
	}
	inline void	get_par_num(unsigned short &i, unsigned short &j)
	{
		i = p0, j = p1;
	}
	inline void set_label(unsigned int val)
	{
		label = val;
	}
	inline unsigned int get_label()
	{
		return label;
	}
	inline double get_a0()
	{
		return a0;
	}
	inline double get_a1()
	{
		return a1;
	}
	inline double get_ro()
	{
		return ro;
	}
	//======= relative position/velocity  ========//
	inline double get_r()
	{
		return r;
	}
	inline double get_reduced_gap()
	{
		return reduced_gap;
	}
	inline double get_gap()
	{
		return r-ro;
	}
	inline vec3d get_nvec()
	{
		return nvec;
	}
	double getNormalVelocity();
	double getRelativeVelocity()
	{
		return relative_velocity.norm();
	}
	double getContactVelocity();
};
#endif /* defined(__LF_DEM__Interaction__) */
