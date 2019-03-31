//
//  ContactDashpot.h
//  LF_DEM
//
//  Copyright (c) 2016 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class ContactDashpot
 \brief ContactDashpot interaction, to be called from a Contact object

	Contains both normal and tangential dashpots.
	As for the force, it formally acts as a distance independant Lubrication object.
	The stresses, however, have to be computed as -xF, from the Contact object.

 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__ContactDashpot__
#define __LF_DEM__ContactDashpot__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Sym2Tensor.h"
#include "MatrixBlocks.h"

class System;
class Interaction;

class ContactDashpot{
private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	Interaction *interaction;
	bool _active;
	//======= particles data  ====================//
	unsigned int p0;
	unsigned int p1;
	unsigned int p0_6;
	unsigned int p1_6;
	double a0;
	double a1;
	double ro;
	double ro_12; // = ro/2
	vec3d *nvec;
	double normal_coeff;
	double tangential_coeff;
	double XA[4]; // ii ij ji jj
	double YA[4]; // ii ij ji jj
	double YB[4]; // ii ij ji jj
	double YC[4]; // ii ij ji jj
	void calcDashpotResistances();

 public:
	ContactDashpot();

	void init(System *sys_, Interaction *int_);
	bool is_active() const {return _active;};
	void activate();
	void deactivate();
	void setParticleData();
	//===== forces/stresses  ==========================
	std::tuple<vec3d, vec3d, vec3d, vec3d> getRFU_Uinf(const vec3d &u_inf_p0,
													   const vec3d &u_inf_p1,
													   const vec3d &omega_inf) const;
	std::tuple<vec3d, vec3d, vec3d, vec3d> getRFU_Ulocal(const vec3d &u_local_p0,
														 const vec3d &u_local_p1,
														 const vec3d &omega_local_p0,
														 const vec3d &omega_local_p1) const;
	vec3d getForceOnP0(const vec3d &vel_p0,
					   const vec3d &vel_p1,
					   const vec3d &ang_vel_p0,
					   const vec3d &ang_vel_p1) const;
	vec3d getForceOnP0_nonaffine(const vec3d &na_vel_p0,
								 const vec3d &na_vel_p1,
								 const vec3d &na_ang_vel_p0,
								 const vec3d &na_ang_vel_p1) const;
	void setDashpotResistanceCoeffs(double kn, double kt,
									double rtime_normal, double rtime_tan);
	//=============  Resistance Matrices ====================/
	struct ODBlock RFU_ODBlock() const;
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks() const;
};
#endif /* defined(__LF_DEM__ContactDashpot__) */
