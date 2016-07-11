//
//  ContactDashpot.h
//  LF_DEM
//
//  Copyright (c) 2016 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class ContactDashpot
 \brief ContactDashpot interaction, to be called from an Interaction object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__ContactDashpot__
#define __LF_DEM__ContactDashpot__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "StressTensor.h"
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
	double range;

	//======= particles data  ====================//
	unsigned int p0;
	unsigned int p1;
	unsigned int p0_6;
	unsigned int p1_6;
	vec3d *nvec;
	double normal_coeff;
	double tangential_coeff;
	double a0;
	double a1;
	double ro;
	double XA[4]; // ii ij ji jj
	double YA[4]; // ii ij ji jj
	double YB[4]; // ii ij ji jj
	double YC[4]; // ii ij ji jj
	double ro_12; // = ro/2
	void calcDashpotResistances();

 public:
	 ContactDashpot();
	void init(System *sys_, Interaction *int_);
	inline bool is_active() {return _active;};
	void activate();
	void deactivate();
	void setParticleData();
	//===== forces/stresses  ========================== //
    vec3d lubforce_p0; // lubforce_p1 = - lubforce_p0
	void calcPairwiseForce();
	double get_lubforce_normal()
	{
		// positive for compression
        //lubforce_p0.cerr();
		return -dot(lubforce_p0, nvec);
	}
	vec3d get_lubforce_tan()
	{
		return lubforce_p0-dot(lubforce_p0, nvec)*(*nvec);
	}
	vec3d get_lubforce()
	{
		return lubforce_p0;
	}
	void setDashpotResistanceCoeffs(double normal_rc, double tangent_rc);
//void setResistanceCoeffTang(double tangent_rc);
	//=============  Resistance Matrices ====================/
	struct ODBlock RFU_ODBlock();
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks();
};
#endif /* defined(__LF_DEM__ContactDashpot__) */
