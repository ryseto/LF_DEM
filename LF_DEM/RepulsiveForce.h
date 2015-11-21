//
//  RepulsiveForce.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/11/15.
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class RepulsiveForce
 \brief Exponentially decaying repulsion, to be called from an Interaction object
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__RepulsiveForce__
#define __LF_DEM__RepulsiveForce__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "StressTensor.h"

class System;
class Interaction;

class RepulsiveForce{
private:
	System* sys;
	Interaction* interaction;
	unsigned short p0;
	unsigned short p1;
	//===== forces and stresses ==================== //
	double geometric_factor;
	double screening_length;
	double max_length;
	double cutoff_roundlength;
	vec3d force_vector; // normal contact force
	double force_norm;
	double reduced_force_norm;
	StressTensor stresslet_XF;
	void calcReducedForceNorm();
	void calcScaledForce();
public:
	RepulsiveForce(): force_norm(0) {};
	~RepulsiveForce(){};
	void init(System* sys_, Interaction* int_);
	void activate();
	//===== forces/stresses  ========================== //
	void calcForce();
	void addUpForce();
	inline double getForceNorm()
	{
		return force_norm;
	}
	vec3d getForceVector()
	{
		return force_vector;
	}
	void calcStressXF();
	StressTensor getStressXF()
	{
		return stresslet_XF;
	}
	double calcEnergy();

};
#endif /* defined(__LF_DEM__RepulsiveForce__) */
