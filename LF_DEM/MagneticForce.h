//
//  MagneticForce.h
//  LF_DEM
//
//  Created by Ryohei Seto on 4/17/15.
//  Copyright (c) 2015 Ryohei Seto. All rights reserved.
//

/**
 \class MagneticForce
 \brief Magnetic dipole interaction
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__MagneticForce__
#define __LF_DEM__MagneticForce__

//#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "StressTensor.h"
using namespace std;
class System;
class Interaction;

class MagneticForce{
private:
	System *sys;
	Interaction *interaction;
	unsigned short p0;
	unsigned short p1;
	//===== forces and stresses ==================== //
	double geometric_factor;
	double length;
	vec3d force_vector; // normal contact force
	double force_norm;
	double reduced_force_norm;
	StressTensor stresslet_XF;
	void calcReducedForceNorm();
	void calcScaledForce();
public:
	MagneticForce(): force_norm(0) {};
	~MagneticForce(){};
	void init(System *sys_, Interaction *int_);
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
};
#endif /* defined(__LF_DEM__MagneticForce__) */
