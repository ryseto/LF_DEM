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
	vec3d force_vector0; // normal contact force
	vec3d torque0;
	vec3d torque1;
	
	double force_norm;
	StressTensor stresslet_XF;
public:
	MagneticForce(): force_norm(0) {};
	~MagneticForce(){};
	void init(System *sys_, Interaction *int_);
	void activate();
	//===== forces/stresses  ========================== //
	void calcForceToruqe();
	void addUpForceTorque();
	inline double getForceNorm()
	{
		return force_norm;
	}
	vec3d getForceVector()
	{
		return force_vector0;
	}
	void calcStressXF();
	StressTensor getStressXF()
	{
		return stresslet_XF;
	}
};
#endif /* defined(__LF_DEM__MagneticForce__) */
