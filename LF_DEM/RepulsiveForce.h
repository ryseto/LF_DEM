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
#include "Sym2Tensor.h"

class System;
class Interaction;

class RepulsiveForce{
private:
	System* sys;
	Interaction* interaction;
	unsigned int p0;
	unsigned int p1;
	bool vdW;
	//===== forces and stresses ==================== //
	double geometric_factor;
	double screening_length;
	double max_length;
	double cutoff_roundlength;
	vec3d force_vector; // normal contact force
	double force_norm;
	double reduced_force_norm;
	void calcReducedForceNorm();
	void calcScaledForce();
	void calcForce_NottBrady();
public:
	RepulsiveForce():
	p0(0),
	p1(0),
	geometric_factor(0),
	screening_length(0),
	max_length(0),
	cutoff_roundlength(0),
	force_vector(0),
	force_norm(0),
	reduced_force_norm(0)
	{};
	void init(System* sys_, Interaction* int_);
	~RepulsiveForce(){};
	void activate();
	//===== forces/stresses  ========================== //
	void calcForce();
	void addUpForce(std::vector<vec3d> &force) const;
	inline double getForceNorm() const
	{
		return force_norm;
	}
	vec3d getForceVector() const
	{
		return force_vector;
	}
	void addUpStressXF(Sym2Tensor &stress_p0, Sym2Tensor &stress_p1);
	double calcEnergy() const;

};
#endif /* defined(__LF_DEM__RepulsiveForce__) */
