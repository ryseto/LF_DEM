//
//  BrownianForce.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__BrownianForce__
#define __LF_DEM__BrownianForce__
#include "System.h"
#include "cholmod.h"
#include "MersenneTwister.h"
using namespace std;
class System;

class BrownianForce{
private:
	System *sys;
	double *forces;
	double kb_T, kb_T2;
	MTRand r_gen;
	double **pair_resistance_matrix;
	double **L_factor;
	double *ran_vector;
	void factorize();
	void addToDiagBlock_RFU(const vec3d &nvec, int ii, double alpha);
	void addToOffDiagBlock_RFU(const vec3d &nvec, double alpha);
protected:
public:
	
	// 	BrownianForce(double kb_T, cholmod_dense *brownian_force_vector, cholmod_common *c, System *sys_);
 	BrownianForce(System *sys_);
	~BrownianForce();
	void init(); // once SuiteSparse algebra cholmod_common object is allocated
	void add_to(cholmod_dense*);
	//	cholmod_dense* generate();
	double* generate();
	double* generate_invLFb();
	double* generate_new_2();
	void generate(double*);
	
};

#endif /* defined(__LF_DEM__BrownianForce__) */
