//
//  Interaction.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__BrownianForce__
#define __LF_DEM__BrownianForce__
#include <iostream>
#include <iomanip>
#include <fstream>
//#include "vec3d.h"
//#include "Interaction.h"
#include "System.h"
#include "cholmod.h"
#include "MersenneTwister.h"


using namespace std;
class System;

class BrownianForce{
 private:
        System *sys;
	cholmod_dense *forces;
        double kb_T, kb_T2;
	MTRand r_gen;

        cholmod_dense *rand_vec;
        cholmod_dense *random_vector();

	cholmod_sparse *L_sparse;
	cholmod_common *c;
	void generate();
 protected:
 public:
        
	// 	BrownianForce(double kb_T, cholmod_dense *brownian_force_vector, cholmod_common *c, System *sys_);
 	BrownianForce(System *sys_);
	~BrownianForce();

	void init(); // once SuiteSparse algebra cholmod_common object is allocated

	//	cholmod_dense* generate();
	void add_to(cholmod_dense*);
	void generate(double*);

};

#endif /* defined(__LF_DEM__BrownianForce__) */
