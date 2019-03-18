//
//  SolventFlow.h
//  LF_DEM
//
//  Created by Ryohei Seto on 2019/03/16.
//  Copyright © 2019 Ryohei Seto. All rights reserved.
//

#ifndef SolventFlow_h
#define SolventFlow_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
#include "cholmod.h"
#include "cholmod_check.h"

#include "vec3d.h"
#include "Contact.h"
#include "Lubrication.h"
#include "RepulsiveForce.h"
#include "TimeActivatedAdhesion.h"
#include "Sym2Tensor.h"
#define DELETE(x) if(x){delete [] x; x = NULL;}
// uncomment below to use long integer cholmod (necessary for GPU)
#ifndef USE_GPU
#define CHOL_FUNC(NAME) cholmod_ ## NAME
typedef int chol_int;
#else
#define CHOL_FUNC(NAME) cholmod_l_ ## NAME
typedef long chol_int;
#endif

using namespace std;

class System;
class SolventFlow {
	
private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	int nx;
	int nz;
	int n;
	double dx;
	double dz;
	vector<double> Pressure;
	vector<double> div_u_particle;
	vector<double> u_solvent_x;
	vector<double> u_solvent_z;
	vector<double> u_particle_x;
	vector<double> u_particle_z;
	vector<vec3d> pos;

	// Staggered grid stores
	// - the pressure at the cell center
	// - the velocities at the cell faces.
	
	// Cholmod variables
	cholmod_factor* chol_L ;
	cholmod_common chol_c ;
	cholmod_dense* chol_rhs;
	cholmod_sparse* chol_matrix;
	cholmod_dense* chol_solution;
	cholmod_dense* chol_solve_workspace;
	// cholmod_dense* chol_PTsolution;
	//	cholmod_dense* chol_Psolution;
	// resistance matrix building
	chol_int dblocks_size;
	chol_int current_column;
	
	
	
public:
	SolventFlow();
	~SolventFlow();
	void init(System* sys_);
	void updateParticleVelocity();
	void calcVelocityDivergence();
	
	
	void outputYaplot(std::ofstream &fout_flow);

};
#endif /* SolventFlow_hpp */
