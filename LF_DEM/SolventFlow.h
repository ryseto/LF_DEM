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
#include "vec3d.h"
#include "Eigen/Sparse"
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

//using namespace Eigen;

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
	// Staggered grid stores
	// - the pressure at the cell center
	// - the velocities at the cell faces.
	std::vector<double> pressure;
	std::vector<double> div_u_particle;
	std::vector<double> u_solvent_x;
	std::vector<double> u_solvent_z;
	std::vector<double> u_particle_x;
	std::vector<double> u_particle_z;
	std::vector<vec3d> pos;
	int q(int xi, int zi){
		int i = xi+zi*nx;
		return i;
	}
	SpMat lap_mat;
	//lap_mat;
	Eigen::VectorXd b;
	Eigen::VectorXd x;
	Eigen::SimplicialLDLT <SpMat> *psolver;
	
public:
	SolventFlow();
	~SolventFlow();
	void init(System* sys_);
	void updateParticleVelocity();
	void calcVelocityDivergence();
	void initPoissonSolver();
	void solvePressure();
	
	void outputYaplot(std::ofstream &fout_flow);

};
#endif /* SolventFlow_hpp */
