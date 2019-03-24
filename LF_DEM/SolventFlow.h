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
	double d_tau;
	double pressure_difference;
	// Staggered grid stores
	// - the pressure at the cell center
	// - the velocities at the cell faces.
	std::vector<double> pressure;
	std::vector<double> div_u_sol_ast;
	std::vector<double> u_sol_x;
	std::vector<double> u_sol_z;
	std::vector<double> u_sol_ast_x;
	std::vector<double> u_sol_ast_z;
	std::vector<double> u_particle_x;
	std::vector<double> u_particle_z;
	std::vector<double> u_x;
	std::vector<double> u_z;
	std::vector<double> phi;

	
	std::vector<vec3d> pos;
	int q(int xi, int zi){
		if (xi >= nx) {
			xi -= nx;
		} else if (xi <= -1) {
			xi += nx;
		}
		if (zi >= nz) {
			zi -= nz ;
		} else if (zi <= -1) {
			zi += nz;
		}
		return xi+zi*nx;
	}
	SpMat lap_mat;
	//lap_mat;
	Eigen::VectorXd b;
	Eigen::VectorXd x;
	Eigen::SimplicialLDLT <SpMat> *psolver;
	
	void calcMeshVelocity();
	
public:
	SolventFlow();
	~SolventFlow();
	void init(System* sys_);
	void update(double pressure_difference);
	void calcVelocityDivergence();
	void initPoissonSolver();
	void solvePressure();
	void updateSolventFlow();
	vec3d localFlow(const vec3d &p);
	double meanVelocity();
	void outputYaplot(std::ofstream &fout_flow);
	void velocityProfile(std::ofstream &fout_fp);
};
#endif /* SolventFlow_hpp */
