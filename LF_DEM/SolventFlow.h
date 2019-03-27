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
	double smooth_length;
	double sq_smooth_length;
	double cell_area;
	double numerical_Re;
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
	std::vector<double> u_diff_x;
	std::vector<double> u_diff_z;
	std::vector<double> phi;
	std::vector<vec3d> pos;
	std::vector <int> u_mesh_nb;
	std::vector <int> phi_mesh_nb;
	std::vector <double> udx_values;
	std::vector <double> udz_values;
	std::vector <double> phi_values;
	
	
	int meshNb(int xi, int zi);
	SpMat lap_mat;
	//lap_mat;
	Eigen::VectorXd b;
	Eigen::VectorXd x;
	Eigen::SimplicialLDLT <SpMat> *psolver;
	
	void particleVelocityDiffToMesh();
	double weightFunc(double r_sq);
	void predictorStep();
	void calcVelocityDivergence();
	void solvePressure();
	void correctorStep();
	double porousResistance(double volume_fraction);
public:
	SolventFlow();
	~SolventFlow();
	double pressure_difference;
	void init(System* sys_);
	void update(double pressure_difference);
	void initPoissonSolver();
	vec3d localFlow(const vec3d &p);
	double meanVelocity();
	void outputYaplot(std::ofstream &fout_flow);
	void velocityProfile(std::ofstream &fout_fp);
	double calcFlux();
};
#endif /* SolventFlow_hpp */
