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
#include <stdexcept>
#include "vec3d.h"
#include "Eigen/Sparse"
#include "Eigen/Core"

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
	int jmax_uz;
	double dx;
	double dz;
	double ux_bot;
	double ux_top;
	double d_tau;
	double smooth_length;
	double sq_smooth_length;
	double cell_area;
	double numerical_Re;
	double target_flux;
	double viscosity;

	bool settling;
	bool channel_flow;
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
	std::vector<double> omega;
	std::vector<double> strain_rate_xx;
	std::vector<double> strain_rate_xz;
	std::vector<double> strain_rate_zz;
	std::vector<double> phi_ux;
	std::vector<double> phi_uz;
	std::vector<vec3d> pos;
	std::vector <int> mesh_nb;
	std::vector <double> udx_values;
	std::vector <double> udz_values;
	std::vector <double> phi_ux_values;
	std::vector <double> phi_uz_values;
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
	void calcOmega();
public:
	SolventFlow();
	~SolventFlow();
	double pressure_difference;
	void init(System* sys_, std::string simulation_type);
	void update(double pressure_difference);
	void initPoissonSolver();
	void localFlow(const vec3d &p, vec3d &u_local, vec3d &omega_local,
				   std::vector<double> &e_local);
	double meanVelocity();
	void outputYaplot(std::ofstream &fout_flow);
	void velocityProfile(std::ofstream &fout_fp);
	double calcFlux();
};
#endif /* SolventFlow_hpp */
