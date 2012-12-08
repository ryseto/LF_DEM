//
//  System.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__System__
#define __LF_DEM__System__
#define CHOLMOD 1
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include <string>
#include <Accelerate/Accelerate.h>
#ifdef CHOLMOD
#include "cholmod.h"
#endif
#include "vec3d.h"
#include "ContactForce.h"

using namespace std;
class Interaction;

class System{
private:
	int n;
	int n3;
	double dx;
	double dy;
	double dz;
#ifdef CHOLMOD
	cholmod_sparse *sparse_res;
	cholmod_dense *v, *rhs_b;
	cholmod_common c ;
	int max_lub_int;
	int stype;
	int sorted;
	int packed;
	int xtype;
	vector <int> rows;
	double *diag_values;
	vector <double> *off_diag_values;
	int *ploc;
	void fillSparseResmatrix();
	void addToDiag(double *nvec, int ii, double alpha);
	void appendToColumn(double *nvec, int jj, double alpha);
#else
	double *res;
	int nrhs;
	int *ipiv;
	int lda;
	int ldb;
	int info;
	double *rhs_b;
	int lwork;
	double *work;
	char UPLO;
#endif
protected:
public:
    /* For DEMsystem
     */
	System(){};
	~System();
	int ts; // time steps
	int dimension;
	vec3d *position;
	int **i_position;
	double *angle; // for 2D visualization
	vec3d *velocity;
	vec3d *ang_velocity;
	vec3d *force;
	vec3d *torque;
	double kn;
	double kt;
	double eta;
	double sq_lub_max;
	double mu_static; // static friction coefficient.
	double mu_dynamic;// dynamic friction coefficient.
	double dynamic_friction_critical_velocity;
	double sq_critical_velocity;
	bool friction;
	bool lub;
	double lubcore;  // = 2a'
	/*************************************************************/
	double lx;
	double ly;
	double lz;
	double lx2; // =lx/2
	double ly2; // =ly/2
	double lz2; // =lz/2
	double shear_disp;
	double shear_rate;
	double volume_fraction;
	double vel_difference;
	double dt;
	vector <int> lubparticle;
	vector <double> lubparticle_vec[3];
#ifdef CHOLMOD
	cholmod_factor *L ;
#endif
	string simu_name;
	/*************************************************************/
	void prepareSimulation(unsigned long number_of_particles);
	void init();
	double sq_norm();
	double sq_distance(vec3d &pos , int i);
	double sq_distance(int i, int j);
	double sq_distanceToCheckContact(int i, int j);
	double sq_neardistance(int i, int j);
	double lubricationForceFactor(int i, int j);
	void displacement(int i, const double &dx_, const double &dy, const double &dz);

	double distance(int i, int j);
	void updateVelocity();
	void updateVelocityLubrication();
	void deltaTimeEvolution();
	void forceReset();
	void torqueReset();
	bool noOverlap();
		
};
#endif /* defined(__LF_DEM__State__) */
