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
#include "BrownianForce.h"

using namespace std;
class Interaction;
class BrownianForce;

class System{
private:
	int n3;
	double dx;
	double dy;
	double dz;
#ifdef CHOLMOD
	cholmod_sparse *sparse_res;
	cholmod_dense *v, *rhs_b;
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
	void buildLubricationTerms();
	void buildBrownianTerms();
	void buildContactTerms();

protected:
public:
    /* For DEMsystem
     */
	System(){};
	~System();
	int n; // number of particles
	int ts; // time steps
	int dimension;
	vec3d *position;
	double *angle; // for 2D visualization
	vec3d *velocity;
	vec3d *ang_velocity;
	vec3d *force;
	vec3d *torque;
	double **stress; // S_xx S_xy S_xz S_yz S_yy
	double mean_stress[5];
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
	bool brownian;
	/*
	 * Leading term of lubrication force is 1/(r-2a).
	 * This can be weakened by using a'<a.
	 * lubcore = 2a'.
	 * lubcore = 2 means full lubrication force.
	 * lubcore < 2 gives weaker lubriaction force that allows particle contact.
	 */
	double lubcore;
	BrownianForce *fb;
	/*************************************************************/

	double lx;
	double ly;
	double lz;
	double lx2; // =lx/2
	double ly2; // =ly/2
	double lz2; // =lz/2
	double shear_disp;
	double shear_rate;
	double kb_T;
	double volume_fraction;
	double vel_difference;
	double dt;
	bool draw_rotation_2d;
	vector <int> lubparticle;
	vector <double> lubparticle_vec[3];

	string simu_name;
	void prepareSimulationName();
	void prepareSimulation();
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
	void stressReset();
	bool noOverlap();
	void calcHydrodynamicStress();
	void calcStressAverage();

	int numpart(){
	  return n;
	}


#ifdef CHOLMOD
	cholmod_factor *L ;
	cholmod_common c ;
	
#endif
	
	void lubricationStress(int i, int j);


};
#endif /* defined(__LF_DEM__State__) */
