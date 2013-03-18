//
//  Interaction.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__Interaction__
#define __LF_DEM__Interaction__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "System.h"
#include "common.h"
using namespace std;
class System;

class Interaction{
private:


	/*********************************
	 *        Members                *
	 *********************************/

	System *sys;

	//======= relative position/velocity data  =========//
	double _r; // center-center distance
	int zshift;
	double _gap_nondim; // gap between particles (dimensionless gap = s - 2, s = 2r/(a1+a2) )
	double lub_reduce_parameter; // small cut-off for ksi: lubrication breakdown
	double lub_coeff; // = 1/(gap + lub_reduce_parameter)
	double lub_coeff_max; // = 1/lub_reduce_parameter
	vec3d r_vec; // normal vector
	vec3d contact_velocity;
	vec3d unit_contact_velocity_tan;
	double sqnorm_contact_velocity;
	vec3d disp_tan; // tangential displacement
	
	//===== forces and stresses ==================== //
	double kn; // spring constant for contact force
	double kt; // spring constant for contact force
	double mu; // friction coeffient
	double r_lub_max;  // max distance for lubrication
	vec3d lubforce_i; // lubforce_j = - lubforce_i
	vec3d contact_force_i; // lubforce_j = - lubforce_i
	stresslet lubstresslet;
	stresslet contactstressletXF;
	stresslet contactstresslet2;

	//===== observables  ========================== //
	double init_nearing_time;
	double init_contact_time; 
	bool nearing_on;
	double nearing_gapnd_cutoff;

	/*********************************
	 *       Private Methods         *
	 *********************************/

	//======= particles data  ====================//
	double lambda, invlambda;  // a1/a0 , a0/a1

	//======= relative position/velocity  ========//
	void r(double new_r);
	void calcDistanceNormalVector();
	void calcContactVelocity();
	void incrementContactTangentialDisplacement();
	void checkBreakupStaticFriction();

	//======= internal state switches  ===========//
	bool updateState();
	void activate_contact();
	void deactivate_contact();

	//===== forces and stresses computations =====//
	double Fc_normal_norm; // normal contact force
	vec3d Fc_normal; // normal contact force
	vec3d Fc_tan; // tangential contact force
	vec3d Tc_0; // contact torque on p0
	vec3d Tc_1; // contact torque on p1
	void calcContactInteraction();
	void calcContactStressTermXF();

protected:

public:

	/*********************************
	 *       Public Methods          *
	 *********************************/

 	Interaction(){};
	~Interaction(){};

	void init(System *sys_);

	//======= state updates  ====================//
	/* Update the follow items:
	 * - r_vec, zshift, _r, and nr_vec
	 * - contact_velocity_tan
	 * - disp_tan
	 * - Fc_normal and Fc_tan
	 * - check breakup of static friction
	 * - State (deactivation, contact)
	 */
	bool updateStatesForceTorque();
	void activate(int i, int j, const vec3d &pos_diff, double distance, int zshift);
	void deactivate();
	bool active;

	//======= particles data  ====================//
	int particle_num[2];
	int partner(int);
	double a0, a1; // radii
	double ro; // ro = a0 + a1

	//======= relative position/velocity  ========//
	vec3d nr_vec; // vector center to center
	inline double r(){return _r;}
	inline double gap_nondim(){return _gap_nondim;}


	//======= internal state =====================//
	bool contact;
	bool friction;

	//=============  Resistance Matrices ====================/
	double XA[4]; // ii ij ji jj
	double XG[4]; // ii ij ji jj
	double XM[4]; // ii ij ji jj
	void GE(double *GEi, double *GEj);
	void calcXA();
	void calcXG();
	void calcXM();

	//===== forces/stresses  ========================== //
	void addUpContactForceTorque();

	double normal_force(){return Fc_normal_norm;};
	vec3d tangential_force(){return Fc_tan;};
	void evaluateLubricationForce();
	double valLubForce();
	double lubStresslet(int i){return lubstresslet.elm[i];}	
//	void addLubricationStress();
	void addHydroStress();
	void addContactStress();
//	void addContactStress2();
	void pairVelocityStresslet(const vec3d &vi, const vec3d &vj, stresslet &stresslet_i, stresslet &stresslet_j);
	void pairVelocityStresslet(double* &vel_array, stresslet &stresslet_i, stresslet &stresslet_j);
	void pairStrainStresslet(stresslet &stresslet_i, stresslet &stresslet_j);

	//=========== observables ===============================//
	double nearing_time();
	double contact_time();
	
};


#endif /* defined(__LF_DEM__Interaction__) */
