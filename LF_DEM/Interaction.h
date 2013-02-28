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
//	double gap_nondim_eff;
	double lub_coeff; // = 1/(gap + lub_reduce_parameter)
	double lub_coeff_max; // = 1/lub_reduce_parameter
	vec3d r_vec; // normal vector
	vec3d contact_velocity;
	vec3d contact_velocity_tan;
	vec3d unit_contact_velocity_tan;
	double sqnorm_contact_velocity;
	vec3d disp_tan; // tangential displacement

	//===== forces and stresses ==================== //
	double kn; // spring constant for contact force
	double kt; // spring constant for contact force
	double r_lub_max;  // max distance for lubrication
	vec3d lubforce_i; // lubforce_j = - lubforce_i
	vec3d contact_force_i; // lubforce_j = - lubforce_i
	stresslet lubstresslet;
	stresslet contactstresslet;

	//===== observables  ========================== //
	double strain_lub_generated; // The strain when this interaction starts.
	double strain_near_contact; // The strain when this interaction starts.
	double lub_time;
	//	double nearing_time; // unused

	/*********************************
	 *       Private Methods         *
	 *********************************/

	//======= particles data  ====================//
	double lambda, invlambda;  // a1/a0 , a0/a1

	//======= relative position/velocity  ========//
	void r(double new_r);
//	void calcNormalVector();
	void calcDistanceNormalVector();
//	void assignDistanceNormalVector(const vec3d &, double, int);
	void calcContactVelocity();
	void incrementContactTangentialDisplacement();

	//======= internal state switches  ===========//
	bool updateState();
	void deactivate();
	void activate_contact();
	void deactivate_contact();

	//===== forces and stresses computations =====//
	double Fc_normal; // normal contact force
	vec3d Fc_tan; // tangential contact force
	vec3d Tc_0; // contact torque on p0
	vec3d Tc_1; // contact torque on p1
	void calcStaticFriction();
	void calcDynamicFriction();
	void calcContactInteraction();

protected:

public:

	/*********************************
	 *       Public Methods          *
	 *********************************/

 	Interaction(){};
	~Interaction(){};

	void init(System *sys_);

	//======= state updates  ====================//
	bool update();
	void activate(int i, int j, const vec3d &pos_diff, double distance, int zshift);
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
	bool static_friction;

	//=============  Resistance Matrices ====================/
	void XA(double &XAii, double &XAij, double &XAji, double &XAjj);
	void XG(double &XGii, double &XGij, double &XGji, double &XGjj);
	void XM(double &XMii, double &XMij, double &XMji, double &XMjj);
	void GE(double GEi[], double GEj[]);


	//===== forces/stresses  ========================== //
//	void addUpContactForce(vec3d &force0, vec3d &force1);
//	void addUpContactTorque(vec3d &torque0, vec3d &torque1);
	void addUpContactForceTorque();
	double normal_force(){return Fc_normal;};
	vec3d tangential_force(){return Fc_tan;};
	void evaluateLubricationForce();
	double valLubForce();
	double lubStresslet(int i){return lubstresslet.elm[i];}	
	void addLubricationStress();
	void addContactStress();
	void pairVelocityStresslet(const vec3d &vi, const vec3d &vj, stresslet &stresslet_i, stresslet &stresslet_j);
	void pairStrainStresslet(stresslet &stresslet_i, stresslet &stresslet_j);

	//===== other observables  ========================== //
	double nearingTime();
	bool near;
	vector <vec3d> trajectory;
	vector <double> gap_history;
	void recordTrajectory();
	void outputTrajectory();
};


#endif /* defined(__LF_DEM__Interaction__) */
