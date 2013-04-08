//
//  Interaction.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__Interaction__
#define __LF_DEM__Interaction__
//#define RECORD_HISTORY 1

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
	double lub_coeff; // = 1/(gap + lub_reduce_parameter)
	vec3d r_vec; // normal vector
	vec3d contact_velocity;
	vec3d unit_contact_velocity_tan;
	double sqnorm_contact_velocity;
	vec3d disp_tan; // tangential displacement
	vec3d disp_tan_predictor; // tangential displacement
	//	vec3d nr_vec_before_predictor;
	//===== forces and stresses ==================== //
	double r_lub_max;  // max distance for lubrication
	vec3d lubforce_i; // lubforce_j = - lubforce_i
	stresslet lubstresslet;
	double colloidal_force_amplitude;
	//===== observables  ========================== //
	double strain_lub_start; // the strain when lubrication object starts.
	double strain_contact_start; // the strain at h=0.
	double duration; // entire lifetime
	double duration_contact; // enture duraction for h < 0
	double max_stress; // Maximum value of stress in the all history of this object.
	int cnt_sliding;  // to count the number of slips.
#ifdef RECORD_HISTORY
	vector <double> gap_history;
	vector <double> overlap_history;
	vector <double> disp_tan_sq_history;
	void outputHistory();
#endif
	/*********************************
	 *       Private Methods         *
	 *********************************/

	//======= particles data  ====================//
	double lambda; // a1/a0
	double invlambda; // a0/a1

	//======= relative position/velocity  ========//
	void r(const double &new_r);
	void calcDistanceNormalVector();
	void calcContactVelocity();
	void checkBreakupStaticFriction();

	//======= internal state switches  ===========//
	void activate_contact();
	void deactivate_contact();
	
	//=======   ===========//
	void outputSummary();

	//===== forces and stresses computations =====//
	double Fc_normal_norm; // normal contact force
	double F_colloidal_norm;
	vec3d Fc_normal; // normal contact force
	vec3d Fc_tan; // tangential contact force
	vec3d Tc_0; // contact torque on p0
	vec3d Tc_1; // contact torque on p1
	vec3d F_colloidal;
	void calcContactInteraction();
	void calcStressTermXF(stresslet &stresslet_,
						  const vec3d &force);
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
	void updateState(bool &deactivated);
	void activate(int i, int j);
	void deactivate();

	//======= particles data  ====================//
	int par_num[2];
	int partner(int);
	double a0; // radii
	double a1; // second raddi > a0
	double ro; // ro = a0 + a1
	double ro_2; // = ro/2
	int label;

	//======= relative position/velocity  ========//
	vec3d nr_vec; // vector center to center
	inline double r(){return _r;}
	inline double gap_nondim(){return _gap_nondim;}
	inline double overlap(){
		if (contact){
			return ro-_r;
		} else {
			return 0;
		}
	}
	//======= internal state =====================//
	bool active;
	bool contact;

	//======= Data ===============================//
	double total_stress_xz;
	double stress_xz_integration;
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
	void addUpColloidalForce();
	double normal_force(){return Fc_normal_norm;}
	double colloidal_force(){return F_colloidal_norm;}
	vec3d tangential_force(){return Fc_tan;}
	double disp_tan_norm(){return disp_tan.norm();}
	void evaluateLubricationForce();
	double valLubForce();
	double lubStresslet(int i){return lubstresslet.elm[i];}
	double getContactVelocity();
	double getNormalVelocity();
	double calcPotentialEnergy();
	void addHydroStress();
	void addContactStress();
	void addColloidalStress();
	void pairVelocityStresslet(const vec3d &vi, const vec3d &vj,
							   stresslet &stresslet_i, stresslet &stresslet_j);
	void pairVelocityStresslet(double* &vel_array, stresslet &stresslet_i, stresslet &stresslet_j);
	void pairStrainStresslet(stresslet &stresslet_i, stresslet &stresslet_j);
	
	void integrateStress();

	//=========== observables ===============================//

};
#endif /* defined(__LF_DEM__Interaction__) */
