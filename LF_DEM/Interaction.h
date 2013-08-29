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
#include "StressTensor.h"

using namespace std;
class System;

class Interaction{
private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	double a0; // radii
	double a1; // second raddi > a0
	double ro; // ro = a0+a1;
	double ro_half; // = ro/2
	//======= internal state =====================//
	bool active;
	unsigned int label;
	unsigned int par_num[2];
	bool contact;
	//======= relative position/velocity data  =========//
	double r; // center-center distance
	int zshift;
	double gap_nondim; // gap between particles (dimensionless gap = s - 2, s = 2r/(a1+a2) )
	double lub_coeff; // = 1/(gap + lub_reduce_parameter)
	double log_lub_coeff; // = log(lub_coeff);
	double lub_coeff_contact; //
	double tangential_dashpot_coeff; //
	vec3d r_vec; // normal vector
	vec3d nr_vec; // vector center to center
	vec3d contact_velocity;
	double normal_relative_velocity;
	vec3d disp_tan; // tangential displacement
	vec3d disp_tan_predictor; // tangential displacement
	vec3d disp_tan_previous;
	//===== forces and stresses ==================== //
	double lub_max_scaled;  // max distance for lubrication
	vec3d lubforce_i; // lubforce_j = - lubforce_i
	double kn_scaled;
	double kt_scaled;
	double colloidalforce_amplitude;
	double colloidalforce_length;
	double XA[4]; // ii ij ji jj
	double YA[4]; // ii ij ji jj
	double YB[4]; // ii ij ji jj
	double YC[4]; // ii ij ji jj
	double XG[4]; // ii ij ji jj
	double YG[4]; // ii ij ji jj
	double YH[4]; // ii ij ji jj
	double XM[4]; // ii ij ji jj
	double YM[4]; // ii ij ji jj
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

	//======= internal state switches  ===========//
	void activate_contact();
	void deactivate_contact();
	
	//=======   ===========//
	void outputSummary();
	
	//===== forces and stresses computations =====//
	double f_contact_normal_norm; // normal contact force
	double f_colloidal_norm;
	vec3d f_contact_normal; // normal contact force
	vec3d f_contact_tan; // tangential contact force
	vec3d f_colloidal;
	StressTensor colloidal_stresslet_XF; //stress tensor of colloidal force
	StressTensor contact_stresslet_XF_normal; //stress tensor of normal contact force
	StressTensor contact_stresslet_XF_tan; //stress tensor of frictional contact force
	void calcContactInteraction();
	void checkBreakupStaticFriction();
	//==========================================================================//
	
	void calcLubConstants();
	double lambda_square;
	double lambda_cubic;
	double lambda_p_1;
	double lambda_p_1_square;
	double lambda_p_1_cubic;
	double cXA[4];
	double cYA[4];
	double cYB[4];
	double cYC[4];
	double cXG[4];
	double cYG[4];
	double cYH[4];
	double cXM[4];
	double cYM[4];
	double a0a0_23;
	double a1a1_23;
	double roro_6;
	double a0a0a0_53;
	double a1a1a1_53;
	double rororo_524;
	double a0a0a0_43;
	double a1a1a1_43;
	double rororo_16;
	double g1_XA;
	double g1_inv_XA;
	double g2_YA;
	double g2_inv_YA;
	double g2_YB;
	double g2_inv_YB;
	double g2_YC;
	double g2_inv_YC;
	double g4_YC;
	double g1_XG;
	double g1_inv_XG;
	double g2_YG;
	double g2_inv_YG;
	double g2_YH;
	double g2_inv_YH;
	double g5_YH;
	double g5_inv_YH;
	double g1_XM;
	double g1_inv_XM;
	double g4_XM;
	double g2_YM;
	double g2_inv_YM;
	double g5_YM;

protected:
public:
	/*********************************
	 *       Public Methods          *
	 *********************************/
	// 	Interaction(){};
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
	void updateStateRelax(bool &deactivated);
	void activate(int i, int j);
	void deactivate();
	inline bool is_overlap(){return r<ro;}
	inline bool is_contact(){return contact;}
	inline bool is_active(){return active;}
	void updateContactModel();
	void calcNormalVectorDistanceGap();
	//======= particles data  ====================//
	inline int
	partner(unsigned int i){
		return (i == par_num[0] ? par_num[1] : par_num[0]);
	}
	inline void
	get_par_num(unsigned int &i, unsigned int &j){
		i = par_num[0], j = par_num[1];
	}
	inline void set_label(unsigned int val){label = val;}
	inline unsigned int get_label(){return label;}
	inline double get_a0(){return a0;}
	inline double get_a1(){return a1;}
	inline void set_ro(double val){
		ro = val;
		ro_half = 0.5*ro;
	}; // ro = a0 + a1
	inline double get_ro(){return ro;}
	//======= relative position/velocity  ========//
	inline double get_r(){return r;}
	inline double get_gap_nondim(){return gap_nondim;}
	inline vec3d get_nr_vec(){return nr_vec;}

	//=============  Resistance Matrices ====================/

	void GE(double *GEi, double *GEj);
	void HE(double *HEi, double *HEj);
	void calcResistanceFunctions();
	void calcXA();
	void calcYA();
	void calcYB();
	void calcYC();
	void calcXG();
	void calcYG();
	void calcXM();
	void calcYM();
	void calcYH();
	inline double get_scaled_XA0(){return a0*XA[0];}
	inline double get_scaled_XA1(){return ro_half*XA[1];}
	inline double get_scaled_XA2(){return ro_half*XA[2];}
	inline double get_scaled_XA3(){return a1*XA[3];}
	inline double get_scaled_YA0(){return a0*YA[0];}
	inline double get_scaled_YA1(){return ro_half*YA[1];}
	inline double get_scaled_YA2(){return ro_half*YA[2];}
	inline double get_scaled_YA3(){return a1*YA[3];}
	inline double get_scaled_YB0(){return a0a0_23*YB[0];}
	inline double get_scaled_YB1(){return roro_6*YB[1];}
	inline double get_scaled_YB2(){return roro_6*YB[2];}
	inline double get_scaled_YB3(){return a1a1_23*YB[3];}
	inline double get_scaled_YC0(){return a0a0a0_43*YC[0];}
	inline double get_scaled_YC3(){return a1a1a1_43*YC[3];}
	inline double get_scaled_YC1(){return rororo_16*YC[1];}
	inline double get_scaled_YC2(){return rororo_16*YC[2];}
	inline double get_scaled_XG0(){return a0a0_23*XG[0];}
	inline double get_scaled_XG1(){return roro_6*XG[1];}
	inline double get_scaled_XG2(){return roro_6*XG[2];}
	inline double get_scaled_XG3(){return a1a1_23*XG[3];}
	inline double get_scaled_YG0(){return a0a0_23*YG[0];}
	inline double get_scaled_YG1(){return roro_6*YG[1];}
	inline double get_scaled_YG2(){return roro_6*YG[2];}
	inline double get_scaled_YG3(){return a1a1_23*YG[3];}
	
	
	
	
	
	//===== forces/stresses  ========================== //
	void calcRelativeVelocities();
	void addUpContactForceTorque();
	void addUpColloidalForce();
	void evaluateLubricationForce();
	double getContactVelocity();
	double getNormalVelocity();
	double getPotentialEnergy();
	inline double get_f_contact_normal_norm(){return f_contact_normal_norm;}
	inline double get_f_contact_tan_norm(){return f_contact_tan.norm();}
	inline double get_f_colloidal_norm(){return f_colloidal_norm;}
	inline double disp_tan_norm(){return disp_tan.norm();}
	inline double getLubForce(){return -dot(lubforce_i, nr_vec);}
	void addHydroStress();
	void addContactStress();
	void addColloidalStress();
	StressTensor getColloidalStressXF(){return colloidal_stresslet_XF;}
	StressTensor getContactStressXF(){return contact_stresslet_XF_normal+contact_stresslet_XF_tan;}
	StressTensor getContactStressXF_normal(){return contact_stresslet_XF_normal;}
	StressTensor getContactStressXF_tan(){return contact_stresslet_XF_tan;}
	void pairVelocityStresslet(const vec3d &vi, const vec3d &vj,
							   StressTensor &stresslet_i, StressTensor &stresslet_j);
	void pairStrainStresslet(StressTensor &stresslet_i, StressTensor &stresslet_j);
	void integrateStress();
	void info(){
		cerr << "contact " << contact << endl;
		cerr << "kn " << kn_scaled << endl;
		cerr << "kn " << colloidalforce_amplitude << endl;
		cerr << "kn " << colloidalforce_length << endl;
	}
	//=========== observables ===============================//
};
#endif /* defined(__LF_DEM__Interaction__) */
