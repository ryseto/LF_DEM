//
//  Contact.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__Contact__
#define __LF_DEM__Contact__
//#define RECORD_HISTORY 1

#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "StressTensor.h"

using namespace std;
class System;
class Interaction;

class Contact{
private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	Interaction *interaction;

	unsigned int i,j;
	//======= internal state =====================//
	vec3d disp_tan; // tangential displacement
	vec3d disp_tan_previous;
	//===== forces and stresses ==================== //
	double kn_scaled;
	double kt_scaled;
	double mu;
	//===== observables  ========================== //
	double strain_contact_start; // the strain at h=0.
	double duration_contact; // entire duration for h < 0
	int cnt_sliding;  // to count the number of slips.
#ifdef RECORD_HISTORY
	vector <double> disp_tan_sq_history;
	void outputHistory();
#endif
	/*********************************
	 *       Private Methods         *
	 *********************************/
	//===== forces and stresses computations =====//
	StressTensor contact_stresslet_XF_normal; //stress tensor of normal contact force
	StressTensor contact_stresslet_XF_tan; //stress tensor of frictional contact force

	double f_contact_normal_norm; // normal contact force
	vec3d f_contact_normal; // normal contact force
	vec3d f_contact_tan; // tangential contact force
	vec3d tvec;
	double supportable_tanforce;
	vec3d resforce_tan;
//	vec3d lubforce_tan;
	vec3d f_test_vec;
	double previous_f_test;
	double previous_supportable_tanforce;
	int just_switched = -1;

protected:
public:
	/*********************************
	 *       Public Methods          *
	 *********************************/
	Contact(){};
	void init(System *sys_, Interaction *int_);
	void getInteractionData();
	void activate();
	void deactivate();
	bool active;
	bool staticfriction;
	void updateContactModel();
	void resetObservables();
	void frictionlaw();
	void frictionlaw_legacy();
	//===== forces/stresses  ========================== //

	void incrementTangentialDisplacement();
	void calcContactInteraction();
	void calcContactInteractionRelax();
	void addUpContactForceTorque();
	double getContactVelocity();
	inline double get_f_contact_normal_norm(){return f_contact_normal_norm;}
	inline double get_f_contact_tan_norm(){return f_contact_tan.norm();}
	inline vec3d get_f_test(){return f_test_vec;}
	
	inline double disp_tan_norm(){return disp_tan.norm();}
	void addContactStress();
	StressTensor getContactStressXF(){return contact_stresslet_XF_normal+contact_stresslet_XF_tan;}
	StressTensor getContactStressXF_normal(){return contact_stresslet_XF_normal;}
	StressTensor getContactStressXF_tan(){return contact_stresslet_XF_tan;}
	void info(){
		cerr << "kn " << kn_scaled << endl;
	}
	inline double get_duration(){
		return duration_contact;
	}
	inline double get_cnt_sliding(){
		return cnt_sliding;
	}
	vec3d get_disp_tan(){return disp_tan;}

};
#endif /* defined(__LF_DEM__Contact__) */
