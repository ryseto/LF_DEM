//
//  Contact.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__Contact__
#define __LF_DEM__Contact__

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
	
	
	void (Contact::*frictionlaw)();

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
	/*********************************
	 *       Private Methods         *
	 *********************************/
	//===== forces and stresses computations =====//
	StressTensor contact_stresslet_XF_normal; //stress tensor of normal contact force
	StressTensor contact_stresslet_XF_tan; //stress tensor of frictional contact force

	double f_contact_normal_norm; // normal contact force
	vec3d f_contact_normal; // normal contact force
	vec3d f_contact_tan; // tangential contact force
	vec3d old_relative_velocity;
	vec3d old_f_test_vec;
	vec3d old_lubforce_tan;
	vec3d old_dashpot;
	vec3d old_spring;
	bool old_state;
	vec3d tvec;
//	double supportable_tanforce;

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
	void frictionlaw_criticalload();
	void frictionlaw_coulomb();
	void frictionlaw_null();

//	void frictionlaw_legacy();
	//===== forces/stresses  ========================== //

	void incrementTangentialDisplacement();
	void calcContactInteraction();
	void calcContactInteractionRelax();
	void addUpContactForceTorque();
	double getContactVelocity();
	inline double get_f_contact_normal_norm(){return f_contact_normal_norm;}
	inline double get_f_contact_tan_norm(){return f_contact_tan.norm();}
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

	vec3d get_disp_tan(){return disp_tan;}

};
#endif /* defined(__LF_DEM__Contact__) */
