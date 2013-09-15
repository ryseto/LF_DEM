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

	void applyFrictionLaw_spring();
	void applyFrictionLaw_spring_dashpot();
	void imposeFrictionLaw_spring();
	double f_contact_normal_norm; // normal contact force
	vec3d f_contact_normal; // normal contact force
	vec3d f_contact_tan; // tangential contact force
	vec3d tvec;
	double supportable_tanforce;
	vec3d dashpot;
	vec3d lubforce_tan;
	
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
	inline double get_cnt_sliding(){
		return cnt_sliding;
	}

};
#endif /* defined(__LF_DEM__Contact__) */
