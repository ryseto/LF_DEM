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
	unsigned int i; // <-- p0?
	unsigned int j; // <-- p1?

	//===== forces and stresses ==================== //
	double kn_scaled;
	double kt_scaled;
	double mu;
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
protected:
public:
	/*********************************
	 *       Public Methods          *
	 *********************************/
	//======= internal state =====================//

	vec3d disp_tan; // tangential displacement
	vec3d prev_disp_tan; // useful for predictor-corrector method: disp_tan in the previous time step

	Contact(){};
	Contact(const Contact& obj){
		disp_tan = obj.disp_tan;
	}
	void init(System *sys_, Interaction *int_);
	void getInteractionData();
	void activate();
	void deactivate();
	bool active;
	bool staticfriction;
	void updateContactModel();
	void frictionlaw_criticalload();
	void frictionlaw_criticalload_mu_inf();
	void frictionlaw_coulomb();
	void frictionlaw_null();
	//===== forces/stresses  ========================== //

	void incrementTangentialDisplacement();
	void calcContactInteraction();
	void calcContactInteractionRelax();
	void addUpContactForceTorque();
	//double getContactVelocity();
	inline double get_f_contact_normal_norm(){return f_contact_normal_norm;}
	inline double get_f_contact_tan_norm(){return f_contact_tan.norm();}
	inline double disp_tan_norm(){return disp_tan.norm();}
	void addContactStress();
	StressTensor getContactStressXF(){return contact_stresslet_XF_normal+contact_stresslet_XF_tan;}
	StressTensor getContactStressXF_normal(){return contact_stresslet_XF_normal;}
	StressTensor getContactStressXF_tan(){return contact_stresslet_XF_tan;}
	bool is_activated_friction() {
		if (disp_tan.is_not_zero()){
			return true;
		} else {
			return false;
		}
	}
	void info(){
		cerr << "kn " << kn_scaled << endl;
	}
	vec3d get_disp_tan(){return disp_tan;}

};
#endif /* defined(__LF_DEM__Contact__) */
