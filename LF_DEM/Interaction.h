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
#include "Contact.h"
#include "Lubrication.h"
#include "StressTensor.h"

using namespace std;
class System;
class Lubrication;
class Contact;

class Interaction{
	friend class Contact;
	friend class Lubrication;
 private:
	/*********************************
	 *        Members                *
	 *********************************/
	System *sys;
	double a0; // radii
	double a1; // second raddi > a0
	double ro, ro_12; // ro = a0+a1;
	//======= internal state =====================//
	bool active;
	unsigned int label;
	unsigned int par_num[2]; // <<--- We may change this to p0 and p1 (?)
	//======= relative position/velocity data  =========//
	double r; // center-center distance
	int zshift;
	double gap_nondim; // gap between particles (dimensionless gap = s - 2, s = 2r/(a1+a2) )
	vec3d rvec; // vector center to center
	vec3d nvec; // normal vector
	vec3d relative_surface_velocity;
	double nxnx;
	double nxny;
	double nxnz;
	double nynz;
	double nyny;
	double nznz;
	double a0_dash; // radius - 0.5*overlap
	double a1_dash; // second radius > a0
	//===== forces and stresses ==================== //
	double interaction_range_scaled;  // max distance for lubrication
	double colloidalforce_amplitude;
	double colloidalforce_length;
	//===== observables  ========================== //
	double strain_lub_start; // the strain when lubrication object starts.
	double duration; // entire lifetime
	double max_stress; // Maximum value of stress in the all history of this object.
	/*********************************
	 *       Private Methods         *
	 *********************************/

	//=======   ===========//
	void outputSummary();
	//===== forces and stresses computations =====//
	double f_colloidal_norm;
	vec3d f_colloidal;
	StressTensor colloidal_stresslet_XF; //stress tensor of colloidal force
protected:
public:
	Contact contact;
	Lubrication lubrication;
	/*********************************
	 *       Public Methods          *
	 *********************************/
	Interaction(): contact(), lubrication(Lubrication(this)) {;}
	Interaction(const Interaction& obj): contact(), lubrication(Lubrication(this)){
		contact = obj.contact;
		//active = obj.is_active();
	}
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
	inline bool is_contact(){return contact.active;}
	inline bool is_friccontact(){
		if (contact.disp_tan.is_not_zero()){
			return true;
		} else {
			return false;
		}
	}
	inline bool is_active(){return active;}
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
		ro_12 = ro/2;
	}; // ro = a0 + a1
	inline double get_ro(){return ro;}
	//======= relative position/velocity  ========//
	inline double get_r(){return r;}
	inline double get_gap_nondim(){return gap_nondim;}
	inline vec3d get_nvec(){return nvec;}
	double getContactVelocity();

	//===== forces/stresses  ========================== //
	void setResistanceCoeff(double, double);
	void calcRelativeVelocities();
	void addUpColloidalForce();
	double getNormalVelocity();
	double getPotentialEnergy();
	inline double get_f_colloidal_norm(){return f_colloidal_norm;}
	void addColloidalStress();
	void calcTestStress();
	StressTensor getColloidalStressXF(){return colloidal_stresslet_XF;}
	void integrateStress();
	void info(){
		cerr << "particles " << par_num[0] << " " << par_num[1] << endl;
		cerr << "contact " << contact.active << endl;
		cerr << "colloidal force amp " << colloidalforce_amplitude << endl;
		cerr << "colloidal force length " << colloidalforce_length << endl;
		contact.info();
	}
	//=========== observables ===============================//
};
#endif /* defined(__LF_DEM__Interaction__) */
