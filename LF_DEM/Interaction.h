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
using namespace std;
class System;

class Interaction{
private:
	System *sys;
	vec3d contact_velocity;
	vec3d contact_velocity_tan;
	vec3d unit_contact_velocity_tan;
	vec3d xi; // tangential displacement

	void calcDistanceNormalVector();
	void assignDistanceNormalVector(const vec3d &, double, int);
	void calcContactVelocity();

	void incrementContactTangentialDisplacement();
	//	vec3d t_tangent;
	int pd_z;
	double sqnorm_contact_velocity;
	
	// state switch
	void deactivate();
	void activate_contact();
	void deactivate_contact();

	double _r; // center-center distance  
	double ksi; // gap between particles
//	double iksi; // inverse gap
	double ksi_cutoff; // small cut-off for ksi: lubrication breakdown
	double ksi_eff;  // max(ksi, ksi_cutoff)
	double iksi_eff;
	double r_lub_max;  // max distance for lubrication
	double strain_0; // 
	//	double twothird, onesixth; // used in lubrication computations;

protected:
	void calcStaticFriction();
	void calcDynamicFriction();
	void calcNormalVector();
public:
 	Interaction(){};
	~Interaction(){};

	void init(System *sys_);

	void activate(int i, int j, const vec3d &pos_diff, double distance, int zshift);
	bool active;

	int particle_num[2];

	inline double r(){
	  return _r;
	}
	void r(double new_r);
	inline double gap(){
		return ksi;
	}
	
	bool update(); // after particles dispacement
	double age();
	
	double a0, a1;
	double ro; // ro = a0 + a1
	double lambda, invlambda;  // a1/a0 , a0/a1
	vec3d r_vec; // vector center to center
	vec3d nr_vec; // normal vector
	int partner(int);

	bool contact;
	bool static_friction;
	
	bool near;
	vector <vec3d> trajectory;
	vector <double> gap_history;
	void recordTrajectory();
	void outputTrajectory();
	
	/* Fc_normal: normal contact force
	 * Fc_tangent: tangneital contact force
	 */
	double Fc_normal;
	vec3d Fc_tangent;

	double valNormalForce();
	void calcContactInteraction();
	void calcContactInteractionNoFriction();
	void calcContactStress();
	void addLubricationStress();
	void addContactStress();
	void pairStresslet(double vi[], double vj[], double stresslet_i[], double stresslet_j[]);


	void XA(double &XAii, double &XAij, double &XAji, double &XAjj);
	void XG(double &XGii, double &XGij, double &XGji, double &XGjj);
	void XM(double &XMii, double &XMij, double &XMji, double &XMjj);
	void GE(double GEi[], double GEj[]);

};


#endif /* defined(__LF_DEM__Interaction__) */
