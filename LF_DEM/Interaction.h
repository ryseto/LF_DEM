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
//#include "ContactForce.h"
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "System.h"
//#include "ContactForce.h"
using namespace std;
class System;
//class ContactForce;

class Interaction{
private:
	System *sys;
//	ContactForce fcontact;
	vec3d contact_velocity;
	vec3d contact_velocity_tan;
	vec3d unit_contact_velocity_tan;
protected:
	void calcStaticFriction();
	void calcDynamicFriction();
	
public:
 	Interaction(){};
	~Interaction(){};
	void init(System *sys_);
	void newContact();
	void create(int i, int j);
	
	
	void calcContactInteraction();
	void calcContactInteractionNoFriction();
	void calcContactStress();
	void normalElement();
	void makeNormalVector();
	void incrementContactTangentialDisplacement();
	

	bool active;
	bool static_friction;
	double sqnorm_contact_velocity;
	int particle_num[2];
	double r; // center-center distance // done
	double f_normal;
	vec3d f_tangent;
	vec3d t_tangent;
	vec3d r_vec; // vector center to center
	vec3d nr_vec; // normal vector
	vec3d xi; // tangential displacement
	int pd_z;

	bool contact;
};

#endif /* defined(__LF_DEM__Interaction__) */
