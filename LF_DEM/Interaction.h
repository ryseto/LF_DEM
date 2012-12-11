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
protected:
	void calcStaticFriction();
	void calcDynamicFriction();
	void calcNormalVector();
public:
 	Interaction(){};
	~Interaction(){};
	void init(System *sys_);
	void newContact();
	void create(int i, int j);
	void calcContactInteraction();
	void calcContactInteractionNoFriction();
	void calcContactStress();
	void calcDistanceNormalVector();
	void calcContactVelocity();
	void addLubricationStress();
	void addContactStress();
	double valNormalForce();
	void incrementContactTangentialDisplacement();
	bool active;
	bool contact;
	bool static_friction;
	double sqnorm_contact_velocity;
	int particle_num[2];
	double r; // center-center distance // done
	vec3d r_vec; // vector center to center
	vec3d nr_vec; // normal vector
	double f_normal;
	vec3d f_tangent;
//	vec3d t_tangent;
	int pd_z;
	
};

#endif /* defined(__LF_DEM__Interaction__) */
