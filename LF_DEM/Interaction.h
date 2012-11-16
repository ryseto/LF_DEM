//
//  Interaction.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__Interaction__
#define __LF_DEM__Interaction__
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Interaction.h"
#include "System.h"
using namespace std;
class System;

class Interaction{
private:
	System *sys;
	vec3d r_vec; // vector center to center
	vec3d rs_velocity; // relative surface velocity
	vec3d vel_ij_t;
	vec3d unit_tang_vel_contact;
protected:
	
public:
 	Interaction(){};
	~Interaction(){};
	void init(System *sys_);
	void create(int i, int j);
	void calcInteraction();
	void incrementTangentialDisplacement(double dt);

/*********************************************/
	bool active;
	bool static_friction;
	int particle_num[2];
	double r; // center-center distance
	double f_normal;
	vec3d f_tangent;
	vec3d t_tangent;
	
	vec3d nr_vec; // normal vector
	vec3d xi; // tangential displacement
	int pd_x;
	int pd_y;
	int pd_z;
};

#endif /* defined(__LF_DEM__Interaction__) */
