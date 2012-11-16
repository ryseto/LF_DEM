//
//  System.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__System__
#define __LF_DEM__System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Interaction.h"
#include <queue>

using namespace std;
class Interaction;

class System{
private:
	int num_particle;
protected:
public:
    /* For DEMsystem
     */
	System(){
	
	};
	~System(){
		delete [] position;
		delete [] velocity;
		delete [] ang_velocity;
	};
	vec3d *position;
	vec3d *velocity;
	vec3d *ang_velocity;
	vec3d *force;
	vec3d *torque;
	double kn;
	double kt;
	double mu_static; // static friction coefficient.
	double mu_dynamic;// dynamic friction coefficient.
	
	double lx;
	double ly;
	double lz;
	double lx2;
	double ly2;
	double lz2;
	double x_shift;
	double shear_rate;

	/*************************************************************/
	void setNumberParticle(int num_particle_);
	void init();
	double sq_norm(double &dx, double &dy, double &dz);
	double sq_distance(vec3d &pos , int i);
	double sq_distance(int i, int j);
	double distance(int i, int j);
	void setRandomPosition(int dimension);
};


#endif /* defined(__LF_DEM__State__) */

