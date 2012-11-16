//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "System.h"

void System::init(){
	lx2 = 0.5*lx;
	ly2 = 0.5*ly;
	lz2 = 0.5*lz;
}


/* Set number of particles.
 * Allocate vectors for the state.
 */
void System::setNumberParticle(int num_particle_){
	num_particle = num_particle_;
	position = new vec3d [num_particle];
	velocity = new vec3d [num_particle];
	ang_velocity = new vec3d [num_particle];
	force = new vec3d [num_particle];
	torque = new vec3d [num_particle];
}

/* Set positions of particles randomly.
 */
void System::setRandomPosition(int dimension){
	for (int i=0; i < num_particle; i++){
		vec3d trial_pos;
		bool overlap;
		do{
			trial_pos.x = lx*drand48();
			if (dimension == 2){
				trial_pos.y = ly2;
			} else {
				trial_pos.y = ly*drand48();
			}
			trial_pos.z = lz*drand48();
			overlap = false;
			for (int j = 0; j < i ; j++){
				if (sq_distance(trial_pos, j) < 4){
					overlap = true;
					break;
				}
			}
			if ( overlap == false)
				break;
			
		} while (true);
		position[i] = trial_pos;
	}
}

/* Distance between particle i and particle j
 */
double System::distance(int i, int j){
	return sqrt(sq_distance(i,j));
}

/* Square norm of vector (dx, dy, dz)
 */
double System::sq_norm(double &dx, double &dy, double &dz){
	if (dz > lz2 ){
		dz -= lz;
		dx -= x_shift;
	} else if (dz < -lz2){
		dz += lz;
		dx += x_shift;
	}
	if (dx > lx2 ){
		dx -= lx;
	} else if (dx < -lx2){
		dx += lx;
	}
	if (dy > ly2 ){
		dy -= ly;
	} else if (dy < -ly2){
		dy += ly;
	}
	return dx*dx + dy*dy + dz*dz;
}

/* Square distance between particle i and particle j
 */
double System::sq_distance(int i, int j){
	double dx = position[i].x - position[j].x;
	double dy = position[i].y - position[j].y;
	double dz = position[i].z - position[j].z;
	return sq_norm(dx, dy, dz);
}

/* Square distance between position pos and particle i
 */
double System::sq_distance(vec3d &pos , int i){
	double dx = position[i].x - pos.x;
	double dy = position[i].y - pos.y;
	double dz = position[i].z - pos.z;
	return sq_norm(dx, dy, dz);
}
