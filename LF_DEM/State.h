//
//  State.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__State__
#define __LF_DEM__State__

#include <iostream>
#include <iomanip>
#include <fstream>
#include "vec3d.h"
#include "Interaction.h"
#include <queue>


using namespace std;
class Interaction;

class State{
private:
	int num_particle;
protected:
public:
    /* For DEMsystem
     */
	State(){
		kn = 100;
		lx = 10;
		ly = 10;
		lz = 10;
		lx2 =lx/2;
		ly2 =ly/2;
		lz2 =lz/2;
		
	};
	~State(){
		delete [] position;
		delete [] velocity;
		delete [] ang_velocity;
	};
	vec3d *position;
	vec3d *velocity;
	vec3d *ang_velocity;
	vec3d *force;
	double kn;
	double kt;
	double lx;
	double ly;
	double lz;
	double lx2;
	double ly2;
	double lz2;
	double x_shift;
//	void hello();

	double sq_distance(vec3d &pos , int i){
		double dx = position[i].x - pos.x;
		double dy = position[i].y - pos.y;
		double dz = position[i].z - pos.z;
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
		//		return sq_dist(position[i], position[j]);
	}
	double  sq_distance(int i, int j){
		double dx = position[i].x - position[j].x;
		double dy = position[i].y - position[j].y;
		double dz = position[i].z - position[j].z;
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
//		return sq_dist(position[i], position[j]);
	}
	
	
	double  distance(int i, int j){
		return sqrt(sq_distance(i,j));
	}
	
	
	void setNumberParticle(int num_particle_){
		num_particle = num_particle_;
		position = new vec3d [num_particle];
		velocity = new vec3d [num_particle];
		ang_velocity = new vec3d [num_particle];
		force = new vec3d [num_particle];
			}
	void setRandomPosition();
		

//st->hello();



};






#endif /* defined(__LF_DEM__State__) */

