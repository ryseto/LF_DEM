//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "Interaction.h"

void Interaction::init(System *sys_){
	active = false;
	sys = sys_;
}

void  Interaction::create(int i, int j){
	active = true;
	particle_num[0] = i;
	particle_num[1] = j;
	return;
}

void Interaction::calcInteraction(){
	if (active){
		// 
		// p0 ---> p1
		
		r_vec = sys->position[ particle_num[1] ] - sys->position[ particle_num[0] ];

		if (abs(r_vec.z) > sys->lz2){
			if (r_vec.z > 0){
				r_vec.z -= sys->lz;
				r_vec.x -= sys->x_shift;
			}else{
				r_vec.z += sys->lz;
				r_vec.x += sys->x_shift;
			}
		}
		
		if (abs(r_vec.x) > sys->lx2){
			if (r_vec.x > 0)
				r_vec.x -= sys->lx;
			else
				r_vec.x += sys->lx;
		}
		if (abs(r_vec.y) > sys->ly2){
			if (r_vec.y > 0)
				r_vec.y -= sys->ly;
			else
				r_vec.y += sys->ly;
		}
		
		r = r_vec.norm();
		nr_vec = r_vec / r;

		double force;
		if ( r < 2){
			force = sys->kn*(r - 2); // force < 0
		}
		sys->force[  particle_num[0] ] += force * nr_vec;
		sys->force[  particle_num[1] ] -= force * nr_vec;
	}
}








