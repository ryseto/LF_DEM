//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "Interaction.h"

void Interaction::init(State *state_){
	active = false;
	state = state_;
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
		
		r_vec = state->position[ particle_num[1] ] - state->position[ particle_num[0] ];

		if (abs(r_vec.z) > state->lz2){
			if (r_vec.z > 0){
				r_vec.z -= state->lz;
				r_vec.x -= state->x_shift;
			}else{
				r_vec.z += state->lz;
				r_vec.x += state->x_shift;
			}
		}
		
		if (abs(r_vec.x) > state->lx2){
			if (r_vec.x > 0)
				r_vec.x -= state->lx;
			else
				r_vec.x += state->lx;
		}
		if (abs(r_vec.y) > state->ly2){
			if (r_vec.y > 0)
				r_vec.y -= state->ly;
			else
				r_vec.y += state->ly;
		}
		
		r = r_vec.norm();
		nr_vec = r_vec / r;

		double force;
		if ( r < 2){
			force = state->kn*(r - 2); // force < 0
		}
		state->force[  particle_num[0] ] += force * nr_vec;
		state->force[  particle_num[1] ] -= force * nr_vec;
	}
}








