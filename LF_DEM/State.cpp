//
//  State.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "State.h"

void State::setRandomPosition(){
	for (int i=0; i < num_particle; i++){
		vec3d trial_pos;
		bool overlap;
		do{
			double x = lx*drand48();
			double y = ly*drand48();
			//				double y = ly2;
			double z = lz*drand48();
			trial_pos.set(x,y,z);
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
		double vx = 0;
		double vy = 0;
		double vz = 0;
		velocity[i].set(vx,vy,vz);
		ang_velocity[i].set(0,0,0);
	}
}
