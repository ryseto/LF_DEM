//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "System.h"
#include <sstream>


void System::init(){
	lx2 = 0.5*lx;
	ly2 = 0.5*ly;
	lz2 = 0.5*lz;
	ostringstream ss_simu_name;
	
	ss_simu_name << "D" << dimension << "L" << lz << "vf" << volume_fraction <<  "ms" << mu_static << "md" << mu_dynamic ;
	
	simu_name = ss_simu_name.str();
	cerr << simu_name << endl;
	sq_critical_velocity = dynamic_friction_critical_velocity * dynamic_friction_critical_velocity;
}
/* Set number of particles.
 * Allocate vectors for the state.
 */
void System::setNumberParticle(int num_particle_){
	num_particle = num_particle_;
	position = new vec3d [num_particle];
	
	i_position = new int * [num_particle];
	for (int i = 0; i<num_particle; i++ ){
		i_position[i] = new int [3];
	}

	
	angle = new double [num_particle];
	velocity = new vec3d [num_particle];
	ang_velocity = new vec3d [num_particle];
	force = new vec3d [num_particle];
	torque = new vec3d [num_particle];
}

/* Set positions of particles randomly.
 */
void System::setRandomPosition(){
	vec3d trial_pos;
	double lx_2 = lx - 2;
	double ly_2 = ly - 2;
	double lz_2 = lz - 2;
	for (int i=0; i < num_particle; i++){
		bool overlap;
		do{
			trial_pos.x = lx*drand48();
			trial_pos.z = lz*drand48();
			if (dimension == 2){
				trial_pos.y = ly2;
			} else {
				trial_pos.y = ly*drand48();
			}
			overlap = false;
			/* Before dr^2 < 4
			 * Check |dz| < 2, |dx| < 2 and |dy| < 2.
			 *
			 */
			
			for (int j = 0; j < i ; j++){
				double dx = position[j].x - trial_pos.x;
				if (dx > lx_2 ){
					dx -= lx;
				} else if (dx < -lx_2){
					dx += lx;
				}
				if (abs(dx) < 2){
					double dz = position[j].z - trial_pos.z;
					if (dz > lz_2 ){
						dz -= lz;
					} else if (dz < -lz_2){
						dz += lz;
					}
					if (abs(dz) < 2){
						double dy = position[j].y - trial_pos.y;
						if (dy > ly_2 ){
							dy -= ly;
						} else if (dy < -ly_2){
							dy += ly;
						}
						if (abs(dy) < 2){
							//cerr << dx*dx +  dy*dy + dz* dz << endl;
							if ( dx*dx +  dy*dy + dz* dz < 4){
								overlap = true;
								break;
							}
						}
					}
				}
			}
		} while (overlap);
		position[i] = trial_pos;
		angle[i] = 0;
		cerr << i << " / " << num_particle << endl;
	}
}


void System::updateVelocity(){
	vec3d O_inf(0, 0.5*shear_rate, 0);
	vec3d U_inf(0, 0, 0);
	for (int i=0; i < num_particle; i++){
		U_inf.x = shear_rate*position[i].z;
		velocity[i] = (1.0/eta)*force[i] + U_inf;
		ang_velocity[i] = (1.33333/eta)*torque[i] + O_inf;
	}
}

void System::deltaTimeEvolution(){
	x_shift += shear_rate*lz*dt;
	if (x_shift > lx){
		x_shift -= lx;
	}
	for (int i=0; i < num_particle; i++){
		position[i] += velocity[i]*dt;
	}
	for (int i=0; i < num_particle; i++){
		if (position[i].z > lz ){
			position[i].z -= lz;
			position[i].x -= x_shift;
		} else if ( position[i].z < 0 ){
			position[i].z += lz;
			position[i].x += x_shift;
		}
		if ( position[i].x > lx ){
			position[i].x -= lx;
		} else if (position[i].x < 0 ){
			position[i].x += lx;
		}
		if ( position[i].y > ly ){
			position[i].y -= ly;
		} else if (position[i].y < 0 ){
			position[i].y += ly;
		}
	}
//	i_position[i][0] = position[i].x / 2.0;
//	i_position[i][1] = position[i].y / 2.0;
//	i_position[i][2] = position[i].z / 2.0;
	if (dimension == 2){
		for (int i=0; i < num_particle; i++){
			angle[i] += ang_velocity[i].y*dt;
		}
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

double System::checkContact(int i, int j){
	double dx = position[i].x - position[j].x;
	double dz = position[i].z - position[j].z;
	if (dz > lz2 ){
		dz -= lz;
		dx -= x_shift;
	} else if (dz < -lz2){
		dz += lz;
		dx += x_shift;
	}
	if (abs(dz) < 2){
		if (dx > lx2 ){
			dx -= lx;
		} else if (dx < -lx2){
			dx += lx;
		}
		if (abs(dx) < 2){
			double dy = position[i].y - position[j].y;
			if (dy > ly2 ){
				dy -= ly;
			} else if (dy < -ly2){
				dy += ly;
			}
			if (abs(dy) < 2){
				return dx*dx + dy*dy + dz*dz;
			}
		}
	}
	return 100;
}








