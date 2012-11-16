//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "Interaction.h"
/* Initialize interaction object.
 */
void Interaction::init(System *sys_){
	active = false;
	sys = sys_;
}

/* Activate interaction between particles i and j.
 */
void  Interaction::create(int i, int j){
	active = true;
	static_friction = true;

	xi.reset();
	particle_num[0] = i;
	particle_num[1] = j;
	return;
}

/* Calculate interaction.
 * Force acts on particle 0 from particle 1.
 *
 */
void Interaction::calcInteraction(){
	if (active){
		r_vec = sys->position[ particle_num[1] ] - sys->position[ particle_num[0] ];
		if (abs(r_vec.z) > sys->lz2){
			if (r_vec.z > 0){
				pd_z = 1; //  p1 (z = lz), p0 (z = 0)
				r_vec.z -= sys->lz;
				r_vec.x -= sys->x_shift;
			}else{
				pd_z = -1; //  p1 (z = 0), p0 (z = lz)
				r_vec.z += sys->lz;
				r_vec.x += sys->x_shift;
			}
		} else{
			pd_z = 0;
		}
		if (abs(r_vec.x) > sys->lx2){
			if (r_vec.x > 0){
				pd_x = 1;  // p1 (x = lx), p0 (x = 0)
				r_vec.x -= sys->lx;
			}else{
				pd_x = -1;  // p1 (x = 0), p0 (x = lx)
				r_vec.x += sys->lx;
			}
		} else {
			pd_x = 0;
		}
		if ( abs(r_vec.y) > sys->ly2 ){
			if ( r_vec.y > 0 ){
				pd_y = 1;  // p1 (y = ly), p0 (y = 0)
				r_vec.y -= sys->ly;
			} else {
				pd_y = -1;  // p1 (y = ly), p0 (y = 0)
				r_vec.y += sys->ly;
			}
		} else {
			pd_y = 0;
		}
		r = r_vec.norm();
		nr_vec = r_vec / r;
		if (r < 2){
			f_normal = sys->kn*(r - 2); // force < 0
			if (static_friction){
				double f_static = - sys->mu_static*f_normal;
				double f_spring = sys->kt*xi.norm();
				if (f_spring < f_static){
					f_tangent = -sys->kt*xi;
				} else {
					double f_dynamic = -sys->mu_dynamic*f_normal;
					unit_tang_vel_contact = vel_ij_t / vel_ij_t.norm();
					f_tangent = -f_dynamic*unit_tang_vel_contact;
					static_friction = false;
				}
			} else {
				double f_spring = sys->kt*xi.norm();
				double f_dynamic = - sys->mu_dynamic*f_normal;
				if (f_spring < f_dynamic){
					static_friction = true;
					f_tangent = -sys->kt*xi;
				} else {
					unit_tang_vel_contact = vel_ij_t / vel_ij_t.norm();
					f_tangent = -f_dynamic*unit_tang_vel_contact;
					xi = -(1.0/sys->kt) * f_tangent;
				}
			}
			sys->force[  particle_num[0] ] += f_tangent;
			sys->force[  particle_num[1] ] -= f_tangent;
			/* nr_vec is from particle 0 to particle 1;
			 */
			t_tangent = cross(f_tangent, nr_vec );
			sys->torque[  particle_num[0] ] -= t_tangent;
			sys->torque[  particle_num[1] ] -= t_tangent;
		}
		sys->force[  particle_num[0] ] += f_normal * nr_vec;
		sys->force[  particle_num[1] ] -= f_normal * nr_vec;
	}
}

void Interaction::incrementTangentialDisplacement(double dt){
	// relative velocity particle 0 from particle 1.
	vec3d vel_ij = sys->velocity[particle_num[0]] - sys->velocity[particle_num[1]];
	vel_ij.x +=  pd_z*sys->lz*sys->shear_rate;

	vel_ij += cross(sys->ang_velocity[ particle_num[0]], nr_vec) + cross(sys->ang_velocity[ particle_num[1]], nr_vec);
	
	vel_ij_t = vel_ij - dot(vel_ij,nr_vec)*nr_vec;
	xi += vel_ij_t*dt;
}
