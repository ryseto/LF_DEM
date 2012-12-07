//
//  ContactForce.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "ContactForce.h"

void ContactForce::init(System *sys_){
	active = false;
	sys = sys_;
}

/* Activate interaction between particles i and j.
 */
void  ContactForce::create(int i, int j){
	active = true;
	static_friction = true;
	xi.reset();
	particle_num[0] = i;
	particle_num[1] = j;
	return;
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1.
 * pd_x, pd_y, pd_z : Periodic boundary condition
 */
void ContactForce::makeNormalVector(){
	r_vec = sys->position[ particle_num[1] ] - sys->position[ particle_num[0] ];

	if (r_vec.z > sys->lz2){
		pd_z = 1; //  p1 (z = lz), p0 (z = 0)
		r_vec.z -= sys->lz;
		r_vec.x -= sys->x_shift;
	} else if (r_vec.z < - sys->lz2){
		pd_z = -1; //  p1 (z = 0), p0 (z = lz)
		r_vec.z += sys->lz;
		r_vec.x += sys->x_shift;
	} else{
		pd_z = 0;
	}

	while (r_vec.x > sys->lx2){
		r_vec.x -= sys->lx;
	}
	while (r_vec.x < -sys->lx2){
		r_vec.x += sys->lx;
	}
	
	if (sys->dimension == 3){
		if ( abs(r_vec.y) > sys->ly2 ){
			if ( r_vec.y > 0 ){
				r_vec.y -= sys->ly;
			} else {
				r_vec.y += sys->ly;
			}
		}
	}
}

void ContactForce::calcStaticFriction(){
	double f_static = -sys->mu_static*f_normal;
	double f_spring = sys->kt*xi.norm();
	if (f_spring < f_static){
		f_tangent = -sys->kt*xi; //
	} else {
		/* switch to dynamic friction */
		static_friction = false;
		calcDynamicFriction();
	}
}

void ContactForce::calcDynamicFriction(){
	double f_dynamic = -sys->mu_dynamic*f_normal;
	/* Use the velocity of one time step before as approximation. */
	unit_contact_velocity_tan = contact_velocity_tan/contact_velocity_tan.norm();
	f_tangent = -f_dynamic*unit_contact_velocity_tan;
}

/* 
 * Calculate interaction.
 * Force acts on particle 0 from particle 1.
 *
 */
void ContactForce::calcInteraction(){
	if (active){
		makeNormalVector();
		r = r_vec.norm();
		if (r < 2){
			nr_vec = r_vec / r;
			f_normal = sys->kn*(r - 2);
			if (static_friction){
				calcStaticFriction();
			} else {
				calcDynamicFriction();
			}
			sys->force[particle_num[0]] += f_tangent;
			sys->force[particle_num[1]] -= f_tangent;
			t_tangent = cross(nr_vec, f_tangent);
			sys->torque[particle_num[0]] += t_tangent;
			sys->torque[particle_num[1]] += t_tangent;
			sys->force[particle_num[0]] += f_normal * nr_vec;
			sys->force[particle_num[1]] -= f_normal * nr_vec;
		} else {
			static_friction = true;
		}
	}
}

void ContactForce::calcInteractionNoFriction(){
	if (active){
		makeNormalVector();
		r = r_vec.norm();
		if (r < 2){
			nr_vec = r_vec / r;
			f_normal = sys->kn*(r - 2);
			sys->force[ particle_num[0] ] += f_normal * nr_vec;
			sys->force[ particle_num[1] ] -= f_normal * nr_vec;
		} else {
			static_friction = true;
		}
	}
}

void ContactForce::incrementTangentialDisplacement(){
	// relative velocity particle 0 from particle 1.
	contact_velocity = sys->velocity[particle_num[0]] - sys->velocity[particle_num[1]];
	if (pd_z != 0){
		contact_velocity.x += pd_z * sys->vel_difference;
	}

	contact_velocity  += cross(sys->ang_velocity[particle_num[0]] + sys->ang_velocity[particle_num[1]], nr_vec);
	contact_velocity_tan = contact_velocity - dot(contact_velocity,nr_vec)*nr_vec;
	if (static_friction){
		xi += contact_velocity_tan*sys->dt;
		// projection
		xi -= dot(xi,nr_vec)*nr_vec;
	} else {
		sqnorm_contact_velocity = contact_velocity_tan.sq_norm();
		if ( sqnorm_contact_velocity < sys->sq_critical_velocity){
			cerr << "contact_velocity.norm()  = " << contact_velocity.norm()  << endl;
			static_friction = true;
			xi = -(1.0/sys->kt) * f_tangent;
		}
	}
}
