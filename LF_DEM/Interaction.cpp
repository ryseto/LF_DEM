//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "Interaction.h"

void Interaction::init(System *sys_){
	contact = false;
	sys = sys_;
	active = false;
	//xi.reset();
	//static_friction = true;
}


/* Activate interaction between particles i and j.
 */
void  Interaction::create(int i, int j){
	active = true;
	particle_num[0] = i;
	particle_num[1] = j;
	return;
}

void Interaction::newContact(){
	contact = true;
	xi.reset();
	static_friction = true;
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 1 to particle 0. ( i --> j)
 * pd_x, pd_y, pd_z : Periodic boundary condition
 */
void Interaction::makeNormalVector(){
	r_vec = sys->position[ particle_num[0]] - sys->position[ particle_num[1]  ];
	if (r_vec.z > sys->lz2){
		pd_z = 1; //  p1 (z = lz), p0 (z = 0)
		r_vec.z -= sys->lz;
		r_vec.x -= sys->shear_disp;
	} else if (r_vec.z < - sys->lz2){
		pd_z = -1; //  p1 (z = 0), p0 (z = lz)
		r_vec.z += sys->lz;
		r_vec.x += sys->shear_disp;
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

void Interaction::calcContactStress(){
	// S_xx
	vec3d f_contact = f_tangent + f_normal* nr_vec;
	double Sxx = f_contact.x * nr_vec.x + f_contact.x * nr_vec.x ;;
	double Sxy = f_contact.x * nr_vec.y + f_contact.y * nr_vec.x ;
	double Sxz = f_contact.x * nr_vec.z + f_contact.z * nr_vec.x ;
	double Syz = f_contact.y * nr_vec.z + f_contact.z * nr_vec.y ;
	double Syy = f_contact.y * nr_vec.y + f_contact.y * nr_vec.y ;
	sys->stress[particle_num[0]][0] += Sxx;
	sys->stress[particle_num[1]][0] += Sxx;
	
	sys->stress[particle_num[0]][1] += Sxy;
	sys->stress[particle_num[1]][1] += Sxy;
	
	sys->stress[particle_num[0]][2] += Sxz;
	sys->stress[particle_num[1]][2] += Sxz;
	
	sys->stress[particle_num[0]][3] += Syz;
	sys->stress[particle_num[1]][3] += Syz;
	
	sys->stress[particle_num[0]][4] += Syy;
	sys->stress[particle_num[1]][4] += Syy;
}

void Interaction::calcStaticFriction(){
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


void Interaction::calcDynamicFriction(){
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
void Interaction::calcContactInteraction(){
	if (active){
		if (r < 2){
			contact = true;
//			nr_vec = r_vec / r;
//			cerr << "r " << r << endl;
//			exit(1);
			f_normal = sys->kn*(r - 2);
			if (static_friction){
				calcStaticFriction();
			} else {
				calcDynamicFriction();
			}
			sys->force[particle_num[0]] += f_tangent;
			sys->force[particle_num[1]] -= f_tangent;
			t_tangent = cross(nr_vec, f_tangent);
			sys->torque[particle_num[0]] -= t_tangent;
			sys->torque[particle_num[1]] -= t_tangent;
			sys->force[particle_num[0]] -= f_normal * nr_vec;
			sys->force[particle_num[1]] += f_normal * nr_vec;
		} else {
			static_friction = true;
			contact = false;
		}
	}
}

void Interaction::calcContactInteractionNoFriction(){
	if (active){
		if (r < 2){
			//nr_vec = r_vec / r;
			f_normal = sys->kn*(r - 2);
			sys->force[ particle_num[0] ] -= f_normal * nr_vec;
			sys->force[ particle_num[1] ] += f_normal * nr_vec;
		} else {
			static_friction = true;
		}
	}
}

void Interaction::normalElement(){
	makeNormalVector();
	r = r_vec.norm();
	nr_vec = r_vec / r;
}

void Interaction::incrementContactTangentialDisplacement(){
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
