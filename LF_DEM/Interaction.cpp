//
//  Interaction.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "Interaction.h"

void
Interaction::init(System *sys_){
	sys = sys_;
	contact = false;
	active = false;
	f_normal = 0;
	
}

/* Activate interaction between particles i and j.
 */
void
Interaction::create(int i, int j){
	if(j>i){
		particle_num[0] = i;
		particle_num[1] = j;
	}		
	else{
		particle_num[0] = j;
		particle_num[1] = i;
	}		
	active = true;
	contact = false;
	ro = 2; // for polydispesity, we will rewrite this to a1+a2
	return;
}

void
Interaction::newContact(){
	contact = true;
	static_friction = true;
	xi.reset();
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 1 to particle 0. ( i --> j)
 * pd_x, pd_y, pd_z : Periodic boundary condition
 */
void
Interaction::calcNormalVector(){
	r_vec = sys->position[particle_num[0]] - sys->position[particle_num[1]];

	sys->periodize_diff(&r_vec, &pd_z);

	//	cout << "p0 " <<  particle_num[0] << " p1 " <<  particle_num[1] << " " << r_vec.x <<" " << r_vec.y <<" " << r_vec.z << endl;

}

void
Interaction::calcDistanceNormalVector(){
	if (active){
		calcNormalVector();
		r = r_vec.norm();
		nr_vec = r_vec / r;
	}
}

void
Interaction::assignDistanceNormalVector(vec3d pos_diff, double distance, int zshift){
	r_vec = pos_diff;
	r = distance;
	nr_vec = r_vec / r;
	pd_z = zshift;
	//	cout << "p0 " <<  particle_num[0] << " p1 " <<  particle_num[1] << " " << r_vec.x <<" " << r_vec.y <<" " << r_vec.z << endl;
}


void
Interaction::calcStaticFriction(){
	double f_static = sys->mu_static*f_normal;
	double f_spring = sys->kt*xi.norm();
	if ( xi.x != 0){
		if (f_spring < f_static){
			/*
			 * f_tangent is force acting on particle 0 from particle 1
			 * xi = r0' - r1'
			 */
			f_tangent = -sys->kt*xi;
		} else {
			/* switch to dynamic friction */
			static_friction = false;
			calcDynamicFriction();
		}
	}
}

void
Interaction::calcDynamicFriction(){
	double f_dynamic = sys->mu_dynamic*f_normal;
	/* Use the velocity of one time step before as approximation. */
	unit_contact_velocity_tan = contact_velocity_tan/contact_velocity_tan.norm();

	f_tangent = -f_dynamic*unit_contact_velocity_tan;
}


/*
 * Calculate interaction.
 * Force acts on particle 0 from particle 1.
 * r_vec = p[0] - p[1]
 * f_normal is positive (by overlapping particles r < 2)
 */
void
Interaction::calcContactInteraction(){
	if (contact){
		f_normal = sys->kn*(ro - r);
		if (static_friction){
			calcStaticFriction();
		} else {
			calcDynamicFriction();
		}
		vec3d f_ij = f_normal * nr_vec + f_tangent; // acting on p0
		vec3d t_ij = cross(-nr_vec, f_tangent); // acting on p0
		sys->force[particle_num[0]] += f_ij;
		sys->force[particle_num[1]] -= f_ij;
		sys->torque[particle_num[0]] = t_ij;
		sys->torque[particle_num[1]] = t_ij;
	}
}

void
Interaction::calcContactInteractionNoFriction(){
	if (contact){
		f_normal = sys->kn*(ro-r);
		sys->force[ particle_num[0] ] += f_normal * nr_vec;
		sys->force[ particle_num[1] ] -= f_normal * nr_vec;
	}
}

/* Relative velocity of particle 0 from particle 1.
 *
 */
void
Interaction::calcContactVelocity(){
	// relative velocity particle 0 from particle 1.
	contact_velocity = sys->velocity[particle_num[0]] - sys->velocity[particle_num[1]];
	if (pd_z != 0){
		// v0 - v1
		//	pd_z = 1; //  p1 (z = lz), p0 (z = 0)
		// v0 - v1
		contact_velocity.x += pd_z * sys->vel_difference;
	}
	contact_velocity += cross(nr_vec, sys->ang_velocity[particle_num[0]] + sys->ang_velocity[particle_num[1]]);
	contact_velocity_tan = contact_velocity - dot(contact_velocity,nr_vec)*nr_vec;
}

void
Interaction::incrementContactTangentialDisplacement(){
	if (active){
		calcContactVelocity();
		if (static_friction){
			xi += contact_velocity_tan*sys->dt;
			// projection
			xi -= dot(xi,nr_vec)*nr_vec;
		} else {
			sqnorm_contact_velocity = contact_velocity_tan.sq_norm();
			if ( sqnorm_contact_velocity < sys->sq_critical_velocity){
				static_friction = true;
				xi = -(1.0/sys->kt) * f_tangent;
			}
		}
	}
}

double
Interaction::valNormalForce(){
	double h = r  - ro;
	double f_normal_total = 0;
	if ( h > 0){
		int i = particle_num[0];
		int j = particle_num[1];
		vec3d rel_vel = sys->velocity[i] - sys->velocity[j];
		rel_vel.x += pd_z * sys->vel_difference;
		if ( h < sys->h_cutoff){
			h = sys->h_cutoff;
		}
		double alpha = 1.0/(4*h);
		f_normal_total += abs(alpha*dot(rel_vel, nr_vec));
	}
	if (contact){
		f_normal_total += f_normal;
	}
	return f_normal_total;
}

void
Interaction::addLubricationStress(){
	double h = r  - ro;
	int i = particle_num[0];
	int j = particle_num[1];
	vec3d rel_vel = sys->velocity[i] - sys->velocity[j];
	double alpha;
	if ( h < sys->h_cutoff){
		h = sys->h_cutoff;
	}
	if (h > 0){
		alpha = 1.0/(4*h);
		rel_vel.x += pd_z * sys->vel_difference;
		double alpha = 1.0/(4*h);
		vec3d force = alpha*dot(rel_vel, nr_vec)*nr_vec;
		double Sxx = 2*(force.x * nr_vec.x);
		double Sxy = force.x * nr_vec.y + force.y * nr_vec.x ;
		double Sxz = force.x * nr_vec.z + force.z * nr_vec.x ;
		double Syz = force.y * nr_vec.z + force.z * nr_vec.y ;
		double Syy = 2*(force.y * nr_vec.y);
		sys->lubstress[i][0] += Sxx;
		sys->lubstress[j][0] += Sxx;
		
		sys->lubstress[i][1] += Sxy;
		sys->lubstress[j][1] += Sxy;
		
		sys->lubstress[i][2] += Sxz;
		sys->lubstress[j][2] += Sxz;
		
		sys->lubstress[i][3] += Syz;
		sys->lubstress[j][3] += Syz;
		
		sys->lubstress[i][4] += Syy;
		sys->lubstress[j][4] += Syy;
	}
}

void
Interaction::addContactStress(){
	int i = particle_num[0];
	int j = particle_num[1];
	if (contact){
		vec3d force = - (f_normal * nr_vec + f_tangent);
		double Sxx = 2*(force.x * nr_vec.x);
		double Sxy = force.x * nr_vec.y + force.y * nr_vec.x ;
		double Sxz = force.x * nr_vec.z + force.z * nr_vec.x ;
		double Syz = force.y * nr_vec.z + force.z * nr_vec.y ;
		double Syy = 2*(force.y * nr_vec.y);
		sys->contactstress[i][0] += Sxx;
		sys->contactstress[j][0] += Sxx;
		
		sys->contactstress[i][1] += Sxy;
		sys->contactstress[j][1] += Sxy;
		
		sys->contactstress[i][2] += Sxz;
		sys->contactstress[j][2] += Sxz;
		
		sys->contactstress[i][3] += Syz;
		sys->contactstress[j][3] += Syz;
		
		sys->contactstress[i][4] += Syy;
		sys->contactstress[j][4] += Syy;
	}
}

int
Interaction::partner(int i){
	if( i == particle_num[0] )
		return particle_num[1];
	else
		return particle_num[0];
}



