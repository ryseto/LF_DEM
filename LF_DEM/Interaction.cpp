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
	
	//	twothird = 2./3.;
	//	onesixth = 1./6.;
}


void 
Interaction::r(double new_r){
	_r = new_r;
	ksi = 2 * _r / ro - 2.;
	iksi = 1./ksi;
	if( ksi < ksi_cutoff )
		ksi_eff = ksi_cutoff;
	else
		ksi_eff = ksi;
	
	iksi_eff = 1./ksi;
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 1 to particle 0. ( i --> j)
 * pd_x, pd_y, pd_z : Periodic boundary condition
 */
void
Interaction::calcNormalVector(){
	r_vec = sys->position[particle_num[1]] - sys->position[particle_num[0]];

	sys->periodize_diff(&r_vec, &pd_z);

	//	cout << "p0 " <<  particle_num[0] << " p1 " <<  particle_num[1] << " " << r_vec.x <<" " << r_vec.y <<" " << r_vec.z << endl;

}

void
Interaction::calcDistanceNormalVector(){
	if (active){
		calcNormalVector();
		r(r_vec.norm());
		nr_vec = r_vec / r();
	}
}

void
Interaction::assignDistanceNormalVector(vec3d pos_diff, double distance, int zshift){
	r_vec = pos_diff;
	r(distance);
	nr_vec = r_vec / r();
	pd_z = zshift;
	//	cout << "p0 " <<  particle_num[0] << " p1 " <<  particle_num[1] << " " << r_vec.x <<" " << r_vec.y <<" " << r_vec.z << endl;
}

/* Lubrication force and contact normal force
 * This function needs to be rewritten for poly disperse expression. 
 *
 *
 */
// @@@@@@@@@@@ THIS NEED TO BE UPDATED FOR POLYDISPERSE EXTENSION @@@@@@@@@@@@@@@@@@@
double
Interaction::valNormalForce(){
	double f_normal_total = 0;
	if ( ksi_eff > 0){
		int i = particle_num[0];
		int j = particle_num[1];
		vec3d rel_vel = sys->velocity[j] - sys->velocity[i];
		rel_vel.x += pd_z * sys->vel_difference;
		
		double alpha = 1.0/(4*ksi_eff);
		f_normal_total += abs(alpha*dot(rel_vel, nr_vec));
	}
	if (contact){
		f_normal_total += f_normal;
	}
	return f_normal_total;
}


/*********************************
*                                *
*	   Contact Forces Methods    *
*                                *
*********************************/
void
Interaction::calcStaticFriction(){
	double f_static = sys->mu_static*f_normal;
	double f_spring = sys->kt*xi.norm();
	if ( xi.x != 0){
		if (f_spring < f_static){
			/*
			 * f_tangent is force acting on particle 0 from particle 1
			 * xi = r1' - r0'
			 */
			f_tangent = sys->kt*xi;
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

	f_tangent = f_dynamic*unit_contact_velocity_tan; // on p0
}


/*
 * Calculate interaction.
 * Force acts on particle 0 from particle 1.
 * r_vec = p[1] - p[0]
 * f_normal is positive (by overlapping particles r < 2)
 */
void
Interaction::calcContactInteraction(){
	if (contact){
		f_normal = sys->kn*ksi;
		if (static_friction){
			calcStaticFriction();
		} else {
			calcDynamicFriction();
		}
		vec3d f_ij = - f_normal * nr_vec + f_tangent; // acting on p0
		vec3d t_ij = cross(nr_vec, f_tangent); // acting on p0
		sys->force[particle_num[0]] += f_ij;
		sys->force[particle_num[1]] -= f_ij;
		sys->torque[particle_num[0]] = a0*t_ij;
		sys->torque[particle_num[1]] = a1*t_ij;
	}
}

void
Interaction::calcContactInteractionNoFriction(){
	if (contact){
		f_normal = sys->kn*ksi;
		sys->force[ particle_num[0] ] -= f_normal * nr_vec;
		sys->force[ particle_num[1] ] += f_normal * nr_vec;
	}
}

/* Relative velocity of particle 1 from particle 0.
 *
 */
void
Interaction::calcContactVelocity(){
	// relative velocity particle 1 from particle 0.
	contact_velocity = sys->velocity[particle_num[1]] - sys->velocity[particle_num[0]];
	if (pd_z != 0){
		//	pd_z = -1; //  p1 (z = lz), p0 (z = 0)
		// v1 - v0
		contact_velocity.x += pd_z * sys->vel_difference;
	}
	contact_velocity -= cross(a0*sys->ang_velocity[particle_num[0]] + a1*sys->ang_velocity[particle_num[1]], nr_vec);
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
				xi = (1.0/sys->kt) * f_tangent;
			}
		}
	}
}


/*********************************
*                                *
*  Lubrication Forces Methods    *
*                                *
*********************************/

// Resistance functions
void
Interaction::XA(double &XAii, double &XAij, double &XAji, double &XAjj){
	double g1_l, g1_il;
	double l1, l13, il1, il13;

	l1 = 1.0 + lambda;
	l13 = l1 * l1 * l1;

	il1 = 1.0 + invlambda;
	il13 = il1 * il1 * il1;
	
	g1_l = 2.0 * lambda * lambda / l13;
	g1_il = 2.0 * invlambda * invlambda / il13;
	
	XAii = g1_l * iksi_eff;
	XAij = - 2 * XAii / l1;
	XAjj = g1_il * iksi_eff;
	XAji = - 2 * XAjj / il1;
}


void
Interaction::XG(double &XGii, double &XGij, double &XGji, double &XGjj){
	double g1_l, g1_il;
	double l1, l13, il1, il13;

	l1 = 1.0 + lambda;
	l13 = l1 * l1 * l1;

	il1 = 1.0 + invlambda;
	il13 = il1 * il1 * il1;
	
	g1_l = 2.0 * lambda * lambda / l13;
	g1_il = 2.0 * invlambda * invlambda / il13;
	
	XGii = 1.5 * g1_l * iksi_eff;
	XGij = - 4 * XGii / l1 / l1 ;
	XGjj = - 1.5 * g1_il * iksi_eff;
	XGji = - 4 * XGjj / il1 / il1 ;
	
}

void
Interaction::XM(double &XMii, double &XMij, double &XMji, double &XMjj){


	// under contruction
	double g1_l, g1_il;
	double l1, l13, il1, il13;

	l1 = 1.0 + lambda;
	l13 = l1 * l1 * l1;

	il1 = 1.0 + invlambda;
	il13 = il1 * il1 * il1;
	
	g1_l = 2.0 * lambda * lambda / l13;
	g1_il = 2.0 * invlambda * invlambda / il13;
	
	XMii = 0.6 * g1_l * iksi_eff;
	XMij = 40 * XMii * lambda / ( 3 * l13 * l13 * l13);
	XMjj = 0.6 * g1_il * iksi_eff;
	XMji = XMij;
	
}

void
Interaction::GE(double GEi[], double GEj[]){
	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	
	double nxnz = nr_vec.x * nr_vec.z;
	
	for(int u=0; u<3; u++){
		GEi[u] = a0 * a0 * XGii * 2. / 3.;
		GEi[u] += ro * ro * XGji / 6.;
	}
	GEi[0] *= sys->shear_rate * nxnz * nr_vec.x;
	GEi[1] *= sys->shear_rate * nxnz * nr_vec.y;
	GEi[2] *= sys->shear_rate * nxnz * nr_vec.z;

	for(int u=0; u<3; u++){
		GEj[u] = a1 * a1 * XGjj * 2. / 3.;
		GEj[u] += ro * ro * XGij / 6.;
	}
	GEj[0] *= sys->shear_rate * nxnz * nr_vec.x;
	GEj[1] *= sys->shear_rate * nxnz * nr_vec.y;
	GEj[2] *= sys->shear_rate * nxnz * nr_vec.z;
}

void
Interaction::addLubricationStress(){

	int i = particle_num[0];
	int j = particle_num[1];

	double n [3];
	n[0] = nr_vec.x;
	n[1] = nr_vec.y;
	n[2] = nr_vec.z;

	double Sixx = 0.;
	double Sixy = 0.;
	double Sixz = 0.;
	double Siyz = 0.;
	double Siyy = 0.;
	double Sjxx = 0.;
	double Sjxy = 0.;
	double Sjxz = 0.;
	double Sjyz = 0.;
	double Sjyy = 0.;

	// First G*(U-Uing) term
	double vi [3];
	double vj [3];

	vi[0] = sys->relative_velocity[i].x;
	vi[1] = sys->relative_velocity[i].y;
	vi[2] = sys->relative_velocity[i].z;

	vj[0] = sys->relative_velocity[j].x;
	vj[1] = sys->relative_velocity[j].y;
	vj[2] = sys->relative_velocity[j].z;

	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double n0n0_13 = ( n[0] * n[0] - 1./3. );
	double n1n1_13 = ( n[1] * n[1] - 1./3. );
	double n0n1 = n[0] * n[1];
	double n0n2 = n[0] * n[2];
	double n1n2 = n[1] * n[2];

	double twothird = 2./3.;
	double onesixth = 1./6.;
	double common_factor_i = 0.;
	double common_factor_j = 0.;
	for(int u=0; u<3; u++){
		common_factor_i += n[u] * ( twothird * a0 * a0 * XGii * vi[u] + onesixth * ro * ro * XGij * vj[u] );
		common_factor_j += n[u] * ( twothird * a1 * a1 * XGjj * vj[u] + onesixth * ro * ro * XGji * vi[u] );
	}

	Sixx += n0n0_13 * common_factor_i;
	Sixy += n0n1 * common_factor_i;
	Sixz += n0n2 * common_factor_i;
	Siyy += n1n1_13 * common_factor_i;
	Siyz += n1n2 * common_factor_i;
	Sjxx += n0n0_13 * common_factor_j;
	Sjxy += n0n1 * common_factor_j;
	Sjxz += n0n2 * common_factor_j;
	Sjyy += n1n1_13 * common_factor_j;
	Sjyz += n1n2 * common_factor_j;


	// Second: MEinf term
	double XMii, XMjj, XMij, XMji;
	XM(XMii, XMij, XMji, XMjj);

	// to be constructed

	sys->lubstress[i][0] += Sixx;
	sys->lubstress[j][0] += Sjxx;
	
	sys->lubstress[i][1] += Sixy;
	sys->lubstress[j][1] += Sjxy;
	
	sys->lubstress[i][2] += Sixz;
	sys->lubstress[j][2] += Sjxz;
	
	sys->lubstress[i][3] += Siyz;
	sys->lubstress[j][3] += Sjyz;
	
	sys->lubstress[i][4] += Siyy;
	sys->lubstress[j][4] += Sjyy;

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



/* Activate interaction between particles i and j.
 */
void
Interaction::activate(int i, int j, vec3d pos_diff, double distance, int zshift){

	active = true;

	if(j>i){
		particle_num[0] = i;
		particle_num[1] = j;
	}		
	else{
		particle_num[0] = j;
		particle_num[1] = i;
	}		

	// tell it to particles i and j
	sys->interaction_list[i].insert(this);
	sys->interaction_list[j].insert(this);

	// tell them their new partner
	sys->interaction_partners[i].insert(j);
	sys->interaction_partners[j].insert(i);

	contact = false;
	a0 = sys->radius[particle_num[0]];
	a1 = sys->radius[particle_num[1]];
	ro = a0+a1; // for polydispesity, we will rewrite this to a1+a2
	lambda = a1 / a0;
	invlambda = 1. / lambda;
	ksi_cutoff = sys->gap_cutoff;
	r_lub_max = 0.5*sys->lub_max*ro;

	assignDistanceNormalVector(pos_diff, distance, zshift); 

	return;
}

void
Interaction::deactivate(){
	// r > lub_max
	active = false;
	int i=particle_num[0];
	int j=particle_num[1];
	sys->interaction_list[i].erase(this);
	sys->interaction_list[j].erase(this);
	sys->interaction_partners[i].erase(j);
	sys->interaction_partners[j].erase(i);
}

void
Interaction::activate_contact(){
	// r < a0 + a1
	contact = true;
	static_friction = true;
	xi.reset();

	calcContactVelocity();
}

void
Interaction::deactivate_contact(){
	// r > a0 + a1
	contact = false;
}

bool
Interaction::update(){
	// update tangential displacement: we do it before updating nr_vec
	if (sys->friction) {
		incrementContactTangentialDisplacement();
	}

	// compute new r_vec and distance
	calcDistanceNormalVector();

	// check new state of the interaction
	if (active){
		if(r() > r_lub_max){
			deactivate();
			return true;
		} else {
			if (contact){
				if (r() > ro ){			
					deactivate_contact();
				}
			} else {
					// contact false:
				if (r() < ro){
					activate_contact();
				}
			}
		}
	}
	return false;
}
