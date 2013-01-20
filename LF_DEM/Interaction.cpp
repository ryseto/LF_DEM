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
	Fc_normal = 0;
	
	//	twothird = 2./3.;
	//	onesixth = 1./6.;
}


void 
Interaction::r(double new_r){
	_r = new_r;
	ksi = 2 * _r / ro - 2.;
	//	iksi = 1./ksi;
	if( ksi < ksi_cutoff )
		ksi_eff = ksi_cutoff;
	else
		ksi_eff = ksi;
	
	iksi_eff = 1./ksi_eff;
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 1 to particle 0. ( i --> j)
 * pd_x, pd_y, pd_z : Periodic boundary condition
 */
void
Interaction::calcNormalVector(){
	r_vec = sys->position[particle_num[1]] - sys->position[particle_num[0]];
	sys->periodize_diff(r_vec, pd_z);
}

void
Interaction::calcDistanceNormalVector(){
	calcNormalVector();
	r(r_vec.norm());
	nr_vec = r_vec / r();
}

void
Interaction::assignDistanceNormalVector(const vec3d &pos_diff, double distance, int zshift){
	r_vec = pos_diff;
	r(distance);
	nr_vec = r_vec / r();
	pd_z = zshift;
	//	cout << "p0 " <<  particle_num[0] << " p1 " <<  particle_num[1] << " " << r_vec.x <<" " << r_vec.y <<" " << r_vec.z << endl;
}


/*********************************
*                                *
*	   Contact Forces Methods    *
*                                *
*********************************/
void
Interaction::calcStaticFriction(){
	double f_static = sys->mu_static*Fc_normal;
	double f_spring = sys->kt*xi.norm();
	if ( xi.x != 0){
		if (f_spring < f_static){
			/*
			 * f_tangent is force acting on particle 0 from particle 1
			 * xi = r1' - r0'
			 */

			Fc_tangent = sys->kt*xi;
		} else {
			/* switch to dynamic friction */
			static_friction = false;
			calcDynamicFriction();
		}
	}
}

void
Interaction::calcDynamicFriction(){
	double f_dynamic = sys->mu_dynamic*Fc_normal;
	/* Use the velocity of one time step before as approximation. */
	unit_contact_velocity_tan = contact_velocity_tan/contact_velocity_tan.norm();

	Fc_tangent = f_dynamic*unit_contact_velocity_tan; // on p0
}


/*
 * Calculate interaction.
 * Force acts on particle 0 from particle 1.
 * r_vec = p[1] - p[0]
 * Fc_normal is positive (by overlapping particles r < ro)
 */
void
Interaction::calcContactInteraction(){
	if (contact){
		Fc_normal = -sys->kn*ksi;
		if (static_friction){
			calcStaticFriction();
		} else {
			calcDynamicFriction();
		}
		
		vec3d f_ij = - Fc_normal * nr_vec + Fc_tangent; // acting on p0 //@@TO BE CHECKED.
		vec3d t_ij = cross(nr_vec, Fc_tangent); // acting on p0
		sys->contact_force[particle_num[0]] += f_ij;
		sys->contact_force[particle_num[1]] -= f_ij;
		sys->total_force[particle_num[0]] += f_ij;
		sys->total_force[particle_num[1]] -= f_ij;
		sys->torque[particle_num[0]] = a0*t_ij;
		sys->torque[particle_num[1]] = a1*t_ij;
//		cerr << Fc_tangent.x << ' ' << Fc_tangent.z << endl;
	}
}

void
Interaction::calcContactInteractionNoFriction(){
	if (contact){
		Fc_normal = sys->kn*ksi;
		sys->contact_force[ particle_num[0] ] -= Fc_normal * nr_vec;
		sys->contact_force[ particle_num[1] ] += Fc_normal * nr_vec;
		sys->total_force[ particle_num[0] ] -= Fc_normal * nr_vec;
		sys->total_force[ particle_num[1] ] += Fc_normal * nr_vec;
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
	calcContactVelocity();
	if (static_friction){
		xi += contact_velocity_tan*sys->dt;
		// projection
		xi -= dot(xi,nr_vec)*nr_vec;
	} else {
		sqnorm_contact_velocity = contact_velocity_tan.sq_norm();
		if ( sqnorm_contact_velocity < sys->sq_critical_velocity){
			static_friction = true;
			xi = (1.0/sys->kt) * Fc_tangent;
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
	double g1_l; //, g1_il;
	double l1, l13; //, il1; //, il13;

	l1 = 1.0 + lambda;
	l13 = l1 * l1 * l1;

//	il1 = 1.0 + invlambda;
	//	il13 = il1 * il1 * il1;
	
	g1_l = 2.0 * lambda * lambda / l13;
	//g1_il = 2.0 * invlambda * invlambda / il13;
//	g1_il = g1_l/lambda;
	
	XAii = g1_l * iksi_eff;
	XAij = - 2 * XAii / l1;
	XAji = XAij;
	XAjj = XAii / lambda;
//	XAjj = g1_il * iksi_eff;
//	XAji = - 2 * XAjj / il1;

}


void
Interaction::XG(double &XGii, double &XGij, double &XGji, double &XGjj){
	double g1_l, g1_il;
	double l1, l13, il1; //, il13;
	l1 = 1.0 + lambda;
	l13 = l1 * l1 * l1;
	il1 = 1.0 + invlambda;
	//	il13 = il1 * il1 * il1;
	
	g1_l = 2.0 * lambda * lambda / l13;
	//g1_il = 2.0 * invlambda * invlambda / il13;
	g1_il = g1_l/lambda;
	
	XGii = 1.5 * g1_l * iksi_eff;
	XGij = - 4 * XGii / l1 / l1 ;
	//	XGjj = - 1.5 * g1_il * iksi_eff;
	XGjj = - XGii / lambda;
	XGji = - 4 * XGjj / il1 / il1 ;
}

void
Interaction::XM(double &XMii, double &XMij, double &XMji, double &XMjj){
	double g1_l; //, g1_il;
	double l1, l13;// , il1; //, il13;
	l1 = 1.0 + lambda;
	l13 = l1 * l1 * l1;
	//	il1 = 1.0 + invlambda;
	//	il13 = il1 * il1 * il1;
	g1_l = 2.0 * lambda * lambda / l13;
	//	g1_il = 2.0 * invlambda * invlambda / il13;
	//g1_il = g1_l/lambda;
	XMii = 0.6 * g1_l * iksi_eff;
	XMij = 40*lambda / (3 * l13) * XMii ;
	XMji = XMij;
	XMjj = XMii/lambda;
}

void
Interaction::GE(double GEi[], double GEj[]){
	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double nxnz = nr_vec.x * nr_vec.z;
	
	for(int u=0; u<3; u++){
		GEi[u] = a0 * a0 * XGii * 2 / 3;
		GEi[u] += ro * ro * XGji / 6;
	}
	GEi[0] *= sys->shear_rate * nxnz * nr_vec.x;
	GEi[1] *= sys->shear_rate * nxnz * nr_vec.y;
	GEi[2] *= sys->shear_rate * nxnz * nr_vec.z;

	for(int u=0; u<3; u++){
		GEj[u] = a1 * a1 * XGjj * 2 / 3;
		GEj[u] += ro * ro * XGij / 6;
	}
	GEj[0] *= sys->shear_rate * nxnz * nr_vec.x;
	GEj[1] *= sys->shear_rate * nxnz * nr_vec.y;
	GEj[2] *= sys->shear_rate * nxnz * nr_vec.z;
}


// computes the contribution to S = R_SU * V (in Brady's notations) [ S = G V in Jefrrey's ones ]
// from pair (i,j).
// ie fills :
// stresslet_i = R_SU^{ii} * vi + R_SU^{ij} * vj
// stresslet_j = R_SU^{ji} * vi + R_SU^{jj} * vj
 
void
Interaction::pairStresslet(double vi[], double vj[], double stresslet_i[], double stresslet_j[]){
	double n [3];
	n[0] = nr_vec.x;
	n[1] = nr_vec.y;
	n[2] = nr_vec.z;

	for (int k=0; k < 5; k ++){
		stresslet_i[k] = 0.;
		stresslet_j[k] = 0.;
	}
	
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

	stresslet_i[0] += n0n0_13 * common_factor_i;
	stresslet_i[1] += n0n1    * common_factor_i;
	stresslet_i[2] += n0n2    * common_factor_i;
	stresslet_i[3] += n1n2    * common_factor_i;
	stresslet_i[4] += n1n1_13 * common_factor_i;

	stresslet_j[0] += n0n0_13 * common_factor_j;
	stresslet_j[1] += n0n1 * common_factor_j;
	stresslet_j[2] += n0n2 * common_factor_j;
	stresslet_j[3] += n1n2 * common_factor_j;
	stresslet_j[4] += n1n1_13 * common_factor_j;

}


/*
 *
 * Stresslet stresslet_[5]
 * 0 = Sxx, 1 = Sxy, 2 = Sxz, 3 = Syz, 4 = Syy
 * (NOTE: I follow the way of RYUON)
 *
 * Unit: 
 *   F0 = 6 pi \eta a U0
 *   U0 = a \dot{\gamma}
 *   L0 = a
 */
void
Interaction::addLubricationStress(){
	int i = particle_num[0];
	int j = particle_num[1];

	double n [3];
	n[0] = nr_vec.x;
	n[1] = nr_vec.y;
	n[2] = nr_vec.z;

	double stresslet_i[5];
	double stresslet_j[5];
//	for (int k=0; k < 5; k ++){
//		stresslet_i[k] = 0.;
//		stresslet_j[k] = 0.;
//	}
	
	// First -G*(U-Uinf) term
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
	double n0n0_13 = ( n[0] * n[0] - 1./3 );
	double n1n1_13 = ( n[1] * n[1] - 1./3 );
	double n0n1 = n[0] * n[1];
	double n0n2 = n[0] * n[2];
	double n1n2 = n[1] * n[2];

	double twothird = 2./3;
	double onesixth = 1./6;
	double common_factor_i = 0.;
	double common_factor_j = 0.;
	for(int u=0; u<3; u++){
		common_factor_i += n[u] * ( twothird * a0 * a0 * XGii * vi[u] + onesixth * ro * ro * XGij * vj[u] );
		common_factor_j += n[u] * ( twothird * a1 * a1 * XGjj * vj[u] + onesixth * ro * ro * XGji * vi[u] );
	}

	stresslet_i[0] = n0n0_13 * common_factor_i;
	stresslet_i[1] = n0n1 * common_factor_i;
	stresslet_i[2] = n0n2 * common_factor_i;
	stresslet_i[3] = n1n2 * common_factor_i;
	stresslet_i[4] = n1n1_13 * common_factor_i;

	stresslet_j[0] = n0n0_13 * common_factor_j;
	stresslet_j[1] = n0n1 * common_factor_j;
	stresslet_j[2] = n0n2 * common_factor_j;
	stresslet_j[3] = n1n2 * common_factor_j;
	stresslet_j[4] = n1n1_13 * common_factor_j;

	// Second: +M*Einf term
	double XMii, XMjj, XMij, XMji;
	XM(XMii, XMij, XMji, XMjj);

	double five24 = 5./24.;
	double fivethird = 5./3.;

	common_factor_i = n0n2 * ( fivethird * a0 * a0 * XMii + five24 * ro * ro * XMij ) * sys->shear_rate;
	common_factor_j = n0n2 * ( fivethird * a1 * a1 * XMjj + five24 * ro * ro * XMji ) * sys->shear_rate;

	stresslet_i[0] += n0n0_13 * common_factor_i;
	stresslet_i[1] += n0n1 * common_factor_i;
	stresslet_i[2] += n0n2 * common_factor_i;
	stresslet_i[3] += n1n2 * common_factor_i;
	stresslet_i[4] += n1n1_13 * common_factor_i;
	stresslet_j[0] += n0n0_13 * common_factor_j;
	stresslet_j[1] += n0n1 * common_factor_j;
	stresslet_j[2] += n0n2 * common_factor_j;
	stresslet_j[3] += n1n2 * common_factor_j;
	stresslet_j[4] += n1n1_13 * common_factor_j;

	for (int k=0; k < 5; k++){
		sys->lubstress[i][k] += stresslet_i[k];
		sys->lubstress[j][k] += stresslet_j[k];
	}
}

void
Interaction::evaluateLubricationForce(){
	int i = particle_num[0];
	int j = particle_num[1];
	
	double n [3];
	n[0] = nr_vec.x;
	n[1] = nr_vec.y;
	n[2] = nr_vec.z;
		
	//double stresslet_j[5];
	
	lubforce_i.reset();
	lubforce_j.reset();
	
	
	// First -G*(U-Uinf) term
	double vi [3];
	double vj [3];
	
	vi[0] = sys->relative_velocity[i].x;
	vi[1] = sys->relative_velocity[i].y;
	vi[2] = sys->relative_velocity[i].z;
	
	vj[0] = sys->relative_velocity[j].x;
	vj[1] = sys->relative_velocity[j].y;
	vj[2] = sys->relative_velocity[j].z;
	
	double XAii, XAij, XAji, XAjj;
	XA(XAii, XAij, XAji, XAjj);
	double common_factor_i = 0;
	double common_factor_j = 0;
	for(int u=0; u<3; u++){
		common_factor_i += a0*XAii*n[u]*vi[u] + 0.5*ro*XAij*n[u]*vj[u];
		common_factor_j += a1*XAjj*n[u]*vj[u] * 0.5*ro*XAji*n[u]*vi[u];
	}


	lubforce_i.x = - n[0]*common_factor_i;
	lubforce_i.y = - n[1]*common_factor_i;
	lubforce_i.z = - n[2]*common_factor_i;
	lubforce_j.x = - n[0]*common_factor_j;
	lubforce_j.y = - n[1]*common_factor_j;
	lubforce_j.z = - n[2]*common_factor_j;
	

//	stresslet_i[0] += n0n0_13 * common_factor_i;

	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double n0n2 = n[0] * n[2];
	
	double twothird = 2./3;
	double onesixth = 1./6;
	
	common_factor_i = (twothird*a0*a0*XGii + onesixth*ro*ro*XGji)* n0n2* sys->shear_rate;
	common_factor_j = (twothird*a1*a1*XGjj + onesixth*ro*ro*XGij)* n0n2* sys->shear_rate;
		

	lubforce_i.x +=  n[0]*common_factor_i;
	lubforce_i.y +=  n[1]*common_factor_i;
	lubforce_i.z +=  n[2]*common_factor_i;
	lubforce_j.x +=  n[0]*common_factor_j;
	lubforce_j.y +=  n[1]*common_factor_j;
	lubforce_j.z +=  n[2]*common_factor_j;

}

double
Interaction::valLubForce(){
	return 	-dot(lubforce_i , nr_vec);
}


void
Interaction::addContactStress(){
	int i = particle_num[0];
	int j = particle_num[1];
	if (contact){
		vec3d force = - Fc_normal * nr_vec + Fc_tangent; //@@TO BE CHECKED.
        
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
Interaction::activate(int i, int j, const vec3d &pos_diff, double distance, int zshift){

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


	a0 = sys->radius[particle_num[0]];
	a1 = sys->radius[particle_num[1]];
	ro = a0+a1; // for polydispesity, we will rewrite this to a1+a2
	if (distance > ro)
		contact = false;
	else
		contact = true;
	lambda = a1 / a0;
	invlambda = 1. / lambda;
	ksi_cutoff = sys->gap_cutoff;
	r_lub_max = 0.5*ro*sys->lub_max;
	assignDistanceNormalVector(pos_diff, distance, zshift);

	strain_0 = sys->shear_strain;
	if ( distance - ro < sys->dist_near*0.5*ro){
		near = true;
	} else {
		near = false;
	}
	return;
}

void
Interaction::deactivate(){
	// r > lub_max
	if (sys->output_trajectory)
		outputTrajectory();
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

/*
 * update()
 *  return `true' if r > r_lub_max 
 *   ---> deactivate
 */
bool
Interaction::update(){
	if (active){
		// update tangential displacement: we do it before updating nr_vec
		if (sys->friction && contact) {
			incrementContactTangentialDisplacement();
		}
		// compute new r_vec and distance
		calcDistanceNormalVector();
		// check new state of the interaction
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
				if (r() <= ro){
					activate_contact();
				} 
			}
		}
		
		if ( near == false ){
			if (ksi < sys->dist_near*0.5*ro){
				near = true;
			}
		} else {
			if (ksi > sys->dist_near*0.5*ro){
				near = false;
			}
		}

	}
	return false;
}

/*
 * Just provisional function for future analysis
 *
 */
void
Interaction::outputTrajectory(){
	for (int k=0; k < trajectory.size(); k++){
		sys->fout_trajectory << trajectory[k].x << ' ' <<  trajectory[k].y << ' '<<  trajectory[k].z << ' ';
		sys->fout_trajectory << gap_history[k] << endl;
	}
	sys->fout_trajectory << endl ;
	trajectory.clear();
	gap_history.clear();
}

double
Interaction::age(){
	return sys->shear_strain - strain_0;
}

void
Interaction::recordTrajectory(){
	trajectory.push_back(r_vec);
	gap_history.push_back(ksi);
}
