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
	prev_iksi_eff = iksi_eff;
}

void
Interaction::prev_XA(double &prev_XAii, double &prev_XAij, double &prev_XAji, double &prev_XAjj){
	double g1_l;
	double l1, l13;

	l1 = 1.0 + lambda;
	l13 = l1 * l1 * l1;

	
	g1_l = 2.0 * lambda * lambda / l13;
	
	prev_XAii = g1_l * prev_iksi_eff;
	prev_XAij = - 2 * prev_XAii / l1;
	prev_XAji = prev_XAij;
	prev_XAjj = prev_XAii / lambda;

}


void
Interaction::XG(double &XGii, double &XGij, double &XGji, double &XGjj){
	double l1 = 1.0 + lambda;
	double l13 = l1 * l1 * l1;
	double il1 = 1.0 + invlambda;
	//	il13 = il1 * il1 * il1;
	double g1_l = 2.0 * lambda * lambda / l13;
	//g1_il = 2.0 * invlambda * invlambda / il13;
	//g1_il = g1_l/lambda;
	
	XGii = 1.5 * g1_l * iksi_eff;
	XGij = - 4 * XGii / l1 / l1 ;
	//	XGjj = - 1.5 * g1_il * iksi_eff;
	XGjj = - XGii / lambda;
	XGji = - 4 * XGjj / il1 / il1 ;
}

void
Interaction::XM(double &XMii, double &XMij, double &XMji, double &XMjj){
	double l1 = 1.0 + lambda;
	double l13 = l1 * l1 * l1;
	//	il1 = 1.0 + invlambda;
	//	il13 = il1 * il1 * il1;
	double g1_l = 2.0 * lambda * lambda / l13;
	//g1_il = 2.0 * invlambda * invlambda / il13;
	//g1_il = g1_l/lambda;
	XMii = 0.6 * g1_l * iksi_eff;
	XMij = 40*lambda / 3 / l13 * XMii ;
	XMji = XMij;
	XMjj = XMii/lambda;
}

void
Interaction::GE(double GEi[], double GEj[]){
	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double nxnz_sr = nr_vec.x * nr_vec.z * sys->shear_rate;
	double common_factor_1 = (a0*a0*XGii*4 + ro*ro*XGji)/6;
	GEi[0] = common_factor_1 * nxnz_sr * nr_vec.x;
	GEi[1] = common_factor_1 * nxnz_sr * nr_vec.y;
	GEi[2] = common_factor_1 * nxnz_sr * nr_vec.z;
	double common_factor_2 = (a1*a1*XGjj*4 + ro*ro*XGij)/6;
	GEj[0] = common_factor_2 * nxnz_sr * nr_vec.x;
	GEj[1] = common_factor_2 * nxnz_sr * nr_vec.y;
	GEj[2] = common_factor_2 * nxnz_sr * nr_vec.z;
}


// computes the contribution to S = R_SU * V (in Brady's notations) [ S = G V in Jefrrey's ones ]
// from pair (i,j).
// ie fills :
// stresslet_i = R_SU^{ii} * vi + R_SU^{ij} * vj
// stresslet_j = R_SU^{ji} * vi + R_SU^{jj} * vj
 
void
Interaction::pairStresslet(const vec3d &vi, const vec3d &vj, stresslet &stresslet_i, stresslet &stresslet_j){
	for (int u=0; u < 5; u ++){
		stresslet_i.elm[u] = 0.;
		stresslet_j.elm[u] = 0.;
	}
	
	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double n0n0_13 = ( nr_vec.x * nr_vec.x - 1./3 );
	double n1n1_13 = ( nr_vec.y * nr_vec.y - 1./3 );
	double n0n1 = nr_vec.x * nr_vec.y;
	double n0n2 = nr_vec.x * nr_vec.z;
	double n1n2 = nr_vec.y * nr_vec.z;

	double twothird = 2./3;
	double onesixth = 1./6;
	double common_factor_i = 0;
	double common_factor_j = 0;

	common_factor_i += dot(nr_vec, ( twothird * a0 * a0 * XGii * vi + onesixth * ro * ro * XGij * vj));
	common_factor_j += dot(nr_vec, ( twothird * a1 * a1 * XGjj * vj + onesixth * ro * ro * XGji * vi));
	

	stresslet_i.elm[0] += n0n0_13 * common_factor_i;
	stresslet_i.elm[1] += n0n1    * common_factor_i;
	stresslet_i.elm[2] += n0n2    * common_factor_i;
	stresslet_i.elm[3] += n1n2    * common_factor_i;
	stresslet_i.elm[4] += n1n1_13 * common_factor_i;

	stresslet_j.elm[0] += n0n0_13 * common_factor_j;
	stresslet_j.elm[1] += n0n1 * common_factor_j;
	stresslet_j.elm[2] += n0n2 * common_factor_j;
	stresslet_j.elm[3] += n1n2 * common_factor_j;
	stresslet_j.elm[4] += n1n1_13 * common_factor_j;

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
	stresslet stresslet_i;
	stresslet stresslet_j;
	vec3d &vi = sys->relative_velocity[i];
	vec3d &vj = sys->relative_velocity[j];
	/*
	 *  First: -G*(U-Uinf) term
	 */
	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double n0n0_13 = nr_vec.x*nr_vec.x - 1./3;
	double n1n1_13 = nr_vec.y*nr_vec.y - 1./3;
	double n0n1 = nr_vec.x*nr_vec.y;
	double n0n2 = nr_vec.x*nr_vec.z;
	double n1n2 = nr_vec.y*nr_vec.z;
	double roro = ro*ro;
	double a0a0 = a0*a0;
	double a1a1 = a1*a1;

	double cf_GU_i = -dot(nr_vec, (2*a0a0*XGii/3)*vi + (roro*XGij/6)*vj);
	double cf_GU_j = -dot(nr_vec, (2*a1a1*XGjj/3)*vj + (roro*XGji/6)*vi);
	/*
	 *  Second: +M*Einf term
	 */
	double XMii, XMjj, XMij, XMji;
	XM(XMii, XMij, XMji, XMjj);
	double cf_ME_i = 5*(a0*a0a0*XMii/3 + ro*roro*XMij/24)*n0n2*sys->shear_rate;
	double cf_ME_j = 5*(a1*a1a1*XMjj/3 + ro*roro*XMji/24)*n0n2*sys->shear_rate;
	
	double common_factor_i = cf_GU_i + cf_ME_i;
	double common_factor_j = cf_GU_j + cf_ME_j;
	stresslet_i.elm[0] = n0n0_13 * common_factor_i;
	stresslet_i.elm[1] = n0n1 * common_factor_i;
	stresslet_i.elm[2] = n0n2 * common_factor_i;
	stresslet_i.elm[3] = n1n2 * common_factor_i;
	stresslet_i.elm[4] = n1n1_13 * common_factor_i;
	stresslet_j.elm[0] = n0n0_13 * common_factor_j;
	stresslet_j.elm[1] = n0n1 * common_factor_j;
	stresslet_j.elm[2] = n0n2 * common_factor_j;
	stresslet_j.elm[3] = n1n2 * common_factor_j;
	stresslet_j.elm[4] = n1n1_13 * common_factor_j;

	for (int u=0; u < 5; u++){
		sys->lubstress[i].elm[u] += stresslet_i.elm[u];
		sys->lubstress[j].elm[u] += stresslet_j.elm[u];
		lubstresslet.elm[u] = stresslet_i.elm[u] + stresslet_j.elm[u];
	}
	/*
	 * Calculate the lubrication stresslet by using lubrication forces
	 */
	
	stresslet stresslet2_i;
	stresslet stresslet2_j;
	vec3d dr_i = a0*nr_vec;
	stresslet2_i.elm[0] = 2*dr_i.x*lubforce_i.x;
	stresslet2_i.elm[1] = dr_i.x*lubforce_i.y + dr_i.y*lubforce_i.x;
	stresslet2_i.elm[2] = dr_i.x*lubforce_i.z + dr_i.z*lubforce_i.x;
	stresslet2_i.elm[3] = dr_i.y*lubforce_i.z + dr_i.z*lubforce_i.y;
	stresslet2_i.elm[4] = 2*dr_i.y*lubforce_i.y;
	vec3d lubforce_j = -lubforce_i;
	vec3d dr_j = -a1*nr_vec;
	stresslet2_j.elm[0] = 2*dr_j.x*lubforce_j.x;
	stresslet2_j.elm[1] = dr_j.x*lubforce_j.y + dr_j.y*lubforce_j.x;
	stresslet2_j.elm[2] = dr_j.x*lubforce_j.z + dr_j.z*lubforce_j.x;
	stresslet2_j.elm[3] = dr_j.y*lubforce_j.z + dr_j.z*lubforce_j.y;
	stresslet2_j.elm[4] = 2*dr_j.y*lubforce_j.y;

	for (int u=0; u < 5; u++){
		sys->lubstress2[i].elm[u] += stresslet2_i.elm[u];
		sys->lubstress2[j].elm[u] += stresslet2_j.elm[u];
	}
}

/* Lubriction force between two particles is calcurated.
 * This part is used for ouput data.
 * lubforce_j = -lubforce_i
 */
void
Interaction::evaluateLubricationForce(){
	int i = particle_num[0];
	int j = particle_num[1];
	lubforce_i.reset();
	/*	lubforce_j.reset(); */

	/*
	 *  First: -A*(U-Uinf) term
	 */
	vec3d &vi = sys->relative_velocity[i];
	vec3d &vj = sys->relative_velocity[j];
	double XAii, XAij, XAji, XAjj;
	XA(XAii, XAij, XAji, XAjj);
	double cf_AU_i = -dot(a0*XAii*vi + 0.5*ro*XAij*vj, nr_vec);
	/*	double cf_AU_j = -dot(a1*XAjj*vj + 0.5*ro*XAji*vi, nr_vec); */

	/*
	 *  Second -tildeG*(-Einf)term
	 */
	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double n0n2 = nr_vec.x * nr_vec.z;
	double cf_GE_i = n0n2*(2*a0*a0*XGii/3 + ro*ro*XGji/6)* sys->shear_rate;
	/*	double cf_GE_j = n0n2*(2*a1*a1*XGjj/3 + ro*ro*XGij/6)* sys->shear_rate; */
	lubforce_i = (cf_AU_i + cf_GE_i) * nr_vec;
	/* lubforce_j +=  common_factor_j * nr_vec; */
}

double
Interaction::valLubForce(){
	return 	-dot(lubforce_i, nr_vec);
}

void
Interaction::addContactStress(){
	int i = particle_num[0];
	int j = particle_num[1];
	if (contact){
		vec3d force = - Fc_normal * nr_vec + Fc_tangent; //@@TO BE CHECKED.
		contactstresslet.elm[0] = 2*(force.x * nr_vec.x); //xx
		contactstresslet.elm[1] = force.x * nr_vec.y + force.y * nr_vec.x ; //xy
		contactstresslet.elm[2] = force.x * nr_vec.z + force.z * nr_vec.x ; //yy
		contactstresslet.elm[3] = force.y * nr_vec.z + force.z * nr_vec.y ; //xz
		contactstresslet.elm[4] = 2*(force.y * nr_vec.y);
		for (int u=0; u < 5; u++){
			sys->contactstress[i].elm[u] += contactstresslet.elm[u];
			sys->contactstress[j].elm[u] += contactstresslet.elm[u];
		}
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
	invlambda = 1 / lambda;
	ksi_cutoff = sys->gap_cutoff;
	r_lub_max = 0.5*ro*sys->lub_max;
	assignDistanceNormalVector(pos_diff, distance, zshift);

	/*
	 * We may consider the particle size
	 *
	 */
	kn = sys->kn/(0.5*ro);
	kt = sys->kt/(0.5*ro);
	
	strain_0 = sys->shear_strain;
	if ( distance - ro < sys->dist_near*0.5*ro){
		near = true;
	} else {
		near = false;
	}

	prev_iksi_eff = 0;
	return;
}

bool 
Interaction::just_off(){
  if(_just_off){
	_just_off = false;
	return true;
  }
  else
	return false;
}

void
Interaction::deactivate(){
	// r > lub_max
	if (sys->out_pairtrajectory)
		outputTrajectory();
	active = false;
	int i=particle_num[0];
	int j=particle_num[1];
	sys->interaction_list[i].erase(this);
	sys->interaction_list[j].erase(this);
	sys->interaction_partners[i].erase(j);
	sys->interaction_partners[j].erase(i);
	_just_off = true; // for resistance matrix update
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
	Fc_normal = 0;
	Fc_tangent.reset();
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
