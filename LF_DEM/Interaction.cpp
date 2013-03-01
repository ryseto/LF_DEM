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
	Fc_normal_norm = 0;
	static_friction = false;
}


void 
Interaction::r(double new_r){
	_r = new_r;
	_gap_nondim = 2 * _r / ro - 2.;
	
	if( _gap_nondim > 0 ){
		lub_coeff = 1/(_gap_nondim + lub_reduce_parameter);
	}else{
		lub_coeff = lub_coeff_max;
	}
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1. ( i --> j)
 * pd_z : Periodic boundary condition
 */

void
Interaction::calcDistanceNormalVector(){
	r_vec = sys->position[particle_num[1]] - sys->position[particle_num[0]];
	sys->periodize_diff(r_vec, zshift);
	r(r_vec.norm());
	nr_vec = r_vec / r();
}

//void
//Interaction::assignDistanceNormalVector(const vec3d &pos_diff, double distance, int zshift_){
//	r_vec = pos_diff;
//	r(distance);
//	nr_vec = r_vec / r();
//	zshift = zshift_;
//}


/*********************************
*                                *
*	   Contact Forces Methods    *
*                                *
*********************************/
/*
 * Calculate interaction.
 * Force acts on particle 0 from particle 1.
 * r_vec = p[1] - p[0]
 * Fc_normal_norm is positive (for overlapping particles r < ro)
 */
void
Interaction::calcContactInteraction(){

	int normal_potential_type = 1;
	switch (normal_potential_type) {
	case 0:
		Fc_normal_norm = sys->kn*(ro - _r);
		break;
	case 1:
		Fc_normal_norm = exp(100.*(ro - _r)) - 1. ;
	
	}
	Fc_normal_norm /= sys->shear_rate; // to have it nondimensionalized
	Fc_normal = - Fc_normal_norm*nr_vec;
	
	if(static_friction)
		Fc_tan = sys_kt*disp_tan;

}

void
Interaction::addUpContactForceTorque(){
	if (contact){
		sys->contact_force[particle_num[0]] += Fc_normal;
		sys->contact_force[particle_num[1]] -= Fc_normal;

		if(static_friction){
			sys->contact_force[particle_num[0]] += Fc_tan;
			sys->contact_force[particle_num[1]] -= Fc_tan;

			vec3d t_ij = cross(nr_vec, Fc_tan);
			sys->contact_torque[particle_num[0]] += t_ij;
			sys->contact_torque[particle_num[1]] += t_ij;
		}
	}
}

/* Relative velocity of particle 1 from particle 0.
 *
 * Use: 
 *  sys->velocity and ang_velocity
 *
 */
void
Interaction::calcContactVelocity(){
	// relative velocity particle 1 from particle 0.
	contact_velocity = sys->velocity[particle_num[1]] - sys->velocity[particle_num[0]];
	if (zshift != 0){
		/*
		 * v1' = v1 - Lz = v1 - zshift*lz;
		 */
		/**** NOTE ********************************************
		 * In the Corrector, this contact_velocity
		 * is also the correcting velocity.
		 * This correcting velocity should not involve the
		 * velocity diffrence due to crossing the z boundary.
		 * fix_interaction_status = true : in the Predictor
		 * fix_interaction_status = false : in the Corrector
		 *
		 * if p1 is upper, zshift = -1.
		 * zshift = -1; //  p1 (z ~ lz), p0 (z ~ 0)
		 *
		 ******************************************************/
		if (sys->integration_method == 0 ){
			/* In the Eular Method
			 */
			contact_velocity.x += zshift * sys->vel_difference;
		} else {
			/* In the Predictor-Corrector Method
			 */
			if (sys->in_predictor){
				contact_velocity.x += zshift * sys->vel_difference;
			}
		}
	}
	contact_velocity -= cross(a0*sys->ang_velocity[particle_num[0]]
							  + a1*sys->ang_velocity[particle_num[1]],
							  nr_vec);

}

void
Interaction::incrementContactTangentialDisplacement(){
	disp_tan += contact_velocity*sys->dt;
	disp_tan -= dot(disp_tan, nr_vec)*nr_vec;// projection
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

	g1_l = 2.0 * lambda * lambda / l13;
	
	XAii = g1_l * lub_coeff;
	XAij = - 2 * XAii / l1;
	XAji = XAij;
	XAjj = XAii / lambda;

}

void
Interaction::XG(double &XGii, double &XGij, double &XGji, double &XGjj){
	double l1 = 1.0 + lambda;
	double l13 = l1 * l1 * l1;
	double il1 = 1.0 + invlambda;
	double g1_l = 2.0 * lambda * lambda / l13;
	XGii = 1.5 * g1_l * lub_coeff;
	XGij = - 4 * XGii / l1 / l1 ;
	XGjj = - XGii / lambda;
	XGji = - 4 * XGjj / il1 / il1 ;
}

void
Interaction::XM(double &XMii, double &XMij, double &XMji, double &XMjj){
	double l1 = 1.0 + lambda;
	double l13 = l1 * l1 * l1;
	double g1_l = 2.0 * lambda * lambda / l13;

	XMii = 0.6 * g1_l * lub_coeff;
	XMij = 8*lambda / l13 * XMii ;
	XMji = XMij;
	XMjj = XMii/lambda;
}

void
Interaction::GE(double GEi[], double GEj[]){
	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double nxnz_sr = nr_vec.x * nr_vec.z;
	double common_factor_1 = (a0*a0*XGii*4 + ro*ro*XGji)/6;
	GEi[0] = common_factor_1 * nxnz_sr * nr_vec.x;
	GEi[1] = common_factor_1 * nxnz_sr * nr_vec.y;
	GEi[2] = common_factor_1 * nxnz_sr * nr_vec.z;
	double common_factor_2 = (a1*a1*XGjj*4 + ro*ro*XGij)/6;
	GEj[0] = common_factor_2 * nxnz_sr * nr_vec.x;
	GEj[1] = common_factor_2 * nxnz_sr * nr_vec.y;
	GEj[2] = common_factor_2 * nxnz_sr * nr_vec.z;
}

// computes the contribution to S = R_SU * V (in Brady's notations) [ S = G V in Jeffrey's ones ]
// from pair (i,j).
// ie fills :
// stresslet_i = R_SU^{ii} * vi + R_SU^{ij} * vj
// stresslet_j = R_SU^{ji} * vi + R_SU^{jj} * vj
 
void
Interaction::pairVelocityStresslet(const vec3d &vi, const vec3d &vj,
								   stresslet &stresslet_i, stresslet &stresslet_j){
	
	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double n0n0_13 = ( nr_vec.x * nr_vec.x - 1./3 );
	double n1n1_13 = ( nr_vec.y * nr_vec.y - 1./3 );
	double n0n1 = nr_vec.x * nr_vec.y;
	double n0n2 = nr_vec.x * nr_vec.z;
	double n1n2 = nr_vec.y * nr_vec.z;

	double common_factor_i = - dot(nr_vec, ( 2. * a0 * a0 * XGii * vi / 3. +  ro * ro * XGij * vj / 6. ));
	double common_factor_j = - dot(nr_vec, ( 2. * a1 * a1 * XGjj * vj / 3. +  ro * ro * XGji * vi / 6. ));

	stresslet_i.elm[0] = n0n0_13 * common_factor_i;
	stresslet_i.elm[1] = n0n1    * common_factor_i;
	stresslet_i.elm[2] = n0n2    * common_factor_i;
	stresslet_i.elm[3] = n1n2    * common_factor_i;
	stresslet_i.elm[4] = n1n1_13 * common_factor_i;

	stresslet_j.elm[0] = n0n0_13 * common_factor_j;
	stresslet_j.elm[1] = n0n1 * common_factor_j;
	stresslet_j.elm[2] = n0n2 * common_factor_j;
	stresslet_j.elm[3] = n1n2 * common_factor_j;
	stresslet_j.elm[4] = n1n1_13 * common_factor_j;

}


// convenient interface for pairVelocityStresslet(const vec3d &vi, const vec3d &vj, stresslet &stresslet_i, stresslet &stresslet_j)
void
Interaction::pairVelocityStresslet(double* &vel_array, stresslet &stresslet_i, stresslet &stresslet_j){
	
	vec3d vi, vj;

	int i = particle_num[0];
	int j = particle_num[1];
	int i3 = 3*i;
	int j3 = 3*j;
	
	vi.x = vel_array[i3  ];
	vi.y = vel_array[i3+1];
	vi.z = vel_array[i3+2];
	
	vj.x = vel_array[j3  ];
	vj.y = vel_array[j3+1];
	vj.z = vel_array[j3+2];
	
	pairVelocityStresslet(vi, vj, stresslet_i, stresslet_j);
}


void
Interaction::pairStrainStresslet(stresslet &stresslet_i, stresslet &stresslet_j){
	double n0n0_13 = nr_vec.x*nr_vec.x - 1./3;
	double n1n1_13 = nr_vec.y*nr_vec.y - 1./3;
	double n0n1 = nr_vec.x*nr_vec.y;
	double n0n2 = nr_vec.x*nr_vec.z;
	double n1n2 = nr_vec.y*nr_vec.z;
	double roro = ro*ro;
	double a0a0 = a0*a0;
	double a1a1 = a1*a1;

	double XMii, XMjj, XMij, XMji;
	XM(XMii, XMij, XMji, XMjj);
	double common_factor_i = 5*(a0*a0a0*XMii/3 + ro*roro*XMij/24)*n0n2;
	double common_factor_j = 5*(a1*a1a1*XMjj/3 + ro*roro*XMji/24)*n0n2;
	
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

	/*
	 * Calculate the lubrication stresslet by using lubrication forces
	 */
	
	stresslet stresslet2_i;
	stresslet stresslet2_j;
	vec3d dr_i = a0*nr_vec;
	//	vec3d dr_i = r_vec;
	stresslet2_i.elm[0] = dr_i.x*lubforce_i.x;
	stresslet2_i.elm[1] = 0.5*(dr_i.x*lubforce_i.y + dr_i.y*lubforce_i.x);
	stresslet2_i.elm[2] = 0.5*(dr_i.x*lubforce_i.z + dr_i.z*lubforce_i.x);
	stresslet2_i.elm[3] = 0.5*(dr_i.y*lubforce_i.z + dr_i.z*lubforce_i.y);
	stresslet2_i.elm[4] = dr_i.y*lubforce_i.y;
	vec3d lubforce_j = -lubforce_i;
	vec3d dr_j = -a1*nr_vec;
	//	vec3d dr_j = -r_vec;
	stresslet2_j.elm[0] = dr_j.x*lubforce_j.x;
	stresslet2_j.elm[1] = 0.5*(dr_j.x*lubforce_j.y + dr_j.y*lubforce_j.x);
	stresslet2_j.elm[2] = 0.5*(dr_j.x*lubforce_j.z + dr_j.z*lubforce_j.x);
	stresslet2_j.elm[3] = 0.5*(dr_j.y*lubforce_j.z + dr_j.z*lubforce_j.y);
	stresslet2_j.elm[4] = dr_j.y*lubforce_j.y;

	for (int u=0; u < 5; u++){
		sys->lubstress2[i].elm[u] += stresslet2_i.elm[u];
		sys->lubstress2[j].elm[u] += stresslet2_j.elm[u];
	}
}


void
Interaction::addHydroStress(){
	int i = particle_num[0];
	int j = particle_num[1];
	int i3=3*i;
	int j3=3*j;

	stresslet stresslet_GU_i;
	stresslet stresslet_GU_j;
	stresslet stresslet_ME_i;
	stresslet stresslet_ME_j;

	vec3d *vi = new vec3d(sys->v_hydro[i3], sys->v_hydro[i3+1], sys->v_hydro[i3+2]);
	vec3d *vj = new vec3d(sys->v_hydro[j3], sys->v_hydro[j3+1], sys->v_hydro[j3+2]);
	// vec3d *vi = new vec3d(sys->v_hydro[i3]+sys->v_cont[i3], sys->v_hydro[i3+1]+sys->v_cont[i3+1], sys->v_hydro[i3+2]+sys->v_cont[i3+2]);
	// vec3d *vj = new vec3d(sys->v_hydro[j3]+sys->v_cont[j3], sys->v_hydro[j3+1]+sys->v_cont[j3+1], sys->v_hydro[j3+2]+sys->v_cont[j3+2]);
	/*
	 *  First: -G*(U-Uinf) term
	 */

	pairVelocityStresslet(*vi, *vj, stresslet_GU_i, stresslet_GU_j);

	/*
	 *  Second: +M*Einf term
	 */
	pairStrainStresslet(stresslet_ME_i, stresslet_ME_j);

	for (int u=0; u < 5; u++){
		sys->lubstress[i].elm[u] += stresslet_GU_i.elm[u] + stresslet_ME_i.elm[u];
		sys->lubstress[j].elm[u] += stresslet_GU_j.elm[u] + stresslet_ME_j.elm[u];
		lubstresslet.elm[u] = stresslet_GU_i.elm[u] + stresslet_ME_i.elm[u] + stresslet_GU_j.elm[u] + stresslet_ME_j.elm[u];
	}

	delete vi;
	delete vj;
}


/* Lubriction force between two particles is calculated. 
 * Note that only the Brownian component of the velocity is NOT included here.
 * This part is used for ouput data.
 * lubforce_j = -lubforce_i
 */
void
Interaction::evaluateLubricationForce(){
	int i = particle_num[0];
	int j = particle_num[1];
	lubforce_i.reset();

	/*
	 *  First: -A*(U-Uinf) term
	 */
	vec3d &vi = sys->relative_velocity_lub_cont[i];
	vec3d &vj = sys->relative_velocity_lub_cont[j];
	double XAii, XAij, XAji, XAjj;
	XA(XAii, XAij, XAji, XAjj);

	double cf_AU_i = -dot(a0*XAii*vi + 0.5*ro*XAij*vj, nr_vec);

	/*
	 *  Second -tildeG*(-Einf)term
	 */
	double XGii, XGjj, XGij, XGji;
	XG(XGii, XGij, XGji, XGjj);
	double n0n2 = nr_vec.x * nr_vec.z;
	double cf_GE_i = n0n2*(2*a0*a0*XGii/3 + ro*ro*XGji/6);

	lubforce_i = (cf_AU_i + cf_GE_i) * nr_vec;

}

double
Interaction::valLubForce(){
	return 	-dot(lubforce_i, nr_vec);
}



// term nr_vec*F
void
Interaction::calcContactStressTermXF(){

	if (contact){
		vec3d force = Fc_normal + Fc_tan;
		contactstresslet.elm[0] = force.x * nr_vec.x; //xx
		contactstresslet.elm[1] = 0.5*(force.x * nr_vec.y + force.y * nr_vec.x) ; //xy
		contactstresslet.elm[2] = 0.5*(force.x * nr_vec.z + force.z * nr_vec.x) ; //yy
		contactstresslet.elm[3] = 0.5*(force.y * nr_vec.z + force.z * nr_vec.y) ; //xz
		contactstresslet.elm[4] = force.y * nr_vec.y;
	}

}

void
Interaction::addContactStress(){

	int i = particle_num[0];
	int j = particle_num[1];
	int i3=3*i;
	int j3=3*j;
	if (contact){
		calcContactStressTermXF();
		for (int u=0; u < 5; u++){
			sys->contactstress[i].elm[u] += contactstresslet.elm[u];
			sys->contactstress[j].elm[u] += contactstresslet.elm[u];
		}
		

		// Add term G*V_cont
		stresslet stresslet_GU_i;
		stresslet stresslet_GU_j;
		vec3d *vi = new vec3d(sys->v_cont[i3], sys->v_cont[i3+1], sys->v_cont[i3+2]);
		vec3d *vj = new vec3d(sys->v_cont[j3], sys->v_cont[j3+1], sys->v_cont[j3+2]);
		pairVelocityStresslet(*vi, *vj, stresslet_GU_i, stresslet_GU_j);
		
		delete vi;
		delete vj;
		for (int u=0; u < 5; u++){
			sys->contactstress[i].elm[u] += stresslet_GU_i.elm[u];
			sys->contactstress[j].elm[u] += stresslet_GU_j.elm[u];
		}
	}

}

void
Interaction::addContactStress2(){
	int i = particle_num[0];
	int j = particle_num[1];

	if (contact){
		calcContactStressTermXF();
		for (int u=0; u < 5; u++){
			sys->contactstress2[i].elm[u] += contactstresslet.elm[u];
			sys->contactstress2[j].elm[u] += contactstresslet.elm[u];
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
Interaction::activate(int i, int j,
					  const vec3d &pos_diff,
					  double distance, int _zshift){
	active = true;

	Fc_normal_norm = 0.;	
	Fc_normal.reset();
	Fc_tan.reset();

	r_vec = pos_diff;
	r(distance);
	nr_vec = r_vec / r();
	zshift = _zshift;

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
	ro = a0 + a1;
	if (distance > ro)
		contact = false;
	else
		contact = true;
	lambda = a1 / a0;
	invlambda = 1 / lambda;
	lub_reduce_parameter = sys->lub_reduce_parameter;
	lub_coeff_max = 1/lub_reduce_parameter;
	
	r_lub_max = 0.5*ro*sys->lub_max;

	//assignDistanceNormalVector(pos_diff, distance, zshift);

	/*
	 * We may consider the particle size
	 *
	 */
	//	kn = sys->kn/(0.5*ro);
	//	kt = sys->kt/(0.5*ro);
	kn = sys->kn;
	kt = sys->kt;

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
	if (sys->friction){
		static_friction = true;
	}
	disp_tan.reset();
	calcContactVelocity();
}

void
Interaction::deactivate_contact(){
	// r > a0 + a1
	contact = false;
	static_friction = false;
	disp_tan.reset();
	Fc_normal_norm = 0;
	Fc_normal.reset();
	Fc_tan.reset();
}

/* 
 * updateState 
 * return value is true, this interaction will be broken up.
 */
bool
Interaction::updateState(){
	if(r() > r_lub_max ){
		deactivate();
		return true; // breakup
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
	return false;
}

/*
 * update()
 *  return `true' if r > r_lub_max 
 *   ---> deactivate
 */
bool
Interaction::updateStatesForceTorque(){
	if (active){
		// compute new r_vec and distance
		calcDistanceNormalVector();
		// update tangential displacement: we do it before updating nr_vec
		if (contact) {
			calcContactVelocity();
			incrementContactTangentialDisplacement();
			calcContactInteraction();
			if ( static_friction && (!sys->in_predictor) ){
				checkBreakupStaticFriction();
			}
		}
		// compute new contact forces if needed
		// check new state of the interaction
		if (!sys->in_predictor){
			return updateState();
		}
	}
	
	if(!static_friction && Fc_tan.sq_norm() > 0.){
		cerr << " Error : Interaction " << particle_num[0] << " " << particle_num[1] << " : " << endl;
		cerr << " No friction required but non-zero tangential force ( " << Fc_tan << " ) after update " << endl;
		exit(1);
	}
		
	return false;

}

void
Interaction::checkBreakupStaticFriction(){
	double f_static = sys->mu_static*Fc_normal_norm;
	if (Fc_tan.sq_norm() > f_static*f_static){
		/**
		 ** switch to dynamic friction
		 **
		 ** A simple imprementation is used temporary.
		 */
		disp_tan.reset();
		Fc_tan.reset();
	}
}
