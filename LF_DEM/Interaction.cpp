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
	colloidal_force_norm = 0;
	nearing_gapnd_cutoff = 0.01;
}

void 
Interaction::r(double new_r){
	_r = new_r;
	_gap_nondim = 2*_r/ro-2;
	if(_gap_nondim > 0){
		lub_coeff = 1/(_gap_nondim+sys->lub_reduce_parameter);
	}else{
		lub_coeff = sys->lub_coeff_contact;
	}
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1. ( i --> j)
 * pd_z : Periodic boundary condition
 */

void
Interaction::calcDistanceNormalVector(){
	r_vec = sys->position[par_num[1]]-sys->position[par_num[0]];
	sys->periodize_diff(r_vec, zshift);
	r(r_vec.norm());
	nr_vec = r_vec/r();
}

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
	if (sys->shearrate_scale_Fc_normal){
		//Fc_normal_norm = exp(kn*(ro-_r))-1;
		Fc_normal_norm = sys->kn*(ro-_r);
	} else {
		Fc_normal_norm = (exp(sys->kn*(ro-_r))-1)/sys->shear_rate;
	}
	Fc_normal = -Fc_normal_norm*nr_vec;
	if (sys->friction){
		Fc_tan = sys->kt*disp_tan;
	}
}

void
Interaction::addUpContactForceTorque(){
	if (contact){
		sys->contact_force[par_num[0]] += Fc_normal;
		sys->contact_force[par_num[1]] -= Fc_normal;
		if(sys->friction){
			sys->contact_force[par_num[0]] += Fc_tan;
			sys->contact_force[par_num[1]] -= Fc_tan;
			vec3d t_ij = cross(nr_vec, Fc_tan);
			sys->contact_torque[par_num[0]] += t_ij;
			sys->contact_torque[par_num[1]] += t_ij;
		}
	}
}

/*
 *
 *
 *
 */
void
Interaction::addUpColloidalForce(){
	if (contact == false){
		colloidal_force_norm = -sys->cf_amp_dl*ro_2*exp(-(_r-ro)/sys->cf_range);
		vec3d colloidal_force = colloidal_force_norm*nr_vec;
		sys->colloidal_force[par_num[0]] += colloidal_force;
		sys->colloidal_force[par_num[1]] -= colloidal_force;
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
	contact_velocity = sys->velocity[par_num[1]]-sys->velocity[par_num[0]];
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
		if (sys->integration_method == 0){
			/* In the Euler Method
			 */
			contact_velocity.x += zshift*sys->vel_difference;
		}else{
			/* In the Predictor-Corrector Method
			 */
			if (sys->in_predictor){
				contact_velocity.x += zshift*sys->vel_difference;
			}
		}
	}
	contact_velocity -= \
	cross(a0*sys->ang_velocity[par_num[0]]+a1*sys->ang_velocity[par_num[1]], \
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
Interaction::calcXA(){
	double l1 = 1+lambda;
	double l13 = l1*l1*l1;
	double g1_l = 2*lambda*lambda/l13;
	XA[0] = g1_l*lub_coeff;
	XA[1] = -2*XA[0]/l1;
	XA[2] = XA[1];
	XA[3] = XA[0]/lambda;
}

void
Interaction::calcXG(){
	double l1 = 1+lambda;// defined in XA
	double l13 = l1*l1*l1; // defined in XA
	double g1_l = 2*lambda*lambda/l13; // defined in XA
	double il1 = invlambda+1;
	XG[0] = 1.5*g1_l*lub_coeff;
	XG[1] = -4*XG[0]/l1/l1;
	XG[3] = -XG[0]/lambda;
	XG[2] = -4*XG[3]/il1/il1;
}

void
Interaction::calcXM(){
	double l1 = 1+lambda;
	double l13 = l1*l1*l1;
	double g1_l = 2*lambda*lambda/l13;
	XM[0] = 0.6*g1_l*lub_coeff;
	XM[1] = 8*lambda/l13*XM[0];
	XM[2] = XM[1];
	XM[3] = XM[0]/lambda;
}

void
Interaction::GE(double *GEi, double *GEj){
	calcXG();
	double nxnz_sr = nr_vec.x*nr_vec.z;
	double common_factor_1 = (4*a0*a0*XG[0]+ro*ro*XG[2])/6;
	GEi[0] = common_factor_1*nxnz_sr*nr_vec.x;
	GEi[1] = common_factor_1*nxnz_sr*nr_vec.y;
	GEi[2] = common_factor_1*nxnz_sr*nr_vec.z;
	double common_factor_2 = (4*a1*a1*XG[3]+ro*ro*XG[1])/6;
	GEj[0] = common_factor_2*nxnz_sr*nr_vec.x;
	GEj[1] = common_factor_2*nxnz_sr*nr_vec.y;
	GEj[2] = common_factor_2*nxnz_sr*nr_vec.z;
}

// computes the contribution to S = R_SU * V (in Brady's notations) [ S = G V in Jeffrey's ones ]
// from pair (i,j).
// ie fills :
// stresslet_i = R_SU^{ii} * vi + R_SU^{ij} * vj
// stresslet_j = R_SU^{ji} * vi + R_SU^{jj} * vj
void
Interaction::pairVelocityStresslet(const vec3d &vi, const vec3d &vj,
								   stresslet &stresslet_i, stresslet &stresslet_j){
	calcXG();
	double n0n0_13 = nr_vec.x*nr_vec.x-1./3;
	double n1n1_13 = nr_vec.y*nr_vec.y-1./3;
	double n0n1 = nr_vec.x*nr_vec.y;
	double n0n2 = nr_vec.x*nr_vec.z;
	double n1n2 = nr_vec.y*nr_vec.z;
	double common_factor_i = -dot(nr_vec, (4*a0*a0*XG[0]*vi+ro*ro*XG[1]*vj)/6);
	double common_factor_j = -dot(nr_vec, (4*a1*a1*XG[3]*vj+ro*ro*XG[2]*vi)/6);
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
	int i3 = 3*par_num[0];
	int j3 = 3*par_num[1];
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

	calcXM();
	double common_factor_i = 5*(a0*a0a0*XM[0]/3+ro*roro*XM[1]/24)*n0n2;
	double common_factor_j = 5*(a1*a1a1*XM[3]/3+ro*roro*XM[2]/24)*n0n2;
	
	stresslet_i.elm[0] = n0n0_13*common_factor_i;
	stresslet_i.elm[1] = n0n1*common_factor_i;
	stresslet_i.elm[2] = n0n2*common_factor_i;
	stresslet_i.elm[3] = n1n2*common_factor_i;
	stresslet_i.elm[4] = n1n1_13*common_factor_i;

	stresslet_j.elm[0] = n0n0_13*common_factor_j;
	stresslet_j.elm[1] = n0n1*common_factor_j;
	stresslet_j.elm[2] = n0n2*common_factor_j;
	stresslet_j.elm[3] = n1n2*common_factor_j;
	stresslet_j.elm[4] = n1n1_13*common_factor_j;
}

void
Interaction::addHydroStress(){
	int i3 = 3*par_num[0];
	int j3 = 3*par_num[1];
	stresslet stresslet_GU_i;
	stresslet stresslet_GU_j;
	stresslet stresslet_ME_i;
	stresslet stresslet_ME_j;
	vec3d vi(sys->v_hydro[i3], sys->v_hydro[i3+1], sys->v_hydro[i3+2]);
	vec3d vj(sys->v_hydro[j3], sys->v_hydro[j3+1], sys->v_hydro[j3+2]);
	/*
	 *  First: -G*(U-Uinf) term
	 */
	pairVelocityStresslet(vi, vj, stresslet_GU_i, stresslet_GU_j);
	/*
	 *  Second: +M*Einf term
	 */
	pairStrainStresslet(stresslet_ME_i, stresslet_ME_j);
	for (int u=0; u<5; u++){
		sys->lubstress[par_num[0]].elm[u] += stresslet_GU_i.elm[u]+stresslet_ME_i.elm[u];
		sys->lubstress[par_num[1]].elm[u] += stresslet_GU_j.elm[u]+stresslet_ME_j.elm[u];
		lubstresslet.elm[u] = \
		stresslet_GU_i.elm[u]+stresslet_ME_i.elm[u]+stresslet_GU_j.elm[u]+stresslet_ME_j.elm[u];
	}
}



/* Lubriction force between two particles is calculated. 
 * Note that only the Brownian component of the velocity is NOT included here.
 * This part is used for ouput data.
 * lubforce_j = -lubforce_i
 */
void
Interaction::evaluateLubricationForce(){
	/*
	 *  First: -A*(U-Uinf) term
	 */
	vec3d vi = sys->velocity[par_num[0]];
	vec3d vj = sys->velocity[par_num[1]];
	vi.x -= sys->position[par_num[0]].z;
	vj.x -= sys->position[par_num[1]].z;
	calcXA();
	double cf_AU_i = -dot(a0*XA[0]*vi+0.5*ro*XA[1]*vj, nr_vec);
	/*
	 *  Second -tildeG*(-Einf)term
	 */
	calcXG();
	double cf_GE_i = nr_vec.x*nr_vec.z*(a0*a0*XG[0]*4+ro*ro*XG[2])/6;
	lubforce_i = (cf_AU_i+cf_GE_i)*nr_vec;
}

double
Interaction::valLubForce(){
	return -dot(lubforce_i, nr_vec);
}

// term nr_vec*F
void
Interaction::calcContactStressTermXF(){
	if (contact){
		vec3d force = Fc_normal+Fc_tan;
		contactstressletXF.elm[0] = force.x*nr_vec.x; //xx
		contactstressletXF.elm[1] = 0.5*(force.x*nr_vec.y+force.y*nr_vec.x); //xy
		contactstressletXF.elm[2] = 0.5*(force.x*nr_vec.z+force.z*nr_vec.x); //yy
		contactstressletXF.elm[3] = 0.5*(force.y*nr_vec.z+force.z*nr_vec.y); //xz
		contactstressletXF.elm[4] = force.y*nr_vec.y;
	}else{
		for (int k=0; k<5; k++){
			contactstressletXF.elm[k] = 0;
		}
	}
}

// term nr_vec*F
void
Interaction::calcColloidalStressTermXF(){
	if (contact == false){
		vec3d force = colloidal_force_norm*nr_vec;
		colloidalstressletXF.elm[0] = force.x*nr_vec.x; //xx
		colloidalstressletXF.elm[1] = 0.5*(force.x*nr_vec.y+force.y*nr_vec.x); //xy
		colloidalstressletXF.elm[2] = 0.5*(force.x*nr_vec.z+force.z*nr_vec.x); //yy
		colloidalstressletXF.elm[3] = 0.5*(force.y*nr_vec.z+force.z*nr_vec.y); //xz
		colloidalstressletXF.elm[4] = force.y*nr_vec.y;
	}else{
		for (int k=0; k<5; k++){
			colloidalstressletXF.elm[k] = 0;
		}
	}
}

void
Interaction::addContactStress(){
	int i3 = 3*par_num[0];
	int j3 = 3*par_num[1];
	if (contact){
		calcContactStressTermXF();
		for (int u=0; u<5; u++){
			sys->contactstressXF[par_num[0]].elm[u] += a0*contactstressletXF.elm[u];
			sys->contactstressXF[par_num[1]].elm[u] += a1*contactstressletXF.elm[u];
		}
		/* I think:
		 * In hard sphere model,
		 * the following stress contribution has no physical meaning.
		 * We should not add this term for the total viscosity.
		 */
		// Add term G*V_cont
		stresslet stresslet_GU_i;
		stresslet stresslet_GU_j;
		vec3d vi(sys->v_cont[i3], sys->v_cont[i3+1], sys->v_cont[i3+2]);
		vec3d vj(sys->v_cont[j3], sys->v_cont[j3+1], sys->v_cont[j3+2]);
		pairVelocityStresslet(vi, vj, stresslet_GU_i, stresslet_GU_j);
		for (int u=0; u<5; u++){
			sys->contactstressGU[par_num[0]].elm[u] += stresslet_GU_i.elm[u];
			sys->contactstressGU[par_num[1]].elm[u] += stresslet_GU_j.elm[u];
		}
	}
}

void
Interaction::addColloidalStress(){
	int i3 = 3*par_num[0];
	int j3 = 3*par_num[1];
	if (!contact){
		calcColloidalStressTermXF();
		for (int u=0; u<5; u++){
			sys->colloidalstressXF[par_num[0]].elm[u] += a0*colloidalstressletXF.elm[u];
			sys->colloidalstressXF[par_num[1]].elm[u] += a1*colloidalstressletXF.elm[u];
		}
		// Add term G*V_cont
		stresslet stresslet_colloid_GU_i;
		stresslet stresslet_colloid_GU_j;
		vec3d vi(sys->v_colloidal[i3], sys->v_colloidal[i3+1], sys->v_colloidal[i3+2]);
		vec3d vj(sys->v_colloidal[j3], sys->v_colloidal[j3+1], sys->v_colloidal[j3+2]);
		pairVelocityStresslet(vi, vj, stresslet_colloid_GU_i, stresslet_colloid_GU_j);
		for (int u=0; u<5; u++){
			sys->colloidalstressGU[par_num[0]].elm[u] += stresslet_colloid_GU_i.elm[u];
			sys->colloidalstressGU[par_num[1]].elm[u] += stresslet_colloid_GU_j.elm[u];
		}
	}
}

int
Interaction::partner(int i){
	if(i == par_num[0]){
		return par_num[1];
	}else{
		return par_num[0];
	}
}

/* Activate interaction between particles i and j.
 * Always j>i is satisfied.
 */
void
Interaction::activate(int i, int j,
					  const vec3d &pos_diff,
					  double distance, int _zshift){
	active = true;
	Fc_normal_norm = 0;
	Fc_normal.reset();
	Fc_tan.reset();
	r_vec = pos_diff;
	r(distance);
	nr_vec = r_vec/r();
	zshift = _zshift;
	if(j > i){
		par_num[0] = i;
		par_num[1] = j;
	}else{
		par_num[0] = j;
		par_num[1] = i;
	}		
	// tell it to particles i and j
	sys->interaction_list[i].insert(this);
	sys->interaction_list[j].insert(this);
	// tell them their new partner
	sys->interaction_partners[i].insert(j);
	sys->interaction_partners[j].insert(i);
	a0 = sys->radius[par_num[0]];
	a1 = sys->radius[par_num[1]];
	ro = a0+a1;
	ro_2 = ro/2;
	if (distance > ro){
		contact = false;
	}else{
		contact = true;
	}
	lambda = a1/a0;
	invlambda = 1/lambda;
	r_lub_max = 0.5*ro*sys->lub_max;
	return;
}

void
Interaction::deactivate(){
	// r > lub_max
	active = false;
	sys->interaction_list[par_num[0]].erase(this);
	sys->interaction_list[par_num[1]].erase(this);
	sys->interaction_partners[par_num[0]].erase(par_num[1]);
	sys->interaction_partners[par_num[1]].erase(par_num[0]);
}

void
Interaction::activate_contact(){
	// r < a0 + a1
	contact = true;
	disp_tan.reset();
	init_contact_time = sys->time();
	colloidal_force_norm = 0;
}

void
Interaction::deactivate_contact(){
	// r > a0 + a1
	contact = false;
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
	if(r() > r_lub_max){
		deactivate();
		return true; // breakup
	}else{
		if (contact){
			if (r() > ro){
				deactivate_contact();
			}
		}else{
			// contact false:
			if (r() <= ro){
				activate_contact();
			}
		}
		// nearing observable
		if(!nearing_on
		   && gap_nondim()<nearing_gapnd_cutoff){
			nearing_on = true;
			init_nearing_time = sys->time();
		}
		if(nearing_on
		   && gap_nondim()>nearing_gapnd_cutoff){
			nearing_on = false;
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
		// update tangential displacement: we do it before updating nr_vec
		// as it should be along the tangential vector defined in the previous time step
		if (contact){
			calcContactVelocity();
			incrementContactTangentialDisplacement();
		}
		// compute new r_vec and distance
		// z_shift is updated
		calcDistanceNormalVector();
		if (contact){
			calcContactInteraction();
			if (!sys->in_predictor){
				checkBreakupStaticFriction();
			}
		}
		// compute new contact forces if needed
		// check new state of the interaction
		if (!sys->in_predictor){
			return updateState();
		}
	}
	return false;
}

void
Interaction::checkBreakupStaticFriction(){
	/* Do we need to add the lubrication force
	 * for the friction law?
	 * -->
	 * It is better to add. 
	 * But, when beta is not huge,
	 * the difference may be neglegible.
	 */
	//	evaluateLubricationForce();
	//	double f_lub = -dot(lubforce_i, nr_vec);
	//	double f_static = sys->mu_static*(Fc_normal_norm+f_lub);
	double f_static = sys->mu_static*(Fc_normal_norm);
	if (Fc_tan.sq_norm() > f_static*f_static){
		/**
		 ** switch to dynamic friction
		 **
		 ** A simple imprementation is used temporary.
		 */
		disp_tan.reset();
	}
}

double 
Interaction::nearing_time(){
	double time = 0;
	if(nearing_on)
		time = sys->time()-init_nearing_time;
	return time;
}

double
Interaction::contact_time(){
	double time = 0;
	if(contact)
		time = sys->time()-init_contact_time;
	return time;
}
