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
	friction = sys->friction;

}


void 
Interaction::r(double new_r){
	_r = new_r;
	_gap = 2 * _r / ro - 2.;
	//	iksi = 1./ksi;
	if( _gap < gap_cutoff )
		gap_eff = gap_cutoff;
	else
		gap_eff = _gap;
	
	inv_gap_eff = 1./gap_eff;
}

/* Make a normal vector
 * Periodic boundaries are checked for all partices.
 * vector from particle 0 to particle 1. ( i --> j)
 * pd_z : Periodic boundary condition
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
		Fc_normal = sys->kn*(_r - ro);
		if(friction){
			if(static_friction){
				calcStaticFriction();
			} else {
				calcDynamicFriction();
			}
			vec3d t_ij = cross(nr_vec, Fc_tangent); // acting on p0
			Tc_0 = a0*t_ij;
			Tc_1 = a1*t_ij;
			//		cerr << Fc_tangent.x << ' ' << Fc_tangent.z << endl;
		}
		vec3d f_ij = Fc_normal * nr_vec + Fc_tangent; // acting on p0 //@@TO BE CHECKED.
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
Interaction::addUpContactForce(vec3d &force0, vec3d &force1){
	vec3d f_ij;
	if (contact){
		
		if(friction)
			f_ij = Fc_normal * nr_vec + Fc_tangent; // acting on p0 //@@TO BE CHECKED.
		else
			f_ij = Fc_normal * nr_vec;
		
		force0 += f_ij;
		force1 -= f_ij;
	}
}

void
Interaction::addUpContactTorque(vec3d &torque0, vec3d &torque1){
  if (contact&&friction){
	  torque0 += Tc_0;
	  torque1 += Tc_1;
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

	g1_l = 2.0 * lambda * lambda / l13;
	
	XAii = g1_l * inv_gap_eff;
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

	
	XGii = 1.5 * g1_l * inv_gap_eff;
	XGij = - 4 * XGii / l1 / l1 ;
	XGjj = - XGii / lambda;
	XGji = - 4 * XGjj / il1 / il1 ;
}

void
Interaction::XM(double &XMii, double &XMij, double &XMji, double &XMjj){
	double l1 = 1.0 + lambda;
	double l13 = l1 * l1 * l1;
	double g1_l = 2.0 * lambda * lambda / l13;

	XMii = 0.6 * g1_l * inv_gap_eff;
	XMij = 8*lambda / l13 * XMii ;
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


// computes the contribution to S = R_SU * V (in Brady's notations) [ S = G V in Jeffrey's ones ]
// from pair (i,j).
// ie fills :
// stresslet_i = R_SU^{ii} * vi + R_SU^{ij} * vj
// stresslet_j = R_SU^{ji} * vi + R_SU^{jj} * vj
 
void
Interaction::pairVelocityStresslet(const vec3d &vi, const vec3d &vj, stresslet &stresslet_i, stresslet &stresslet_j){
	
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
	double common_factor_i = 5*(a0*a0a0*XMii/3 + ro*roro*XMij/24)*n0n2*sys->shear_rate;
	double common_factor_j = 5*(a1*a1a1*XMjj/3 + ro*roro*XMji/24)*n0n2*sys->shear_rate;
	
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
	stresslet2_i.elm[0] = 2*dr_i.x*lubforce_i.x;
	stresslet2_i.elm[1] = dr_i.x*lubforce_i.y + dr_i.y*lubforce_i.x;
	stresslet2_i.elm[2] = dr_i.x*lubforce_i.z + dr_i.z*lubforce_i.x;
	stresslet2_i.elm[3] = dr_i.y*lubforce_i.z + dr_i.z*lubforce_i.y;
	stresslet2_i.elm[4] = 2*dr_i.y*lubforce_i.y;
	vec3d lubforce_j = -lubforce_i;
	vec3d dr_j = -a1*nr_vec;
	//	vec3d dr_j = -r_vec;
	stresslet2_j.elm[0] = 2*dr_j.x*lubforce_j.x;
	stresslet2_j.elm[1] = dr_j.x*lubforce_j.y + dr_j.y*lubforce_j.x;
	stresslet2_j.elm[2] = dr_j.x*lubforce_j.z + dr_j.z*lubforce_j.x;
	stresslet2_j.elm[3] = dr_j.y*lubforce_j.z + dr_j.z*lubforce_j.y;
	stresslet2_j.elm[4] = 2*dr_j.y*lubforce_j.y;

	for (int u=0; u < 5; u++){
		sys->lubstress2[i].elm[u] += stresslet2_i.elm[u];
		sys->lubstress2[j].elm[u] += stresslet2_j.elm[u];
	}

	//  vec3d &vi = sys->relative_velocity_lub_cont[i];
	//  vec3d &vj = sys->relative_velocity_lub_cont[j];
	// if(i==98){
	//   	cout << " version2 " <<i << " " << j << " " << stresslet2_i.elm[2] << " " <<0.25*a0*a0*dot(nr_vec, vi-vj)*inv_gap_eff*nr_vec.x*nr_vec.z<< " " << vi.x <<" " << vj.x <<  endl;
	// 	//		cout << " version2 " <<i << " " << j << " " << " " << vi.x <<" " << vj.x <<  endl;
	// }
	// if(j==98)
	//   	cout << " version2 " << i << " " << j << " " << stresslet2_j.elm[2] <<" " <<0.25*a0*a0*dot(nr_vec, vi-vj)*inv_gap_eff*nr_vec.x*nr_vec.z << " " << vi.x <<" " << vj.x << endl;

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

	// if(i==98){
	// // 	//	  	cout << i << " " << j << " "<< "stresslet " << stresslet_GU_i.elm[2] << " " << (3.*inv_gap_eff/8.)*(2/3.)*a0*a0*dot(nr_vec, *vi-*vj)*nr_vec.x*nr_vec.z   << endl;
	//  	cout << i << " " << j << " " <<"stresslet " << " " << stresslet_GU_i.elm[2]  << " " << 2*(stresslet_GU_i.elm[2]+stresslet_ME_i.elm[2]) << " " << vi->x <<" " << vj->x << " " << sys->relative_velocity_lub_cont[i].x << " " << sys->v_hydro[3*i] << " " << sys->v_cont[3*i] << " " << sys->v_hydro[3*i]+sys->v_cont[3*i] << " " << _gap<< " " << r() << endl;
	// 	//		cout << i << " " << j << " " << vi->x <<" " << vj->x <<endl;
	// }


	// if(j==98){
	// 	//	  	cout << "stresslet " << stresslet_GU_i.elm[2] << " " << (3.*inv_gap_eff/8.)*(2/3.)*a0*a0*dot(nr_vec, (*vi-*vj) )*nr_vec.x*nr_vec.z << endl;
	// 	cout << i << " " << j << " " << "stresslet " << stresslet_GU_j.elm[2] << " " << 2*(stresslet_GU_i.elm[2]+stresslet_ME_i.elm[2]) << " " << vi->x <<" " << vj->x << endl;
	// }

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
	/*	lubforce_j.reset(); */

	/*
	 *  First: -A*(U-Uinf) term
	 */
	vec3d &vi = sys->relative_velocity_lub_cont[i];
	vec3d &vj = sys->relative_velocity_lub_cont[j];
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

	lubforce_i = (cf_AU_i + cf_GE_i) * nr_vec; //---> original

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
	int i3=3*i;
	int j3=3*j;

	stresslet stresslet_GU_i;
	stresslet stresslet_GU_j;

	if (contact){
		vec3d force =  - Fc_normal * nr_vec + Fc_tangent; //@@TO BE CHECKED.
		contactstresslet.elm[0] = 2*(force.x * nr_vec.x); //xx
		contactstresslet.elm[1] = force.x * nr_vec.y + force.y * nr_vec.x ; //xy
		contactstresslet.elm[2] = force.x * nr_vec.z + force.z * nr_vec.x ; //yy
		contactstresslet.elm[3] = force.y * nr_vec.z + force.z * nr_vec.y ; //xz
		contactstresslet.elm[4] = 2*(force.y * nr_vec.y);
		for (int u=0; u < 5; u++){
			sys->contactstress[i].elm[u] += 0.5*contactstresslet.elm[u];
			sys->contactstress[j].elm[u] += 0.5*contactstresslet.elm[u];
		}

	}

	// Add term G*V_cont
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

void
Interaction::addContactStress2(){
	int i = particle_num[0];
	int j = particle_num[1];
	int i3=3*i;
	int j3=3*j;

	stresslet stresslet_GU_i;
	stresslet stresslet_GU_j;

	if (contact){
		vec3d force = - Fc_normal * nr_vec + Fc_tangent; //@@TO BE CHECKED.
		contactstresslet2.elm[0] = 2*(force.x * nr_vec.x); //xx
		contactstresslet2.elm[1] = force.x * nr_vec.y + force.y * nr_vec.x ; //xy
		contactstresslet2.elm[2] = force.x * nr_vec.z + force.z * nr_vec.x ; //yy
		contactstresslet2.elm[3] = force.y * nr_vec.z + force.z * nr_vec.y ; //xz
		contactstresslet2.elm[4] = 2*(force.y * nr_vec.y);
		for (int u=0; u < 5; u++){
			sys->contactstress2[i].elm[u] += contactstresslet2.elm[u];
			sys->contactstress2[j].elm[u] += contactstresslet2.elm[u];
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
	ro = a0+a1;
	if (distance > ro)
		contact = false;
	else
		contact = true;
	lambda = a1 / a0;
	invlambda = 1 / lambda;
	gap_cutoff = sys->gap_cutoff;
	r_lub_max = 0.5*ro*sys->lub_max;
	assignDistanceNormalVector(pos_diff, distance, zshift);

	/*
	 * We may consider the particle size
	 *
	 */
	//	kn = sys->kn/(0.5*ro);
	//	kt = sys->kt/(0.5*ro);
	kn = sys->kn;
	kt = sys->kt;
	/*
	 * Record the strain when this lub interaction starts.
	 */
	strain_lub_generated = sys->shear_strain;
	if ( distance - ro < sys->dist_near*0.5*ro){
		strain_near_contact = sys->shear_strain;
		near = true;
	} else {
		near = false;
	}
	return;
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


bool
Interaction::updateState(bool switch_off_allowed){
	
	if(r() > r_lub_max && switch_off_allowed){
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
	return false;
}

/*
 * update()
 *  return `true' if r > r_lub_max 
 *   ---> deactivate
 */
bool
Interaction::update(const bool switch_off_allowed){
	bool switched_off = false;
	if (active){
		// update tangential displacement: we do it before updating nr_vec
		if (sys->friction && contact) {
			incrementContactTangentialDisplacement();
		}
		// compute new r_vec and distance
		calcDistanceNormalVector();
		
		// check new state of the interaction
		switched_off = updateState(switch_off_allowed);

		// compute new contact forces if needed
		calcContactInteraction();

		if ( near == false ){
			if (_gap < sys->dist_near*0.5*ro){
				near = true;
				strain_near_contact = sys->shear_strain;
			}
		} else {
			if (_gap > sys->dist_near*0.5*ro){
				/*
				 * Separate
				 */
				double age_approaching = sys->shear_strain - strain_near_contact;
				sys->nearing_time_record.push_back(age_approaching);

				near = false;
			}
		}
	}
	return switched_off;
}

/*
 * Just provisional function for future analysis
 *
 */
void
Interaction::outputTrajectory(){
	for (unsigned int k=0; k < trajectory.size(); k++){
		sys->fout_trajectory << trajectory[k].x << ' ' <<  trajectory[k].y << ' '<<  trajectory[k].z << ' ';
		sys->fout_trajectory << gap_history[k] << endl;
	}
	sys->fout_trajectory << endl ;
	trajectory.clear();
	gap_history.clear();
}

double
Interaction::nearingTime(){
	double nearing_time = -1;
	if (near == true){
		nearing_time = sys->shear_strain - strain_near_contact;
	}
	return nearing_time;
}

void
Interaction::recordTrajectory(){
	trajectory.push_back(r_vec);
	gap_history.push_back(_gap);
}
