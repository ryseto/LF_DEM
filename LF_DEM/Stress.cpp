//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/21/13.
//  Copyright (c) 2013 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"

void
System::calcStressPerParticle(){
	/////////////////////////////////////////////////
	// from the velocities V_H, V_C, V_Rep, V_B,
	// compute stresses R_SV * V
	// and then add rF stress for F_C and F_Rep
	stressReset();
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			if (lubrication_model == 1) {
				interaction[k].lubrication.calcXFunctionsStress();
			} else if (lubrication_model == 2) {
				interaction[k].lubrication.calcXYFunctionsStress();
			} else {
				cerr << "lubrication_model = 3 is not implemented" << endl;
				exit(1);
			}
			interaction[k].lubrication.addHydroStress(); // R_SE:Einf-R_SU*v 
			interaction[k].contact.calcContactStress(); // - rF_cont
			if (repulsiveforce) {
				interaction[k].calcRepulsiveStress(); // - rF_rep
			}
		}
	}
	if (brownian) {
		if (in_predictor) {
			for (int i=0; i<np; i++) {
				brownianstressGU_predictor[i] = brownianstressGU[i];
			}
		} else {
			for (int i=0; i<np; i++) {
				/*
				 * [ Banchio & Brady 2003 ] [ Ball & Melrose 1997 ]
				 */
				brownianstressGU[i] = 0.5*(brownianstressGU[i]-brownianstressGU_predictor[i]);
			}
		}
	}
}

void
System::calcStress(){
	//////////////////////////////////////////////////////////////
	total_hydro_stress.reset();
	for (int i=0; i<np; i++) {
		total_hydro_stress += lubstress[i];
	}
	if (dimensionless_shear_rate > 0) {
		total_hydro_stress /= system_volume;
	} else {
		total_hydro_stress /= -system_volume;
	}
	//////////////////////////////////////////////////////////////
	total_contact_stressGU.reset();
	for (int i=0; i<np; i++) {
		total_contact_stressGU += contactstressGU[i];
	}
	total_contact_stressGU /= system_volume;
	//////////////////////////////////////////////////////////////
	total_contact_stressXF_normal.reset();
	total_contact_stressXF_tan.reset();
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()) {
			total_contact_stressXF_normal += interaction[k].contact.getContactStressXF_normal();
			total_contact_stressXF_tan += interaction[k].contact.getContactStressXF_tan();
		}
	}
	total_contact_stressXF_normal /= system_volume;
	total_contact_stressXF_tan /= system_volume;
	if (stress_controlled && p.unscaled_contactmodel) {
		total_contact_stressXF_normal /= abs(dimensionless_shear_rate);
		total_contact_stressXF_tan /= abs(dimensionless_shear_rate);
	}
	//////////////////////////////////////////////////////////////
	if (repulsiveforce) {
		total_repulsive_stressXF.reset();
		for (int k=0; k<nb_interaction; k++) {
			total_repulsive_stressXF += interaction[k].getRepulsiveStressXF();
		}
		total_repulsive_stressXF /= system_volume;
		if (stress_controlled) {
			total_repulsive_stressXF /= abs(dimensionless_shear_rate);
		}
		//////////////////////////////////////////////////////////
		total_repulsive_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_repulsive_stressGU += repulsivestressGU[i];
		}
		total_repulsive_stressGU /= system_volume;
	}
	//////////////////////////////////////////////////////////////
	if (brownian) {
		total_brownian_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_brownian_stressGU += brownianstressGU[i];
		}
		total_brownian_stressGU /= system_volume;
	}
}
