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
			interaction[k].lubrication.addHydroStress(); // - R_SU * v
			interaction[k].contact.addContactStress(); //  - rF_cont
			if (repulsiveforce) {
				interaction[k].addRepulsiveStress(); //  - rF_rep
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

// void
// System::calcStress(){


// 	//////////////////////////////////////////////////////////////
// 	total_hydro_stress.reset();
// 	for (int i=0; i<np; i++) {
// 		total_hydro_stress += lubstress[i];
// 	}
// 	total_hydro_stress /= System_volume();
// 	//////////////////////////////////////////////////////////////
// 	total_contact_stressGU.reset();
// 	for (int i=0; i<np; i++) {
// 		total_contact_stressGU += contactstressGU[i];
// 	}
// 	total_contact_stressGU /= System_volume();
// 	//////////////////////////////////////////////////////////////
// 	total_contact_stressXF_normal.reset();
// 	total_contact_stressXF_tan.reset();
// 	for (int k=0; k<nb_interaction; k++) {
// 		if (interaction[k].is_contact()) {
// 			total_contact_stressXF_normal += interaction[k].contact.getContactStressXF_normal();
// 			total_contact_stressXF_tan += interaction[k].contact.getContactStressXF_tan();
// 		}
// 	}
// 	total_contact_stressXF_normal /= System_volume();
// 	total_contact_stressXF_tan /= System_volume();
// 	//////////////////////////////////////////////////////////////
// 	if (repulsiveforce) {
// 		total_repulsive_stressGU.reset();
// 		for (int i=0; i<np; i++) {
// 			total_repulsive_stressGU += repulsivestressGU[i];
// 		}
// 		total_repulsive_stressGU /= System_volume();
// 		//////////////////////////////////////////////////////////////
// 		total_repulsive_stressXF.reset();
// 		for (int k=0; k<nb_interaction; k++) {
// 			if (interaction[k].is_active()) {
// 				total_repulsive_stressXF += interaction[k].getRepulsiveStressXF();
// 			}
// 		}
// 		total_repulsive_stressXF /= System_volume();
// 	}
// 	//////////////////////////////////////////////////////////////
// 	if (brownian) {
// 		total_brownian_stressGU.reset();
// 		for (int i=0; i<np; i++) {
// 			total_brownian_stressGU += brownianstressGU[i];
// 		}
// 		total_brownian_stressGU /= System_volume();
// 	}
// }


void
System::calcStress(){
	//////////////////////////////////////////////////////////////
	total_hydro_stress.reset();
	for (int i=0; i<np; i++) {
		total_hydro_stress += avg_lubstress[i];
	}
	total_hydro_stress /= (avg_stress_nb*System_volume());
	//////////////////////////////////////////////////////////////
	total_contact_stressGU.reset();
	for (int i=0; i<np; i++) {
		total_contact_stressGU += avg_contactstressGU[i];
	}
	total_contact_stressGU /= (avg_stress_nb*System_volume());
	//////////////////////////////////////////////////////////////
	total_contact_stressXF_normal = avg_contactstressXF_normal/(avg_stress_nb*System_volume());
	total_contact_stressXF_tan = avg_contactstressXF_tan/(avg_stress_nb*System_volume());
	//////////////////////////////////////////////////////////////
	if (repulsiveforce) {
		total_repulsive_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_repulsive_stressGU += avg_repulsivestressGU[i];
		}
		total_repulsive_stressGU /= (avg_stress_nb*System_volume());
		//////////////////////////////////////////////////////////
		total_repulsive_stressXF = avg_repulsivestressXF/(avg_stress_nb*System_volume());
	}
	//////////////////////////////////////////////////////////////
	if (brownian) {
		total_brownian_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_brownian_stressGU += avg_brownianstressGU[i];
		}
		total_brownian_stressGU /= (avg_stress_nb*System_volume());
	}
}
