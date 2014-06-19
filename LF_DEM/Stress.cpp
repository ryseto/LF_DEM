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
	// from the velocities V_H, V_C, V_Coll, V_B, 
	// compute stresses R_SV * V
	// and then add rF stress for F_C and F_Coll
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			if (lubrication_model == 1) {
				interaction[k].lubrication.calcXFunctionsStress();
			} else if (lubrication_model == 2) {
				interaction[k].lubrication.calcXYFunctionsStress();
			} else if (lubrication_model == 3) {
				cerr << "lubrication_model = 3 is not implemented" << endl;
				exit(1);
				//	if (interaction[k].is_contact()){
				//		interaction[k].lubrication.calcXYFunctionsStress();
				//	} else {
				//		interaction[k].lubrication.calcXFunctionsStress();
				//	}
			}
			interaction[k].lubrication.addHydroStress(); // - R_SU * v
			interaction[k].contact.addContactStress(); //  - rF_cont
			interaction[k].addColloidalStress(); //  - rF_colloid
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
// 	if (colloidalforce) {
// 		total_colloidal_stressGU.reset();
// 		for (int i=0; i<np; i++) {
// 			total_colloidal_stressGU += colloidalstressGU[i];
// 		}
// 		total_colloidal_stressGU /= System_volume();
// 		//////////////////////////////////////////////////////////////
// 		total_colloidal_stressXF.reset();
// 		for (int k=0; k<nb_interaction; k++) {
// 			if (interaction[k].is_active()) {
// 				total_colloidal_stressXF += interaction[k].getColloidalStressXF();
// 			}
// 		}
// 		total_colloidal_stressXF /= System_volume();
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
		total_hydro_stress += avg_lubstress[i]/avg_stress_nb;
	}
	total_hydro_stress /= System_volume();
	//////////////////////////////////////////////////////////////
	total_contact_stressGU.reset();
	for (int i=0; i<np; i++) {
		total_contact_stressGU += avg_contactstressGU[i]/avg_stress_nb;
	}
	total_contact_stressGU /= System_volume();
	//////////////////////////////////////////////////////////////
	total_contact_stressXF_normal = avg_contactstressXF_normal/avg_stress_nb;
	total_contact_stressXF_tan = avg_contactstressXF_tan/avg_stress_nb;
	total_contact_stressXF_normal /= System_volume();
	total_contact_stressXF_tan /= System_volume();
	//////////////////////////////////////////////////////////////
	if (colloidalforce) {
		total_colloidal_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_colloidal_stressGU += avg_colloidalstressGU[i]/avg_stress_nb;
		}
		total_colloidal_stressGU /= System_volume();
		//////////////////////////////////////////////////////////////
		total_colloidal_stressXF = avg_colloidalstressXF/avg_stress_nb;
		total_colloidal_stressXF /= System_volume();
	}
	//////////////////////////////////////////////////////////////
	if (brownian) {
		total_brownian_stressGU.reset();
		for (int i=0; i<np; i++) {
			total_brownian_stressGU += avg_brownianstressGU[i]/avg_stress_nb;
		}
		total_brownian_stressGU /= System_volume();
	}
}
