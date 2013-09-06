//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/21/13.
//  Copyright (c) 2013 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"

//void
//System::calcStressesHydroContactBrownian(){
//	int zero_2Dsimu;
//	if (dimension == 2){
//		zero_2Dsimu = 0;
//	}else{
//		zero_2Dsimu = 1;
//	}
//	/**************************************************
//	 2. and 3.: Stresses from
//	 2-body lubrication and contacts  **/
//	
//	// first obtain hydrodynamic part of velocity
//    stokes_solver.resetRHS();
//    stokes_solver.resetResistanceMatrix("direct");
//    addStokesDrag();
//    buildLubricationTerms();
//    stokes_solver.completeResistanceMatrix();
//    stokes_solver.solve(v_hydro);
//	
//	// then obtain contact forces, adn contact part of velocity
//    stokes_solver.resetRHS();
//	setContactForceToParticle();
//    buildContactTerms();
//    stokes_solver.solve(v_cont);
//    stokes_solver.solvingIsDone();
//	
//	// from that, compute stresses
//	for (int k = 0; k < nb_interaction; k++) {
//		if (interaction[k].is_active()) {
//			interaction[k].addHydroStress(); // - R_SU * v_hydro
//			interaction[k].addContactStress(); // - R_SU * v_cont - rF_cont
//		}
//	}
//	
//	/**************************************************
//	 4. : Stresses from Brownian forces
//	 **/
//	
//	// This needs to be done with a mid-point algorithm,
//	// due to the drift term of Brownian motion.
//	// The steps are similar to the ones done for the
//	// actual dynamics, although the motions are reverted
//	// at the very end, to let the system back in the initial
//	// state.
//	StressTensor stresslet_i_init;
//    StressTensor stresslet_j_init;
//    StressTensor stresslet_i_mid;
//    StressTensor stresslet_j_mid;
//    StressTensor *step_stresslet = new StressTensor [np];
//    for (int i=0; i < np; i++) {
//		step_stresslet[i].reset();
//    }
//    vec3d vi;
//    vec3d vj;
//	/*********************************************************/
//	/*                    First Step                         */
//	/*********************************************************/
//    stokes_solver.resetRHS();
//	
//    stokes_solver.resetResistanceMatrix("direct");
//	
//    addStokesDrag();
//    buildLubricationTerms(true);
//	
//    stokes_solver.completeResistanceMatrix();
//    buildContactTerms();
//    stokes_solver.solve(v_lub_cont);
//	
//    // now the Brownian part of the velocity:
//    // predictor-corrector algortithm (see Melrose & Ball, 1997)
//    //
//    // we do not call solvingIsDone() before new solve(), because
//    // R_FU has not changed, so same factorization is safely used
//	
//	stokes_solver.setRHS( fb->generate_invLFb() );
//	stokes_solver.solve_CholTrans( v_Brownian_init );
//	stokes_solver.solvingIsDone();
//	
//    for (int i=0; i < np; i++) {
//		v_Brownian_init[3*i+1] *= zero_2Dsimu;
//    }
//	
//	/****** Brownian Stress: term R_SU * v_Brownian_init ****/
//	
//    for (int k = 0; k < nb_interaction; k++) {
//		if (interaction[k].is_active()) {
//			interaction[k].pairVelocityStresslet(v_Brownian_init, stresslet_i_init, stresslet_j_init);
//			unsigned int i, j;
//			interaction[k].get_par_num(i, j);
//			step_stresslet[i] -= 0.5*stresslet_i_init;
//			step_stresslet[j] -= 0.5*stresslet_j_init;
//		}
//    }
//    // move particles to intermediate point
//    for (int i=0; i < np; i++) {
//		int i3 = 3*i;
//		vec3d dr(v_Brownian_init[i3]*dt, v_Brownian_init[i3+1]*dt, v_Brownian_init[i3+2]*dt);
//		displacement(i, dr);
//    }
//    updateInteractions();
//	/*********************************************************/
//	/*                   Second Step                         */
//	/*********************************************************/
//	
//    // build new Resistance matrix after move
//    stokes_solver.resetResistanceMatrix("direct");
//    addStokesDrag();
//    buildLubricationTerms(false); // false: don't modify rhs
//    stokes_solver.completeResistanceMatrix();
//	
//    // get the intermediate brownian velocity
//	stokes_solver.solve_CholTrans( v_Brownian_mid );
//	
//    stokes_solver.solvingIsDone();
//	
//	/**** Brownian Stress: term  -R_SU_mid * v_Brownian_mid **/
//    for (int i=0; i < np; i++) {
//		v_Brownian_mid[3*i+1] *= zero_2Dsimu;
//	}
//	
//    for (int k = 0; k < nb_interaction; k++) {
//		if(interaction[k].is_active()) {
//			interaction[k].pairVelocityStresslet(v_Brownian_mid, stresslet_i_mid, stresslet_j_mid);
//			unsigned int i, j;
//			interaction[k].get_par_num(i, j);
//			step_stresslet[i] += 0.5*stresslet_i_mid;
//			step_stresslet[j] += 0.5*stresslet_j_mid;
//		}
//    }
//	/**********************************************************/
//	/*  Finishing stress computation                          */
//	/*  and leaving particle positions as they initially were */
//	/**********************************************************/
//    for (int i=0; i < np; i++) {
//		brownianstress[i] += step_stresslet[i];
//    }
//	// move particles back to initial point, and update interactions
//    for (int i=0; i < np; i++) {
//		int i3 = 3*i;
//		vec3d dr(-v_Brownian_init[i3]*dt, -v_Brownian_init[i3+1]*dt, -v_Brownian_init[i3+2]*dt);
//		displacement(i, dr);
//    }
//    updateInteractions();
//	delete [] step_stresslet;
//}

void
System::calcStressesHydroContact(){
	/**************************************************
	 1. Stress from background flow
	 **/
 
	/**************************************************
	 2. and 3.: Stresses from
	 2-body lubrication and contacts  **/

	// first obtain hydrodynamic part of velocity
	stokes_solver.resetRHS();
	int nb_of_active_interactions = nb_interaction-deactivated_interaction.size();
    stokes_solver.resetResistanceMatrix("direct",nb_of_active_interactions);
	addStokesDrag();
	buildLubricationTerms();
	stokes_solver.completeResistanceMatrix();
    stokes_solver.solve(v_hydro);
	// then obtain contact forces, and contact part of velocity
	stokes_solver.resetRHS();
	setContactForceToParticle();
    buildContactTerms();
	stokes_solver.solve(v_cont);
	// then obtain colloidal forces
    stokes_solver.resetRHS();
	setColloidalForceToParticle();
    buildColloidalForceTerms();
    stokes_solver.solve(v_colloidal);
	/////////////////////////////////////////////////
	// from that, compute stresses
//	cout << nb_interaction << endl;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			if (lubrication_model == 1){
				interaction[k].calcXFunctionsStress();
			} else if (lubrication_model == 2){
				interaction[k].calcXYFunctionsStress();
			}
			interaction[k].addHydroStress(); // - R_SU * v_hydro
			interaction[k].addContactStress(); //  - R_SU * v_cont - rF_cont
			interaction[k].addColloidalStress(); //  - R_SU * v_colloid - rF_colloid
		}
	}
    stokes_solver.solvingIsDone();

		// >>>>  testing : compare with stress computation from forces
	// Note that the definition of Hydrodynamic Stress and Contact Stress
	// are different: the part coming from v_cont is included in Hydro stress.
	// There is also a factor 2 difference coming from the way the two methods
	// define the stress. This is known, and is not a bug.
	// The total stress should the same though (well, more precisely 2*total_stress=total_stress2).
	
	//	for (int k = 0; k < num_interaction; k++){
	//		if (interaction[k].active){
	//			interaction[k].evaluateLubricationForce();
	//			interaction[k].addLubricationStress();    // - R_SU * (v_hydro + v_cont)
	//			interaction[k].addContactStress2();       // - rF_cont
	//		}
	//	}
	// <<<< end of testing
}

void
System::calcStress(){
	static double previous_strain = 0;
	previous_strain = shear_strain;
	stressReset();
	if (brownian) {
	 	//calcStressesHydroContactBrownian();
	} else {
		calcStressesHydroContact();
	}
	total_hydro_stress.reset();
	total_contact_stressGU.reset();
	total_colloidal_stressGU.reset();
	total_contact_stressXF_normal.reset();
	total_contact_stressXF_tan.reset();
	total_colloidal_stressXF.reset();
	total_brownian_stress.reset();
	for (int i=0; i<np; i++) {
		total_hydro_stress += lubstress[i];
		total_contact_stressGU += contactstressGU[i];
		total_colloidal_stressGU += colloidalstressGU[i];
		if (brownian) {
			total_brownian_stress += brownianstress[i];
		}
	}
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			total_colloidal_stressXF += interaction[k].getColloidalStressXF();
		}
		if (interaction[k].is_contact()) {
			total_contact_stressXF_normal += interaction[k].getContactStressXF_normal();
			total_contact_stressXF_tan += interaction[k].getContactStressXF_tan();
		}
	}
	total_hydro_stress /= System_volume();
	total_contact_stressGU /= System_volume();
	total_contact_stressXF_normal /= System_volume();
	total_contact_stressXF_tan /= System_volume();
	total_colloidal_stressGU /= System_volume();
	total_colloidal_stressXF /= System_volume();
	total_brownian_stress /= System_volume();
	stressBrownianReset();
}
