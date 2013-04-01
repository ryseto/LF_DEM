//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 02/21/13.
//  Copyright (c) 2013 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"

void System::calcStressesHydroContactBrownian(){
	int zero_2Dsimu;
	if (dimension == 2){
		zero_2Dsimu = 0;
	}else{
		zero_2Dsimu = 1;
	}
	/**************************************************
              1. Stress from background flow 
	                                                **/
    for (int i = 0; i < _np; i++){
		double a = radius[i];
		bgfstress[i].elm[2] = (5.0/9)*bgf_factor*a*a*a;
	}

	/**************************************************
       2. and 3.: Stresses from 
			       2-body lubrication and contacts  **/

	// first obtain hydrodynamic part of velocity
    stokes_solver.resetRHS();
    stokes_solver.prepareNewBuild_RFU("direct");
    addStokesDrag();
    buildLubricationTerms();
    stokes_solver.complete_RFU();
    stokes_solver.solve(v_hydro);

	// then obtain contact forces, adn contact part of velocity
    stokes_solver.resetRHS();
	setContactForceToParticle();
    buildContactTerms();
    stokes_solver.solve(v_cont);
    stokes_solver.solvingIsDone();

	// from that, compute stresses
	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			interaction[k].addHydroStress(); // - R_SU * v_hydro
			interaction[k].addContactStress(); // - R_SU * v_cont - rF_cont
		}
	}

	/**************************************************
       4. : Stresses from Brownian forces 
			                                        **/
	
	// This needs to be done with a mid-point algorithm,
	// due to the drift term of Brownian motion.
	// The steps are similar to the ones done for the  
	// actual dynamics, although the motions are reverted
	// at the very end, to let the system back in the initial
	// state.
	stresslet stresslet_i_init;
    stresslet stresslet_j_init;
    stresslet stresslet_i_mid;
    stresslet stresslet_j_mid;
    stresslet *step_stresslet = new stresslet [_np];
    for (int i=0; i < _np; i++){
	  for (int u=0; u < 5; u++){
	    step_stresslet[i].elm[u] = 0.;
	  }
    }
    vec3d vi;
    vec3d vj;
	/*********************************************************/
	/*                    First Step                         */
	/*********************************************************/
    stokes_solver.resetRHS();

    stokes_solver.prepareNewBuild_RFU("direct");

    addStokesDrag();
    buildLubricationTerms(true);

    stokes_solver.complete_RFU();
    buildContactTerms();
    stokes_solver.solve(v_lub_cont);

    // now the Brownian part of the velocity:
    // predictor-corrector algortithm (see Melrose & Ball, 1997)
    // 
    // we do not call solvingIsDone() before new solve(), because 
    // R_FU has not changed, so same factorization is safely used

	stokes_solver.setRHS( fb->generate_invLFb() );
	stokes_solver.solve_CholTrans( v_Brownian_init );
	stokes_solver.solvingIsDone();

    for (int i=0; i < _np; i++){
		v_Brownian_init[3*i+1] *= zero_2Dsimu;
    }

	/****** Brownian Stress: term R_SU * v_Brownian_init ****/

    for (int k = 0; k < num_interaction; k++) {
		if (interaction[k].active) {
			interaction[k].pairVelocityStresslet(v_Brownian_init, stresslet_i_init, stresslet_j_init);
			int i = interaction[k].par_num[0];
			int j = interaction[k].par_num[1];
			for (int u=0; u < 5; u++) {
				step_stresslet[i].elm[u] -= 0.5*stresslet_i_init.elm[u];
				step_stresslet[j].elm[u] -= 0.5*stresslet_j_init.elm[u];
			}
		}
    }
    // move particles to intermediate point
    for (int i=0; i < _np; i++) {
		int i3 = 3*i;
		vec3d dr(v_Brownian_init[i3]*dt, v_Brownian_init[i3+1]*dt, v_Brownian_init[i3+2]*dt);
		displacement(i, dr);
    }
    updateInteractions();
	/*********************************************************/
	/*                   Second Step                         */
	/*********************************************************/

    // build new Resistance matrix after move
    stokes_solver.prepareNewBuild_RFU("direct");
    addStokesDrag();
    buildLubricationTerms(false); // false: don't modify rhs
    stokes_solver.complete_RFU();

    // get the intermediate brownian velocity
	stokes_solver.solve_CholTrans( v_Brownian_mid );

    stokes_solver.solvingIsDone();

	/**** Brownian Stress: term  -R_SU_mid * v_Brownian_mid **/
    for (int i=0; i < _np; i++){
		v_Brownian_mid[3*i+1] *= zero_2Dsimu;
	}

    for (int k = 0; k < num_interaction; k++){
		if(interaction[k].active){
			interaction[k].pairVelocityStresslet(v_Brownian_mid, stresslet_i_mid, stresslet_j_mid);
			
			int i = interaction[k].par_num[0];
			int j = interaction[k].par_num[1];
			for (int u=0; u < 5; u++){
				step_stresslet[i].elm[u] += 0.5*stresslet_i_mid.elm[u];
				step_stresslet[j].elm[u] += 0.5*stresslet_j_mid.elm[u];
			}
		}
    }
	/**********************************************************/
	/*  Finishing stress computation                          */
	/*  and leaving particle positions as they initially were */
	/**********************************************************/
    for (int i=0; i < _np; i++){
	  for (int u=0; u < 5; u++){
		brownianstress[i].elm[u] += step_stresslet[i].elm[u];
	  }
    }
	// move particles back to initial point, and update interactions
    for (int i=0; i < _np; i++){
		int i3 = 3*i;
		vec3d dr(-v_Brownian_init[i3]*dt, -v_Brownian_init[i3+1]*dt, -v_Brownian_init[i3+2]*dt);
		displacement(i, dr);
    }
    updateInteractions();
	delete [] step_stresslet;
}

void
System::calcStressesHydroContact(){
	/**************************************************
	 1. Stress from background flow
	 **/
    for (int i=0; i<_np; i++){
		bgfstress[i].elm[2] = (5./9)*bgf_factor*radius_cubic[i];
	}
	/**************************************************
	 2. and 3.: Stresses from
	 2-body lubrication and contacts  **/
	// first obtain hydrodynamic part of velocity
    stokes_solver.resetRHS();
    stokes_solver.prepareNewBuild_RFU("direct");
    addStokesDrag();
    buildLubricationTerms();
    stokes_solver.complete_RFU();
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
    stokes_solver.solvingIsDone();
	// from that, compute stresses
	for (int k=0; k<num_interaction; k++){
		if (interaction[k].active){
			interaction[k].addHydroStress(); // - R_SU * v_hydro
			interaction[k].addContactStress(); //  - R_SU * v_cont - rF_cont
			interaction[k].addColloidalStress(); //  - R_SU * v_colloid - rF_colloid
		}
	}
	/*
	 * Calculate lubrication force to output
	 */
	for (int k=0; k<num_interaction; k++){
		if (interaction[k].active){
			interaction[k].evaluateLubricationForce();
		}
	}
	
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
	stressReset();
	if(brownian){
	 	calcStressesHydroContactBrownian();
	}else{
		calcStressesHydroContact();
	}
	for (int u=0; u < 5; u++){
		total_hydro_stress[u] = 0;
		total_contact_stressXF[u] = 0;
		total_contact_stressGU[u] = 0;
		total_colloidal_stressXF[u] = 0;
		total_colloidal_stressGU[u] = 0;
		total_brownian_stress[u] = 0;
	}
	for (int i=0; i < _np; i++){
		for (int u=0; u < 5; u++){
			total_hydro_stress[u] += lubstress[i].elm[u]+bgfstress[i].elm[u];
			total_contact_stressXF[u] += contactstressXF[i].elm[u];
			total_contact_stressGU[u] += contactstressGU[i].elm[u];
			total_colloidal_stressXF[u] += colloidalstressXF[i].elm[u];
			total_colloidal_stressGU[u] += colloidalstressGU[i].elm[u];
			total_brownian_stress[u] += brownianstress[i].elm[u];
		}
	}
	for (int u=0; u < 5; u++){
		total_hydro_stress[u] /= valSystemVolume();
		total_contact_stressXF[u] /= valSystemVolume();
		total_contact_stressGU[u] /= valSystemVolume();
		total_colloidal_stressXF[u] /= valSystemVolume();
		total_colloidal_stressGU[u] /= valSystemVolume();
		total_brownian_stress[u] /= valSystemVolume();
	}
	stressBrownianReset();
}
