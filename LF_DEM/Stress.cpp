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
	

	stresslet stresslet_i_init;
    stresslet stresslet_j_init;
    stresslet stresslet_i_mid;
    stresslet stresslet_j_mid;
    stresslet *step_stresslet = new stresslet [np];
    stresslet total_step_stresslet;
    for (int i=0; i < np; i++){
	  for (int u=0; u < 5; u++){
	    step_stresslet[i].elm[u] = 0.;
	  }
    }

	for (int u=0; u < 5; u++){
	  total_step_stresslet.elm[u] = 0.;
	}


    vec3d vi;
    vec3d vj;


	/*********************************************************/
	/*                    First Step                         */
	/*********************************************************/
    stokes_solver->resetRHS();

    stokes_solver->prepareNewBuild_RFU("direct");

    addStokesDrag();
    buildLubricationTerms(true);

    stokes_solver->complete_RFU();
    buildContactTerms();
    stokes_solver->solve(v_lub_cont);

    // now the Brownian part of the velocity:
    // predictor-corrector algortithm (see Melrose & Ball, 1997)
    // 
    // we do not call solvingIsDone() before new solve(), because 
    // R_FU has not changed, so same factorization is safely used

	stokes_solver->setRHS( fb->generate_invLFb() );
	stokes_solver->solve_CholTrans( v_Brownian_init );
	stokes_solver->solvingIsDone();


	/****** Brownian Stress: term R_SU * v_Brownian_init ****/
    for (int k = 0; k < num_interaction; k++){
	  for (int u=0; u < 5; u ++){
 		stresslet_i_init.elm[u] = 0.;
		stresslet_j_init.elm[u] = 0.;
	  }

	  int i = interaction[k].particle_num[0];
	  int j = interaction[k].particle_num[1];
	  int i3 = 3*i;
	  int j3 = 3*j;

	  vi.x = v_Brownian_init[i3  ];
	  vi.y = v_Brownian_init[i3+1];
	  vi.z = v_Brownian_init[i3+2];

	  vj.x = v_Brownian_init[j3  ];
	  vj.y = v_Brownian_init[j3+1];
	  vj.z = v_Brownian_init[j3+2];

	  interaction[k].pairVelocityStresslet(vi, vj, stresslet_i_init, stresslet_j_init);

	  for (int u=0; u < 5; u++){
	    step_stresslet[i].elm[u] -= stresslet_i_init.elm[u];
	    step_stresslet[j].elm[u] -= stresslet_j_init.elm[u];
	  }
    }

	// for (int u=0; u < 5; u++){
	//   for (int i=0; i < np; i++){
	// 	total_step_stresslet.elm[u] += step_stresslet[i].elm[u];
	//   }
	//   //cout << total_step_stresslet.elm[u]/np << " " ;
	// }

    // move particles to intermediate point
    for (int i=0; i < np; i++){
		int i3 = 3*i;
		vec3d dr(v_Brownian_init[i3]*dt,
				 v_Brownian_init[i3+1]*dt*zero_2Dsimu,
				 v_Brownian_init[i3+2]*dt);
		displacement(i, dr);
    }
    updateInteractions(false);
	


	/*********************************************************/
	/*                   Second Step                         */
	/*********************************************************/

    // build new Resistance matrix after move
    stokes_solver->prepareNewBuild_RFU("direct");
    addStokesDrag();
    buildLubricationTerms(false); // false: don't modify rhs
    stokes_solver->complete_RFU();

    // get the intermediate brownian velocity
	stokes_solver->solve_CholTrans( v_Brownian_mid );

    stokes_solver->solvingIsDone();

	/**** Brownian Stress: term  -R_SU_mid * v_Brownian_mid **/
    for (int k = 0; k < num_interaction; k++){
	  for (int u=0; u < 5; u ++){
		stresslet_i_mid.elm[u] = 0.;
		stresslet_j_mid.elm[u] = 0.;
	  }

	  int i = interaction[k].particle_num[0];
	  int j = interaction[k].particle_num[1];
	  int i3 = 3*i;
	  int j3 = 3*j;

	  vi.x = v_Brownian_mid[i3  ];
	  vi.y = v_Brownian_mid[i3+1];
	  vi.z = v_Brownian_mid[i3+2];

	  vj.x = v_Brownian_mid[j3  ];
	  vj.y = v_Brownian_mid[j3+1];
	  vj.z = v_Brownian_mid[j3+2];

	  interaction[k].pairVelocityStresslet(vi, vj, stresslet_i_mid, stresslet_j_mid);

	  for (int u=0; u < 5; u++){
	    step_stresslet[i].elm[u] += stresslet_i_mid.elm[u];
	    step_stresslet[j].elm[u] += stresslet_j_mid.elm[u];
		//		total_step_stresslet.elm[u] -= stresslet_i_mid.elm[u] + stresslet_j_mid.elm[u];
	  }
    }



	/**********************************************************/
	/*  Finishing stress computation                          */
	/*  and leaving particle positions as they initially were */
	/**********************************************************/

	/** Overall Brownian Stress: multiply all by  0.5 **/
    for (int i=0; i < np; i++){
	  for (int u=0; u < 5; u++){
	    step_stresslet[i].elm[u] *= 0.5;
		brownianstress[i].elm[u] += step_stresslet[i].elm[u];
	  }
    }
	//cout << " total ";
	for (int u=0; u < 5; u++){
	  total_step_stresslet.elm[u] *= 0.5;
	  //cout << total_step_stresslet.elm[u]/np << " " ;	  
	}
	//cout << endl;


	// move particles back to initial point, and update interactions

    for (int i=0; i < np; i++){
		int i3 = 3*i;
		vec3d dr(-v_Brownian_init[i3]*dt,
				 -v_Brownian_init[i3+1]*dt*zero_2Dsimu,
				 -v_Brownian_init[i3+2]*dt);
		displacement(i, dr);
    }
    updateInteractions(false);

}


void
System::calcStressesHydroContact(){

	setContactForceToParticle();

    stokes_solver->resetRHS();

    stokes_solver->prepareNewBuild_RFU("direct");

    addStokesDrag();
    buildLubricationTerms();

    stokes_solver->complete_RFU();

    stokes_solver->solve(v_hydro);


    stokes_solver->resetRHS();
    buildContactTerms();
	
    stokes_solver->solve(v_cont);
    stokes_solver->solvingIsDone();

	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			interaction[k].addHydroStress();
			interaction[k].addContactStress();
		}
	}

	
	// testing : compare with stress computation from forces
    for (int i = 0; i < np; i++){
		int i3 = 3*i;
		relative_velocity_lub_cont[i].x = v_hydro[i3  ] + v_cont[i3  ];
		relative_velocity_lub_cont[i].y = v_hydro[i3+1] + v_cont[i3+1];
		relative_velocity_lub_cont[i].z = v_hydro[i3+2] + v_cont[i3+2];
		
		velocity[i].x = relative_velocity_lub_cont[i].x + shear_rate*position[i].z;
		velocity[i].y = relative_velocity_lub_cont[i].y;
		velocity[i].z = relative_velocity_lub_cont[i].z;
		
    }

	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			interaction[k].evaluateLubricationForce();
			interaction[k].addLubricationStress();
			interaction[k].addContactStress2();
		}
	}


}
