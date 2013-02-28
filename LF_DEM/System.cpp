//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
#include <sstream>

System::~System(){

	if (!position)
		delete [] position;
	if (!radius)
		delete [] radius;
	if (!angle)
		delete [] angle;

	if (!velocity)
		delete [] velocity;
	if (!relative_velocity)
		delete [] relative_velocity;
	if (!relative_velocity_lub_cont)
		delete [] relative_velocity_lub_cont;
	if (!relative_velocity_brownian)
		delete [] relative_velocity_brownian;

	if (!ang_velocity)
		delete [] ang_velocity;

//	if (!total_force)
//		delete [] total_force;
	if (!lubrication_force)
		delete [] lubrication_force;
	if (!contact_force)
		delete [] contact_force;

	if (!brownian_force)
	  delete [] brownian_force;

//	if (!torque)
//		delete [] torque;

	if (!interaction)
		delete [] interaction;

	for(int i=0; i<np; i++){
	  interaction_list[i].clear();
	}
	if (!interaction_list)
		delete [] interaction_list;
	if (!interaction_partners)
		delete [] interaction_partners;
	if (!fb){
		delete [] fb;
	}
	
	if (!v_lub_cont)
		delete [] v_lub_cont;

	if(brownian){
	    if (!v_Brownian_init)
		delete [] v_Brownian_init;

	    if (!v_Brownian_mid)
		delete [] v_Brownian_mid;

		if (!v_lub_cont_mid)
			delete [] v_lub_cont_mid;

		if (!lub_cont_forces_init)
			delete [] lub_cont_forces_init;
	}

	delete stokes_solver;
};


void
System::allocateRessources(){
	position = new vec3d [np];
	radius = new double [np];
	angle = new double [np];
	velocity = new vec3d [np];
	velocity_predictor.resize(np);
	relative_velocity = new vec3d [np];
	relative_velocity_lub_cont = new vec3d [np];
	relative_velocity_brownian = new vec3d [np];
	ang_velocity = new vec3d [np];
	ang_velocity_predictor.resize(np);
//	total_force = new vec3d [np];
	lubrication_force = new vec3d [np];
	contact_force = new vec3d [np];
	contact_torque = new vec3d [np];
	brownian_force = new vec3d [np];
//	torque = new vec3d [np];
	lub_force = new vec3d [np];

	lubstress.resize(np);
	contactstress.resize(np);
	brownianstress.resize(np);


	int maxnum_interactionpair_per_particle = 15;
	maxnum_interactionpair = (int)(maxnum_interactionpair_per_particle*np);
	interaction = new Interaction [maxnum_interactionpair];
	interaction_list = new set <Interaction*> [np];
	interaction_partners = new set <int> [np];
	
	lubstress2.resize(np);
	
	dof = 3;
	linalg_size = dof*np;
	
	v_lub_cont = new double [linalg_size];

	if(brownian){
	    v_Brownian_init = new double [linalg_size];
	    v_Brownian_mid = new double [linalg_size];
		v_lub_cont_mid = new double [linalg_size];
		lub_cont_forces_init = new double [linalg_size];
	}


	stokes_solver = new StokesSolver(np, brownian);
	fb = new BrownianForce(this); 
}

void
System::setupSystem(const vector<vec3d> &initial_positions,
					const vector <double> &radii){

	if (kb_T == 0){
		brownian = false;
	} else {
		brownian = true;
		integration_method = 2; // > force Euler
	}

	allocateRessources();
	radius_cubic.resize(np);
	for (int i=0; i < np; i++){
		position[i] = initial_positions[i];
		radius[i] = radii[i];
		radius_cubic[i] = radius[i]*radius[i]*radius[i];
		angle[i] = 0;
	}
	
	for (int k=0; k < maxnum_interactionpair ; k++){
		interaction[k].init(this);
	}
	for (int i=0; i < np; i++){
		velocity[i].x=0.;
		velocity[i].y=0.;
		velocity[i].z=0.;
		relative_velocity[i].x=0.;
		relative_velocity[i].y=0.;
		relative_velocity[i].z=0.;
		relative_velocity_lub_cont[i].x=0.;
		relative_velocity_lub_cont[i].y=0.;
		relative_velocity_lub_cont[i].z=0.;
		relative_velocity_brownian[i].x=0.;
		relative_velocity_brownian[i].y=0.;
		relative_velocity_brownian[i].z=0.;
	}
	initializeBoxing();
	checkNewInteraction();

	dt = dt * 1.0/radius_max;
	shear_strain = 0;
	num_interaction = 0;
	
	sq_critical_velocity = \
	dynamic_friction_critical_velocity*dynamic_friction_critical_velocity;
	sq_lub_max = lub_max*lub_max; // square of lubrication cutoff length.
	ts = 0;
	shear_disp = 0;
	vel_difference = _lz;
	/*
	 * dt_mid: the intermediate time step for the mid-point
	 * algortithm. dt/dt_mid = dt_ratio
	 * Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 * ASD code from Brady has dt_ratio=150
	 */
	dt_mid = dt/dt_ratio;// UNUSED NOW: ALGO EQUIVALENT TO dt_ratio = 1 (Melrose & Ball 1997)

	stokes_solver->initialize();
	
	// initialize the brownian force after the solver, as it assumes
	// the cholmod_common of the solver is already initialized
	if(brownian){
		fb->init();
	}

	brownianstress_calc_nb = 0;
	fix_interaction_status = false;
}

void
System::initializeBoxing(){// need to know radii first
	
	double max_radius=0.;
	for (int i=0; i < np; i++){
		if(radius[i]>max_radius){
			max_radius=radius[i];
		}
	}

	boxset = new BoxSet(lub_max*max_radius, this);
	for (int i=0; i < np; i++){
		boxset->box(i);
	}
	boxset->update();
}

void
System::timeEvolutionEulersMethod(){
	setContactForceToParticle();
	updateVelocityLubrication();
	deltaTimeEvolution();
}

void
System::deltaTimeEvolution(){
	// evolve PBC
	shear_disp += vel_difference*dt;
	if (shear_disp > lx()){
		shear_disp -= lx();
	}
	// move particles
	for (int i=0; i < np; i++){
		displacement(i, velocity[i]*dt);
	}
	if (dimension == 2){
		for (int i=0; i < np; i++){
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	// update boxing system
	boxset->update();
	checkNewInteraction();
	//
	fix_interaction_status = false;
	updateInteractions();
}

void
System::timeEvolutionPredictorCorrectorMethod(){
	/* Ball&Melrose 1997
	 * x'(t+dt) = x(t) + V^{-}
	 */
	setContactForceToParticle();
	updateVelocityLubrication();
	deltaTimeEvolutionPredictor();
	/*
	 * Keep V^{-} to use them in the corrector.
	 */
	for (int i=0; i < np; i++){
		velocity_predictor[i] = velocity[i];
		ang_velocity_predictor[i] = ang_velocity[i];
	}
	/*
	 * x(t + dt) = x(t) + 0.5*(V^{+}+V^{-})
	 *           = x'(t+dt) + 0.5*(V^{+}-V^{-})
	 */
	setContactForceToParticle();
	updateVelocityLubrication();
	deltaTimeEvolutionCorrector();
}

void
System::deltaTimeEvolutionPredictor(){
	shear_disp += vel_difference*dt;
	if (shear_disp > lx()){
		shear_disp -= lx();
	}
	// evolve PBC
	// move particles
	for (int i=0; i < np; i++){
		displacement(i, velocity[i]*dt);
	}
	if (dimension == 2){
		for (int i=0; i < np; i++){
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	fix_interaction_status = true;
	// update boxing system
	//boxset->update();
	//checkNewInteraction();
	updateInteractions();
}

void
System::deltaTimeEvolutionCorrector(){
	/*
	 *
	 * x'(t+dt) = x(t) + V^{-}dt
	 * x(t)     = x'(t+dt) - V^{-}dt
	 * x(t + dt) = x(t)     + 0.5*(V^{+}+V^{-})*dt
	 *           = x'(t+dt) + 0.5*(V^{+}-V^{-})*dt
	 */
	for (int i=0; i < np; i++){
		velocity[i] = 0.5*(velocity[i] - velocity_predictor[i]);
		ang_velocity[i] = 0.5*(ang_velocity[i] - ang_velocity_predictor[i]);
	}

	for (int i=0; i < np; i++){
		displacement(i, velocity[i]*dt);
	}
	if (dimension == 2){
		for (int i=0; i < np; i++){
			angle[i] += ang_velocity[i].y*dt;
		}
	}

	// update boxing system
	fix_interaction_status = false;
	boxset->update();
	checkNewInteraction();
	updateInteractions();
	
	revertRealVelocity();
}

void
System::revertRealVelocity(){
	/* In deltaTimeEvolutionCorrector,
	 * velocity[] and ang_velocity[]
	 * are virtual velocities to correct the predictor.
	 * This function reverts them
	 */
	for (int i=0; i < np; i++){
		velocity[i] += velocity_predictor[i];
		ang_velocity[i] += ang_velocity_predictor[i];
	}
}


void System::timeEvolutionBrownian(){
	int zero_2Dsimu;
	if (dimension == 2){
		zero_2Dsimu = 0;
	}else{
		zero_2Dsimu = 1;
	}
	/***************************************************************************
	 *  This routine implements a predictor-corrector algorithm                 *
	 *  for the dynamics with Brownian motion.                                  *
	 *  The		algorithm is the one of Melrose & Ball 1997.                        *
	 *                                                                          *
	 *  X_pred = X(t) + V(t)*dt                                                 *
	 *  X(t+dt) = X_pred + 0.5*(V(t+dt)-V(t))*dt                                *
	 *                                                                          *
	 *  Parameters                                                              *
	 *  They essentially controls what enter V(t) and V(t+dt):                  *
	 *   * displubcont : - if true  V(t) = V^{B}(t) + V^{H}(t) + V^{C}(t)       *
	 *                   - if false V(t) = V^{B}(t)                             *
	 *   * flubcont_update : - requires displubcont                             *
	 *                       - allows dt^2 scheme for contact and hydro         *
	 *                         velocities                                       *
	 *                       - if true :                                        *
	 *                          R_FU(t+dt) V^{*}(t+dt) = F^{*}(t+dt)            *
	 *                          with * = C or H                                 *
	 *                       - if false :                                       *
	 *                          R_FU(t+dt) V^{*}(t+dt) = F^{*}(t)               *
	 *                                                                          *
	 ****************************************************************************/
	bool displubcont = true;
	bool flubcont_update = true;
	/*********************************************************/
	/*                    Predictor                          */
	/*********************************************************/
	setContactForceToParticle();
	stokes_solver->resetRHS();
    stokes_solver->prepareNewBuild_RFU("direct");
	
    addStokesDrag();
    buildLubricationTerms(true);
	
    stokes_solver->complete_RFU();
    buildContactTerms();
    stokes_solver->solve(v_lub_cont);
	
	if(displubcont){
		stokes_solver->getRHS(lub_cont_forces_init);
	}
    // now the Brownian part of the velocity:
    // predictor-corrector algortithm (see Melrose & Ball, 1997)
    //
    // we do not call solvingIsDone() before new solve(), because
    // R_FU has not changed, so same factorization is safely used
	
	stokes_solver->setRHS( fb->generate_invLFb() );
	stokes_solver->solve_CholTrans( v_Brownian_init );
	stokes_solver->solvingIsDone();
    // move particles to intermediate point

    for (int i=0; i < np; i++){
		int i3 = 3*i;
		velocity[i].set(v_Brownian_init[i3],
						v_Brownian_init[i3+1]*zero_2Dsimu,
						v_Brownian_init[i3+2]);
		velocity[i].add(v_lub_cont[i3],
						v_lub_cont[i3+1],
						v_lub_cont[i3+2]);
		velocity[i].x += position[i].z;
		displacement(i, velocity[i]*dt);
    }
	if(friction){
		double O_inf_y = 0.5;
		for (int i=0; i < np; i++){
			ang_velocity[i] = 0.75*contact_torque[i]/radius_cubic[i];
			ang_velocity[i].y += O_inf_y;
		}
    }
	// evolve PBC
	shear_disp += vel_difference*dt;
	if (shear_disp > lx()){
		shear_disp -= lx();
	}
    updateInteractions();
	
	for (int i=0; i < np; i++){
		velocity_predictor[i] = velocity[i];
		ang_velocity_predictor[i] = ang_velocity[i];
	}
	
	/*********************************************************/
	/*                   Corrector                           */
	/*********************************************************/
	setContactForceToParticle();
    // build new Resistance matrix after move
    stokes_solver->prepareNewBuild_RFU("direct");
    addStokesDrag();
    buildLubricationTerms(false); // false: don't modify rhs, as we want to keep same Brownian force
    stokes_solver->complete_RFU();
	
    // get the intermediate brownian velocity
	stokes_solver->solve_CholTrans( v_Brownian_mid );

	if(flubcont_update){  // rebuild rhs
		stokes_solver->resetRHS();
		buildLubricationRHS();
	}
	else{  // don't rebuild rhs
		stokes_solver->setRHS(lub_cont_forces_init);
	}
	stokes_solver->solve(v_lub_cont_mid);
	
    stokes_solver->solvingIsDone();
    // update total velocity
    // first term is hydrodynamic + contact velocities
    // second term is Brownian velocities
    // third term is Brownian drift
    // fourth term for vx is the shear rate
    for (int i = 0; i < np; i++){
		int i3 = 3*i;
		velocity[i].set(0.5*(v_lub_cont_mid[i]    + v_Brownian_mid[i3] + position[i].z),
						0.5*(v_lub_cont_mid[i3+1] + v_Brownian_mid[i3+1]*zero_2Dsimu),
						0.5*(v_lub_cont_mid[i3+2] + v_Brownian_mid[i3+2]));

		velocity[i] -= 0.5*velocity_predictor[i];
    }
	if(friction){
		//
//		double O_inf_y = 0.5;
//		for (int i=0; i < np; i++){
//			ang_velocity[i] = 0.75*contact_torque[i]/radius_cubic[i];
//			ang_velocity[i].y += 0.5*O_inf_y;
//			ang_velocity[i] -= 0.5*ang_velocity_predictor[i];
//		}
//		//	ang_velocity[i] = 0.5*(ang_velocity_mp1st[i] + ang_velocity_mp2nd[i]);
	}
	for (int i=0; i < np; i++){
		displacement(i, velocity[i]*dt);
	}
	if (dimension == 2 ){
		for (int i=0; i < np; i++){
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	// update boxing system
	boxset->update();
	//	ksi_min=1.;
	checkNewInteraction();
	updateInteractions();
}
	
void
System::timeEvolution(int time_step){
	int ts_next = ts + time_step;
	checkNewInteraction();
	while (ts < ts_next){
		switch (integration_method) {
			case 0:
				timeEvolutionEulersMethod();
				break;
			case 1:
				timeEvolutionPredictorCorrectorMethod();
				break;
			case 2:
				timeEvolutionBrownian();
		}
		ts ++;
		shear_strain += dt;
	}
}

void
System::checkNewInteraction(){
	vector<int>::iterator it;
	vector<int>::iterator it_beg;
	vector<int>::iterator it_end;
	vec3d pos_diff;
	int zshift;
	double sq_dist;
	for (int i=0; i < np-1; i++){
		it_beg = boxset->neighborhood_begin(i);
		it_end = boxset->neighborhood_end(i);
		for (it = it_beg; it != it_end; it++){
			int j=*it;
			if(j>i){
				if ( interaction_partners[i].find(j) == interaction_partners[i].end() ){
					// distance is done in 3 steps because we need each information for Interaction creation
					pos_diff = position[j] - position[i];
					periodize_diff(pos_diff, zshift);
					sq_dist = pos_diff.sq_norm();

					double ri_rj_2 = 0.5*(radius[i] + radius[j]);
					double sq_dist_lim = sq_lub_max * ri_rj_2 * ri_rj_2;
					if ( sq_dist < sq_dist_lim){
						int interaction_new;
						if (deactivated_interaction.empty()){
							// add an interaction object.
							interaction_new = num_interaction;
							num_interaction ++;
						} else {
							// fill a deactivated interaction object.
							interaction_new = deactivated_interaction.front();
							deactivated_interaction.pop();
						}
						// new interaction
						interaction[interaction_new].activate(i, j, +pos_diff, sqrt(sq_dist), zshift);
					}
				}
			}
		}
	}
}

/* Check the distance between separating particles.
 * i < j
 *
 * A patch-up prescription to aboid
 * contact_pair[i][j] < 0 indicates separating particles to be checked.
 * contact_pair[i][j] = -1, the particles are near contact. So every time step, distance should be checked.a
 * contact_pair[i][j] < -1, the particles have some distance.
 */
void
System::updateInteractions(){
	//	int active_int=0;
	for (int k = 0; k < num_interaction; k++){
		bool switch_off = interaction[k].update();
		if(fix_interaction_status == false){
			if ( switch_off ){
				deactivated_interaction.push(k);
			}
		}
	}
}

void
System::stressReset(){
	for (int i=0; i < np; i++){
		for (int u=0; u < 5; u++){
			lubstress[i].elm[u]=0;
			contactstress[i].elm[u]=0;
			lubstress2[i].elm[u]=0;
		}
	}
}
void
System::stressBrownianReset(){
	for (int i=0; i < np; i++){
		for (int u=0; u < 5; u++){
			brownianstress[i].elm[u]=0;
		}
	}
	brownianstress_calc_nb = 0;
}

///*
// * Free-draining approximation
// */
//void
//System::updateVelocity(){
//	vec3d U_inf(0, 0, 0);
//	for (int i=0; i < np; i++){
//		U_inf.x = shear_rate*position[i].z;
//		relative_velocity[i] = (1.0/eta)*total_force[i];
//		velocity[i] = relative_velocity[i] + U_inf;
//	}
//	if(friction){
//		double O_inf_y = 0.5*shear_rate;
//		for (int i=0; i < np; i++){
//			ang_velocity[i] = 0.75*torque[i];
//			ang_velocity[i].y += O_inf_y;
//		}
//	}
//}

void
System::addStokesDrag(){
    for (int i = 0; i < np; i ++){
		stokes_solver->addToDiag_RFU(i, bgf_factor*radius[i]);
    }
}


// We solve A*(U-Uinf) = Gtilde*Einf ( in Jeffrey's notations )
// This method computes:
//  - elements of matrix A 
//  - vector Gtilde*Einf if rhs is true (default behavior)
void
System::buildLubricationTerms(bool rhs=true){
    
    double XAii, XAjj, XAij, XAji;
	
    double GEi[3];
    double GEj[3];
    
    set<Interaction*>::iterator it;
    int j;
    Interaction *inter;
    for (int i = 0; i < np - 1; i ++){
		for (it = interaction_list[i].begin() ; it != interaction_list[i].end(); it ++){
			inter=*it;
			j=inter->partner(i);
			if(j>i){
				inter->XA(XAii, XAij, XAji, XAjj);
				
				stokes_solver->addToDiagBlock_RFU(inter->nr_vec, i, inter->a0 * XAii);
				stokes_solver->addToDiagBlock_RFU(inter->nr_vec, j, inter->a1 * XAjj);
				stokes_solver->appendToOffDiagBlock_RFU(inter->nr_vec, i, j, 0.5 * inter->ro * XAji);
				
				if(rhs){
					inter->GE(GEi, GEj);  // G*E_\infty term
					for(int u=0; u<3; u++){
						stokes_solver->addToRHS( 3*i + u, GEi[u] );
						stokes_solver->addToRHS( 3*j + u, GEj[u] );
					}
				}
			}
		}
		stokes_solver->doneBlocks(i);
    }
}

void
System::buildLubricationRHS(){
    
    double GEi[3];
    double GEj[3];
    
    set<Interaction*>::iterator it;
    int j;
    Interaction *inter;
    for (int i = 0; i < np - 1; i ++){
		for (it = interaction_list[i].begin() ; it != interaction_list[i].end(); it ++){
			inter=*it;
			j=inter->partner(i);
			if(j>i){
				inter->GE(GEi, GEj);  // G*E_\infty term
				for(int u=0; u<3; u++){
					stokes_solver->addToRHS( 3*i + u, GEi[u] );
					stokes_solver->addToRHS( 3*j + u, GEj[u] );
				}
			}
		}
    }
}


void
System::setContactForceToParticle(){
	for (int i = 0; i < np; i++){
		contact_force[i].reset();
		contact_torque[i].reset();
	}
	for (int k=0; k < num_interaction; k++){
		interaction[k].addUpContactForceTorque();
	}
}

void
System::buildContactTerms(){
	//	calcContactForces();
    // add contact force
    for (int i = 0; i < np; i++){
		int i3 = 3*i;
		stokes_solver->addToRHS(i3  , contact_force[i].x);
		stokes_solver->addToRHS(i3+1, contact_force[i].y);
		stokes_solver->addToRHS(i3+2, contact_force[i].z);
    }
}

/*
 *
 *
 * This function deteremins velocity[] and ang_velocity[].
 *
 */
void
System::updateVelocityLubrication(){
    stokes_solver->resetRHS();
    stokes_solver->prepareNewBuild_RFU("direct");
//	stokes_solver->prepareNewBuild_RFU("iterative");

    addStokesDrag();
    buildLubricationTerms();
    stokes_solver->complete_RFU();
    buildContactTerms();
    stokes_solver->solve(v_lub_cont);
	//stokes_solver->print_RFU();
	/* TEST IMPLEMENTATION
	 * SDFF : Stokes drag force factor:
	 * SDFF = 1.0 : full drag forces from the undisturbed background flow.
	 * SDFF = 0.0 : no drag force from the undisturbed background flow.
	 */
    for (int i = 0; i < np; i++){
		int i3 = 3*i;
		relative_velocity_lub_cont[i].x = v_lub_cont[i3];
		relative_velocity_lub_cont[i].y = v_lub_cont[i3+1];
		relative_velocity_lub_cont[i].z = v_lub_cont[i3+2];
		velocity[i].x = relative_velocity_lub_cont[i].x + position[i].z;
		velocity[i].y = relative_velocity_lub_cont[i].y;
		velocity[i].z = relative_velocity_lub_cont[i].z;
    }

	// Tc - 8 pi eta a^3 (omega - omega_inf) = 0
	// U0 = a gammadot
	// F0 = 6 pi eta a U0 = 6 pi eta a^2 gammadot
	// Tc/a F0 - 8 pi eta a'^3 (omega - omega_inf) / a F0 = 0
	// \hat{Tc} - (4/3) (a'/a)^3(omega - omega_inf) / gammadot = 0
	//  omega = (3/4) (a/a')^3\hat{Tc} gammadot + omega_inf;
	/*
	 * This must be modified for polidisperse system.
	 *
	 *
	 */
	double O_inf_y = 0.5;
	for (int i=0; i < np; i++){
		ang_velocity[i] = 0.75*contact_torque[i]/radius_cubic[i];
		ang_velocity[i].y += O_inf_y;
	}
    
    stokes_solver->solvingIsDone();
}

void System::updateVelocityLubricationBrownian(){
	/***************************************************************************
    *  This routine implements a predictor-corrector algorithm                 *
    *  for the dynamics with Brownian motion.                                  *
    *  The algorithm is the one of Melrose & Ball 1997.                        * 
    *                                                                          *
	*  X_pred = X(t) + V(t)*dt                                                 *
    *  X(t+dt) = X_pred + 0.5*(V(t+dt)-V(t))*dt                                *
    *                                                                          *
    *  Parameters                                                              *
	*  They essentially controls what enter V(t) and V(t+dt):                  *
	*   * displubcont : - if true  V(t) = V^{B}(t) + V^{H}(t) + V^{C}(t)       *
    *                   - if false V(t) = V^{B}(t)                             *
	*   * flubcont_update : - requires displubcont                             *
	*                       - allows dt^2 scheme for contact and hydro         *
    *                         velocities                                       *
	*                       - if true :                                        *
    *                          R_FU(t+dt) V^{*}(t+dt) = F^{*}(t+dt)            *
    *                          with * = C or H                                 *
	*                       - if false :                                       *
    *                          R_FU(t+dt) V^{*}(t+dt) = F^{*}(t)               *
    *                                                                          *
	****************************************************************************/
	bool displubcont = true;
	bool flubcont_update = true;

	int zero_2Dsimu;
	if (dimension == 2){
		zero_2Dsimu = 0;
	}else{
		zero_2Dsimu = 1;
	}
	
	/*********************************************************/
	/*                    First Step                         */
	/*********************************************************/

	stokes_solver->resetRHS();

    stokes_solver->prepareNewBuild_RFU("direct");

    addStokesDrag();
    buildLubricationTerms();

    stokes_solver->complete_RFU();
    buildContactTerms();
    stokes_solver->solve(v_lub_cont);
	
	if(displubcont){
		stokes_solver->getRHS(lub_cont_forces_init);
	}
    // now the Brownian part of the velocity:
    // predictor-corrector algortithm (see Melrose & Ball, 1997)
    // 
    // we do not call solvingIsDone() before new solve(), because 
    // R_FU has not changed, so same factorization is safely used

	stokes_solver->setRHS( fb->generate_invLFb() );
	stokes_solver->solve_CholTrans( v_Brownian_init );

	stokes_solver->solvingIsDone();



    // move particles to intermediate point
    for (int i=0; i < np; i++){
		int i3 = 3*i;
		
		vec3d dr(v_Brownian_init[i3]*dt,
				 v_Brownian_init[i3+1]*dt*zero_2Dsimu,
				 v_Brownian_init[i3+2]*dt);
		
		if(displubcont){
			vec3d dr_lub_cont(v_lub_cont[i3]*dt,
							  v_lub_cont[i3+1]*dt,
							  v_lub_cont[i3+2]*dt);
			dr += dr_lub_cont;
		}

		displacement(i, dr);
    }
    updateInteractions();

	/*********************************************************/
	/*                   Second Step                         */
	/*********************************************************/

    // build new Resistance matrix after move
    stokes_solver->prepareNewBuild_RFU("direct");
    addStokesDrag();
    buildLubricationTerms(false); // false: don't modify rhs, as we want to keep same Brownian force
    stokes_solver->complete_RFU();

    // get the intermediate brownian velocity
	stokes_solver->solve_CholTrans( v_Brownian_mid );

	if(displubcont){
		if(flubcont_update){  // rebuild rhs
			stokes_solver->resetRHS();
			buildLubricationRHS();
		}
		else{  // don't rebuild rhs
			stokes_solver->setRHS(lub_cont_forces_init);
		}
		stokes_solver->solve(v_lub_cont_mid);
	}

    stokes_solver->solvingIsDone();
    // update total velocity
    // first term is hydrodynamic + contact velocities
    // second term is Brownian velocities
    // third term is Brownian drift
    // fourth term for vx is the shear rate
	vec3d relative_velocity_lub_cont;
	vec3d relative_velocity_brownian;
    for (int i = 0; i < np; i++){
		int i3 = 3*i;
		if(displubcont){
			relative_velocity_lub_cont.x = 0.5*( v_lub_cont_mid[i3  ] - v_lub_cont[i3  ] );
			relative_velocity_lub_cont.y = 0.5*( v_lub_cont_mid[i3+1] - v_lub_cont[i3+1] );
			relative_velocity_lub_cont.z = 0.5*( v_lub_cont_mid[i3+2] - v_lub_cont[i3+2] );
		}
		else{
			relative_velocity_lub_cont.x = v_lub_cont[i3];
			relative_velocity_lub_cont.y = v_lub_cont[i3+1];
			relative_velocity_lub_cont.z = v_lub_cont[i3+2];
		}

		relative_velocity_brownian.x = 0.5*( v_Brownian_mid[i3  ] - v_Brownian_init[i3  ] );
		relative_velocity_brownian.y = 0.5*( v_Brownian_mid[i3+1] - v_Brownian_init[i3+1] )*zero_2Dsimu;
		relative_velocity_brownian.z = 0.5*( v_Brownian_mid[i3+2] - v_Brownian_init[i3+2] );
		
		velocity[i] = relative_velocity_lub_cont + relative_velocity_brownian;
		velocity[i].x += position[i].z;
    }
	
    if(friction){
		double O_inf_y = 0.5;
		for (int i=0; i < np; i++){
			ang_velocity[i] = 0.75*contact_torque[i]/radius_cubic[i];
			ang_velocity[i].y += O_inf_y;
		}
    }
	
	static int cnt = 0;
	if (cnt ++ % 10 == 0){
		for (int k = 0; k < num_interaction; k++){
			if (interaction[k].gap_nondim()  < 0.1){
				int i = interaction[k].particle_num[0];
				int j = interaction[k].particle_num[1];
				int i3 = i*3;
				vec3d v_brownian1(v_Brownian_mid[i3  ] - v_Brownian_init[i3],
								  v_Brownian_mid[i3 +1 ] - v_Brownian_init[i3 +1],
								  v_Brownian_mid[i3 +2 ] - v_Brownian_init[i3 +2]);
				i3 = j*3;
				vec3d v_brownian2(v_Brownian_mid[i3  ] - v_Brownian_init[i3],
								  v_Brownian_mid[i3 +1 ] - v_Brownian_init[i3 +1],
								  v_Brownian_mid[i3 +2 ] - v_Brownian_init[i3 +2]);
				
				double dr_normal = dot(v_brownian1-v_brownian2, interaction[k].nr_vec)*dt;

				cout << interaction[k].gap_nondim() << ' ' << abs(dr_normal) << endl;
			}
		}
	}
}




void
System::displacement(int i, const vec3d &dr){
	position[i] += dr;
	periodize(position[i]);
	boxset->box(i);
}

// [0,l]
void
System::periodize(vec3d &pos){
	if (pos.z > lz() ){
		pos.z -= lz();
		pos.x -= shear_disp;
	} else if ( pos.z < 0 ){
		pos.z += lz();
		pos.x += shear_disp;
	}
	while ( pos.x > lx() ){
		pos.x -= lx();
	}
	while (pos.x < 0 ){
		pos.x += lx();
	}
	if (dimension == 3){
		if ( pos.y > ly() ){
			pos.y -= ly();
		} else if (pos.y < 0 ){
			pos.y += ly();
		}
	}
}

// [-l/2,l/2]
void
System::periodize_diff(vec3d &pos_diff){
	if (pos_diff.z > lz2() ){
		pos_diff.z -= lz();
		pos_diff.x -= shear_disp;
	} else if ( pos_diff.z < -lz2() ){
		pos_diff.z += lz();
		pos_diff.x += shear_disp;
	}
	while ( pos_diff.x > lx2() ){
		pos_diff.x -= lx();
	}
	while (pos_diff.x < -lx2() ){
		pos_diff.x += lx();
	}
	if (dimension == 3){
		if ( pos_diff.y > ly2() ){
			pos_diff.y -= ly();
		} else if (pos_diff.y < -ly2() ){
			pos_diff.y += ly();
		}
	}
}

// periodize + give z_shift= number of boundaries crossed in z-direction
void
System::periodize_diff(vec3d &pos_diff, int &zshift){
	/*
	 * The displacement of the second particle along z direction
	 * is zshift * lz;
	 */
	if (pos_diff.z > lz2() ){
		pos_diff.z -= lz();
		pos_diff.x -= shear_disp;
		zshift = -1;
	} else if ( pos_diff.z < -lz2() ){
		pos_diff.z += lz();
		pos_diff.x += shear_disp;
		zshift = +1;
	} else{
		zshift = 0;
	}
	while ( pos_diff.x > lx2() ){
		pos_diff.x -= lx();
	}
	while (pos_diff.x < -lx2() ){
		pos_diff.x += lx();
	}
	if (dimension == 3){
		if ( pos_diff.y > ly2() ){
			pos_diff.y -= ly();
		} else if (pos_diff.y < -ly2() ){
			pos_diff.y += ly();
		}
	}
}

/*
 * Distance between particle i and particle j
 */
double
System::distance(int i, int j){
	return sqrt(sq_distance(i,j));
}

/*
 * Square distance between particle i and particle j
 */
double
System::sq_distance(int i, int j){
	vec3d pos_diff = position[j] - position[i];
	periodize_diff(pos_diff);
	if (dimension == 3){
		return pos_diff.sq_norm();
	} else {
	 	return pos_diff.sq_norm_xz();
	}
}

void
System::calcStress(){
	stressReset();

	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			interaction[k].evaluateLubricationForce();
			interaction[k].addLubricationStress();
			if (interaction[k].contact){
				interaction[k].addContactStress();
			}
		}
	}
	if(brownian)
	 	calcBrownianStress();

	for (int u=0; u < 5; u++){
		total_lub_stress[u] = 0;
		total_contact_stress[u] = 0;
		total_brownian_stress[u] = 0;
	}
	
	for (int i=0; i < np; i++){
		for (int u=0; u < 5; u++){
			total_lub_stress[u] += lubstress[i].elm[u];
			total_contact_stress[u] += contactstress[i].elm[u];
			total_brownian_stress[u] += brownianstress[i].elm[u];
		}
	}
	
	/*
	 * The term 5.0/9 is the one-body part
	 *
	 */
	total_stress_bgf = 0;
	for (int i=0; i < np; i++){
		double a = radius[i];
		total_stress_bgf += (5.0/9)*bgf_factor*a*a*a;
	}

	stressBrownianReset();
}

void
System::analyzeState(){
	minvalue_gap_nondim = lz();
	// for analysis
	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			if (interaction[k].gap_nondim() < minvalue_gap_nondim){
				minvalue_gap_nondim = interaction[k].gap_nondim();
			}
			interaction[k].recordTrajectory();
		}
	}
}

void
System::setSystemVolume(){
	if (dimension == 2){
		system_volume = _lx*_lz*2*radius_max;
	} else {
		system_volume = _lx*_ly*_lz;
	}
}

