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

	if (!total_force)
		delete [] total_force;
	if (!lubrication_force)
		delete [] lubrication_force;
	if (!contact_force)
		delete [] contact_force;

	if (!brownian_force)
	  delete [] brownian_force;

	if (!torque)
		delete [] torque;

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
	velocity_mp1st.resize(np);
	velocity_mp2nd.resize(np);
	relative_velocity = new vec3d [np];
	relative_velocity_lub_cont = new vec3d [np];
	relative_velocity_brownian = new vec3d [np];
	ang_velocity = new vec3d [np];
	ang_velocity_mp1st.resize(np);
	ang_velocity_mp2nd.resize(np);
	total_force = new vec3d [np];
	lubrication_force = new vec3d [np];
	contact_force = new vec3d [np];
	brownian_force = new vec3d [np];
	torque = new vec3d [np];
	lub_force = new vec3d [np];
	nearing_number.resize(np);
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
		integration_method = 0; // > force Euler
	}

	allocateRessources();

	for (int i=0; i < np; i++){
		position[i] = initial_positions[i];
		radius[i] = radii[i];
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
	double O_inf_y = 0.5*shear_rate/2.0;
	for (int i=0; i < np; i++){
		ang_velocity[i].set(0, O_inf_y, 0);
		torque[i].reset();
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
	vel_difference = shear_rate*_lz;
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
	gap_min = 1;
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
	if (lubrication){
		// Lubrication dynamics
		if(brownian){
			updateVelocityLubricationBrownian();
		} else{
			updateVelocityLubrication();
		}
	} else {
		// Free-draining approximation
		updateVelocity();
	}
	deltaTimeEvolution();
}

void
System::timeEvolutionPredictorCorrectorMethod(){
	/* Ball&Melrose 1997
	 * x'(t+dt) = x(t) + V^{-}
	 */
	setContactForceToParticle();
	
	
	updateVelocityLubrication(velocity_mp1st, ang_velocity_mp1st);
	deltaTimeEvolution_firststep();
	/*
	 * x(t + dt) = x(t) + 0.5*(V^{+}+V^{-})
	 *           = x'(t+dt) + 0.5*(V^{+}-V^{-})
	 */
	setContactForceToParticle();
	updateVelocityLubrication(velocity_mp2nd, ang_velocity_mp2nd);
	deltaTimeEvolution_secondstep();
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
		displacement(i, (velocity[i]*dt));
	}
	if (draw_rotation_2d){
		for (int i=0; i < np; i++){
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	
	// update boxing system
	boxset->update();
	
	//	ksi_min=1.;
	checkNewInteraction();
	updateInteractions();
	//	cout << " ksi_min " << ksi_min << endl;
}

void
System::deltaTimeEvolution_firststep(){
	shear_disp += vel_difference*dt;
	if (shear_disp > lx()){
		shear_disp -= lx();
	}
	// evolve PBC
	// move particles
	for (int i=0; i < np; i++){
		displacement(i, velocity_mp1st[i]*dt);
		ang_velocity[i] = ang_velocity_mp1st[i];
	}
	// update boxing system
	boxset->update();
	checkNewInteraction();
	updateInteractions();
}


void
System::deltaTimeEvolution_secondstep(){
	/*
	 *
	 * x'(t+dt) = x(t) + V^{-}
	 * x(t) = x'(t+dt)- V^{-}
	 * x(t + dt) = x(t) + 0.5*(V^{+}+V^{-})*dt
	 *           = x'(t+dt) + 0.5*(V^{+}-V^{-})*dt
	 */
	// evolve PBC
	// move particles
	double dt_2 = 0.5*dt;
	for (int i=0; i < np; i++){
		displacement(i, (velocity_mp2nd[i] - velocity_mp1st[i])*dt_2);
		ang_velocity[i] = 0.5*(ang_velocity_mp1st[i]+ang_velocity_mp2nd[i]);
	}
	// update boxing system
	boxset->update();
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
		}
		ts ++;
		shear_strain += shear_rate*dt;
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
System::updateInteractions(const bool switch_off_allowed){
	//	int active_int=0;
	for (int k = 0; k < num_interaction; k++){
		bool switch_off = interaction[k].update(switch_off_allowed);
		if(switch_off)
			deactivated_interaction.push(k);
	}
}

//void
//System::forceReset(){
//	for (int i=0; i < np; i++){
//		total_force[i].reset();
//		lubrication_force[i].reset();
//		contact_force[i].reset();
//		brownian_force[i].reset();
//	}
//}
//
//void
//System::torqueReset(){
//	for (int i=0; i < np; i++){
//		torque[i].reset();
//	}
//}

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

/*
 * Free-draining approximation
 */
void
System::updateVelocity(){
	vec3d U_inf(0, 0, 0);
	for (int i=0; i < np; i++){
		U_inf.x = shear_rate*position[i].z;
		relative_velocity[i] = (1.0/eta)*total_force[i];
		velocity[i] = relative_velocity[i] + U_inf;
	}
	if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < np; i++){
			ang_velocity[i] = (1.33333/eta)*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
	}
}

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
		torque[i].reset();
	}
	for (int k=0; k < num_interaction; k++){
		int i = interaction[k].particle_num[0];
		int j = interaction[k].particle_num[1];
		interaction[k].addUpContactForce(contact_force[i], contact_force[j]);
		//	  interaction[k].addUpContactTorque(contact_torque[i], contact_torque[j]);
		total_force[i] += contact_force[i];
		total_force[j] += contact_force[j];
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
		
		velocity[i].x = relative_velocity_lub_cont[i].x + shear_rate*position[i].z;
		velocity[i].y = relative_velocity_lub_cont[i].y;
		velocity[i].z = relative_velocity_lub_cont[i].z;
    }
	
    if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < np; i++){
			ang_velocity[i] = 1.33333*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
    }
    
    stokes_solver->solvingIsDone();
}

void System::updateVelocityLubrication(vector<vec3d> &velocity_, vector<vec3d> &ang_velocity_){
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
		
		velocity_[i].x = relative_velocity_lub_cont[i].x + shear_rate*position[i].z;
		velocity_[i].y = relative_velocity_lub_cont[i].y;
		velocity_[i].z = relative_velocity_lub_cont[i].z;
    }
	
    if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < np; i++){
			ang_velocity_[i] = 1.33333*torque[i];
			ang_velocity_[i].y += O_inf_y;
		}
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
	bool flubcont_update = false;


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
    updateInteractions(false);
	

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
    for (int i = 0; i < np; i++){
		int i3 = 3*i;
		if(displubcont){
			relative_velocity_lub_cont[i].x = 0.5*( v_lub_cont_mid[i3  ] - v_lub_cont[i3  ] );
			relative_velocity_lub_cont[i].y = 0.5*( v_lub_cont_mid[i3+1] - v_lub_cont[i3+1] );
			relative_velocity_lub_cont[i].z = 0.5*( v_lub_cont_mid[i3+2] - v_lub_cont[i3+2] );
		}
		else{
			relative_velocity_lub_cont[i].x = v_lub_cont[i3];
			relative_velocity_lub_cont[i].y = v_lub_cont[i3+1];
			relative_velocity_lub_cont[i].z = v_lub_cont[i3+2];
		}

		relative_velocity_brownian[i].x = 0.5*( v_Brownian_mid[i3  ] - v_Brownian_init[i3  ] );
		relative_velocity_brownian[i].y = 0.5*( v_Brownian_mid[i3+1] - v_Brownian_init[i3+1] )*zero_2Dsimu;
		relative_velocity_brownian[i].z = 0.5*( v_Brownian_mid[i3+2] - v_Brownian_init[i3+2] );
		
		velocity[i].x = relative_velocity_lub_cont[i].x + relative_velocity_brownian[i].x + shear_rate*position[i].z;
		velocity[i].y = relative_velocity_lub_cont[i].y + relative_velocity_brownian[i].y;
		velocity[i].z = relative_velocity_lub_cont[i].z + relative_velocity_brownian[i].z;
    }


    if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < np; i++){
			ang_velocity[i] = 1.33333*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
    }
}


void System::calcBrownianStress(){
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
    buildLubricationTerms();

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
	    step_stresslet[i].elm[u] += stresslet_i_init.elm[u];
	    step_stresslet[j].elm[u] += stresslet_j_init.elm[u];
	  }
    }

	for (int u=0; u < 5; u++){
	  for (int i=0; i < np; i++){
		total_step_stresslet.elm[u] += step_stresslet[i].elm[u];
	  }
	  //cout << total_step_stresslet.elm[u]/np << " " ;
	}

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
	    step_stresslet[i].elm[u] -= stresslet_i_mid.elm[u];
	    step_stresslet[j].elm[u] -= stresslet_j_mid.elm[u];
		total_step_stresslet.elm[u] -= stresslet_i_mid.elm[u] + stresslet_j_mid.elm[u];
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

	for(int u=0; u < 10 ; u++){
		cnt_nearing_number[u] = 0;
	}
	for (int i =0; i < np; i++){
		cnt_nearing_number[ nearing_number[i] ] ++;
	}
	
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
	gap_min = lz();
	double sum_overlap = 0;
	int cnt_overlap = 0;
	max_nearing_time = 0;
	
	for (int i=0; i < np; i++){
		nearing_number[i] = 0;
	}
	
	num_nearing = 0;
	// for analysis
	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			
			if (interaction[k].gap() < 0){
				sum_overlap +=interaction[k].gap();
				cnt_overlap ++;
			}
			
			if (interaction[k].gap() < gap_min){
				gap_min = interaction[k].gap();
			}

			interaction[k].recordTrajectory();
			
			if (interaction[k].near){
				num_nearing ++;
				nearing_number[interaction[k].particle_num[0]] ++;
				nearing_number[interaction[k].particle_num[1]] ++;
				if ( max_nearing_time < interaction[k].nearingTime()){
					max_nearing_time = interaction[k].nearingTime();
				}
			}
		}
	}
	
	double sum_nearing_time = 0;
	for(int k = 0; k < nearing_time_record.size(); k++){
		sum_nearing_time += nearing_time_record[k];
	}
	ave_overlap = sum_overlap / cnt_overlap;
	if (nearing_time_record.size() > 0){
		ave_nearing_time = sum_nearing_time / nearing_time_record.size();
	} else {
		ave_nearing_time = 0;
	}
	nearing_time_record.clear();
	
}

void
System::setSystemVolume(){
	if (dimension == 2){
		system_volume = _lx*_lz*2*radius_max;
	} else {
		system_volume = _lx*_ly*_lz;
	}
}

