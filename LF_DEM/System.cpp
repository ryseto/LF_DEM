//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
#include <sstream>
#define DELETE(x) if(x){delete [] x; x = NULL;}

System::System(){
	brownian = false;
};

System::~System(){
	DELETE(position);
	DELETE(radius);
	DELETE(radius_cubic);
	DELETE(angle);
	DELETE(velocity);
	DELETE(ang_velocity);
	if (integration_method >= 1) {
		DELETE(velocity_predictor);
		DELETE(ang_velocity_predictor);
	}
	DELETE(contact_force);
	DELETE(contact_torque);
	DELETE(lubstress);
	DELETE(bgfstress);
	DELETE(contactstressGU);
	DELETE(colloidalstressGU);
	DELETE(brownianstress);
	DELETE(interaction);
	DELETE(interaction_list);
	DELETE(interaction_partners);
	DELETE(v_lub_cont);
	DELETE(v_hydro);
	DELETE(v_cont);
	DELETE(v_colloidal);
	if(brownian){
		DELETE(fb);
		DELETE(v_Brownian_init);
		DELETE(v_Brownian_mid);
		DELETE(v_lub_cont_mid);
		DELETE(lub_cont_forces_init);
	}
};

void
System::allocateRessources(){
	position = new vec3d [_np];
	radius = new double [_np];
	radius_cubic = new double [_np];
	angle = new double [_np];
	velocity = new vec3d [_np];
	ang_velocity = new vec3d [_np];
	if (integration_method >= 1) {
		ang_velocity_predictor = new vec3d [_np];
		velocity_predictor = new vec3d [_np];
	}
	contact_force = new vec3d [_np];
	contact_torque = new vec3d [_np];
	colloidal_force = new vec3d [_np];
	lubstress = new stresslet [_np];
	bgfstress = new stresslet [_np];
	contactstressGU = new stresslet [_np];
	colloidalstressGU = new stresslet [_np];
	brownianstress = new stresslet [_np];
	int maxnum_interactionpair_per_particle = 15;
	maxnum_interactionpair = maxnum_interactionpair_per_particle*_np;
	interaction = new Interaction [maxnum_interactionpair];
	interaction_list = new set <Interaction*> [_np];
	interaction_partners = new set <int> [_np];
	dof = 3;
	linalg_size = dof*_np;
	v_lub_cont = new double [linalg_size];
	v_cont = new double [linalg_size];
	v_hydro = new double [linalg_size];
	v_colloidal = new double [linalg_size];
	for (int i=0; i<linalg_size; i++) {
		v_lub_cont[i] = 0;
		v_cont[i] = 0;
		v_hydro[i] = 0;
		v_colloidal[i] = 0;
	}
	if (brownian) {
	    v_Brownian_init = new double [linalg_size];
	    v_Brownian_mid = new double [linalg_size];
		v_lub_cont_mid = new double [linalg_size];
		lub_cont_forces_init = new double [linalg_size];
		for(int i=0;i<linalg_size;i++){
			v_Brownian_init[i] = 0;
			v_Brownian_mid[i] = 0;
			lub_cont_forces_init[i] = 0;
		}
		fb = new BrownianForce(this);
	}
	stokes_solver.init(_np, brownian);
}


void
System::setupSystemForGenerateInit(){
	for (int i=0; i < _np; i++) {
		radius_cubic[i] = radius[i]*radius[i]*radius[i];
		angle[i] = 0;
	}
	for (int k=0; k<maxnum_interactionpair ; k++) {
		interaction[k].init(this);
		interaction[k].label = k;
	}
	for (int i=0; i<_np; i++) {
		velocity[i].reset();
	}
	shear_strain = 0;
	shear_disp = 0;
	num_interaction = 0;
	sq_lub_max = lub_max*lub_max; // square of lubrication cutoff length.
	contact_relaxzation_time = 1e-3;
	kn = 2000;
	kt = 0;
	friction = false;
	colloidalforce = true;
	if (contact_relaxzation_time < 0) {
		// 1/(h+c) --> 1/c
		lub_coeff_contact = 1/lub_reduce_parameter;
	} else {
		/* t = beta/kn
		 *  beta = t*kn
		 * lub_coeff_contact = 4*beta = 4*kn*contact_relaxzation_time
		 */
		lub_coeff_contact = 4*kn*contact_relaxzation_time;
	}
	ts = 0;
	stokes_solver.initialize();
	initializeBoxing();
	checkNewInteraction();
}


void
System::setupSystem(const vector<vec3d> &initial_positions,
					const vector<double> &radii){
	if (kb_T == 0) {
		brownian = false;
	} else {
		brownian = true;
		integration_method = 2; // > force Euler
	}
	if (mu_static > 0) {
		friction = true;
	} else {
		friction = false;
	}
	if (colloidalforce_length > 0) {
		/*
		 * The diemnsionless shear rate is defined as follows:
		 * dimensionless_shear_rate = 6pi*eta*a^2*shear_rate/colloidalforce_amplitude
		 * Under the unit of this simulation
		 * 6pi*eta*a^2*shear_rate is set to 1.
		 */
		colloidalforce_amplitude = 1/dimensionless_shear_rate;
		colloidalforce = true;
		cerr << "Colloidal force" << endl;
	} else {
		colloidalforce = false;
		cerr << "No colloidal force" << endl;
	}
	allocateRessources();
	for (int i=0; i<_np; i++) {
		position[i] = initial_positions[i];
		radius[i] = radii[i];
		radius_cubic[i] = radius[i]*radius[i]*radius[i];
		angle[i] = 0;
	}
	for (int k=0; k<maxnum_interactionpair ; k++) {
		interaction[k].init(this);
		interaction[k].label = k;
	}
	for (int i=0; i<_np; i++) {
		velocity[i].set(position[i].z, 0, 0);
		ang_velocity[i].set(0, 0.5, 0);
	}
	shear_strain = 0;
	shear_disp = 0;
	num_interaction = 0;
	sq_lub_max = lub_max*lub_max; // square of lubrication cutoff length.
	if (contact_relaxzation_time < 0) {
		// 1/(h+c) --> 1/c
		lub_coeff_contact = 1/lub_reduce_parameter;
	} else {
		/* t = beta/kn
		 *  beta = t*kn
		 * lub_coeff_contact = 4*beta = 4*kn*contact_relaxzation_time
		 */
		lub_coeff_contact = 4*kn*contact_relaxzation_time;
	}
	cerr << "lub_coeff_contact = " << lub_coeff_contact << endl;
	ts = 0;
	shear_disp = 0;
	vel_difference = _lz;
	/*
	 * dt_mid: the intermediate time step for the mid-point
	 * algortithm. dt/dt_mid = dt_ratio
	 * Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 * ASD code from Brady has dt_ratio=150
	 */
	//dt_mid = dt/dt_ratio; // UNUSED NOW: ALGO EQUIVALENT TO dt_ratio = 1 (Melrose & Ball 1997)
	stokes_solver.initialize();
	// initialize the brownian force after the solver, as it assumes
	// the cholmod_common of the solver is already initialized
	if (brownian) {
		fb->init();
	}
	brownianstress_calc_nb = 0;
	initializeBoxing();
	checkNewInteraction();
	cnt_monitored_data = 0;
	setSystemVolume();
	for (int i=0; i<_np; i++) {
		bgfstress[i].set(0, 0, (5.0/9)*bgf_factor*radius_cubic[i], 0, 0, 0);
	}
}

void
System::initializeBoxing(){// need to know radii first
	double max_radius = 0;
	for (int i=0; i < _np; i++) {
		if (radius[i] > max_radius) {
			max_radius = radius[i];
		}
	}
	boxset.init(lub_max*max_radius, this);
	for (int i=0; i<_np; i++) {
		boxset.box(i);
	}
	boxset.update();
}

void
System::timeEvolutionEulersMethod(){
	setContactForceToParticle();
	setColloidalForceToParticle();
	updateVelocityLubrication();
	deltaTimeEvolution();
}

void
System::timeEvolutionPredictorCorrectorMethod(){
	/*
	 * Predictor
	 * x'(t+dt) = x(t) + V^{-}dt
	 * x(t)     = x'(t+dt) - V^{-}dt
	 * Corrector
	 * x(t + dt) = x(t)     + 0.5*(V^{+}+V^{-})*dt
	 *           = x'(t+dt) + 0.5*(V^{+}-V^{-})*dt
	 */
	/* predictore */
	setContactForceToParticle();
	setColloidalForceToParticle();
	updateVelocityLubrication();
	deltaTimeEvolutionPredictor();
	/* corrector */
	setContactForceToParticle();
	setColloidalForceToParticle();
	updateVelocityLubrication();
	deltaTimeEvolutionCorrector();
}

void
System::deltaTimeEvolution(){
	// evolve PBC
	shear_disp += vel_difference*dt;
	if (shear_disp > lx()) {
		shear_disp -= lx();
	}
	// move particles
	for (int i=0; i<_np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (dimension == 2) {
		for (int i=0; i<_np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	// update boxing system
	boxset.update();
	checkNewInteraction();
	in_predictor = true;
	in_corrector = true;
	updateInteractions();
}

void
System::deltaTimeEvolutionPredictor(){
	/* The periodic boundary condition is updated in predictor.
	 * It must not be updated in corrector.
	 */
	shear_disp += vel_difference*dt;
	if (shear_disp >= lx()) {
		shear_disp -= lx();
	}
	for (int i=0; i<_np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (dimension == 2) {
		for (int i=0; i<_np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	/* In predictor, the values of interactions is updated,
	 * but the statuses are fixed by using boolean `fix_interaction_status'
	 */
	in_predictor = true;
	in_corrector = false;
	updateInteractions();
	/*
	 * Keep V^{-} to use them in the corrector.
	 */
	for (int i=0; i<_np; i++) {
		velocity_predictor[i] = velocity[i];
		ang_velocity_predictor[i] = ang_velocity[i];
	}
}

void
System::deltaTimeEvolutionCorrector(){
	for (int i=0; i<_np; i++) {
		velocity[i] = 0.5*(velocity[i]-velocity_predictor[i]);
		ang_velocity[i] = 0.5*(ang_velocity[i]-ang_velocity_predictor[i]);
	}
	for (int i=0; i<_np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (dimension == 2) {
		for (int i=0; i<_np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	// update boxing system
	boxset.update();
	checkNewInteraction();
	/*
	 * Interaction
	 *
	 */
	in_predictor = false;
	in_corrector = true;
	updateInteractions(); // false --> in corrector
	/* In deltaTimeEvolutionCorrector,
	 * velocity[] and ang_velocity[]
	 * are virtual velocities to correct the predictor.
	 *
	 * This function reverts them by (2):
	 * V = 0.5*(V^{+}-V^{-})   (1)
	 * V += V^{-}              (2)
	 * V = 0.5*(V^{+}+V^{-})   (3)
	 *
	 * [Note]
	 * Contact velocities in the interaction objects
	 * are not right ones.
	 */
	for (int i=0; i<_np; i++) {
		velocity[i] += velocity_predictor[i];
		ang_velocity[i] += ang_velocity_predictor[i];
	}
}

void System::timeEvolutionBrownian(){
	int zero_2Dsimu;
	if (dimension == 2) {
		zero_2Dsimu = 0;
	} else {
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
	stokes_solver.resetRHS();
    stokes_solver.prepareNewBuild_RFU("direct");
    addStokesDrag();
    buildLubricationTerms();
	
    stokes_solver.complete_RFU();
    buildContactTerms();
    stokes_solver.solve(v_lub_cont);
	
	if (displubcont) {
		stokes_solver.getRHS(lub_cont_forces_init);
	}
    // now the Brownian part of the velocity:
    // predictor-corrector algortithm (see Melrose & Ball, 1997)
    //
    // we do not call solvingIsDone() before new solve(), because
    // R_FU has not changed, so same factorization is safely used
	
	stokes_solver.setRHS( fb->generate_invLFb() );
	stokes_solver.solve_CholTrans( v_Brownian_init );
	stokes_solver.solvingIsDone();
	
    // move particles to intermediate point
	
    for (int i=0; i<_np; i++) {
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
	if (friction) {
		double O_inf_y = 0.5;
		for (int i=0; i<_np; i++) {
			ang_velocity[i] = 0.75*contact_torque[i]/radius_cubic[i];
			ang_velocity[i].y += O_inf_y;
		}
    }
	// evolve PBC
	shear_disp += vel_difference*dt;
	if (shear_disp >= lx()) {
		shear_disp -= lx();
	}
    updateInteractions();
	for (int i=0; i<_np; i++) {
		velocity_predictor[i] = velocity[i];
		ang_velocity_predictor[i] = ang_velocity[i];
	}
	
	/*********************************************************/
	/*                   Corrector                           */
	/*********************************************************/
	setContactForceToParticle();
    // build new Resistance matrix after move
    stokes_solver.prepareNewBuild_RFU("direct");
    addStokesDrag();
    buildLubricationTerms(false); // false: don't modify rhs, as we want to keep same Brownian force
    stokes_solver.complete_RFU();
    // get the intermediate brownian velocity
	stokes_solver.solve_CholTrans( v_Brownian_mid );
	if (flubcont_update) {  // rebuild rhs
		stokes_solver.resetRHS();
		buildLubricationRHS();
	} else {  // don't rebuild rhs
		stokes_solver.setRHS(lub_cont_forces_init);
	}
	stokes_solver.solve(v_lub_cont_mid);
    stokes_solver.solvingIsDone();
    // update total velocity
    // first term is hydrodynamic + contact velocities
    // second term is Brownian velocities
    // third term is Brownian drift
    // fourth term for vx is the shear rate
    for (int i = 0; i<_np; i++) {
		int i3 = 3*i;
		velocity[i].set(0.5*(v_lub_cont_mid[i]+v_Brownian_mid[i3]+position[i].z),
						0.5*(v_lub_cont_mid[i3+1]+v_Brownian_mid[i3+1]*zero_2Dsimu),
						0.5*(v_lub_cont_mid[i3+2]+v_Brownian_mid[i3+2]));
		velocity[i] -= 0.5*velocity_predictor[i];
    }
	if (friction) {
		//		double O_inf_y = 0.5;
		//		for (int i=0; i < _np; i++){
		//			ang_velocity[i] = 0.75*contact_torque[i]/radius_cubic[i];
		//			ang_velocity[i].y += 0.5*O_inf_y;
		//			ang_velocity[i] -= 0.5*ang_velocity_predictor[i];
		//		}
		//		//	ang_velocity[i] = 0.5*(ang_velocity_mp1st[i] + ang_velocity_mp2nd[i]);
	}
	for (int i=0; i < _np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (dimension == 2) {
		for (int i=0; i < _np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	// update boxing system
	boxset.update();
	checkNewInteraction();
	updateInteractions();
}

void
System::timeEvolution(double strain_interval){
	double shear_strain_next = shear_strain+strain_interval-1e-6;
	if (shear_strain == 0) {
		checkNewInteraction();
	}
	do {
		switch (integration_method) {
			case 0:
				timeEvolutionEulersMethod();
				break;
			case 1:
				timeEvolutionPredictorCorrectorMethod();
				break;
			case 2:
				timeEvolutionBrownian();
				break;
		}
		ts++;
		shear_strain += dt; //
	} while (shear_strain < shear_strain_next);
}

void
System::timeEvolutionRelax(int time_step){
	int ts_next = ts+time_step;
	checkNewInteraction();
	while (ts < ts_next) {
		setContactForceToParticle();
		setColloidalForceToParticle();
		in_predictor = true;
		in_corrector = true;
		updateVelocityRestingFluid();
		deltaTimeEvolution();
		ts++;
	}
}

void
System::checkNewInteraction(){
	vec3d pos_diff;
	int zshift;
	double sq_dist;
	for (int i=0; i<_np-1; i++) {
		vector<int>::iterator it_beg = boxset.neighborhood_begin(i);
		vector<int>::iterator it_end = boxset.neighborhood_end(i);
		for (vector<int>::iterator it = it_beg; it != it_end; it++) {
			if (*it > i) {
				if (interaction_partners[i].find(*it) == interaction_partners[i].end()) {
					// distance is done in 3 steps
					// because we need each information for Interaction creation
					pos_diff = position[*it]-position[i];
					periodize_diff(pos_diff, zshift);
					sq_dist = pos_diff.sq_norm();
					double ro_2 = 0.5*(radius[i]+radius[*it]);
					double sq_dist_lim = sq_lub_max*ro_2*ro_2;
					if (sq_dist < sq_dist_lim) {
						int interaction_new;
						if (deactivated_interaction.empty()) {
							// add an interaction object.
							interaction_new = num_interaction;
							num_interaction ++;
						} else {
							// fill a deactivated interaction object.
							interaction_new = deactivated_interaction.front();
							deactivated_interaction.pop();
						}
						// new interaction
						interaction[interaction_new].activate(i, *it);
					}
				}
			}
		}
	}
}

/* Check the distance between separating particles.
 * i < j
 *
 * A patch-up prescription to avoid
 * contact_pair[i][j] < 0 indicates separating particles to be checked.
 * contact_pair[i][j] = -1, the particles are near contact. So every time step, distance should be checked.a
 * contact_pair[i][j] < -1, the particles have some distance.
 */
void
System::updateInteractions(){
	/* default value of `_in_predictor' is false
	 *
	 */
	bool deactivated;
	for (int k=0; k<num_interaction; k++) {
		interaction[k].updateState(deactivated);
		if (deactivated) {
			deactivated_interaction.push(k);
		}
	}
}

void
System::stressReset(){
	for (int i=0; i<_np; i++) {
		lubstress[i].reset();
		bgfstress[i].reset();
		contactstressGU[i].reset();
		colloidalstressGU[i].reset();
	}
}

void
System::stressBrownianReset(){
	for (int i=0; i<_np; i++) {
		brownianstress[i].reset();
	}
	brownianstress_calc_nb = 0;
}

void
System::addStokesDrag(){
    for (int i=0; i<_np; i++) {
		stokes_solver.addToDiag_RFU(i, bgf_factor*radius[i]);
    }
}

/* We solve A*(U-Uinf) = Gtilde*Einf ( in Jeffrey's notations )
 * This method computes:
 *  - elements of matrix A
 *  - vector Gtilde*Einf if rhs is true (default behavior)
 */
void
System::buildLubricationTerms(bool rhs){
    double GEi[3];
    double GEj[3];
	/* interaction_list[i] includes all partners j (j > i and j < i).
	 * This range i < _np - 1 is ok?
	 */
    for (int i=0; i<_np-1; i ++) {
		int i3 = 3*i;
		for (set<Interaction*>::iterator it = interaction_list[i].begin();
			 it != interaction_list[i].end(); it ++) {
			int j = (*it)->partner(i);
			if (j > i) {
				(*it)->calcXA();
				stokes_solver.addToDiagBlock_RFU((*it)->nr_vec, i, (*it)->a0*(*it)->XA[0]);
				stokes_solver.addToDiagBlock_RFU((*it)->nr_vec, j, (*it)->a1*(*it)->XA[3]);
				stokes_solver.appendToOffDiagBlock_RFU((*it)->nr_vec, i, j,
													   0.5*(*it)->ro*(*it)->XA[2]);
				if (rhs) {
					int j3 = 3*j;
					(*it)->GE(GEi, GEj);  // G*E_\infty term
					for (int u=0; u<3; u++) {
						stokes_solver.addToRHS(i3+u, GEi[u]);
						stokes_solver.addToRHS(j3+u, GEj[u]);
					}
				}
			}
		}
		stokes_solver.doneBlocks(i);
    }
}

void
System::buildLubricationRHS(){
	double GEi[3];
    double GEj[3];
    for (int i=0; i<_np-1; i ++) {
		int i3 = 3*i;
		for (set<Interaction*>::iterator it = interaction_list[i].begin();
			 it != interaction_list[i].end(); it ++) {
			int j = (*it)->partner(i);
			if (j > i) {
				int j3 = 3*j;
				(*it)->GE(GEi, GEj);  // G*E_\infty term
				for (int u=0; u<3; u++) {
					stokes_solver.addToRHS(i3+u, GEi[u]);
					stokes_solver.addToRHS(j3+u, GEj[u]);
				}
			}
		}
    }
}

void
System::setContactForceToParticle(){
	for (int i=0; i<_np; i++) {
		contact_force[i].reset();
		contact_torque[i].reset();
	}
	for (int k=0; k<num_interaction; k++) {
		interaction[k].addUpContactForceTorque();
	}
}

void
System::setColloidalForceToParticle(){
	if (colloidalforce) {
		for (int i=0; i<_np; i++) {
			colloidal_force[i].reset();
		}
		for (int k=0; k<num_interaction; k++) {
			if (interaction[k].active) {
				interaction[k].addUpColloidalForce();
			}
		}
	}
}

void
System::buildContactTerms(){
    // add contact force
    for (int i=0; i<_np; i++) {
		int i3 = 3*i;
		stokes_solver.addToRHS(i3  , contact_force[i].x);
		stokes_solver.addToRHS(i3+1, contact_force[i].y);
		stokes_solver.addToRHS(i3+2, contact_force[i].z);
    }
}

void
System::buildColloidalForceTerms(){
	if (colloidalforce) {
		for (int i=0; i<_np; i++) {
			int i3 = 3*i;
			stokes_solver.addToRHS(i3  , colloidal_force[i].x);
			stokes_solver.addToRHS(i3+1, colloidal_force[i].y);
			stokes_solver.addToRHS(i3+2, colloidal_force[i].z);
		}
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
    stokes_solver.resetRHS();
    stokes_solver.prepareNewBuild_RFU("direct");
	//	stokes_solver->prepareNewBuild_RFU("iterative");
    addStokesDrag();
	buildLubricationTerms();
    stokes_solver.complete_RFU();
    buildContactTerms();
	buildColloidalForceTerms();
    stokes_solver.solve(v_lub_cont);
	//stokes_solver->print_RFU();
	/* TEST IMPLEMENTATION
	 * SDFF : Stokes drag force factor:
	 * SDFF = 1.0 : full drag forces from the undisturbed background flow.
	 * SDFF = 0.0 : no drag force from the undisturbed background flow.
	 */
    for (int i=0; i<_np; i++) {
		int i3 = 3*i;
		velocity[i].x = v_lub_cont[i3];
		velocity[i].y = v_lub_cont[i3+1];
		velocity[i].z = v_lub_cont[i3+2];
    }
	for (int i=0; i<_np; i++) {
		ang_velocity[i] = 0.75*contact_torque[i]/radius_cubic[i];
	}
	if (dimensionless_shear_rate != 0) {
		for (int i=0; i<_np; i++) {
			velocity[i].x += position[i].z;
		}
		double O_inf_y = 0.5;
		for (int i=0; i<_np; i++) {
			ang_velocity[i].y += O_inf_y;
		}
	}
    stokes_solver.solvingIsDone();
}

void
System::updateVelocityRestingFluid(){
    for (int i=0; i<_np; i++) {
		velocity[i] = contact_force[i]+colloidal_force[i];
    }
}


void
System::displacement(int i, const vec3d &dr){
	position[i] += dr;
	int z_shift = periodize(position[i]);
	/* Note:
	 * When the position of the particle is periodized,
	 * we need to modify the velocity, which was already evaluated.
	 * The position and velocity will be used to calculate the contact forces.
	 */
	if (z_shift) {
		velocity[i].x += z_shift*lz();
	}
	boxset.box(i);
}

// [0,l]
int
System::periodize(vec3d &pos){
	int z_shift = 0;
	vec3d tmp = pos;
	if (pos.z >= _lz) {
		pos.z -= _lz;
		pos.x -= shear_disp;
		z_shift--;
	} else if (pos.z < 0) {
		pos.z += _lz;
		pos.x += shear_disp;
		z_shift++;
	}
	if (pos.x >= _lx) {
		pos.x -= _lx;
		if (pos.x >= _lx){
			pos.x -= _lx;
		}
	} else if (pos.x < 0) {
		pos.x += _lx;
		if (pos.x < 0){
			pos.x += _lx;
		}
	}
	if (pos.y >= _ly) {
		pos.y -= _ly;
	} else if (pos.y < 0) {
		pos.y += _ly;
	}
	return z_shift;
}

// [-l/2,l/2]
void
System::periodize_diff(vec3d &pos_diff){
	if (abs(pos_diff.z) > _lz_half) {
		if (pos_diff.z > 0) {
			pos_diff.z -= _lz;
			pos_diff.x -= shear_disp;
		} else {
			pos_diff.z += _lz;
			pos_diff.x += shear_disp;
		}
	}
	if (abs(pos_diff.x) > _lx_half) {
		if (pos_diff.x > 0) {
			pos_diff.x -= _lx;
			if (pos_diff.x > _lx_half) {
				pos_diff.x -= _lx;
			}
		} else {
			pos_diff.x += _lx;
			if (pos_diff.x < -_lx_half) {
				pos_diff.x += _lx;
			}
		}
	}
	if (pos_diff.y > _ly_half) {
		pos_diff.y -= _ly;
	} else if (pos_diff.y < -_ly_half) {
		pos_diff.y += _ly;
	}
}

// periodize + give z_shift= number of boundaries crossed in z-direction
void
System::periodize_diff(vec3d &pos_diff, int &zshift){
	/*
	 * The displacement of the second particle along z direction
	 * is zshift * lz;
	 */
	if (abs(pos_diff.z) > _lz_half) {
		if (pos_diff.z > 0) {
			pos_diff.z -= _lz;
			pos_diff.x -= shear_disp;
			zshift = -1;
		} else {
			pos_diff.z += _lz;
			pos_diff.x += shear_disp;
			zshift = +1;
		}
	} else {
		zshift = 0;
	}
	if (abs(pos_diff.x) > _lx_half) {
		if (pos_diff.x > 0) {
			pos_diff.x -= _lx;
			if (pos_diff.x > _lx_half) {
				pos_diff.x -= _lx;
			}
		} else {
			pos_diff.x += _lx;
			if (pos_diff.x < -_lx_half) {
				pos_diff.x += _lx;
			}
		}
	}
	if (pos_diff.y > _ly_half) {
		pos_diff.y -= _ly;
	} else if (pos_diff.y < -_ly_half) {
		pos_diff.y += _ly;
	}
}

/*
 * Distance between particle i and particle j
 */
double
System::distance(int i, int j){
	return sqrt(sq_distance(i, j));
}

/*
 * Square distance between particle i and particle j
 */
double
System::sq_distance(int i, int j){
	vec3d pos_diff = position[j] - position[i];
	periodize_diff(pos_diff);
	if (dimension == 3) {
		return pos_diff.sq_norm();
	} else {
	 	return pos_diff.sq_norm_xz();
	}
}


void
System::evaluateMaxContactVelocity(){
	max_contact_velo_tan = 0;
	max_contact_velo_normal = 0;
	for (int k=0; k<num_interaction; k++) {
		if (interaction[k].active) {
			if (interaction[k].contact) {
				if (interaction[k].getContactVelocity() > max_contact_velo_tan) {
					max_contact_velo_tan = interaction[k].getContactVelocity();
				}
				if (abs(interaction[k].getNormalVelocity()) > max_contact_velo_normal) {
					max_contact_velo_normal = abs(interaction[k].getNormalVelocity());
				}
			}
		}
	}
}

double
System::evaluateMaxVelocity(){
	double _max_velocity = 0;
	for (int i = 0; i < _np; i++) {
		vec3d velocity_deviation = velocity[i];
		velocity_deviation.x -= position[i].z;
		if (velocity_deviation.norm() > _max_velocity) {
			_max_velocity = velocity_deviation.norm();
		}
	}
	return _max_velocity;
}

double
System::evaluateMaxAngVelocity(){
	double _max_ang_velocity = 0;
	for (int i = 0; i < _np; i++) {
		vec3d ang_velocity_deviation = ang_velocity[i];
		ang_velocity_deviation.y -= 0.5;
		if (ang_velocity_deviation.norm() > _max_ang_velocity) {
			_max_ang_velocity = ang_velocity_deviation.norm();
		}
	}
	return _max_ang_velocity;
}

void
System::analyzeState(){
	max_velocity = evaluateMaxVelocity();
	max_ang_velocity = evaluateMaxAngVelocity();
	evaluateMaxContactVelocity();
	contact_nb = 0;
	max_disp_tan = 0;
	min_gap_nondim = lub_max;
	double sum_Fc_normal_norm = 0;
	max_Fc_normal_norm = 0;
	for (int k=0; k<num_interaction; k++) {
		if (interaction[k].active) {
			if (interaction[k].gap_nondim() < min_gap_nondim) {
				min_gap_nondim = interaction[k].gap_nondim();
			}
			if (interaction[k].contact) {
				sum_Fc_normal_norm += interaction[k].getFcNormal();
				contact_nb ++;
				if (interaction[k].getFcNormal() > max_Fc_normal_norm) {
					max_Fc_normal_norm = interaction[k].getFcNormal();
				}
				if (interaction[k].disp_tan_norm() > max_disp_tan) {
					max_disp_tan = interaction[k].disp_tan_norm();
				}
			}
			//			cout << interaction[k].gap_nondim() << ' ' << interaction[k].getNormalVelocity() << endl;
		}
	}
	if (contact_nb > 0) {
		average_Fc_normal_norm = sum_Fc_normal_norm/contact_nb;
	} else {
		average_Fc_normal_norm = 0;
	}
}

void
System::setSystemVolume(){
	if (dimension == 2) {
		system_volume = _lx*_lz*2*radius_max;
	} else {
		system_volume = _lx*_ly*_lz;
	}
}

void
System::openFileInteractionData(){
	string int_daat_filename = "irecord_" + simu_name + ".dat";
	fout_int_data.open(int_daat_filename.c_str());
}

double
System::evaluateMaxOverlap(){
	double _max_overlap = 0;
	for (int k=0; k<num_interaction; k++) {
		if (interaction[k].active &&
			-interaction[k].gap_nondim() > _max_overlap) {
			_max_overlap = -interaction[k].gap_nondim();
		}
	}
	return _max_overlap;
}

double
System::evaluateMaxDispTan(){
	double _max_disp_tan = 0;
	for (int k= 0; k<num_interaction; k++) {
		if (interaction[k].active &&
			interaction[k].disp_tan_norm() > _max_disp_tan) {
			_max_disp_tan = interaction[k].disp_tan_norm();
		}
	}
	return _max_disp_tan;
}

double
averageList(list<double> &_list, bool remove_max_min){
	double sum = 0;
	double max_one = 0;
	double min_one = 999999;
	for(list<double>::iterator j=_list.begin(); j != _list.end(); ++j) {
		if (remove_max_min) {
			if (*j > max_one) {
				max_one = *j;
			}
			if (*j < min_one) {
				min_one = *j;
			}
		}
		sum += *j;
	}
	if (remove_max_min) {
		return (sum-max_one-min_one)/(_list.size()-2);
	} else {
		return sum/_list.size();
	}
}

void
System::adjustContactModelParameters(int nb_average){
	//	double strain_interval_for_average = 5;
	//	int num_average = strain_interval_for_average/strain_interval_output;
	/*
	 * Averaged max Fn, over a strain interval
	 */
	static list<double> max_Fn_list;
	max_overlap = evaluateMaxOverlap();
	if (max_overlap > 0) {
		max_Fn_list.push_back(kn*max_overlap);
	}
	if (max_Fn_list.size() > 10) {
		double ave_max_Fn = averageList(max_Fn_list, true);
		kn = ave_max_Fn/overlap_target;
		lub_coeff_contact = 4*kn*contact_relaxzation_time;
	}
	if (max_Fn_list.size() == nb_average){
		max_Fn_list.pop_front();
	}
	/*
	 * Averaged max Ft, over a strain interval
	 */
	static list<double> max_Ft_list;
	if (mu_static > 0) {
		max_disp_tan = evaluateMaxDispTan();
		if (max_disp_tan > 0){
			max_Ft_list.push_back(kt*max_disp_tan);
		}
		if (max_Ft_list.size() > 10){
			double ave_max_Ft = averageList(max_Ft_list, true);
			kt = ave_max_Ft/disp_tan_target;
		}
		if (max_Ft_list.size() == nb_average){
			max_Ft_list.pop_front();
		}
	}
	cerr << "(kn, kt, lub_coeff_contact) = " << kn << ' ' << kt << ' ' << lub_coeff_contact << endl;
}

void
System::calcTotalPotentialEnergy(){
	total_energy = 0;
	for (int k=0; k<num_interaction; k++) {
		if (interaction[k].active){
			total_energy += interaction[k].getPotentialEnergy();
		}
	}
}







