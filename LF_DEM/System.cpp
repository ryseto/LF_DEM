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
	DELETE(contactstressGU);
	DELETE(colloidalstressGU);
	DELETE(brownianstress);
	DELETE(interaction);
	DELETE(interaction_list);
	DELETE(interaction_partners);
	DELETE(v_total);
	DELETE(v_hydro);
	DELETE(v_cont);
	DELETE(v_colloidal);
	//	if(brownian){
	//		DELETE(fb);
	//		DELETE(v_Brownian_init);
	//		DELETE(v_Brownian_mid);
	//		DELETE(v_lub_cont_mid);
	//		DELETE(lub_cont_forces_init);
	//	}
};

void
System::allocateRessources(){
	radius_cubic = new double [np];
	angle = new double [np];
	velocity = new vec3d [np];
	ang_velocity = new vec3d [np];
	if (integration_method >= 1) {
		ang_velocity_predictor = new vec3d [np];
		velocity_predictor = new vec3d [np];
	}
	contact_force = new vec3d [np];
	contact_torque = new vec3d [np];
	colloidal_force = new vec3d [np];
	lubstress = new StressTensor [np];
	contactstressGU = new StressTensor [np];
	colloidalstressGU = new StressTensor [np];
	brownianstress = new StressTensor [np];
	int maxnb_interactionpair_per_particle = 15;
	maxnb_interactionpair = maxnb_interactionpair_per_particle*np;
	interaction = new Interaction [maxnb_interactionpair];
	interaction_list = new set <Interaction*> [np];
	interaction_partners = new set <int> [np];
	linalg_size = 6*np;
	v_total = new double [linalg_size];
	v_cont = new double [linalg_size];
	v_hydro = new double [linalg_size];
	v_colloidal = new double [linalg_size];
	//	if (brownian) {
	//	    v_Brownian_init = new double [linalg_size];
	//	    v_Brownian_mid = new double [linalg_size];
	//		v_lub_cont_mid = new double [linalg_size];
	//		lub_cont_forces_init = new double [linalg_size];
	//		for(int i=0;i<linalg_size;i++){
	//			v_Brownian_init[i] = 0;
	//			v_Brownian_mid[i] = 0;
	//			lub_cont_forces_init[i] = 0;
	//		}
	//		fb = new BrownianForce(this);
	//	}
	stokes_solver.init(np, false);
}


void
System::setupSystemForGenerateInit(){
	for (int i=0; i < np; i++) {
		radius_cubic[i] = radius[i]*radius[i]*radius[i];
		angle[i] = 0;
	}
	for (int k=0; k<maxnb_interactionpair ; k++) {
		interaction[k].init(this);
		interaction[k].set_label(k);
	}
	for (int i=0; i<np; i++) {
		velocity[i].reset();
	}
	kb_T = 0;
	shear_strain = 0;
	shear_disp = 0;
	nb_interaction = 0;
	sq_lub_max = lub_max*lub_max; // square of lubrication cutoff length.
	contact_relaxzation_time = 1e-3;
	kn = 2000;
	kt = 0;
	friction = false;
	dimensionless_shear_rate = 1;
	colloidalforce = false;
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
System::allocatePositionRadius(){
	position = new vec3d [np];
	radius = new double [np];
}

void
System::setConfiguration(const vector <vec3d> &initial_positions,
						 const vector <double> &radii,
						 double lx_, double ly_, double lz_){
	set_np(initial_positions.size());
	setBoxSize(lx_, ly_, lz_);
	allocatePositionRadius();
	for (int i=0; i<np; i++) {
		position[i] = initial_positions[i];
		radius[i] = radii[i];
	}
	if (ly_ == 0) {
		dimension = 2;
	} else {
		dimension = 3;
	}
}

void
System::setupSystem(){
	//	if (kb_T == 0) {
	//		brownian = false;
	//	} else {
	//		brownian = true;
	//		integration_method = 2; // > force Euler
	//	}
	if (mu_static > 0) {
		friction = true;
	} else {
		friction = false;
	}
	if (colloidalforce_length > 0) {
		/*
		 * The diemnsionless shear rate is defined as follows:
		 * dimensionless_shear_rate = F0/colloidalforce_amplitude
		 * F0 = 6pi*eta*a^2*shear_rate
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
	for (int k=0; k<maxnb_interactionpair ; k++) {
		interaction[k].init(this);
		interaction[k].set_label(k);
	}
	for (int i=0; i<np; i++) {
		radius_cubic[i] = radius[i]*radius[i]*radius[i];
		angle[i] = 0;
		velocity[i].set(position[i].z, 0, 0);
		ang_velocity[i].set(0, 0.5, 0);
	}
	for (int i=0; i<linalg_size; i++) {
		v_total[i] = 0;
		v_cont[i] = 0;
		v_hydro[i] = 0;
		v_colloidal[i] = 0;
	}
	shear_strain = 0;
	shear_disp = 0;
	nb_interaction = 0;
	cnt_static_to_dynamic = 0;

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
	cerr << "1/lub_reduce_parameter = " <<  1/lub_reduce_parameter << endl;
	/* t = beta/kn
	 *  beta = t*kn
	 * lub_coeff_contact = 4*beta = 4*kn*contact_relaxzation_time
	 */
	/* If a contact is in sliding mode,
	 * lubrication and dashpot forces are activated.
	 * `log_lub_coeff_contactlub' is the parameter for lubrication during dynamic friction.
	 *
	 */
	if (lubrication_model == 1) {
		log_lub_coeff_contact_tan_lubrication = 0;
		log_lub_coeff_contact_tan_dashpot = 0;
	} else if (lubrication_model == 2) {
		log_lub_coeff_contact_tan_lubrication = log(1/lub_reduce_parameter);
		log_lub_coeff_contact_tan_dashpot = 6*kt*contact_relaxzation_time_tan;
	}
	log_lub_coeff_contact_tan_total = log_lub_coeff_contact_tan_dashpot+log_lub_coeff_contact_tan_lubrication;
	ratio_dashpot_total = log_lub_coeff_contact_tan_dashpot/log_lub_coeff_contact_tan_total;
	
	if (lubrication_model == 1) {
		ratio_dashpot_total = 0;
	}
	
	cerr << "log_lub_coeff_contact_tan_lubrication = " << log_lub_coeff_contact_tan_total << endl;
	cerr << "log_lub_coeff_contact_tan_dashpot = " << log_lub_coeff_contact_tan_dashpot << endl;
	cerr << "ratio_dashpot_total = " << ratio_dashpot_total << endl;
	
	ts = 0;
	shear_disp = 0;
	vel_difference = lz;
	after_parameter_changed = false;
	/*
	 * dt_mid: the intermediate time step for the mid-point
	 * algortithm. dt/dt_mid = dt_ratio
	 * Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 * ASD code from Brady has dt_ratio=150
	 */
	stokes_solver.initialize();
	// initialize the brownian force after the solver, as it assumes
	// the cholmod_common of the solver is already initialized
	//	if (brownian) {
	//		fb->init();
	//		brownianstress_calc_nb = 0;
	//	}
	dt = dt_max;
	initializeBoxing();
	checkNewInteraction();
	if (dimension == 2) {
		twodimension = true;
		setSystemVolume(2*radius[np-1]);
	} else {
		twodimension = false;
		setSystemVolume();
	}
	cnt_prameter_convergence = 0;
}

void
System::initializeBoxing(){// need to know radii first
	double max_radius = 0;
	for (int i=0; i < np; i++) {
		if (radius[i] > max_radius) {
			max_radius = radius[i];
		}
	}
	boxset.init(lub_max*max_radius, this);
	for (int i=0; i<np; i++) {
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
System::evaluateFrictionalState(){
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()){
			interaction[k].contact.frictionlaw();
		}
	}
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
	/* evolve PBC */
	shear_disp += vel_difference*dt;
	if (shear_disp > lx) {
		shear_disp -= lx;
	}
	/* move particles */
	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (twodimension) {
		for (int i=0; i<np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	/* update boxing system */
	boxset.update();
	checkNewInteraction();
	in_predictor = true;
	in_corrector = true;
	updateInteractions();
}

void
System::deltaTimeEvolutionRelax(){
	/* evolve PBC */
	shear_disp += vel_difference*dt;
	if (shear_disp > lx) {
		shear_disp -= lx;
	}
	/* move particles */
	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (twodimension) {
		for (int i=0; i<np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	/* update boxing system */
	boxset.update();
	checkNewInteraction();
	in_predictor = true;
	in_corrector = true;
	bool deactivated;
	for (int k=0; k<nb_interaction; k++) {
		interaction[k].updateStateRelax(deactivated);
		if (deactivated) {
			deactivated_interaction.push(k);
		}
	}
}

void
System::deltaTimeEvolutionPredictor(){
	/* The periodic boundary condition is updated in predictor.
	 * It must not be updated in corrector.
	 */
	shear_disp += vel_difference*dt;
	if (shear_disp >= lx) {
		shear_disp -= lx;
	}
	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (twodimension) {
		for (int i=0; i<np; i++) {
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
	for (int i=0; i<np; i++) {
		velocity_predictor[i] = velocity[i];
		ang_velocity_predictor[i] = ang_velocity[i];
	}
}

void
System::deltaTimeEvolutionCorrector(){
	for (int i=0; i<np; i++) {
		velocity[i] = 0.5*(velocity[i]-velocity_predictor[i]);
		ang_velocity[i] = 0.5*(ang_velocity[i]-ang_velocity_predictor[i]);
	}
	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (twodimension) {
		for (int i=0; i<np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	/* update boxing system
	 */
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
	for (int i=0; i<np; i++) {
		velocity[i] += velocity_predictor[i];
		ang_velocity[i] += ang_velocity_predictor[i];
	}
}

//void System::timeEvolutionBrownian(){
//	int zero_2Dsimu;
//	if (dimension == 2) {
//		zero_2Dsimu = 0;
//	} else {
//		zero_2Dsimu = 1;
//	}
//	/***************************************************************************
//	 *  This routine implements a predictor-corrector algorithm                 *
//	 *  for the dynamics with Brownian motion.                                  *
//	 *  The		algorithm is the one of Melrose & Ball 1997.                        *
//	 *                                                                          *
//	 *  X_pred = X(t) + V(t)*dt                                                 *
//	 *  X(t+dt) = X_pred + 0.5*(V(t+dt)-V(t))*dt                                *
//	 *                                                                          *
//	 *  Parameters                                                              *
//	 *  They essentially controls what enter V(t) and V(t+dt):                  *
//	 *   * displubcont : - if true  V(t) = V^{B}(t) + V^{H}(t) + V^{C}(t)       *
//	 *                   - if false V(t) = V^{B}(t)                             *
//	 *   * flubcont_update : - requires displubcont                             *
//	 *                       - allows dt^2 scheme for contact and hydro         *
//	 *                         velocities                                       *
//	 *                       - if true :                                        *
//	 *                          R_FU(t+dt) V^{*}(t+dt) = F^{*}(t+dt)            *
//	 *                          with * = C or H                                 *
//	 *                       - if false :                                       *
//	 *                          R_FU(t+dt) V^{*}(t+dt) = F^{*}(t)               *
//	 *                                                                          *
//	 ****************************************************************************/
//	bool displubcont = true;
//	bool flubcont_update = true;
//	/*********************************************************/
//	/*                    Predictor                          */
//	/*********************************************************/
//	setContactForceToParticle();
//	stokes_solver.resetRHS();
//    stokes_solver.resetResistanceMatrix("direct");
//    addStokesDrag();
//    buildLubricationTerms();
//
//    stokes_solver.completeResistanceMatrix();
//    buildContactTerms();
//    stokes_solver.solve(v_lub_cont);
//
//	if (displubcont) {
//		stokes_solver.getRHS(lub_cont_forces_init);
//	}
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
//    // move particles to intermediate point
//    for (int i=0; i<np; i++) {
//		int i3 = 3*i;
//		velocity[i].set(v_Brownian_init[i3],
//						v_Brownian_init[i3+1]*zero_2Dsimu,
//						v_Brownian_init[i3+2]);
//		velocity[i].add(v_lub_cont[i3],
//						v_lub_cont[i3+1],
//						v_lub_cont[i3+2]);
//		velocity[i].x += position[i].z;
//		displacement(i, velocity[i]*dt);
//    }
//	if (friction) {
//		double O_inf_y = 0.5;
//		for (int i=0; i<np; i++) {
//			ang_velocity[i] = 0.75*contact_torque[i]/radius_cubic[i];
//			ang_velocity[i].y += O_inf_y;
//		}
//    }
//	// evolve PBC
//	shear_disp += vel_difference*dt;
//	if (shear_disp >= lx) {
//		shear_disp -= lx;
//	}
//    updateInteractions();
//	for (int i=0; i<np; i++) {
//		velocity_predictor[i] = velocity[i];
//		ang_velocity_predictor[i] = ang_velocity[i];
//	}
//
//	/*********************************************************/
//	/*                   Corrector                           */
//	/*********************************************************/
//	setContactForceToParticle();
//    // build new Resistance matrix after move
//    stokes_solver.resetResistanceMatrix("direct");
//    addStokesDrag();
//    buildLubricationTerms(false); // false: don't modify rhs, as we want to keep same Brownian force
//    stokes_solver.completeResistanceMatrix();
//    // get the intermediate brownian velocity
//	stokes_solver.solve_CholTrans( v_Brownian_mid );
//	if (flubcont_update) {  // rebuild rhs
//		stokes_solver.resetRHS();
//		buildLubricationRHS();
//	} else {  // don't rebuild rhs
//		stokes_solver.setRHS(lub_cont_forces_init);
//	}
//	stokes_solver.solve(v_lub_cont_mid);
//    stokes_solver.solvingIsDone();
//    // update total velocity
//    // first term is hydrodynamic + contact velocities
//    // second term is Brownian velocities
//    // third term is Brownian drift
//    // fourth term for vx is the shear rate
//    for (int i = 0; i<np; i++) {
//		int i3 = 3*i;
//		velocity[i].set(0.5*(v_lub_cont_mid[i]+v_Brownian_mid[i3]+position[i].z),
//						0.5*(v_lub_cont_mid[i3+1]+v_Brownian_mid[i3+1]*zero_2Dsimu),
//						0.5*(v_lub_cont_mid[i3+2]+v_Brownian_mid[i3+2]));
//		velocity[i] -= 0.5*velocity_predictor[i];
//    }
//	if (friction) {
//		//		double O_inf_y = 0.5;
//		//		for (int i=0; i < np; i++){
//		//			ang_velocity[i] = 0.75*contact_torque[i]/radius_cubic[i];
//		//			ang_velocity[i].y += 0.5*O_inf_y;
//		//			ang_velocity[i] -= 0.5*ang_velocity_predictor[i];
//		//		}
//		//		//	ang_velocity[i] = 0.5*(ang_velocity_mp1st[i] + ang_velocity_mp2nd[i]);
//	}
//	for (int i=0; i < np; i++) {
//		displacement(i, velocity[i]*dt);
//	}
//	if (twodimension) {
//		for (int i=0; i < np; i++) {
//			angle[i] += ang_velocity[i].y*dt;
//		}
//	}
//	// update boxing system
//	boxset.update();
//	checkNewInteraction();
//	updateInteractions();
//}

void
System::timeEvolution(double strain_next){
	static bool firsttime = true;
	if (firsttime) {
		checkNewInteraction();
		firsttime = false;
	}
	while (shear_strain < strain_next-1e-8) {
		switch (integration_method) {
			case 0:
				timeEvolutionEulersMethod();
				break;
			case 1:
				timeEvolutionPredictorCorrectorMethod();
				break;
			case 2:
				//				timeEvolutionBrownian();
				break;
		}
		ts++;
		shear_strain += dt; //
	};
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
		deltaTimeEvolutionRelax();
		ts++;
	}
}

void
System::checkNewInteraction(){
	vec3d pos_diff;
	int zshift;
	double sq_dist;
	for (int i=0; i<np-1; i++) {
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
							interaction_new = nb_interaction;
							nb_interaction ++;
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
	 */
	if (friction) {
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_contact()){
				interaction[k].calcRelativeVelocities();
				interaction[k].contact.incrementTangentialDisplacement();
			}
		}
	}
	
	
	for (int k=0; k<nb_interaction; k++) {
		bool deactivated = false;
		if (interaction[k].is_active()){
			interaction[k].updateState(deactivated);
			if (deactivated) {
				deactivated_interaction.push(k);
			}
		}
	}
	//	evaluateFrictionalState();
}

void
System::stressReset(){
	for (int i=0; i<np; i++) {
		lubstress[i].reset();
		contactstressGU[i].reset();
		colloidalstressGU[i].reset();
	}
}

void
System::stressBrownianReset(){
	for (int i=0; i<np; i++) {
		brownianstress[i].reset();
	}
	brownianstress_calc_nb = 0;
}

void
System::addStokesDrag(){
	double torque_factor = 4./3;
    for (int i=0; i<np; i++) {
		stokes_solver.addToDiag(i, bgf_factor*radius[i], bgf_factor*torque_factor*radius[i]*radius[i]*radius[i]);
    }
}

/* We solve A*(U-Uinf) = Gtilde*Einf ( in Jeffrey's notations )
 * This method computes:
 *  - elements of matrix A
 *  - vector Gtilde*Einf if rhs is true (default behavior)
 */
void
System::buildLubricationTerms(bool rhs){
	/* interaction_list[i] includes all partners j (j > i and j < i).
	 * This range i < np - 1 is ok?
	 */
	switch (lubrication_model) {
		case 1:
			for (int i=0; i<np-1; i ++) {
				for (set<Interaction*>::iterator it = interaction_list[i].begin();
					 it != interaction_list[i].end(); it ++) {
					int j = (*it)->partner(i);
					if (j > i) {
						vec3d nr_vec = (*it)->get_nvec();
						(*it)->lubrication.calcXFunctions();						
						stokes_solver.addToDiagBlock(nr_vec, i, (*it)->lubrication.scaledXA0(), 0, 0, 0);
						stokes_solver.addToDiagBlock(nr_vec, j, (*it)->lubrication.scaledXA3(), 0, 0, 0);
						stokes_solver.setOffDiagBlock(nr_vec, i, j, (*it)->lubrication.scaledXA2(), 0, 0, 0, 0);
						if (rhs) {
							double GEi[3];
							double GEj[3];
							(*it)->lubrication.calcGE(GEi, GEj);  // G*E_\infty term
							stokes_solver.addToRHSForce(i, GEi);
							stokes_solver.addToRHSForce(j, GEj);
						}
					}
				}
				stokes_solver.doneBlocks(i);
			}
			break;
		case 2:
			for (int i=0; i<np-1; i ++) {
				for (set<Interaction*>::iterator it = interaction_list[i].begin();
					 it != interaction_list[i].end(); it ++) {
					int j = (*it)->partner(i);
					if (j > i) {
						// if(i==58&&j==97&&(*it)->is_contact()){
						// 	cout << " pause " << endl;
						// 	getchar();
						// }

						vec3d nr_vec = (*it)->get_nvec();
						(*it)->lubrication.calcXYFunctions();
						stokes_solver.addToDiagBlock(nr_vec, i, (*it)->lubrication.scaledXA0(), (*it)->lubrication.scaledYA0(),
													 (*it)->lubrication.scaledYB0(), (*it)->lubrication.scaledYC0());
						stokes_solver.addToDiagBlock(nr_vec, j, (*it)->lubrication.scaledXA3(), (*it)->lubrication.scaledYA3(),
													 (*it)->lubrication.scaledYB3(), (*it)->lubrication.scaledYC3());
						stokes_solver.setOffDiagBlock(nr_vec, i, j, (*it)->lubrication.scaledXA1(), (*it)->lubrication.scaledYA1(),
													  (*it)->lubrication.scaledYB2(), (*it)->lubrication.scaledYB1(), (*it)->lubrication.scaledYC1());
						if (rhs) {
							double GEi[3];
							double GEj[3];
							double HEi[3];
							double HEj[3];
							(*it)->lubrication.calcGEHE(GEi, GEj, HEi, HEj);  // G*E_\infty term
							stokes_solver.addToRHSForce(i, GEi);
							stokes_solver.addToRHSForce(j, GEj);
							stokes_solver.addToRHSTorque(i, HEi);
							stokes_solver.addToRHSTorque(j, HEj);
						}
					}
				}
				stokes_solver.doneBlocks(i);
			}
			break;
		case 3:
			for (int i=0; i<np-1; i ++) {
				for (set<Interaction*>::iterator it = interaction_list[i].begin();
					 it != interaction_list[i].end(); it ++) {
					int j = (*it)->partner(i);
					if (j > i) {
						vec3d nr_vec = (*it)->get_nvec();
						if ((*it)->is_contact()){
							(*it)->lubrication.calcXYFunctions();
							stokes_solver.addToDiagBlock(nr_vec, i, (*it)->lubrication.scaledXA0(), (*it)->lubrication.scaledYA0(),
														 (*it)->lubrication.scaledYB0(), (*it)->lubrication.scaledYC0());
							stokes_solver.addToDiagBlock(nr_vec, j, (*it)->lubrication.scaledXA3(), (*it)->lubrication.scaledYA3(),
														 (*it)->lubrication.scaledYB3(), (*it)->lubrication.scaledYC3());
							stokes_solver.setOffDiagBlock(nr_vec, i, j, (*it)->lubrication.scaledXA1(), (*it)->lubrication.scaledYA1(),
														  (*it)->lubrication.scaledYB2(), (*it)->lubrication.scaledYB1(), (*it)->lubrication.scaledYC1());
							if (rhs) {
								double GEi[3];
								double GEj[3];
								double HEi[3];
								double HEj[3];
								(*it)->lubrication.calcGEHE(GEi, GEj, HEi, HEj);  // G*E_\infty term
								stokes_solver.addToRHSForce(i, GEi);
								stokes_solver.addToRHSForce(j, GEj);
								stokes_solver.addToRHSTorque(i, HEi);
								stokes_solver.addToRHSTorque(j, HEj);
							}
						} else {
							(*it)->lubrication.calcXFunctions();
							stokes_solver.addToDiagBlock(nr_vec, i, (*it)->lubrication.scaledXA0(), 0, 0, 0);
							stokes_solver.addToDiagBlock(nr_vec, j, (*it)->lubrication.scaledXA3(), 0, 0, 0);
							stokes_solver.setOffDiagBlock(nr_vec, i, j, (*it)->lubrication.scaledXA2(), 0, 0, 0, 0);
							if (rhs) {
								double GEi[3];
								double GEj[3];
								(*it)->lubrication.calcGE(GEi, GEj);  // G*E_\infty term
								stokes_solver.addToRHSForce(i, GEi);
								stokes_solver.addToRHSForce(j, GEj);
							}
						}
					}
				}
				stokes_solver.doneBlocks(i);
			}
			break;
		default:
			cerr << "lubrication_model = 0 is not implemented yet.\n";
			exit(1);
			break;
	}
}

//void
//System::buildLubricationRHS(){
//	double GEi[3];
//    double GEj[3];
//    for (int i=0; i<np-1; i ++) {
//		for (set<Interaction*>::iterator it = interaction_list[i].begin();
//			 it != interaction_list[i].end(); it ++) {
//			int j = (*it)->partner(i);
//			if (j > i) {
//				(*it)->GE(GEi, GEj);  // G*E_\infty term
//				stokes_solver.addToRHSForce(i, GEi);
//				stokes_solver.addToRHSForce(j, GEj);
//			}
//		}
//    }
//}

void
System::setContactForceToParticle(){
	for (int i=0; i<np; i++) {
		contact_force[i].reset();
		contact_torque[i].reset();
	}
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			interaction[k].contact.addUpContactForceTorque();
		}
	}
}

void
System::setColloidalForceToParticle(){
	if (colloidalforce) {
		for (int i=0; i<np; i++) {
			colloidal_force[i].reset();
		}
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_active()) {
				interaction[k].addUpColloidalForce();
			}
		}
	}
}

void
System::buildContactTerms(){
    // add contact force
    for (int i=0; i<np; i++) {
		stokes_solver.addToRHSForce(i, contact_force[i]);
		stokes_solver.addToRHSTorque(i, contact_torque[i]);
    }
}

void
System::buildColloidalForceTerms(){
	if (colloidalforce) {
		for (int i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, colloidal_force[i]);
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
	nb_of_active_interactions = nb_interaction-deactivated_interaction.size();
    stokes_solver.resetResistanceMatrix("direct", nb_of_active_interactions);
	//	stokes_solver->resetResistanceMatrix("iterative");
    addStokesDrag();
	buildLubricationTerms();
    stokes_solver.completeResistanceMatrix();
	buildContactTerms();
	buildColloidalForceTerms();
    stokes_solver.solve(v_total);
	stokes_solver.solvingIsDone();
	//stokes_solver->printResistanceMatrix();
	/* TEST IMPLEMENTATION
	 * SDFF : Stokes drag force factor:
	 * SDFF = 1.0 : full drag forces from the undisturbed background flow.
	 * SDFF = 0.0 : no drag force from the undisturbed background flow.
	 */
    for (int i=0; i<np; i++) {
		int i6 = 6*i;
		velocity[i].x = v_total[i6];
		velocity[i].y = v_total[i6+1];
		velocity[i].z = v_total[i6+2];
		ang_velocity[i].x = v_total[i6+3];
		ang_velocity[i].y = v_total[i6+4];
		ang_velocity[i].z = v_total[i6+5];
		//cout << " in System : particle  " << i << " has a velocity " << velocity[i].x << " " << velocity[i].y << " "  << velocity[i].z << " " << ang_velocity[i].x << " " << ang_velocity[i].y << " "  << ang_velocity[i].z << endl;
    }
	if (dimensionless_shear_rate != 0) {
		for (int i=0; i<np; i++) {
			velocity[i].x += position[i].z;
			ang_velocity[i].y += 0.5;
		}
	}
}

void
System::updateVelocityRestingFluid(){
    for (int i=0; i<np; i++) {
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
		velocity[i].x += z_shift*lz;
	}
	boxset.box(i);
}

// [0,l]
int
System::periodize(vec3d &pos){
	int z_shift = 0;
	vec3d tmp = pos;
	if (pos.z >= lz) {
		pos.z -= lz;
		pos.x -= shear_disp;
		z_shift--;
	} else if (pos.z < 0) {
		pos.z += lz;
		pos.x += shear_disp;
		z_shift++;
	}
	if (pos.x >= lx) {
		pos.x -= lx;
		if (pos.x >= lx){
			pos.x -= lx;
		}
	} else if (pos.x < 0) {
		pos.x += lx;
		if (pos.x < 0){
			pos.x += lx;
		}
	}
	if (pos.y >= ly) {
		pos.y -= ly;
	} else if (pos.y < 0) {
		pos.y += ly;
	}
	return z_shift;
}

// [-l/2,l/2]
void
System::periodize_diff(vec3d &pos_diff){
	if (abs(pos_diff.z) > lz_half) {
		if (pos_diff.z > 0) {
			pos_diff.z -= lz;
			pos_diff.x -= shear_disp;
		} else {
			pos_diff.z += lz;
			pos_diff.x += shear_disp;
		}
	}
	if (abs(pos_diff.x) > lx_half) {
		if (pos_diff.x > 0) {
			pos_diff.x -= lx;
			if (pos_diff.x > lx_half) {
				pos_diff.x -= lx;
			}
		} else {
			pos_diff.x += lx;
			if (pos_diff.x < -lx_half) {
				pos_diff.x += lx;
			}
		}
	}
	if (pos_diff.y > ly_half) {
		pos_diff.y -= ly;
	} else if (pos_diff.y < -ly_half) {
		pos_diff.y += ly;
	}
}

// periodize + give z_shift= number of boundaries crossed in z-direction
void
System::periodize_diff(vec3d &pos_diff, int &zshift){
	/*
	 * The displacement of the second particle along z direction
	 * is zshift * lz;
	 */
	if (abs(pos_diff.z) > lz_half) {
		if (pos_diff.z > 0) {
			pos_diff.z -= lz;
			pos_diff.x -= shear_disp;
			zshift = -1;
		} else {
			pos_diff.z += lz;
			pos_diff.x += shear_disp;
			zshift = +1;
		}
	} else {
		zshift = 0;
	}
	if (abs(pos_diff.x) > lx_half) {
		if (pos_diff.x > 0) {
			pos_diff.x -= lx;
			if (pos_diff.x > lx_half) {
				pos_diff.x -= lx;
			}
		} else {
			pos_diff.x += lx;
			if (pos_diff.x < -lx_half) {
				pos_diff.x += lx;
			}
		}
	}
	if (pos_diff.y > ly_half) {
		pos_diff.y -= ly;
	} else if (pos_diff.y < -ly_half) {
		pos_diff.y += ly;
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
	in_predictor = true;
	double sum_contact_velo_tan = 0;
	double sum_contact_velo_normal = 0;
	double sum_sliding_velocity = 0;
	int cnt_contact = 0;
	int cnt_sliding = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()) {
			interaction[k].calcRelativeVelocities();
			cnt_contact++;
			sum_contact_velo_tan += interaction[k].getContactVelocity();
			sum_contact_velo_normal += abs(interaction[k].getNormalVelocity());
			if (interaction[k].contact.staticfriction == false) {
				cnt_sliding++;
				sum_sliding_velocity += interaction[k].getContactVelocity();
			}
			if (interaction[k].getContactVelocity() > max_contact_velo_tan) {
				max_contact_velo_tan = interaction[k].getContactVelocity();
			}
			if (abs(interaction[k].getNormalVelocity()) > max_contact_velo_normal) {
				max_contact_velo_normal = abs(interaction[k].getNormalVelocity());
			}
			
		}
	}
	ave_contact_velo_tan = sum_contact_velo_tan/cnt_contact;
	ave_contact_velo_normal = sum_contact_velo_normal/cnt_contact;
	ave_sliding_velocity = sum_sliding_velocity/cnt_sliding;
}

double
System::evaluateMaxVelocity(){
	double sq_max_velocity = 0;
	for (int i = 0; i < np; i++) {
		vec3d velocity_deviation = velocity[i];
		velocity_deviation.x -= position[i].z;
		if (velocity_deviation.sq_norm() > sq_max_velocity) {
			sq_max_velocity = velocity_deviation.sq_norm();
		}
	}
	return sqrt(sq_max_velocity);
}

double
System::evaluateMaxAngVelocity(){
	double _max_ang_velocity = 0;
	for (int i = 0; i < np; i++) {
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
	static double previous_strain = 0;
	double strain_interval = shear_strain-previous_strain;
	previous_strain = shear_strain;
	
	max_velocity = evaluateMaxVelocity();

	max_ang_velocity = evaluateMaxAngVelocity();
	evaluateMaxContactVelocity();
	sliding_velocity_history.push_back(max_contact_velo_tan);
	contact_nb = 0;
	max_disp_tan = 0;
	min_gap_nondim = lx;
	double sum_fc_normal = 0;
	max_fc_normal = 0;
	max_fc_tan = 0;
	intr_max_fc_normal = -1;
	intr_max_fc_tan = -1;
	int cnt_sliding_contact=0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			if (interaction[k].get_gap_nondim() < min_gap_nondim) {
				min_gap_nondim = interaction[k].get_gap_nondim();
			}
			if (interaction[k].is_contact()) {
				contact_nb ++;
				if (interaction[k].contact.staticfriction == false){
					cnt_sliding_contact++;
				}
				sum_fc_normal += interaction[k].contact.get_f_contact_normal_norm();
				if (interaction[k].contact.get_f_contact_normal_norm() > max_fc_normal) {
					max_fc_normal = interaction[k].contact.get_f_contact_normal_norm();
					intr_max_fc_normal = k;
				}
				if (interaction[k].contact.get_f_contact_tan_norm() > max_fc_tan) {
					max_fc_tan = interaction[k].contact.get_f_contact_tan_norm();
					intr_max_fc_tan = k;
				}
				if (interaction[k].contact.disp_tan_norm() > max_disp_tan) {
					max_disp_tan = interaction[k].contact.disp_tan_norm();
				}
			}
		}
	}
	/*
	 * History is recorded after the relaxation.
	 *
	 */
	static double max_fc_normal_previous = 0;
	if (after_parameter_changed) {
		if (max_fc_normal > max_fc_normal_previous){
			after_parameter_changed = false;
		}
	}
	max_fc_normal_previous = max_fc_normal;
	if (after_parameter_changed == false) {
		max_fc_normal_history.push_back(max_fc_normal);
		max_fc_tan_history.push_back(max_fc_tan);
	}
	if (contact_nb > 0) {
		average_fc_normal = sum_fc_normal/contact_nb;
	} else {
		average_fc_normal = 0;
	}
	rate_static_to_dynamic = cnt_static_to_dynamic/(strain_interval*np);
	ratio_dynamic_friction = (contact_nb-cnt_sliding_contact)*(1./contact_nb);
	cnt_static_to_dynamic = 0;
}

void
System::setSystemVolume(double depth){
	if (dimension == 2) {
		system_volume = lx*lz*depth;
	} else {
		system_volume = lx*ly*lz;
	}
}

void
System::openFileInteractionData(){
	string int_data_filename = "irecord_" + simu_name + ".dat";
	fout_int_data.open(int_data_filename.c_str());
}

double
System::evaluateMaxDispTan(){
	double _max_disp_tan = 0;
	for (int k= 0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].contact.disp_tan_norm() > _max_disp_tan) {
			_max_disp_tan = interaction[k].contact.disp_tan_norm();
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

void calcMean_StdDev(vector<double> history,
					 double &mean,
					 double &std_dev){
	int ne = history.size();
	double sum = 0;
	for (int i=0; i<ne; i++) {
		sum += history[i];
	}
	mean = sum/ne;
	double sum_sq_deviation = 0;
	for (int i=0; i<ne; i++) {
		double tmp = history[i]-mean;
		sum_sq_deviation += tmp*tmp;
	}
	std_dev = sqrt(sum_sq_deviation/(ne-1));
}


int
System::adjustContactModelParameters(){
	/*
	 * kn, kt and dt are determined in one test simulation.
	 * We give small values of kn and kt as the initial values.
	 * With target values of spring strech,
	 * spring constants are determined by the maximum foces in a certain interaval.
	 * In order to avoid unusual large values of forces,
	 * the maximum force in the interaval is defined by mean + std_dev of the maximum values.
	 * Only increases of kn and kt are accepted.
	 */
	double kn_previous = kn;
	double kt_previous = kt;
	double max_increment = 1000;
	/* determination of kn
	 */
	double mean_max_fc_normal, stddev_max_fc_normal;
	calcMean_StdDev(max_fc_normal_history, mean_max_fc_normal, stddev_max_fc_normal);
	double kn_try = mean_max_fc_normal/overlap_target;
	if (kn_try > kn) {
		if (kn_try > kn*max_increment){
			kn_try = kn*max_increment;
		}
		kn = kn_try;
		lub_coeff_contact = 4*kn*contact_relaxzation_time;
	}
	/* determination of kt
	 */
	double mean_max_fc_tan, stddev_max_fc_tan;
	calcMean_StdDev(max_fc_tan_history, mean_max_fc_tan, stddev_max_fc_tan);
	double kt_try = mean_max_fc_tan/disp_tan_target;
	if (kt_try > kt){
		if (kt_try > kt*max_increment){
			kt_try = kt*max_increment;
		}
		kt = kt_try;
	}
	double average_max_tanvelocity = 0;
	for (int j=0; j<sliding_velocity_history.size(); j++){
		average_max_tanvelocity += sliding_velocity_history[j];
	}
	average_max_tanvelocity = average_max_tanvelocity/sliding_velocity_history.size();
	cerr << "average_max_velocity = " << average_max_tanvelocity << endl;
	cerr << "disp_max " << disp_max << endl;
	double dt_try = disp_max/average_max_tanvelocity;
	if (dt_try < dt){
		dt = dt_try;
	}
	cerr << "dt_try = " << dt_try << endl;
	
	for (int k=0; k<nb_interaction; k++) {
		interaction[k].contact.updateContactModel();
	}
	max_fc_normal_history.clear();
	max_fc_tan_history.clear();
	sliding_velocity_history.clear();
	after_parameter_changed = true;
	
	if (kn_previous == kn && kt_previous == kt ){
		cnt_prameter_convergence ++;
	} else {
		cnt_prameter_convergence = 0;
	}
	if (cnt_prameter_convergence == 5 || kn > max_kn){
		return 1;
	}
	return 0;
}

void
System::calcTotalPotentialEnergy(){
	total_energy = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()){
			total_energy += interaction[k].getPotentialEnergy();
		}
	}
}

void
System::calcLubricationForce(){
	/*
	 * Calculate lubrication force to output
	 */
	nb_of_active_interactions = nb_interaction-deactivated_interaction.size();
	stokes_solver.resetResistanceMatrix("direct", nb_of_active_interactions);
    addStokesDrag();
	stokes_solver.resetRHS();
    buildLubricationTerms();
    setContactForceToParticle();
	buildContactTerms();
	setColloidalForceToParticle();
	buildColloidalForceTerms();
	stokes_solver.completeResistanceMatrix();
	stokes_solver.solve(v_total);
	stokes_solver.solvingIsDone();
	
	if (lubrication_model == 1){
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_active()) {
				interaction[k].lubrication.calcXFunctions();
			}
		}
	} else if (lubrication_model == 2){
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_active()) {
				interaction[k].lubrication.calcXYFunctions();
			}
		}
	} else if (lubrication_model == 3){
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_active()) {
				if (interaction[k].is_contact()) {
					interaction[k].lubrication.calcXYFunctions();
				} else {
					interaction[k].lubrication.calcXFunctions();
				}
			}
		}
	}
	
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			interaction[k].lubrication.calcLubricationForce();
		}
	}
}


