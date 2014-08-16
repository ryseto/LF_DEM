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
#define GRANDOM ( r_gen->randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.

System::System() :
maxnb_interactionpair_per_particle(15),
brownian(false),
zero_shear(false),
friction_model(-1),
repulsiveforce(false)
{}

System::~System(){
	DELETE(position);
	DELETE(radius);
	DELETE(radius_cubed);
	DELETE(resistance_matrix_dblock);
	if (twodimension) {
		DELETE(angle);
	}
	DELETE(velocity);
	DELETE(ang_velocity);
	DELETE(na_velocity);
	DELETE(na_ang_velocity);
	if (integration_method == 1) {
		DELETE(velocity_predictor);
		DELETE(ang_velocity_predictor);
	}
	DELETE(vel_contact);
	DELETE(ang_vel_contact);
	DELETE(vel_hydro);
	DELETE(ang_vel_hydro);
	DELETE(contact_force);
	DELETE(contact_torque);
	DELETE(lubstress);
	DELETE(contactstressGU);
	DELETE(avg_lubstress);
	DELETE(avg_contactstressGU);
	DELETE(interaction);
	DELETE(interaction_list);
	DELETE(interaction_partners);
	if (brownian) {
		DELETE(vel_brownian);
		DELETE(ang_vel_brownian);
		DELETE(brownian_force);
		DELETE(brownianstressGU);
		DELETE(brownianstressGU_predictor);
		DELETE(avg_brownianstressGU);
	}
	if (repulsiveforce) {
		DELETE(repulsivestressGU);
		DELETE(vel_repulsive);
		DELETE(ang_vel_repulsive);
		DELETE(avg_repulsivestressGU);
	}
};

void
System::allocateRessources(){
	linalg_size = 6*np;
	maxnb_interactionpair = maxnb_interactionpair_per_particle*np;
	radius_cubed = new double [np];
	resistance_matrix_dblock = new double [18*np];
	// Configuration
	if (twodimension) {
		angle = new double [np];
	}
	// Velocity
	velocity = new vec3d [np];
	ang_velocity = new vec3d [np];
	na_velocity = new vec3d [np];
	na_ang_velocity = new vec3d [np];
	if (integration_method == 1) {
		velocity_predictor = new vec3d [np];
		ang_velocity_predictor = new vec3d [np];
	}
	vel_contact = new vec3d [np];
	ang_vel_contact = new vec3d [np];
	vel_hydro = new vec3d [np];
	ang_vel_hydro = new vec3d [np];
	if (repulsiveforce) {
		vel_repulsive = new vec3d [np];
		ang_vel_repulsive = new vec3d [np];
	}
	if (brownian) {
		vel_brownian = new vec3d [np];
		ang_vel_brownian = new vec3d [np];
	}
	// Forces
	contact_force = new vec3d [np];
	contact_torque = new vec3d [np];
	if (repulsiveforce) {
		repulsive_force = new vec3d [np];
		repulsivestressGU = new StressTensor [np];
	}
	if (brownian) {
		brownian_force = new double [linalg_size];
	}
	// Stress
	lubstress = new StressTensor [np];
	contactstressGU = new StressTensor [np];
	avg_lubstress = new StressTensor [np];
	avg_contactstressGU = new StressTensor [np];
	if (repulsiveforce) {
		avg_repulsivestressGU = new StressTensor [np];
	}
	if (brownian) {
		brownianstressGU = new StressTensor [np];
		brownianstressGU_predictor = new StressTensor [np];
		avg_brownianstressGU = new StressTensor [np];
	}
	interaction = new Interaction [maxnb_interactionpair];
	interaction_list = new set <Interaction*> [np];
	interaction_partners = new set <int> [np];
	stokes_solver.init(np);
}

void
System::setInteractions_GenerateInitConfig(){
	for (int k=0; k<maxnb_interactionpair; k++) {
		interaction[k].init(this);
		interaction[k].set_label(k);
	}
	nb_interaction = 0;
	sq_lub_max = lub_max*lub_max; // square of lubrication cutoff length.
	shear_strain = 0;
	shear_disp = 0;
	vel_difference = 0;
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
		twodimension = true;
	} else {
		twodimension = false;
	}
}

void
System::setupBrownian(){
	if (brownian) {
		/* kb_T is dimensionless.
		 * dimensionless_shear_rate is Peclet number
		 * Dimensional kb_T is L0*F0 / Pe,
		 * where F0 = 6*pi*eta0*a^2*rate and L0 = a
		 *
		 */
		kb_T = 1/dimensionless_shear_rate;
		if (dimensionless_shear_rate < Pe_switch) {
			// scale_factor_SmallPe > 1
			scale_factor_SmallPe = Pe_switch/dimensionless_shear_rate;
			kn = scale_factor_SmallPe*kn_lowPeclet;
			kt = scale_factor_SmallPe*kt_lowPeclet;
			dt_max = dt_lowPeclet/scale_factor_SmallPe;
			contact_relaxation_time = contact_relaxation_time/scale_factor_SmallPe;
			contact_relaxation_time_tan = contact_relaxation_time_tan/scale_factor_SmallPe; // should be zero.
			shear_strain_end /= scale_factor_SmallPe;
			strain_interval_output_data *= 1/scale_factor_SmallPe;
			strain_interval_output *= 1/scale_factor_SmallPe;
			cerr << "[[small Pe mode]]" << endl;
			cerr << "  kn = " << kn << endl;
			cerr << "  kt = " << kt << endl;
			cerr << "  dt_max = " << dt_max << endl;
			cerr << "  strain_interval_output_data = " << strain_interval_output_data << endl;
			cerr << "  strain_interval_output = " << strain_interval_output << endl;
		}
	}
}

void
System::setupSystem(){
	/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	 * @ We have to consider contact_relaxation_time in Brownian case.
	 * @ The resistance coeffient affects Brownian force.
	 * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	 */
	/* Giving a seed for debugging (Brownian)
	 * r_gen = new MTRand(71);
	 */
	if (integration_method == 0) {
		timeEvolutionDt = &System::timeEvolutionEulersMethod;
	} else if (integration_method == 1) {
		timeEvolutionDt = &System::timeEvolutionPredictorCorrectorMethod;
	} else {
		cerr << "integration_method = " << integration_method << endl;
		cerr << "The integration method is not impremented yet." << endl;
		exit(1);
	}
	if (lubrication_model == 1) {
		buildLubricationTerms = &System::buildLubricationTerms_squeeze;
	} else if (lubrication_model == 2) {
		buildLubricationTerms = &System::buildLubricationTerms_squeeze_tangential;
	} else {
		cerr << "lubrication_model = 0 is not implemented yet.\n";
		exit(1);
	}
	if (friction_model == 0) {
		cerr << "friction_model = 0" << endl;
		mu_static = 0;
		friction = false;
	} else if (friction_model == 1) {
		cerr << "friction_model = 1" << endl;
		friction = true;
	} else if (friction_model == 2 || friction_model == 3) {
		cerr << "friction_model " << friction_model << endl;
		/*
		 * The dimensionless shear rate is defined as follows:
		 * dimensionless_shear_rate = F0/F^{*}
		 * F0 = 6pi*eta*a^2*shear_rate
		 * The force unit is changed by the considering shear rate.
		 * In the simulation, the critical force F^{*} = \tilde{F^{*}} F_0.
		 * Thus, \tilde{F^{*}} = F^{*} / F_0 = 1/dimensionless_shear_rate.
		 *
		 */
		friction = true;
		cerr << "critical_normal_force = " << critical_normal_force << endl;
	} else {
		cerr << "friction_model..." << endl;
		exit(1);
	}
	allocateRessources();
	for (int k=0; k<maxnb_interactionpair ; k++) {
		interaction[k].init(this);
		interaction[k].set_label(k);
	}
	for (int i=0; i<np; i++) {
		radius_cubed[i] = radius[i]*radius[i]*radius[i];
		if (twodimension) {
			angle[i] = 0;
		}
		velocity[i].reset();
		na_velocity[i].reset();
		ang_velocity[i].reset();
		na_ang_velocity[i].reset();
		vel_contact[i].reset();
		ang_vel_contact[i].reset();
		vel_hydro[i].reset();
		ang_vel_hydro[i].reset();
		if (repulsiveforce) {
			vel_repulsive[i].reset();
			ang_vel_repulsive[i].reset();
		}
	}
	/* Prepare
	 *
	 */
	for (int i=0; i<18*np; i++) {
		resistance_matrix_dblock[i] = 0;
	}
	double torque_factor = 4.0/3;
	for (int i=0; i<np; i++) {
		int i18 = 18*i;
		double FUvalue = sd_coeff*radius[i];
		double TWvalue = sd_coeff*torque_factor*radius_cubed[i];
		resistance_matrix_dblock[i18   ] = FUvalue;
		resistance_matrix_dblock[i18+6 ] = FUvalue;
		resistance_matrix_dblock[i18+10] = FUvalue;
		resistance_matrix_dblock[i18+12] = TWvalue;
		resistance_matrix_dblock[i18+15] = TWvalue;
		resistance_matrix_dblock[i18+17] = TWvalue;
	}
	shear_strain = 0;
	shear_disp = 0;
	nb_interaction = 0;
	sq_lub_max = lub_max*lub_max; // square of lubrication cutoff length.
	if (contact_relaxation_time < 0) {
		// 1/(h+c) --> 1/c
		lub_coeff_contact = 1/lub_reduce_parameter;
	} else {
		/* t = beta/kn
		 *  beta = t*kn
		 * lub_coeff_contact = 4*beta = 4*kn*contact_relaxation_time
		 *
		 * For Low Peclet mode:
		 * kn = scale_factor_SmallPe*kn_lowPeclet;
		 * contact_relaxation_time = contact_relaxation_time/scale_factor_SmallPe;
		 * This is why the coeffient is not scaled.
		 * scale_factor_SmallPe*kn_lowPeclet * contact_relaxation_time/scale_factor_SmallPe;
		 * = kn_lowPeclet * contact_relaxation_time
		 */
		lub_coeff_contact = 4*kn*contact_relaxation_time;
	}
	cerr << "lub_coeff_contact = " << lub_coeff_contact << endl;
	cerr << "1/lub_reduce_parameter = " <<  1/lub_reduce_parameter << endl;
	/* t = beta/kn
	 *  beta = t*kn
	 * lub_coeff_contact = 4*beta = 4*kn*contact_relaxation_time
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
		/* [Note]
		 * We finally do not introduce a dashpot for the sliding mode.
		 * This is set in the parameter file, i.e. contact_relaxation_time_tan = 0
		 * So log_lub_coeff_contact_tan_dashpot = 0;
		 */
		log_lub_coeff_contact_tan_dashpot = 6*kt*contact_relaxation_time_tan;
	} else {
		cerr << "lubrication_model..." << endl;
		exit(1);
	}
	log_lub_coeff_contact_tan_total = log_lub_coeff_contact_tan_dashpot+log_lub_coeff_contact_tan_lubrication;
	if (brownian) {
		/* In developing and debugging phases,
		 * we give a seed to generate the same series of random number.
		 *
		 * r_gen = new MTRand;
		 */
		r_gen = new MTRand(17);
	}
	cerr << "log_lub_coeff_contact_tan_lubrication = " << log_lub_coeff_contact_tan_total << endl;
	cerr << "log_lub_coeff_contact_tan_dashpot = " << log_lub_coeff_contact_tan_dashpot << endl;
	ts = 0;
	shear_disp = 0;
	/* shear rate is fixed to be 1 in dimensionless simulation
	 */
	vel_difference = lz;
	//after_parameter_changed = false;
	stokes_solver.initialize();
	dt = dt_max;
	initializeBoxing();
	checkNewInteraction();
	if (twodimension) {
		setSystemVolume(2*radius[np-1]);
	} else {
		setSystemVolume();
	}
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
System::timeStepBoxing(){
	// evolve PBC
	if (!zero_shear) {
		shear_disp += vel_difference*dt;
		if (shear_disp > lx) {
			shear_disp -= lx;
		}
	}
	boxset.update();
}

void
System::timeEvolutionEulersMethod(bool calc_stress){
	in_predictor = true;
	setContactForceToParticle();
	setRepulsiveForceToParticle();
	computeVelocities(calc_stress);
	if (calc_stress) {
		calcStressPerParticle();
		avgStressUpdate();

		/****** <exp> ********/
		if(stress_controlled){
			StressTensor total_contact_stress;
			StressTensor total_repulsive_stress;

			avgStressReset();
			calcStressPerParticle();
			avgStressUpdate();
			calcStress();
			total_contact_stress = total_contact_stressXF_normal+total_contact_stressXF_tan + total_contact_stressGU;
			if (repulsiveforce) {
				total_repulsive_stress = total_repulsive_stressXF+total_repulsive_stressGU;
			}
			if (brownian) {
				cerr << " Stress controlled Brownian simulations not implemented ! " << endl;
				exit(1);
			}
			double sr = target_stress - total_repulsive_stress.getStressXZ();
			sr /= total_hydro_stress.getStressXZ() + total_contact_stress.getStressXZ();
			dimensionless_shear_rate = sr/repulsiveforce_amplitude;
			
			for (int i=0; i<np; i++) {
				na_velocity[i] += vel_repulsive[i]*( 1/sr - 1 );
				na_ang_velocity[i] += ang_vel_repulsive[i]*( 1/sr - 1 );
				velocity[i] += vel_repulsive[i]*( 1/sr - 1 );
				ang_velocity[i] += ang_vel_repulsive[i]*( 1/sr - 1 );
			}
			for (int i=0; i<np; i++) {
				vel_repulsive[i] /= sr;
				ang_vel_repulsive[i] /= sr;
			}
			avgStressReset();
			calcStressPerParticle();
			avgStressUpdate();
			calcStress();
			total_contact_stress = total_contact_stressXF_normal + total_contact_stressXF_tan + total_contact_stressGU;
			if (repulsiveforce) {
				total_repulsive_stress = total_repulsive_stressXF+total_repulsive_stressGU;
			}
			StressTensor total_stress = total_hydro_stress + total_contact_stress + total_repulsive_stress;
			//			cout << total_contact_stress.getStressXZ() << " " << total_hydro_stress.getStressXZ() << " " << total_repulsive_stress.getStressXZ() << endl; 
			//			cout << dimensionless_shear_rate << " " << total_stress.getStressXZ()*sr << endl;
			//			getchar();
		}
		/****** </exp> ********/


	}
	timeStepMove();
}

/****************************************************************************************************
 ******************************************** Mid-Point Scheme ***************************************
 ****************************************************************************************************/

/************************************************************************************************
 *                                   non-Brownian Case                                           *
 *************************************************************************************************
 * Simple mid-point method to solve at dt^2 order
 * R.V = F_H + F_C + F_B
 * where R is the resistance, F_H/F_C are hydro/contact forces.
 *
 * 1st step:
 * x'(t+dt) = x(t) + V^{-}dt
 * x(t)     = x'(t+dt) - V^{-}dt
 *
 * 2nd step
 * x(t + dt) = x(t)     + 0.5*(V^{+}+V^{-})*dt
 *           = x'(t+dt) + 0.5*(V^{+}-V^{-})*dt
 */

/*************************************************************************************************
 *                                   Brownian Case                                                *
 **************************************************************************************************
 * This routine implements a two-step algorithm for the dynamics with Brownian motion,
 * initially derived by [ Fixman 1978 ].
 * The basis of this algorithm is exposed in [ Ball & Melrose 1997 ] and [ Banchio & Brady 2003 ].
 *
 * The equation of motion is:
 *
 * R.V = F_H + F_C + F_B
 * where R is the resistance, F_H/F_C/F_B are hydro/contact/brownian forces.
 *
 * Integrating this, we get a first order update of the positions as:
 * X(t+dt) = X(t) + dt*R^{-1}.( F_H + F_C ) + X_B
 * with <X_B> = kT*dt*div R^{-1} and <X_B X_B> - < X_B >^2 = (2kTdt)R^{-1}
 *
 * The divergence term comes from our naive way of taking the inertialess limit in the Langevin
 * equation correlated with the fact that our time step dt >> "Brownian_time" (the typical
 * time between 2 Brownian kicks). The sum of the many displacements due to Brownian kicks
 * happening during dt has a non-zero mean anywhere the mobility is non uniform,
 * and this is taken care of by the div term.
 * This sum has a variance scaling as \sqrt(kT*dt), taken care of by the X_B term, as
 * <X_B X_B> - < X_B >^2 = (2kTdt)R^{-1} (which is nothing but the FD theorem).
 *
 * Reminding that we obtain the Cholesky decomposition R = L L^T in the stokes_solver,
 * we obtain X_B in a 2-step algorithm [ B & M 1997 ]:
 *   - generate a "Brownian force" F_B(t) such that:
 *       < F_B F_B > = (2kT/dt) R(t) and <F_B> = 0
 *   - (1) solve R(X(t)).V_{-} = F_H(t) + F_C(t) + F_B(t)
 *   - move to X_pred = X(t) + V_{-}dt
 *   - (2) solve R(X_pred).V_{+} = F_H(t+dt) + F_C(t+dt) + F_B(t)  (note: F_B at time t)
 *   - set X(t+dt) = X(t) + 0.5*( V_{-} + V_{+} )*dt = X_pred + 0.5*( V_{+} - V_{-} )*dt
 *
 * One can check that the velocity U = 0.5*( V_{-} + V_{+} ) has the properties:
 * <U> = R^{-1}.( F_H + F_C ) + kT*div R^{-1}
 * <UU> - <U>^2 = (2kT/dt)R^{-1}
 * which gives the correct mean and variance for X_B.
 *
 * F_B is obtained as F_B = \sqrt(2kT/dt) * L.A, where A is a gaussian random vector with
 * <A> = 0 and <AA> = 1.
 *
 *
 ***************************************************************************************************/
void
System::timeEvolutionPredictorCorrectorMethod(bool calc_stress){
	/* predictor */
	in_predictor = true;
	setContactForceToParticle();
	setRepulsiveForceToParticle();
	computeVelocities(calc_stress);
	if (calc_stress) {
		calcStressPerParticle();
	}
	timeStepMovePredictor();
	/* corrector */
	in_predictor = false;
	setContactForceToParticle();
	setRepulsiveForceToParticle();
	computeVelocities(calc_stress);
	if (calc_stress) {
		calcStressPerParticle();
		avgStressUpdate();
	}
	timeStepMoveCorrector();
}

/*
 * timeStepMove is used only for Euler method.
 */
void
System::timeStepMove(){
	/* evolve PBC */
	timeStepBoxing();
	/* move particles */
	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (twodimension) {
		for (int i=0; i<np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	checkNewInteraction();
	updateInteractions();
}

void
System::timeStepMovePredictor(){
	/* The periodic boundary condition is updated in predictor.
	 * It must not be updated in corrector.
	 */
	timeStepBoxing();
	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (twodimension) {
		for (int i=0; i<np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
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
System::timeStepMoveCorrector(){
	for (int i=0; i<np; i++) {
		velocity[i] = 0.5*(velocity[i]+velocity_predictor[i]);  // real velocity, in predictor and in corrector
		ang_velocity[i] = 0.5*(ang_velocity[i]+ang_velocity_predictor[i]);
	}
	for (int i=0; i<np; i++) {
		displacement(i, (velocity[i]-velocity_predictor[i])*dt);
	}
	if (twodimension) {
		for (int i=0; i<np; i++) {
			angle[i] += (ang_velocity[i].y-ang_velocity_predictor[i].y)*dt;
		}
	}
	checkNewInteraction();
	updateInteractions();
}

void
System::timeEvolution(double strain_next){
	static bool firsttime = true;
	if (firsttime) {
		checkNewInteraction();
		firsttime = false;
	}
	avgStressReset();
	if(strain_controlled){
		while (shear_strain < strain_next-dt-dt*0.001) { // integrate until strain_next - 1 time step
			(this->*timeEvolutionDt)(false); // no stress computation
			ts++;
			shear_strain += dt;
		};
		(this->*timeEvolutionDt)(true); // last time step, compute the stress
	}

	if(stress_controlled){
		while (shear_strain < strain_next-dt*0.001) { // integrate until strain_next
			(this->*timeEvolutionDt)(true); // stress computation
			ts++;
			shear_strain += dt;
		};
	}

	ts++;
	shear_strain += dt;
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
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			bool deactivated = false;
			interaction[k].updateState(deactivated);
			if (deactivated) {
				deactivated_interaction.push(k);
			}
		}
	}
}

void
System::stressReset(){
	for (int i=0; i<np; i++) {
		lubstress[i].reset();
		contactstressGU[i].reset();
	}
	if (repulsiveforce) {
		for (int i=0; i<np; i++) {
			repulsivestressGU[i].reset();
		}
	}
	if (brownian) {
		for (int i=0; i<np; i++) {
			brownianstressGU[i].reset();
		}
	}
}

void
System::avgStressReset(){
	avg_stress_nb = 0;
	for (int i=0; i<np; i++) {
		avg_lubstress[i].reset();
		avg_contactstressGU[i].reset();
	}
	avg_contactstressXF_normal.reset();
	avg_contactstressXF_tan.reset();
	if (repulsiveforce) {
		avg_repulsivestressXF.reset();
		for (int i=0; i<np; i++) {
			avg_repulsivestressGU[i].reset();
		}
	}
	if (brownian) {
		for (int i=0; i<np; i++) {
			avg_brownianstressGU[i].reset();
		}
	}
}

void
System::avgStressUpdate(){
	avg_stress_nb++;
	for (int i=0; i<np; i++) {
		avg_lubstress[i] += lubstress[i];
		avg_contactstressGU[i] += contactstressGU[i];
	}
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()) {
			avg_contactstressXF_normal += interaction[k].contact.getContactStressXF_normal();
			avg_contactstressXF_tan += interaction[k].contact.getContactStressXF_tan();
		}
	}
	if (repulsiveforce) {
		for (int i=0; i<np; i++) {
			avg_repulsivestressGU[i] += repulsivestressGU[i];
		}
		for (int k=0; k<nb_interaction; k++) {
			avg_repulsivestressXF += interaction[k].getRepulsiveStressXF();
		}
	}
	if (brownian) {
		for (int i=0; i<np; i++) {
			avg_brownianstressGU[i] += brownianstressGU[i];
		}
	}
}

void
System::buildHydroTerms(bool build_res_mat, bool build_force_GE){
	// Builds the following terms, according to the value of 'build_res_mat' and 'build_force_GE':
	//  - elements of the resistance matrix if 'build_res_mat' is true
	//       (only terms diverging as 1/h if lubrication_model == 1, terms in 1/h and log(1/h) for lubrication_model>1 )
	//  - vector Gtilde*Einf if 'build_force_GE' is true (default behavior)
	//
	// Note that it ADDS the rhs of the solver as rhs += GE. You need to call stokes_solver.resetRHS() before this routine
	// if you want GE to be the only rhs.
	//static int nb_of_active_interactions_in_predictor;
	if (build_res_mat) {
		// create a new resistance matrix in stokes_solver
		nb_of_active_interactions = nb_interaction-deactivated_interaction.size();
		stokes_solver.resetResistanceMatrix("direct", nb_of_active_interactions, resistance_matrix_dblock);
		/* [note]
		 * The resistance matrix is reset with resistance_matrix_dblock,
		 * which is calculated at the begining.
		 * addStokesDrag() is no more used.
		 */
		// add GE in the rhs and lubrication terms in the resistance matrix
		(this->*buildLubricationTerms)(true, build_force_GE); // false: don't modify rhs, as we want rhs=F_B
		stokes_solver.completeResistanceMatrix();
	} else {
		// add GE in the rhs
		(this->*buildLubricationTerms)(false, build_force_GE); // false: don't modify rhs, as we want rhs=F_B
	}
}

/* We solve A*(U-Uinf) = Gtilde*Einf ( in Jeffrey's notations )
 * This method computes:
 *  - elements of the resistance matrix if 'mat' is true
 *       (only terms diverging as 1/h if lubrication_model == 1, terms in 1/h and log(1/h) for lubrication_model>1 )
 *  - vector Gtilde*Einf if 'rhs' is true (default behavior)
 */
void
System::buildLubricationTerms_squeeze(bool mat, bool rhs){
	for (int i=0; i<np-1; i ++) {
		for (set<Interaction*>::iterator it = interaction_list[i].begin();
			 it != interaction_list[i].end(); it ++) {
			int j = (*it)->partner(i);
			if (j > i) {
				if (mat) {
					vec3d nr_vec = (*it)->get_nvec();
					(*it)->lubrication.calcXFunctions();
					stokes_solver.addToDiagBlock(nr_vec, i,
												 (*it)->lubrication.scaledXA0(), 0, 0, 0);
					stokes_solver.addToDiagBlock(nr_vec, j,
												 (*it)->lubrication.scaledXA3(), 0, 0, 0);
					stokes_solver.setOffDiagBlock(nr_vec, i, j,
												  (*it)->lubrication.scaledXA2(), 0, 0, 0, 0);
				}
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
}

void
System::buildLubricationTerms_squeeze_tangential(bool mat, bool rhs){
	for (int i=0; i<np-1; i ++) {
		for (set<Interaction*>::iterator it = interaction_list[i].begin();
			 it != interaction_list[i].end(); it ++) {
			int j = (*it)->partner(i);
			if (j > i) {
				if (mat) {
					vec3d nr_vec = (*it)->get_nvec();
					(*it)->lubrication.calcXYFunctions();
					stokes_solver.addToDiagBlock(nr_vec, i,
												 (*it)->lubrication.scaledXA0(),
												 (*it)->lubrication.scaledYA0(),
												 (*it)->lubrication.scaledYB0(),
												 (*it)->lubrication.scaledYC0());
					stokes_solver.addToDiagBlock(nr_vec, j,
												 (*it)->lubrication.scaledXA3(),
												 (*it)->lubrication.scaledYA3(),
												 (*it)->lubrication.scaledYB3(),
												 (*it)->lubrication.scaledYC3());
					stokes_solver.setOffDiagBlock(nr_vec, i, j,
												  (*it)->lubrication.scaledXA1(),
												  (*it)->lubrication.scaledYA1(),
												  (*it)->lubrication.scaledYB2(),
												  (*it)->lubrication.scaledYB1(),
												  (*it)->lubrication.scaledYC1());
				}
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
}

void
System::generateBrownianForces(){
	// generates a Brownian force F_B with <F_B> = 0, and <F_B F_B> = (2kT/dt)*R
	// where R is the current resistance matrix stored in the stokes_solver.
	// note that it SETS the rhs of the solver as rhs = F_B
	// F_B is also stored in sys->brownian_force
	//
	// kb_T = 1/Pe and dt = dt'*Pe/Pe0
	// kb_T/dt = (1/Pe)/(dt'*Pe/Pe0) = Pe0/(dt'*Pe*Pe)
	// Fb = sqrt(kb_T/dt) = sqrt(Pe0/dt')*1/Pe
	// U = Fb
	// U*dt = sqrt(Pe0/dt')*1/Pe * dt'*Pe/Pe0 = sqrt(dt'/Pe0)
	double sqrt_kbT2_dt = sqrt(2*kb_T/dt); // proportional to 1/Pe.
	for (int i=0; i<linalg_size; i++) {
		brownian_force[i] = sqrt_kbT2_dt*GRANDOM; // random vector A
	}
	stokes_solver.setRHS(brownian_force);
	stokes_solver.compute_LTRHS(brownian_force); // F_B = \sqrt(2kT/dt) * L^T * A
	if (twodimension) {
		for (int i=0; i<np; i++) {
			brownian_force[6*i+1] = 0; // Fy
			brownian_force[6*i+3] = 0; // Tx
			brownian_force[6*i+5] = 0; // Tz
		}
	}
}

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
System::setRepulsiveForceToParticle(){
	if (repulsiveforce) {
		for (int i=0; i<np; i++) {
			repulsive_force[i].reset();
		}
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_active()) {
				interaction[k].addUpRepulsiveForce();
			}
		}
	}
}

void
System::buildContactTerms(bool set_or_add){
    // sets or adds ( set_or_add = t or f resp) contact forces to the rhs of the stokes_solver.
	if (set_or_add) {
		for (int i=0; i<np; i++) {
			stokes_solver.setRHSForce(i, contact_force[i]);
			stokes_solver.setRHSTorque(i, contact_torque[i]);
		}
	} else {
		for (int i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, contact_force[i]);
			stokes_solver.addToRHSTorque(i, contact_torque[i]);
		}
    }
}

void
System::buildRepulsiveForceTerms(bool set_or_add){
	if (set_or_add) {
		for (int i=0; i<np; i++) {
			stokes_solver.setRHSForce(i, repulsive_force[i]);
		}
	} else {
		for (int i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, repulsive_force[i]);
		}
	}
}

void
System::computeVelocities(bool divided_velocities){
	stokes_solver.resetRHS();
	if (divided_velocities) {
		// in case we want to compute the stress contributions
		if (!zero_shear) {
			buildHydroTerms(true, true); // build matrix and rhs force GE
		} else {
			buildHydroTerms(true, false); // zero shear-rate
		}
		stokes_solver.solve(vel_hydro, ang_vel_hydro); // get V_H
		buildContactTerms(true); // set rhs = F_C
		stokes_solver.solve(vel_contact, ang_vel_contact); // get V_C
		for (int i=0; i<np; i++) {
			na_velocity[i] = vel_hydro[i]+vel_contact[i];
			na_ang_velocity[i] = ang_vel_hydro[i]+ang_vel_contact[i];
		}
		if (repulsiveforce) {
			buildRepulsiveForceTerms(true); // set rhs = F_repulsive
			stokes_solver.solve(vel_repulsive, ang_vel_repulsive); // get V_repulsive
			for (int i=0; i<np; i++) {
				na_velocity[i] += vel_repulsive[i];
				na_ang_velocity[i] += ang_vel_repulsive[i];
			}
		}
	} else {
		// for most of the time evolution
		if (!zero_shear) {
			buildHydroTerms(true, true); // build matrix and rhs force GE
		} else {
			buildHydroTerms(true, false); // zero shear-rate
		}
		buildContactTerms(false); // add rhs += F_C
		if (repulsiveforce) {
			buildRepulsiveForceTerms(false); // add rhs += F_repulsive
		}
		stokes_solver.solve(na_velocity, na_ang_velocity); // get V
	}
	
	if (brownian) {
		if (in_predictor) {
			/* generate new F_B only in predictor
			 * Resistance matrix is used.
			 */
			generateBrownianForces();
		}
		stokes_solver.setRHS(brownian_force); // set rhs = F_B
		stokes_solver.solve(vel_brownian, ang_vel_brownian); // get V_B
		for (int i=0; i<np; i++) {
			na_velocity[i] += vel_brownian[i];
			na_ang_velocity[i] += ang_vel_brownian[i];
		}
	}
	for (int i=0; i<np; i++) {
		velocity[i] = na_velocity[i];
		ang_velocity[i] = na_ang_velocity[i];
		if (!zero_shear) {
			velocity[i].x += position[i].z;
			ang_velocity[i].y += 0.5;
		}
	}
	stokes_solver.solvingIsDone();
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
	if (z_shift != 0) {
		velocity[i].x += z_shift*vel_difference;
	}
	boxset.box(i);
}

// [0,l]
int
System::periodize(vec3d &pos){
	int z_shift;
	if (pos.z >= lz) {
		pos.z -= lz;
		pos.x -= shear_disp;
		z_shift = -1;
	} else if (pos.z < 0) {
		pos.z += lz;
		pos.x += shear_disp;
		z_shift = 1;
	} else {
		z_shift = 0;
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

// periodize + give z_shift= number of boundaries crossed in z-direction
void
System::periodize_diff(vec3d &pos_diff, int &zshift){
	/*
	 * The displacement of the second particle along z direction
	 * is zshift * lz;
	 */
	if (pos_diff.z > lz_half) {
		pos_diff.z -= lz;
		pos_diff.x -= shear_disp;
		zshift = -1;
	} else if (pos_diff.z < -lz_half) {
		pos_diff.z += lz;
		pos_diff.x += shear_disp;
		zshift = 1;
	} else {
		zshift = 0;
	}
	if (pos_diff.x > lx_half) {
		pos_diff.x -= lx;
		if (pos_diff.x > lx_half) {
			pos_diff.x -= lx;
		}
	} else if (pos_diff.x < -lx_half) {
		pos_diff.x += lx;
		if (pos_diff.x < -lx_half) {
			pos_diff.x += lx;
		}
	}
	if (pos_diff.y > ly_half) {
		pos_diff.y -= ly;
	} else if (pos_diff.y < -ly_half) {
		pos_diff.y += ly;
	}
}

void
System::evaluateMaxContactVelocity(){
	max_contact_velo_tan = 0;
	max_contact_velo_normal = 0;
	max_relative_velocity = 0;
	max_sliding_velocity = 0;
	double sum_contact_velo_tan = 0;
	double sum_contact_velo_normal = 0;
	double sum_sliding_velocity = 0;
	int cnt_contact = 0;
	int cnt_sliding = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_contact()) {
			cnt_contact++;
			sum_contact_velo_tan += interaction[k].getContactVelocity();
			sum_contact_velo_normal += abs(interaction[k].getNormalVelocity());
			if (interaction[k].getContactVelocity() > max_contact_velo_tan) {
				// relative_surface_velocity for both static friction and sliding state.
				max_contact_velo_tan = interaction[k].getContactVelocity();
			}
			if (abs(interaction[k].getNormalVelocity()) > max_contact_velo_normal) {
				max_contact_velo_normal = abs(interaction[k].getNormalVelocity());
			}
			if (interaction[k].getRelativeVelocity() > max_relative_velocity) {
				max_relative_velocity = interaction[k].getRelativeVelocity();
			}
			if (interaction[k].contact.state == 3) {
				/*
				 * relative_surface_velocity for only sliding state.
				 */
				cnt_sliding++;
				sum_sliding_velocity += interaction[k].getContactVelocity();
				if (interaction[k].getContactVelocity() > max_sliding_velocity) {
					max_sliding_velocity = interaction[k].getContactVelocity();
				}
			}
		}
	}
	if (cnt_contact > 0) {
		ave_contact_velo_tan = sum_contact_velo_tan/cnt_contact;
		ave_contact_velo_normal = sum_contact_velo_normal/cnt_contact;
	} else {
		ave_contact_velo_tan = 0;
		ave_contact_velo_normal = 0;
	}
	if (cnt_sliding > 0) {
		ave_sliding_velocity = sum_sliding_velocity/cnt_sliding;
	} else {
		ave_sliding_velocity = 0;
	}
}

double
System::evaluateMaxVelocity(){
	double sq_max_velocity = 0;
	for (int i = 0; i < np; i++) {
		vec3d na_velocity_tmp = velocity[i];
		if (zero_shear) {
			na_velocity_tmp.x -= position[i].z;
		}
		if (na_velocity_tmp.sq_norm() > sq_max_velocity) {
			sq_max_velocity = na_velocity_tmp.sq_norm();
		}
	}
	return sqrt(sq_max_velocity);
}

double
System::evaluateMaxAngVelocity(){
	double _max_ang_velocity = 0;
	for (int i = 0; i < np; i++) {
		vec3d na_ang_velocity_tmp = ang_velocity[i];
		na_ang_velocity_tmp.y -= 0.5;
		if (na_ang_velocity_tmp.norm() > _max_ang_velocity) {
			_max_ang_velocity = na_ang_velocity_tmp.norm();
		}
	}
	return _max_ang_velocity;
}

double
System::evaluateMinGap(){
	double _min_gap_nondim = 100000;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].get_gap_nondim() < _min_gap_nondim) {
			_min_gap_nondim = interaction[k].get_gap_nondim();
		}
	}
	return _min_gap_nondim;
}

double
System::evaluateMaxDispTan(){
	double _max_disp_tan = 0;
	for (int k= 0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].contact.disp_tan.norm() > _max_disp_tan) {
			_max_disp_tan = interaction[k].contact.disp_tan.norm();
		}
	}
	return _max_disp_tan;
}

double
System::evaluateMaxFcNormal(){
	double max_fc_normal_ = 0;
 	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].is_contact() &&
			interaction[k].contact.get_f_contact_normal_norm() > max_fc_normal_) {
			max_fc_normal_ = interaction[k].contact.get_f_contact_normal_norm();
		}
	}
	return max_fc_normal_;
}

double
System::evaluateMaxFcTangential(){
	double max_fc_tan_ = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].is_contact() &&
			interaction[k].contact.get_f_contact_tan_norm() > max_fc_tan_) {
			max_fc_tan_ = interaction[k].contact.get_f_contact_tan_norm();
		}
	}
	return max_fc_tan_;
}

void
System::countNumberOfContact(){
	contact_nb = 0;
	fric_contact_nb = 0;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].is_contact()) {
			contact_nb ++;
			if (interaction[k].is_friccontact()) {
				fric_contact_nb ++;
			}
		}
	}
}

void
System::analyzeState(){
	max_velocity = evaluateMaxVelocity();
	max_ang_velocity = evaluateMaxAngVelocity();
	evaluateMaxContactVelocity();
	min_gap_nondim = evaluateMinGap();
	max_disp_tan = evaluateMaxDispTan();
	max_fc_normal = evaluateMaxFcNormal();
	max_fc_tan = evaluateMaxFcTangential();
	countNumberOfContact();
}

void
System::setSystemVolume(double depth){
	if (twodimension) {
		system_volume = lx*lz*depth;
	} else {
		system_volume = lx*ly*lz;
	}
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
	/* determination of kn
	 */
	cerr << "We should make a simplar rule." << endl;
	exit(1);
	//	double mean_max_fc_normal, stddev_max_fc_normal;
	//	calcMean_StdDev(max_fc_normal_history, mean_max_fc_normal, stddev_max_fc_normal);
	//	double kn_try = mean_max_fc_normal/overlap_target;
	//	kn = kn_try;
	//	lub_coeff_contact = 4*kn*contact_relaxation_time;
	/* determination of kt
	 */
	//	double mean_max_fc_tan, stddev_max_fc_tan;
	//	calcMean_StdDev(max_fc_tan_history, mean_max_fc_tan, stddev_max_fc_tan);
	//	double kt_try = mean_max_fc_tan/disp_tan_target;
	//	kt = kt_try;
	//	double average_max_tanvelocity = 0;
	//	double max_max_tanvelocity = 0;
	//for (unsigned int j=0; j<sliding_velocity_history.size(); j++){
	//	average_max_tanvelocity += sliding_velocity_history[j];
	//	if (max_max_tanvelocity < sliding_velocity_history[j]){
	//		max_max_tanvelocity = sliding_velocity_history[j];
	//		}
	//}
	//average_max_tanvelocity = average_max_tanvelocity/sliding_velocity_history.size();
	//	double average_max_relative_velocity = 0;
	//	for (unsigned int j=0; j<relative_velocity_history.size(); j++){
	//		average_max_relative_velocity += relative_velocity_history[j];
	//	}
	//	average_max_relative_velocity = average_max_relative_velocity/relative_velocity_history.size();
	double tmp_max_velocity = 0;
	//	if (average_max_relative_velocity > average_max_tanvelocity){
	//		tmp_max_velocity = average_max_relative_velocity ;
	//	} else {
	//		tmp_max_velocity = average_max_tanvelocity ;
	//	}
	//	if (max_max_tanvelocity > 1000){
	//		cerr << "max_max_tanvelocity = " << max_max_tanvelocity << endl;
	//		return 1;
	//	}
	double dt_try = disp_max/tmp_max_velocity;
	if (dt_try < dt_max){
		dt = dt_try;
	}
	for (int k=0; k<nb_interaction; k++) {
		interaction[k].contact.updateContactModel();
	}
	//	max_fc_normal_history.clear();
	//	max_fc_tan_history.clear();
	//	sliding_velocity_history.clear();
	//after_parameter_changed = true;
	if (kn > max_kn){
		cerr << "kn = " << kn << endl;
		cerr << " kn > max_kn : exit" << endl;
		return 1;
	}
	return 0;
}

void
System::calcLubricationForce(){
	/*
	 * Calculate lubrication force to output
	 */
	stokes_solver.resetRHS();
	if (!zero_shear) {
		buildHydroTerms(true, true);
	} else {
		buildHydroTerms(true, false); // no GE
	}
    setContactForceToParticle();
	buildContactTerms(false);
	if (repulsiveforce) {
		setRepulsiveForceToParticle();
		buildRepulsiveForceTerms(false);
	}
	stokes_solver.solve(na_velocity, na_ang_velocity);
	stokes_solver.solvingIsDone();
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			interaction[k].lubrication.calcLubricationForce();
		}
	}
}

