//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
#include <sstream>
#include <cmath>
#define DELETE(x) if(x){delete [] x; x = NULL;}
#define GRANDOM ( r_gen->randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.

System::System():
maxnb_interactionpair_per_particle(15),
lowPeclet(false),
brownian(false),
zero_shear(false),
friction_model(-1),
repulsiveforce(false),
new_contact_gap(0),
init_strain_shear_rate_limit(0),
init_shear_rate_limit(999)
{}

System::~System()
{
	DELETE(position);
	DELETE(radius);
	DELETE(radius_squared);
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
	DELETE(interaction);
	DELETE(interaction_list);
	DELETE(interaction_partners);
	if (brownian) {
		DELETE(vel_brownian);
		DELETE(ang_vel_brownian);
		DELETE(brownian_force);
		DELETE(brownianstressGU);
		DELETE(brownianstressGU_predictor);
	}
	if (repulsiveforce) {
		DELETE(repulsivestressGU);
		DELETE(vel_repulsive);
		DELETE(ang_vel_repulsive);
	}
};

void System::importParameterSet(ParameterSet &ps)
{
	p = ps;
	friction_model = p.friction_model;
	rolling_friction = p.rolling_friction;
	set_lub_max_gap(p.lub_max_gap);
	lub_reduce_parameter = p.lub_reduce_parameter;
	set_lubrication_model(p.lubrication_model);
	kn = p.kn;
	kt = p.kt;
	kr = p.kr;
	if (repulsiveforce) {
		set_repulsiveforce_length(p.repulsive_length);
	} else {
		set_repulsiveforce_length(0);
	}
	set_sd_coeff(p.sd_coeff);
	set_integration_method(p.integration_method);
	mu_static = p.mu_static;
	if (p.mu_dynamic == -1) {
		mu_dynamic = p.mu_static;
	} else {
		mu_dynamic = p.mu_dynamic;
	}
	mu_rolling = p.mu_rolling;
	set_disp_max(p.disp_max);
}

void System::allocateRessources()
{
	linalg_size = 6*np;
	maxnb_interactionpair = maxnb_interactionpair_per_particle*np;
	radius_cubed = new double [np];
	radius_squared = new double [np];
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
	if (brownian) {
		brownianstressGU = new StressTensor [np];
		brownianstressGU_predictor = new StressTensor [np];
	}
	interaction = new Interaction [maxnb_interactionpair];
	interaction_list = new set <Interaction*> [np];
	interaction_partners = new set <int> [np];
	stokes_solver.init(np);
	//
	if (p.auto_determine_knkt) {
		kn_avg = new Averager<double>(p.memory_strain_avg);
		kt_avg = new Averager<double>(p.memory_strain_avg);
		overlap_avg = new Averager<double>(p.memory_strain_avg);
		max_disp_tan_avg = new Averager<double>(p.memory_strain_avg);
	}
	if (brownian) {
		if (lowPeclet) {
			double stress_avg_relaxation_parameter = 0; // 0 --> no average
			stress_avg = new Averager<StressTensor>(stress_avg_relaxation_parameter);
		}
	}
}

void System::setInteractions_GenerateInitConfig()
{
	calcInteractionRange = &System::calcLubricationRange;
	for (int k=0; k<maxnb_interactionpair; k++) {
		interaction[k].init(this);
		interaction[k].set_label(k);
	}
	nb_interaction = 0;
	shear_strain = 0;
	shear_disp = 0;
	vel_difference = 0;
	initializeBoxing();
	checkNewInteraction();
}

void System::allocatePositionRadius()
{
	position = new vec3d [np];
	radius = new double [np];
}

void System::setConfiguration(const vector <vec3d> &initial_positions,
						 const vector <double> &radii,
						 double lx_, double ly_, double lz_)
{
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
	if (twodimension) {
		/* [note]
		 * The depth of mono-layer is the diameter of the largest particles.e
		 * The sample needs to be labeled from smaller particles to larger particles
		 * in the configuration file.
		 */
		double largest_diameter = 2*radius[np-1];
		setSystemVolume(largest_diameter);
	} else {
		setSystemVolume();
	}
	double particle_volume = 0;
	for (int i=0; i<np; i++) {
		particle_volume += (4*M_PI/3)*pow(radius[i],3);
	}
	volume_fraction = particle_volume/system_volume;
	cerr << "volume_fraction = " << volume_fraction << endl;
}

void System::updateUnscaledContactmodel()
{
	if (abs(target_stress) != 0) {
		/* What is the reasonable way
		 * when the target stress is changed during a simulation?
		 *
		 * [temporally change]
		 * For destressing tests, the spring constants are fixed at the previous values.
		 */
		if (!cohesion) {
			kn = kn_master*abs(target_stress);
			kt = kt_master*abs(target_stress);
			kr = kr_master*abs(target_stress);
		} else {
			kn = kn_master;
			kt = kt_master;
			kr = kr_master;
		}
		cout << " kn " << kn << "  kn_master " << kn_master << " target_stress "  << target_stress << endl;
	}
	lub_coeff_contact = 4*kn*p.contact_relaxation_time;
	if(stress_controlled){
		lub_coeff_contact /= abs(dimensionless_number);
	}
	
	if (lowPeclet) {
		lub_coeff_contact *= p.Pe_switch;
	}

	if (lubrication_model == 1) {
		log_lub_coeff_contact_tan_lubrication = 0;
		log_lub_coeff_contact_tan_dashpot = 0;
	} else if (lubrication_model == 2) {
		log_lub_coeff_contact_tan_lubrication = log(1/lub_reduce_parameter);
		/* [Note]
		 * We finally do not introduce a dashpot for the sliding mode.
		 * This is set in the parameter file, i.e. p.contact_relaxation_time_tan = 0
		 * So log_lub_coeff_contact_tan_dashpot = 0;
		 */
		log_lub_coeff_contact_tan_dashpot = 6*kt*p.contact_relaxation_time_tan;
	} else {
		cerr << "lubrication_model..." << endl;
		exit(1);
	}
	log_lub_coeff_contact_tan_total = log_lub_coeff_contact_tan_dashpot+log_lub_coeff_contact_tan_lubrication;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			interaction[k].contact.getInteractionData();
		}
	}
}

void System::setupBrownian()
{
	if (brownian) {
		if (dimensionless_number < p.Pe_switch) {
			// scale_factor_SmallPe > 1
			lowPeclet = true;
			scale_factor_SmallPe = p.Pe_switch/dimensionless_number;
			p.memory_strain_k /= scale_factor_SmallPe;
			p.memory_strain_avg /= scale_factor_SmallPe;
			p.start_adjust /= scale_factor_SmallPe;
			p.dt *= p.Pe_switch; // to make things continuous at Pe_switch
			cerr << "[[small Pe mode]]" << endl;
			cerr << "  kn = " << kn << endl;
			cerr << "  kt = " << kt << endl;
			cerr << "  dt = " << p.dt << endl;
		} else {
			lowPeclet = false; // <--- Always true?
		}
	}
}

void System::setupSystem(string control)
{
	/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	 * @ We have to consider p.contact_relaxation_time in Brownian case.
	 * @ The resistance coeffient affects Brownian force.
	 * @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
	calcInteractionRange = &System::calcLubricationRange;
	friction = false;
	if (friction_model == 0) {
		cerr << "friction_model = 0" << endl;
		mu_static = 0;
		friction = false;
	} else if (friction_model == 1) {
		cerr << "friction_model = 1" << endl;
		friction = true;
	} else if (friction_model == 2 || friction_model == 3) {
		cerr << "friction_model " << friction_model << endl;
		friction = true;
		cerr << "critical_normal_force = " << critical_normal_force << endl;
	} else if (friction_model == 4) {
		cerr << "friction_model = 4" << endl;
		friction = true;
	} else {
		cerr << "friction_model..." << endl;
		exit(1);
	}
	allocateRessources();
	for (int k=0; k<maxnb_interactionpair; k++) {
		interaction[k].init(this);
		interaction[k].set_label(k);
	}
	for (int i=0; i<np; i++) {
		radius_squared[i] = pow(radius[i],2);
		radius_cubed[i] = pow(radius[i],3);
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
	shear_strain = 0;
	nb_interaction = 0;
	shear_direction = 1;
	if (p.unscaled_contactmodel) {
		kn_master = kn;
		kt_master = kt;
		kr_master = kr;
		updateUnscaledContactmodel();
	}
	if (p.contact_relaxation_time < 0) {
		// 1/(h+c) --> 1/c
		lub_coeff_contact = 1/lub_reduce_parameter;
	} else {
		/* t = beta/kn
		 *  beta = t*kn
		 * lub_coeff_contact = 4*beta = 4*kn*p.contact_relaxation_time
		 *
		 * For Low Peclet mode:
		 * kn = scale_factor_SmallPe*kn_lowPeclet;
		 * p.contact_relaxation_time = p.contact_relaxation_time/scale_factor_SmallPe;
		 * This is why the coeffient is not scaled.
		 * scale_factor_SmallPe*kn_lowPeclet * p.contact_relaxation_time/scale_factor_SmallPe;
		 * = kn_lowPeclet * p.contact_relaxation_time
		 */
		lub_coeff_contact = 4*kn*p.contact_relaxation_time;
	}
	cerr << "lub_coeff_contact = " << lub_coeff_contact << endl;
	cerr << "1/lub_reduce_parameter = " << 1/lub_reduce_parameter << endl;
	/* t = beta/kn
	 *  beta = t*kn
	 * lub_coeff_contact = 4*beta = 4*kn*p.contact_relaxation_time
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
		 * This is set in the parameter file, i.e. p.contact_relaxation_time_tan = 0
		 * So log_lub_coeff_contact_tan_dashpot = 0;
		 */
		log_lub_coeff_contact_tan_dashpot = 6*kt*p.contact_relaxation_time_tan;
	} else {
		cerr << "lubrication_model..." << endl;
		exit(1);
	}
	log_lub_coeff_contact_tan_total = log_lub_coeff_contact_tan_dashpot+log_lub_coeff_contact_tan_lubrication;
	if (brownian) {
#ifdef DEV
		/* In developing and debugging phases,
		 * we give a seed to generate the same series of random number.
		 * DEV is defined as a preprocessor option in the Makefile
		 */
		r_gen = new MTRand(17);	cerr << " WARNING : debug mode: hard coded seed is given to the RNG " << endl;
#endif
#ifndef DEV
		r_gen = new MTRand;
#endif
	}
	cerr << "log_lub_coeff_contact_tan_lubrication = " << log_lub_coeff_contact_tan_total << endl;
	cerr << "log_lub_coeff_contact_tan_dashpot = " << log_lub_coeff_contact_tan_dashpot << endl;
	time = 0;
	/* shear rate is fixed to be 1 in dimensionless simulation
	 */
	vel_difference = lz;
	stokes_solver.initialize();
	dt = p.dt;
	initializeBoxing();
	checkNewInteraction();
	if (control == "rate") {
		rate_controlled = true;
	}
	if (control == "stress") {
		rate_controlled = false;
	}
	stress_controlled = !rate_controlled;
	//	dimensionless_number_averaged = 1;
	/* Pre-calculation
	 */
	/* einstein_stress may be affected by sd_coeff.
	 * However, the reason to set value of sd_coeff is not very certain for the moment.
	 * This is why we limit sd_coeff dependence only the diagonal constant.
	 */
	einstein_stress = (1+2.5*volume_fraction)*shear_rate/(6*M_PI); // 6M_PI because  6\pi eta_0/T_0 = F_0/L_0^2. In System, stresses are in F_0/L_0^2
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
}

void System::initializeBoxing()
{
	/**
		\brief Initialize the boxing system.
		
		Initialize the BoxSet instance using as a minimal Box size the maximal interaction range between any two particles in the System.
	 */
	double range;
	double max_range = 0;
	for (int i=0; i < np-1; i++) { // N^2 init, sorry :(
		for (int j=i+1; j < np; j++) {
			range = (this->*calcInteractionRange)(i,j);
			if (range > max_range) {
				max_range = range;
			}
		}
	}
	boxset.init(max_range, this);
	for (int i=0; i<np; i++) {
		boxset.box(i);
	}
	boxset.update();
}

void System::timeStepBoxing(const double strain_increment)
{
	/**
		\brief Apply a strain step to the boxing system.
	 */
	if (!zero_shear) {
		shear_strain += strain_increment;
		shear_disp += strain_increment*lz;
		int m = (int)(shear_disp/lx);
		if (shear_disp < 0) {
			m--;
		}
		shear_disp = shear_disp-m*lx;
	}
	boxset.update();
}

void System::timeEvolutionEulersMethod(bool calc_stress)
{
	/**
	 \brief One full time step, Euler's method.
	 
	 This method is never used when running a Brownian simulation.
	 */
	in_predictor = true;
	in_corrector = true;
	setContactForceToParticle();
	setRepulsiveForceToParticle();
	computeVelocities(calc_stress);
	if (calc_stress) {
		calcStressPerParticle();
	}
	timeStepMove();
}

/****************************************************************************************************
 ******************************************** Mid-Point Scheme ***************************************
 ****************************************************************************************************/

void System::timeEvolutionPredictorCorrectorMethod(bool calc_stress)
{
	/**
	 \brief One full time step, predictor-corrector method.
	 
	 This method is always used when running a Brownian simulation.
	 
	 ### non-Brownian Case
	 
	 Simple mid-point method to solve at dt^2 order
	 \f$ \bm{A}(\bm{U}-\bm{U}_{\infty}) = \bm{F} \f$
	 where \f$\bm{A} \f$ (in Jeffrey notations \cite
	 jeffrey_calculation_1992, that is \f$\bm{R}_{\mathrm{FU}}\f$ in Bossis and Brady
	 \cite brady_stokesian_1988 notations) is the current resistance matrix
	 stored in the stokes_solver, and \f$\bm{F} \f$ are the forces included by the parameter files.
	 
	 - 1st step:
	 - \f$ \bm{U}^{-} = \bm{A}^{-1}( \bm{X}(t) ) \bm{F} ( \bm{X}(t) )  \f$
	 - \f$ \bm{X}' = \bm{X}(t) + \bm{U}^{-}dt \f$
	 
	 - 2nd step:
	 - \f$ \bm{U}^{+} = \bm{A}^{-1}( \bm{X}^{-} ) \bm{F} ( \bm{X}^{-} )  \f$
	 - \f$ \bm{X}(t + dt) = \bm{X}(t) + \frac{1}{2}(\bm{U}^{+}+\bm{U}^{-})dt =  \bm{X}' + \frac{1}{2}(\bm{U}^{+}-\bm{U}^{-})dt \f$
	 
	 ### Brownian Case
	 
	 This routine implements a two-step algorithm for the dynamics with Brownian motion,
	 initially derived by Fixman \cite fixman_1978.
	 The basis of this algorithm is exposed in \cite Ball_1997 and \cite banchio_accelerated_2003.
	 
	 The equation of motion is:
	 \f$ \bm{A}(\bm{U}-\bm{U}_{\infty}) = \bm{F} + \bm{F}_\mathrm{B} + kT \bm{A} \nabla \bm{A}^{-1} \f$
	 with
	 \f$ \langle \bm{F}_\mathrm{B} \rangle = 0\f$
	 and
	 \f$\langle \bm{F}_\mathrm{B} \bm{F}_\mathrm{B}\rangle = \frac{2kT}{dt} \bm{A}\f$,
	 and \f$\bm{F} \f$ are all the other forces included by the parameter files.
	 
	 Reminding that we obtain the Cholesky decomposition \f$ \bm{A} = \bm{L} \bm{L}^T \f$ in the stokes_solver,
	 we obtain X_B in a 2-step algorithm [ B & M 1997 ]:
		- 1st step:
	 + generate a Brownian force \f$ \bm{F}_\mathrm{B}= \sqrt{\frac{2}{dt}} \bm{L} \psi \f$ with \f$ \langle \psi \rangle = 0 \f$ and \f$ \langle \psi \psi \rangle = 1\f$
	 + \f$ \bm{U}^{-} = \bm{A}^{-1}( \bm{X}(t) )( \bm{F}_\mathrm{B} + \bm{F} ( \bm{X}(t) ) )  \f$
	 + \f$ \bm{X}' = \bm{X}(t) + \bm{U}^{-}dt \f$
		- 2nd step:
	 + \f$ \bm{U}^{+} = \bm{A}^{-1}( \bm{X}^{-} ) ( \bm{F}_\mathrm{B} + \bm{F} ( \bm{X}^{-} ) )  \f$ (\b same \f$\bm{F}_\mathrm{B}\f$ as in the first step)
	 + \f$ \bm{X}(t + dt) = \bm{X}(t) + \frac{1}{2}(\bm{U}^{+}+\bm{U}^{-})dt =  \bm{X}' + \frac{1}{2}(\bm{U}^{+}-\bm{U}^{-})dt \f$
	 
	 */
	
	/* predictor */
	in_predictor = true;
	in_corrector = false;
	setContactForceToParticle();
	setRepulsiveForceToParticle();
	computeVelocities(calc_stress);
	if (calc_stress) {
		calcStressPerParticle();
	}
	timeStepMovePredictor();
	/* corrector */
	in_predictor = false;
	in_corrector = true;
	setContactForceToParticle();
	setRepulsiveForceToParticle();
	computeVelocities(calc_stress);
	if (calc_stress) {
		calcStressPerParticle();
	}
	if (lowPeclet) {
		calcStress();
	}
	timeStepMoveCorrector();

	// try to adapt dt
	if (max_velocity > max_sliding_velocity) {
		dt = disp_max/max_velocity;
	} else {
		dt = disp_max/max_sliding_velocity;
	}
}

void System::timeStepMove()
{
	/**
	 \brief Moves particle positions according to previously computed velocities, Euler method step.
	 */
	
	/* Changing dt for every timestep
	 * So far, this is only in Euler method.
	 */
	if (max_velocity > max_sliding_velocity) {
		dt = disp_max/max_velocity;
	} else {
		dt = disp_max/max_sliding_velocity;
	}
	time += dt;

	/* evolve PBC */
	double strain_increment = 0;
	if(stress_controlled){ // to be fixed: the stress_controlled case should not behave differently
		strain_increment = dt*abs(shear_rate)/abs(dimensionless_number); // @@@ << ???
	} else {
		strain_increment = dt*abs(shear_rate); // the dimensionless_number should not be here I think. \delta\gamma = \tilde{\delta t} * \tilde{\dot\gamma} is all we need here.
	}
	timeStepBoxing(strain_increment);
	
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

void System::timeStepMovePredictor()
{
	/**
	 \brief Moves particle positions according to previously computed velocities, predictor step.
	 */
	if (!brownian) { // adaptative time-step for non-Brownian cases
		//dt = disp_max/max_velocity;
		if (max_velocity > max_sliding_velocity) {
			dt = disp_max/max_velocity;
		} else {
			dt = disp_max/max_sliding_velocity;
		}
	}
	time += dt;

	/* evolve PBC */
	/* The periodic boundary condition is updated in predictor.
	 * It must not be updated in corrector.
	 */
	double strain_increment = dt*abs(shear_rate);
	timeStepBoxing(strain_increment);

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

void System::timeStepMoveCorrector()
{
	/**
	 \brief Moves particle positions according to previously computed velocities, corrector step.
	 */
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

void System::timeEvolution(double time_output_data)
{
	/**
	 \brief Main time evolution routine. Evolves the system until strain_output_data or time_output_data if time_output_data>0.
	 
	 This method essentially loops the appropriate one time step method method, according to the Euler vs predictor-corrector or strain rate vs stress controlled choices.
	 */
	static bool firsttime = true;
	in_predictor = false;
	in_corrector = false;
	if (firsttime) {
		checkNewInteraction();
		updateInteractions();
		firsttime = false;
	}
	bool calc_stress = false;
	if (lowPeclet) {
		calc_stress = true;
	}
	// if (time_output_data == 0) {
	// 	/* integrate until strain_next - 1 time step */
	// 	while (shear_strain < strain_output_data-dt*shear_rate) {
	// 		(this->*timeEvolutionDt)(calc_stress); // no stress computation except at low Peclet
	// 	};
	// 	(this->*timeEvolutionDt)(true); // last time step, compute the stress
	// } else {
	while (time < time_output_data-dt) { // integrate until strain_next
		(this->*timeEvolutionDt)(calc_stress); // no stress computation except at low Peclet
	};
	(this->*timeEvolutionDt)(true); // last time step, compute the stress
	//	}
	if (p.auto_determine_knkt && shear_strain>p.start_adjust){
		adjustContactModelParameters();
	}
}

void System::checkNewInteraction()
{
	/**
	 \brief Checks if there are new pairs of interacting particles. If so, creates and sets up the corresponding Interaction objects.
	 
	 To be called after particle moved.
	 */
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
					double range = (this->*calcInteractionRange)(i, *it);
					double sq_dist_lim = range*range;
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
						interaction[interaction_new].activate(i, *it, range);
					}
				}
			}
		}
	}
}

void System::updateInteractions()
{
	/**
	 \brief Updates the state of active interactions.
	 
	 To be called after particle moved.
	 Note that this routine does not look for new interactions (this is done by System::checkNewInteraction), it only updates already known active interactions.
	 It however desactivate interactions removes interactions that became inactive (ie when the distance between particles gets larger than the interaction range).
	 
	 */
	double sq_max_sliding_velocity = 0;
	/*
	 * In the simulation, the unit of force is proportional to the shear rate.
	 * In the dimensionless simulation,
	 * the cohesive force
	 */
	if (cohesion) {
		dimensionless_cohesive_force = cohesive_force/abs(dimensionless_number);
	}
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active()) {
			bool deactivated = false;
			interaction[k].updateState(deactivated);
			if (deactivated) {
				deactivated_interaction.push(k);
			}
			if (interaction[k].is_contact()) {
				double sq_sliding_velocity = interaction[k].relative_surface_velocity.sq_norm();
				if (sq_sliding_velocity > sq_max_sliding_velocity) {
					sq_max_sliding_velocity = sq_sliding_velocity;
				}
			}
		}
	}
	max_sliding_velocity = sqrt(sq_max_sliding_velocity);
}

void System::stressReset()
{
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

void System::buildHydroTerms(bool build_res_mat, bool build_force_GE)
{
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
void System::buildLubricationTerms_squeeze(bool mat, bool rhs)
{
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

void System::buildLubricationTerms_squeeze_tangential(bool mat, bool rhs)
{
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

void System::generateBrownianForces()
{
	/**
	 \brief Generates a Brownian force realization and sets is as the RHS of the stokes_solver.
	 
	 The generated Brownian force \f$F_B\f$ satisfies
	 \f$ \langle F_\mathrm{B} \rangle = 0\f$
	 and
	 \f$\langle F_\mathrm{B} F_\mathrm{B}\rangle = \frac{2kT}{dt} A\f$
	 where \f$A \f$ (in Jeffrey notations \cite
	 jeffrey_calculation_1992, that is \f$R_{\mathrm{FU}}\f$ in Bossis and Brady
	 \cite brady_stokesian_1988 notations) is the current resistance matrix
	 stored in the stokes_solver.
	 
	 Note that it \b sets the rhs of the solver as \f$ rhs = F_B \f$.
	 
	 \f$ F_B\f$ is also stored in sys->brownian_force.
	 */
	double sqrt_2_dt_amp = sqrt(2/dt)*amplitudes.sqrt_temperature;
	for (int i=0; i<linalg_size; i++) {
		brownian_force[i] = sqrt_2_dt_amp*GRANDOM; // random vector A
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

void System::setContactForceToParticle()
{
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

void System::setRepulsiveForceToParticle()
{
	if (repulsiveforce) {
		for (int i=0; i<np; i++) {
			repulsive_force[i].reset();
		}
		for (int k=0; k<nb_interaction; k++) {
			if (interaction[k].is_active()) {
				interaction[k].repulsion.addUpForce();
			}
		}
	}
}

void System::buildContactTerms(bool set_or_add)
{
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

void System::buildRepulsiveForceTerms(bool set_or_add)
{
	// sets or adds ( set_or_add = t or f resp) repulsive forces to the rhs of the stokes_solver.
	if (set_or_add) {
		for (int i=0; i<np; i++) {
			stokes_solver.setRHSForce(i, repulsive_force[i]);
		}
		stokes_solver.resetRHStorque();
	} else {
		for (int i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, repulsive_force[i]);
		}
	}
}

void System::computeVelocitiesStressControlled()
{
		double shearstress_rep = 0;
		// Compute the stress contributions
		buildHydroTerms(true, true); // build matrix and rhs force GE
		stokes_solver.solve(vel_hydro, ang_vel_hydro); // get V_H
		buildContactTerms(true); // set rhs = F_C
		stokes_solver.solve(vel_contact, ang_vel_contact); // get V_C

		// Back out the shear rate
		dimensionless_number = 1; // To obtain normalized stress from repulsive force.
		amplitudes.repulsion = 1/abs(dimensionless_number);

		if (repulsiveforce) {
			buildRepulsiveForceTerms(true); // set rhs = F_repulsive
			stokes_solver.solve(vel_repulsive, ang_vel_repulsive); // get V_repulsive
		}
		calcStressPerParticle();
		calcStress();
		if (p.unscaled_contactmodel) {
			double shearstress_con = total_contact_stressXF_normal.getStressXZ() \
			+total_contact_stressXF_tan.getStressXZ()+total_contact_stressGU.getStressXZ();
			double shear_rate_numerator = target_stress-shearstress_con;
			if (repulsiveforce) {
				shearstress_rep = total_repulsive_stressXF.getStressXZ()+total_repulsive_stressGU.getStressXZ();
				shear_rate_numerator -= shearstress_rep;
			}
			double shearstress_hyd = einstein_stress+total_hydro_stress.getStressXZ();

			dimensionless_number = shear_rate_numerator/shearstress_hyd;
			if (shear_strain < init_strain_shear_rate_limit) {
				if (dimensionless_number > init_shear_rate_limit) {
					dimensionless_number = init_shear_rate_limit;
				}
			}
			amplitudes.repulsion = 1/abs(dimensionless_number);
			if (repulsiveforce) {
				for (int i=0; i<np; i++) {
					vel_repulsive[i] /= dimensionless_number;
					ang_vel_repulsive[i] /= dimensionless_number;
					vel_contact[i] /= dimensionless_number;
					ang_vel_contact[i] /= dimensionless_number;
				}
				for (int i=0; i<np; i++) {
					na_velocity[i] = vel_hydro[i]+vel_contact[i]+vel_repulsive[i];
					na_ang_velocity[i] = ang_vel_hydro[i]+ang_vel_contact[i]+ang_vel_repulsive[i];
				}
			} else {
				/*
				 * Contact velocity remains the same direction for shear rate < 0.
				 * Velocity from strain should become opposite directino.
				 * To reduce the number of calculation step, the overall sign of velocities
				 * will be invesed later.
				 */
				for (int i=0; i<np; i++) {
					vel_contact[i] /= dimensionless_number;
					ang_vel_contact[i] /= dimensionless_number;
				}
				for (int i=0; i<np; i++) {
					na_velocity[i] = vel_hydro[i]+vel_contact[i];
					na_ang_velocity[i] = ang_vel_hydro[i]+ang_vel_contact[i];
				}
			}
		} else {
			cerr << "not implemented yet" << endl;
			exit(1);
			double shearstress_con = total_contact_stressXF_normal.getStressXZ() \
			+total_contact_stressXF_tan.getStressXZ()+total_contact_stressGU.getStressXZ();
			double shear_rate_numerator = target_stress-shearstress_con;
			if (repulsiveforce) {
				shearstress_rep = total_repulsive_stressXF.getStressXZ()+total_repulsive_stressGU.getStressXZ();
				shear_rate_numerator -= shearstress_rep;
			}
			double shearstress_hyd = einstein_stress+total_hydro_stress.getStressXZ();
			dimensionless_number = (target_stress-shearstress_rep)/(shearstress_hyd+shearstress_con);
			for (int i=0; i<np; i++) {
				vel_repulsive[i] /= dimensionless_number;
				ang_vel_repulsive[i] /= dimensionless_number;
			}
		}
}
	

void System::computeVelocities(bool divided_velocities)
{ // this function is slowly becoming a mess. We should refactor to restore readability.
	stokes_solver.resetRHS();

	if (stress_controlled) {
		computeVelocitiesStressControlled();
	} else {
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
	}
	/*
	 * The max velocity is used to find dt from max displacement
	 * at each time step.
	 */
	if (in_predictor) {
		double sq_max_na_velocity = 0;
		double sq_na_velocity, sq_na_ang_velocity;
		for (int i=0; i<np; i++) {
			sq_na_velocity = na_velocity[i].sq_norm();
			if (sq_max_na_velocity < sq_na_velocity) {
				sq_max_na_velocity = sq_na_velocity;
			}
			sq_na_ang_velocity = na_ang_velocity[i].sq_norm()*radius_squared[i];
			if (sq_max_na_velocity < sq_na_ang_velocity) {
				sq_max_na_velocity = sq_na_ang_velocity;
			}
		}
		max_velocity = sqrt(sq_max_na_velocity);
	}
	if (!zero_shear) {
		/* Background flow should be diffrent direction
		 * if shear rate is negative.
		 * But, the sign of velocity will be invesed later.
		 */
		for (int i=0; i<np; i++) {
			velocity[i] = na_velocity[i];
			ang_velocity[i] = na_ang_velocity[i];
			velocity[i].x += position[i].z;
			ang_velocity[i].y += 0.5;
		}
	} else {
		for (int i=0; i<np; i++) {
			velocity[i] = na_velocity[i];
			ang_velocity[i] = na_ang_velocity[i];
		}
	}
	if (dimensionless_number < 0) {
		shear_direction = -1;
		vel_difference = -lz;
		for (int i=0; i<np; i++) {
			velocity[i] *= -1;
			ang_velocity[i] *= -1;
		}
	} else {
		shear_direction = 1;
		vel_difference = lz;
	}
	stokes_solver.solvingIsDone();
}

void System::displacement(int i, const vec3d &dr)
{
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
int System::periodize(vec3d &pos)
{
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
void System::periodize_diff(vec3d &pos_diff, int &zshift)
{
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

void System::evaluateMaxContactVelocity()
{
	max_contact_velo_tan = 0;
	max_contact_velo_normal = 0;
	max_relative_velocity = 0;
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

double System::evaluateMaxVelocity()
{
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

double System::evaluateMaxAngVelocity()
{
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

double System::evaluateMinGap()
{
	double _min_reduced_gap = 100000;
	for (int k=0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].get_reduced_gap() < _min_reduced_gap) {
			_min_reduced_gap = interaction[k].get_reduced_gap();
			// unsigned short i, j;
			// interaction[k].get_par_num(i,j);
			// cout << i << " " << j << " " << interaction[k].get_a0() << " " << interaction[k].get_a1() << " " << interaction[k].get_reduced_gap() << endl;
		}
	}
	return _min_reduced_gap;
}

double System::evaluateMaxDispTan()
{
	double _max_disp_tan = 0;
	for (int k= 0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].contact.disp_tan.norm() > _max_disp_tan) {
			_max_disp_tan = interaction[k].contact.disp_tan.norm();
		}
	}
	return _max_disp_tan;
}

double System::evaluateMaxDispRolling()
{
	double _max_disp_rolling = 0;
	for (int k= 0; k<nb_interaction; k++) {
		if (interaction[k].is_active() &&
			interaction[k].contact.disp_rolling.norm() > _max_disp_rolling) {
			_max_disp_rolling = interaction[k].contact.disp_rolling.norm();
		}
	}
	return _max_disp_rolling;
	
}

double System::evaluateMaxFcNormal()
{
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

double System::evaluateMaxFcTangential()
{
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

void System::countNumberOfContact()
{
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

void System::analyzeState()
{
	//	max_velocity = evaluateMaxVelocity();
	max_ang_velocity = evaluateMaxAngVelocity();
	evaluateMaxContactVelocity();
	min_reduced_gap = evaluateMinGap();
	max_disp_tan = evaluateMaxDispTan();
	if (rolling_friction) {
		max_disp_rolling = evaluateMaxDispRolling();
	}
	max_fc_normal = evaluateMaxFcNormal();
	max_fc_tan = evaluateMaxFcTangential();
	countNumberOfContact();
}

void System::setSystemVolume(double depth)
{
	if (twodimension) {
		system_volume = lx*lz*depth;
		cerr << "lx = " << lx << " lz = " << lz << " ly = "  << depth << endl;
	} else {
		system_volume = lx*ly*lz;
		cerr << "lx = " << lx << " lz = " << lz << " ly = "  << ly << endl;
	}
}

void System::adjustContactModelParameters()
{
	/**
	 * This method tries to determine
	 * the values of the contact parameters kn, kt
	 * required to reach given maximum normal / tangential spring stretches
	 * (= overlap / tangential displacement ) in steady state.
	 *
	 * The algorithm is a simple adaptative scheme trying to resp. increase or decrease
	 * the spring constants when the stretches are too large or too small.
	 * It works reasonably well away from transitions.
	 * Close to discontinuous shear thickening, one cannot expect satisfying results,
	 * as gigantic stress fluctuations confuse the algorithm.
	 *
	 * The stretch values used are average of maximal values
	 * weighted with an exponential memory kernel:
	 *    \f$k_{avg}(\gamma_n) = C_n \sum_{i=1}^n k(\gamma_i)e^{(\gamma_n-\gamma_i)/\gamma_{avg}} \f$
	 * where \f$\gamma_{avg}\f$ is the user-defined parameter ParameterSet::memory_strain_avg.
	 *
	 * The kn and kt are bounded by user-defined parameters ParameterSet::min_kn, ParameterSet::max_kn, ParameterSet::min_kt, ParameterSet::max_kn.
	 *
	 * The target stretches are given by ParameterSet::overlap_target and ParameterSet::disp_tan_target.
	 *
	 * Additionally, this routine estimates the time step dt.
	 */
	
	analyzeState();

	double overlap = -min_reduced_gap;
	overlap_avg->update(overlap, shear_strain);
	max_disp_tan_avg->update(max_disp_tan, shear_strain);
	kn_avg->update(kn, shear_strain);
	kt_avg->update(kt, shear_strain);
	
	static double previous_shear_strain = 0;
	double deltagamma = (shear_strain-previous_shear_strain);
	double kn_target = kn_avg->get()*overlap_avg->get()/p.overlap_target;
	double dkn = (kn_target-kn)*deltagamma/p.memory_strain_k;
	
	kn += dkn;
	if (kn < p.min_kn) {
		kn = p.min_kn;
	}
	if (kn > p.max_kn) {
		kn = p.max_kn;
	}
	double kt_target = kt_avg->get()*max_disp_tan_avg->get()/p.disp_tan_target;
	double dkt = (kt_target-kt)*deltagamma/p.memory_strain_k;
	kt += dkt;
	if (kt < p.min_kt) {
		kt = p.min_kt;
	}
	if (kt > p.max_kt) {
		kt = p.max_kt;
	}
	if (max_velocity > max_sliding_velocity) {
		dt = disp_max/max_velocity;
	} else {
		dt = disp_max/max_sliding_velocity;
	}
	previous_shear_strain = shear_strain;
}

void System::calcLubricationForce()
{
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

double System::calcLubricationRange(const int& i, const int& j)
{
	double rad_ratio = radius[i]/radius[j];
	if (rad_ratio < 2 && rad_ratio > 0.5) {
		return (2+lub_max_gap)*0.5*(radius[i]+radius[j]);
	} else {
		double minradius = (radius[i]<radius[j] ? radius[i] : radius[j]);
		return radius[i]+radius[j]+lub_max_gap*minradius;
	}
}
