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
#include <stdexcept>
#include <assert.h>
#include "SystemHelperFunctions.h"
#include "global.h"
#include "States.h"
#ifndef USE_DSFMT
#include "MersenneTwister.h"
#endif

extern volatile sig_atomic_t sig_caught;

#ifndef USE_DSFMT
#define GRANDOM ( r_gen->randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.
#endif
#ifdef USE_DSFMT
#define GRANDOM  ( sqrt( -2.0 * log( 1.0 - dsfmt_genrand_open_open(&r_gen) ) ) * cos(2.0 * 3.14159265358979323846264338328 * dsfmt_genrand_close_open(&r_gen) ) ) // RNG gaussian with mean 0. and variance 1.
#endif
using namespace std;

#ifdef USE_DSFMT
inline unsigned long
wagnerhash(time_t t, clock_t c)
{
	/**
		\brief Utility function to start up the DSFMT RNG with a nice seed.

	 From MersenneTwister v1.0 by Richard J. Wagner
	 comments below are from the original code.

	 Get a unsigned long from t and c
	 Better than unsigned long(x) in case x is floating point in [0,1]
	 Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)
	*/

	static unsigned long differ = 0; // guarantee time-based seeds will change
	unsigned long h1 = 0;
	unsigned char *pp = (unsigned char *) &t;
	for (size_t i=0; i<sizeof(t); ++i){
		h1 *= UCHAR_MAX + 2U;
		h1 += pp[i];
	}
	unsigned long h2 = 0;
	pp = (unsigned char *) &c;
	for (size_t j=0; j<sizeof(c); ++j) {
		h2 *= UCHAR_MAX + 2U;
		h2 += pp[j];
	}
	return (h1 + differ++)^h2;
}
#endif


System::System(Parameters::ParameterSet& ps,
	             list <Event>& ev,
					     State::BasicCheckpoint chkp):
pairwise_resistance_changed(true),
clk(chkp.clock),
shear_rate(0),
omega_inf(0),
events(ev),
p(ps),
brownian(false),
friction(false),
rolling_friction(false),
repulsiveforce(false),
delayed_adhesion(false),
cohesion(false),
critical_load(false),
brownian_dominated(false),
twodimension(false),
control(Parameters::ControlVariable::rate),
zero_shear(false),
wall_rheology(false),
mobile_fixed(false),
couette_stress(false),
dt(0),
avg_dt(0),
shear_disp(0),
target_stress(0),
init_strain_shear_rate_limit(0),
init_shear_rate_limit(999),
vel_difference(0),
z_top(-1),
eventLookUp(NULL)
{
	max_sliding_velocity = 0;
	lx = 0;
	ly = 0;
	lz = 0;
	shear_strain = 0;
}

System::~System()
{
	interaction.clear();
	interaction_list.clear();
}

void System::allocateRessources()
{
	if (np <= 0) {
		throw runtime_error("System::allocateRessources() :  np is 0 or negative, cannot allocate this.");
	}
	radius_cubed.resize(np);
	radius_squared.resize(np);
	if (!pairwise_resistance) {
		// Stokes-drag simulation
		stokesdrag_coeff_f.resize(np);
		stokesdrag_coeff_f_sqrt.resize(np);
		stokesdrag_coeff_t.resize(np);
		stokesdrag_coeff_t_sqrt.resize(np);
	}
	// Configuration
	if (twodimension) {
		angle.resize(np, 0);
	}
	// Velocity
	velocity.resize(np);
	for (auto &v: velocity) {
		v.reset();
	}
	ang_velocity.resize(np);
	for (auto &v: ang_velocity) {
		v.reset();
	}
	na_velocity.resize(np);
	na_ang_velocity.resize(np);
	if (p.integration_method == 1) {
		velocity_predictor.resize(np);
		ang_velocity_predictor.resize(np);
	}
	u_inf.resize(np);
	for (auto &v: u_inf) {
		v.reset();
	}
	na_disp.resize(np);
	if (mobile_fixed) {
		non_rate_proportional_wall_force.resize(p.np_fixed);
		non_rate_proportional_wall_torque.resize(p.np_fixed);
		rate_proportional_wall_force.resize(p.np_fixed);
		rate_proportional_wall_torque.resize(p.np_fixed);
	}
	declareForceComponents();
	declareVelocityComponents();
	interaction_list.resize(np);
	interaction_partners.resize(np);
	nb_blocks_ff.resize(p.np_fixed, 0);
	nb_blocks_mm.resize(np-p.np_fixed, 0);
	nb_blocks_mf.resize(np-p.np_fixed, 0);
	// Forces and Stress
	forceResultant.resize(np);
	torqueResultant.resize(np);
	declareStressComponents();
	total_stress_pp.resize(np);
	phi6.resize(np);
}

void System::declareForceComponents()
{
	// Only declare in force components the forces on the rhs of
	// R_FU*(U-U_inf) = RHS
	// These forces will be used to compute the na_velo_components,
	// a set of components that must add up exactly to the total non-affine velocity

	bool torque = true;

	/******* Contact force, spring part ***********/
	if (friction) {
		force_components["contact"] = ForceComponent(np, RATE_INDEPENDENT, torque, &System::setContactForceToParticle);
	} else {
		force_components["contact"] = ForceComponent(np, RATE_INDEPENDENT, !torque, &System::setContactForceToParticle);
	}

	/*********** Hydro force, i.e.  R_FE:E_inf *****************/
	if (!zero_shear) {
		if (p.lubrication_model == "normal") {
			force_components["hydro"] = ForceComponent(np, RATE_PROPORTIONAL, !torque, &System::setHydroForceToParticle_squeeze);
		}
		if (p.lubrication_model == "tangential") {
			force_components["hydro"] = ForceComponent(np, RATE_PROPORTIONAL, torque, &System::setHydroForceToParticle_squeeze_tangential);
		}
		/******* Contact force, dashpot part, U_inf only, i.e. R_FU^{dashpot}.U_inf ***********/
		// note that we keep the torque on, as there is nothing else to prevent tangential motion when in contact
		// (apart from Stokes drag, which can be set arbitrarily small)
		force_components["dashpot"] = ForceComponent(np, RATE_PROPORTIONAL, torque, &System::setDashpotForceToParticle);
	}

	if (repulsiveforce) {
		force_components["repulsion"] = ForceComponent(np, RATE_INDEPENDENT, !torque, &System::setRepulsiveForceToParticle);
	}

	if (brownian) {
		force_components["brownian"] = ForceComponent(np, RATE_INDEPENDENT, torque, &System::setBrownianForceToParticle);
	}

	if (delayed_adhesion) {
		force_components["delayed_adhesion"] = ForceComponent(np, RATE_INDEPENDENT, !torque, &System::setTActAdhesionForceToParticle);
	}

	/********** Force R_FU^{mf}*(U^f-U^f_inf)  *************/
	if (mobile_fixed) {
		// rate proportional with walls, but this can change
		if (p.lubrication_model != "normal") {
			force_components["from_fixed"] = ForceComponent(np, RATE_PROPORTIONAL, !torque, &System::setFixedParticleForceToParticle);
		}
		if (p.lubrication_model != "tangential") {
			force_components["from_fixed"] = ForceComponent(np, RATE_PROPORTIONAL, torque, &System::setFixedParticleForceToParticle);
		}
	}
}

void System::declareVelocityComponents()
{
	// Only declare in na_velo_components a set of components that add up
	// exactly to the non-affine velocity
	for (const auto &fc: force_components) {
		na_velo_components[fc.first] = VelocityComponent(fc.second.force.size(), fc.second.rate_dependence);
	}
}

void System::setInteractions_GenerateInitConfig()
{
	calcInteractionRange = &System::calcLubricationRange;

	shear_strain = {0, 0, 0};
	clk.cumulated_strain = 0;
	shear_disp.reset();
	vel_difference.reset();
	initializeBoxing();
	checkNewInteraction();

	for (unsigned int k=0; k<interaction.size(); k++) {
		interaction[k].label = k;
	}
}

void System::setConfiguration(const vector <vec3d>& initial_positions,
							  const vector <double>& radii)
{
	/**
		\brief Set positions of the particles for initialization.
	 */
	string indent = "  System::\t";
	set_np(initial_positions.size());
	np_mobile = np-p.np_fixed;
	if (np_mobile <= 0) {
		throw runtime_error("np_fixed>=np");
	}
	if (p.np_fixed > 0) {
		mobile_fixed = true;
	}
	position.resize(np);
	radius.resize(np);
	for (int i=0; i<np; i++) {
		position[i] = initial_positions[i];
		radius[i] = radii[i];
	}
	radius_wall_particle = radius[np-1];
	setSystemVolume();

	particle_volume = 0;
	if (twodimension) {
		for (auto r: radius) {
			particle_volume += r*r;
		}
		particle_volume *= M_PI;
	} else {
		for (auto r: radius) {
			particle_volume += r*r*r;
		}
		particle_volume *= 4*M_PI/3.;
	}
}

void System::setFixedVelocities(const vector <vec3d>& vel)
{
	fixed_velocities = vel;
}

void System::setContacts(const vector <struct contact_state>& cs)
{
	/**
		\brief Set a list of contacts with their state variables.

		Used to restart the simulation from a given state.
	 */
	for (const auto& c : cs) {
		bool existing_interaction = false;
		for (auto &inter: interaction) {
			unsigned int p0, p1;
			std::tie(p0, p1) = inter.get_par_num();
			if (p0 == c.p0 && p1 == c.p1) {
				inter.contact.setState(c);
				existing_interaction = true;
			}
		}
		assert(existing_interaction); // otherwise you probably messed up the initialisation
	}
}

vector <struct contact_state> System::getContacts() const
{
	/**
		\brief Get the list of contacts with their state variables.

		Used to output a configuration including contact info. Useful if you want to restart from exact same configuration.
	 */
	vector <struct contact_state> cs;
	for (const auto &inter: interaction) {
		if (inter.contact.is_active()) {
			cs.push_back(inter.contact.getState());
		}
	}
	return cs;
}

void System::setupParametersIntegrator()
{
	string indent = "  System::\t";
	if (p.integration_method == 0) {
		timeEvolutionDt = &System::timeEvolutionEulersMethod;
	} else if (p.integration_method == 1) {
		timeEvolutionDt = &System::timeEvolutionPredictorCorrectorMethod;
	} else {
		ostringstream error_str;
		error_str << indent << "integration_method = " << p.integration_method << endl << indent << "The integration method is not impremented yet." << endl;
		throw runtime_error(error_str.str());
	}
}

void System::setupParametersLubrication()
{
	string indent = "  System::\t";
	lubrication = p.lubrication_model != "none";
	if (p.lub_max_gap < 0) {
		throw runtime_error(indent+"lub_max_gap<0 is forbidden.");
	}
	if (p.lub_reduce_parameter > 1) {
		cout << indent+" p.lub_reduce_parameter>1, log terms in lubrication set to 0." << endl;
	}

	if (p.lubrication_model != "normal" &&
		p.lubrication_model != "none" &&
		p.lubrication_model != "tangential") {
		throw runtime_error(indent+"unknown lubrication_model "+p.lubrication_model+"\n");
	}

	if (p.lubrication_model == "tangential" && p.lub_max_gap >= 1) {
		/* The tangential part of lubrication is approximated as log(1/h).
		 * To keep log(1/h) > 0, h needs to be less than 1.
		 */
		throw runtime_error(indent+"lub_max_gap must be smaller than 1\n");
	}
}

void System::setupParametersContacts()
{
	string indent = "  System::\t";
	if (p.friction_model == 0) {
		cout << indent+"friction model: no friction" << endl;
		p.mu_static = 0;
		friction = false;
	} else if (p.friction_model == 1) {
		cout << indent+"friction model: Coulomb" << endl;
		friction = true;
	} else if (p.friction_model == 2 || p.friction_model == 3) {
		cout << indent+"friction model: Coulomb + Critical Load" << endl;
		friction = true;
	} else if (p.friction_model == 5) {
		cout << indent+"friction_model: Max tangential force" << endl;
		friction = true;
	} else if (p.friction_model == 6) {
		cout << indent+"friction_model: Coulomb law + Max tangential force" << endl;
		friction = true;
	} else {
		throw runtime_error(indent+"Error: unknown friction model\n");
	}
	if (p.mu_dynamic < 0) {
		p.mu_dynamic = p.mu_static;
	}
	if (p.mu_rolling > 0) {
		rolling_friction = true;
		if (friction == false) {
			throw runtime_error(indent+"Error: Rolling friction without sliding friction?\n");
		}
	}
}

void System::setupParameters()
{
	setupParametersIntegrator();
	setupParametersLubrication();
	setupParametersContacts();
	pairwise_resistance = lubrication || p.contact_relaxation_time != 0 || p.contact_relaxation_time_tan != 0;
	if (p.interaction_range == -1) {
		/* If interaction_range is not indicated,
		 * interaction object is created at the lubrication cutoff.
		 */
		calcInteractionRange = &System::calcLubricationRange;
		p.interaction_range = 2+p.lub_max_gap;
	} else {
		calcInteractionRange = &System::calcInteractionRangeDefault;
	}

	if (p.repulsive_length <= 0) {
		repulsiveforce = false;
		p.repulsive_length = 0;
	}
	setShearDirection(p.theta_shear);

	if (p.auto_determine_knkt) {
		kn_avg.setRelaxationTime(p.memory_strain_avg);
		kt_avg.setRelaxationTime(p.memory_strain_avg);
		overlap_avg.setRelaxationTime(p.memory_strain_avg);
		max_disp_tan_avg.setRelaxationTime(p.memory_strain_avg);
	}

	if (brownian) {
		if (brownian_dominated) {
			double stress_avg_relaxation_parameter = 10*p.output.time_interval_output_data.value; // 0 --> no average
			stress_avg.setRelaxationTime(stress_avg_relaxation_parameter);
			rate_prop_shearstress_rate1_ave.setRelaxationTime(p.brownian_relaxation_time);
			rate_indep_shearstress_ave.setRelaxationTime(p.brownian_relaxation_time);
		}
		if (pairwise_resistance && p.integration_method != 1) {
			ostringstream error_str;
			error_str << "Brownian simulation needs to use the Predictor-Corrector method." << endl;
			error_str << "Modify the parameter file." << endl;
			throw runtime_error(error_str.str());
		}
	}
	if (p.simulation_mode == 31) {
		p.sd_coeff = 1e-6;
	}
	if (p.cohesion > 0) {
		cohesion = true;
	}
	if (p.TA_adhesion.adhesion_max_force > 0) {
		delayed_adhesion = true;
	}

}

void System::setupBrownian()
{

#ifdef DEV
	/* In developing and debugging phases,
	 * we give a seed to generate the same series of random number.
	 * DEV is defined as a preprocessor option in the Makefile
	 */
#ifndef USE_DSFMT
	r_gen = new MTRand(17);	cerr << " WARNING : debug mode: hard coded seed is given to the RNG " << endl;
#endif
#ifdef USE_DSFMT
	dsfmt_init_gen_rand(&r_gen, 17);	cerr << " WARNING : debug mode: hard coded seed is given to the RNG " << endl;
#endif
#endif

#ifndef DEV
#ifndef USE_DSFMT
	r_gen = new MTRand;
#endif
#ifdef USE_DSFMT
	dsfmt_init_gen_rand(&r_gen, wagnerhash(std::time(NULL), clock()) ) ; // hash of time and clock trick from MersenneTwister v1.0 by Richard J. Wagner
#endif
#endif
}

template<typename T>
void System::setupGenericConfiguration(T conf, Parameters::ControlVariable control_)
{
	string indent = "  System::\t";
	cout << indent << "Setting up System... " << endl;
	np = conf.position.size();
	np_mobile = np - p.np_fixed;
	control = control_;

	setBoxSize(conf.lx,conf.ly,conf.lz);
	twodimension = ly == 0;

	setupParameters();
	// Memory
	allocateRessources();

	if (brownian) {
		setupBrownian();
	}

	total_num_timesteps = 0;

	angle_output = false;
	if (twodimension) {
		angle_output = true;
	}
	cout << indent << "Setting up System... [ok]" << endl;

	shear_disp = conf.lees_edwards_disp;
	if (p.keep_input_strain) {
		shear_strain = shear_disp/lz;
	} else {
		shear_strain = {0, 0, 0};
	}

	setConfiguration(conf.position, conf.radius);
	if (ext_flow) {
		// extensional flow
		// @@@ Some variables need to be set before initializeBoxing() @@@
		// @@@ This part is very messy. Needs to be rewritten careffuly.
		strain_retrim_interval = 2*log(0.5*(3+sqrt(5))); // every this strain, the main simulation box is retrimmed.
		strain_retrim = strain_retrim_interval; // Setting the first value of strain to retrim.
		double cos_ma = cos(p.magic_angle);
		double sin_ma = sin(p.magic_angle);
		sq_cos_ma = cos_ma*cos_ma;
		sq_sin_ma = sin_ma*sin_ma;
		cos_ma_sin_ma = cos_ma*sin_ma;
		updateH(); // cumulated_strain = 0
	}
	initializeBoxing();
	checkNewInteraction();
	setContacts(conf.contact_states);
	setupSystemPostConfiguration();
}


void System::setupConfiguration(struct base_shear_configuration conf, Parameters::ControlVariable control_)
{
	setupGenericConfiguration(conf, control_);
}

void System::setupConfiguration(struct fixed_velo_configuration conf, Parameters::ControlVariable control_)
{
	p.np_fixed = conf.fixed_velocities.size();
	setupGenericConfiguration(conf, control_);
	setFixedVelocities(conf.fixed_velocities);
}

void System::setupConfiguration(struct circular_couette_configuration conf, Parameters::ControlVariable control_)
{
	p.np_fixed = conf.np_wall1 + conf.np_wall2;
	np_wall1 = conf.np_wall1;
	np_wall2 = conf.np_wall2;
	radius_in = conf.radius_in;
	radius_out = conf.radius_out;

	setupGenericConfiguration(conf, control_);
}

struct base_shear_configuration confConvertBase2Shear(const struct base_configuration &conf,
													  vec3d lees_edwards_disp)
{
	struct base_shear_configuration base_shear;
	base_shear.lx = conf.lx;
	base_shear.ly = conf.ly;
	base_shear.lz = conf.lz;
	base_shear.volume_or_area_fraction = conf.volume_or_area_fraction;
	base_shear.position = conf.position;
	base_shear.radius = conf.radius;
	base_shear.angle = conf.angle;
	base_shear.contact_states = conf.contact_states;

	base_shear.lees_edwards_disp = lees_edwards_disp;
	return base_shear;
}


void System::setupConfiguration(const struct delayed_adhesion_configuration &conf, 
								Parameters::ControlVariable control_)
{
	setupGenericConfiguration(confConvertBase2Shear(conf.base, conf.lees_edwards_disp), control_);
	TActAdhesion::setupInteractions(interaction, conf.adhesion_states, get_time());
}

void System::setupSystemPostConfiguration()
{
	for (int i=0; i<np; i++) {
		radius_squared[i] = pow(radius[i], 2);
		radius_cubed[i] = pow(radius[i], 3);
	}
	if (pairwise_resistance) {
		resistance_matrix_dblock.resize(np);
		for (int i=0; i<np; i++) {
			resetDBlock(resistance_matrix_dblock[i]);
		}
	}
	for (int i=0; i<np; i++) {
		double FUvalue = p.sd_coeff*radius[i];
		double TWvalue = p.sd_coeff*radius_cubed[i]*4.0/3;
		if (!pairwise_resistance) {
			// Stokes drag simulation
			stokesdrag_coeff_f[i] = FUvalue;
			stokesdrag_coeff_f_sqrt[i] = sqrt(FUvalue);
			stokesdrag_coeff_t[i] = TWvalue;
			stokesdrag_coeff_t_sqrt[i] = sqrt(TWvalue);
		} else {
			resistance_matrix_dblock[i].col0[0] = FUvalue;
			resistance_matrix_dblock[i].col1[0] = FUvalue;
			resistance_matrix_dblock[i].col2[0] = FUvalue;
			resistance_matrix_dblock[i].col3[0] = TWvalue;
			resistance_matrix_dblock[i].col4[0] = TWvalue;
			resistance_matrix_dblock[i].col5[0] = TWvalue;
		}
	}
	omega_wheel_in  = 0;
	omega_wheel_out = 0;
	if (p.simulation_mode >= 10 && p.simulation_mode <= 20) {
		origin_of_rotation = {lx_half, 0, lz_half};
		for (int i=np_mobile; i<np; i++) {
			angle[i] = -atan2(position[i].z-origin_of_rotation.z,
							  position[i].x-origin_of_rotation.x);
		}
		double omega_wheel = (radius_out-radius_in)*shear_rate/radius_in;
		if (p.simulation_mode == 11) {
			omega_wheel_in  = 0;
			omega_wheel_out = -omega_wheel;
		} else if (p.simulation_mode == 12) {
			omega_wheel_in  = omega_wheel;
			omega_wheel_out = 0;
		} else if (p.simulation_mode == 13) {
			omega_wheel_in  = 0.5*omega_wheel;
			omega_wheel_out = -0.5*omega_wheel;
		} else if (p.simulation_mode == 10) {
			omega_wheel_in  = 0;
			omega_wheel_out = 0;
		}
		couette_stress = true; // output stress per perticle
	} else if (p.simulation_mode == 51) {
		double omega_wheel = (radius_out-radius_in)*shear_rate/radius_out;
		omega_wheel_out = -omega_wheel;
		omega_wheel_in  = omega_wheel*radius_out/radius_in;
	}
	if (pairwise_resistance) {
		stokes_solver.init(np, np_mobile);
	}
	dt = p.dt;
	if (p.fixed_dt) {
		avg_dt = dt;
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
			range = (this->*calcInteractionRange)(i, j);
			if (range > max_range) {
				max_range = range;
			}
		}
	}
	if (!ext_flow) {
		// simple shear
		boxset.init(max_range, this);
		for (int i=0; i<np; i++) {
			boxset.box(i);
		}
		boxset.update();
	} else {
		// extensional flow
		double dl = max_range;
		int num_x = (int)(lx/dl);
		int num_z = (int)(lz/dl);
		lx_ext_flow = 3*dl*(num_x+1)+2*dl;
		ly_ext_flow = ly;
		lz_ext_flow = 2*dl*(num_z+1)+2*dl;
		vec3d box_origin(dl*(num_x+2), 0, dl*(num_z+2));
		boxset.initExtFlow(max_range, box_origin, this);
		for (int i=0; i<np; i++) {
			boxset.box(i);
		}
		boxset.updateExtFlow();
	}
}

struct base_configuration System::getBaseConfiguration() const
{
	base_configuration conf;
	conf.lx = lx;
	conf.ly = ly;
	conf.lz = lz;
	conf.volume_or_area_fraction = particle_volume/system_volume;

	conf.position = position;
	conf.radius = radius;
	if (twodimension) {
		conf.angle = angle;
	}
	conf.contact_states = getContacts();
	return conf;
}

void System::timeStepBoxing()
{
	/**
		\brief Apply a strain step to the boxing system.
	 */
	if (!zero_shear) {
		double strain_increment = shear_rate*dt;
		clk.cumulated_strain += strain_increment;
		if (!ext_flow) {
			// simple shear flow
			vec3d shear_strain_increment = 2*dot(E_infinity, {0, 0, 1})*dt;
			shear_strain += shear_strain_increment;
			shear_disp += shear_strain_increment*lz;
			int m = (int)(shear_disp.x/lx);
			if (shear_disp.x < 0) {
				m--;
			}
			shear_disp.x = shear_disp.x-m*lx;
			if (!twodimension) {
				m = (int)(shear_disp.y/ly);
				if (shear_disp.y < 0) {
					m--;
				}
				shear_disp.y = shear_disp.y-m*ly;
			}
		}
	} else {
		if (wall_rheology || p.simulation_mode == 31) {
			vec3d strain_increment = 2*dot(E_infinity, {0, 0, 1})*dt;
			clk.cumulated_strain += strain_increment.norm();
			shear_strain += strain_increment;
			angle_wheel += dt*(omega_wheel_in-omega_wheel_out);
		}
	}
	if (ext_flow) {
		updateH(); // update cell axes
	}

	if (!ext_flow) {
		/**** simple shear flow ****/
		boxset.update();
	} else {
		/**** extensional flow ****/
		boxset.updateExtFlow();
	}
}

void System::eventShearJamming()
{
	/**
	 \brief Create an event when the shear rate is negative
	*/
	if (shear_rate < 0) {
		Event ev;
		ev.type = "negative_shear_rate";
		events.push_back(Event(ev));
	}
}

void System::forceResultantInterpaticleForces()
{
	auto &contact_force = force_components["contact"].force;
	for (int i=0; i<np; i++) {
		forceResultant[i] += contact_force[i];
	}
	if (friction) {
		auto &contact_torque = force_components["contact"].torque;
		for (int i=0; i<np; i++) {
			torqueResultant[i] += contact_torque[i];
		}
	}
	if (repulsiveforce) {
		auto &repulsive_force = force_components["repulsion"].force;
		for (int i=0; i<np; i++) {
			forceResultant[i] += repulsive_force[i];
		}
	}
	if (delayed_adhesion) {
		auto &adhesion_force = force_components["repulsion"].force;
		for (int i=0; i<np; i++) {
			forceResultant[i] += adhesion_force[i];
		}
	}
}

void System::wallForces()
{
	if (wall_rheology) {
		double max_total_force = 0;
		double max_total_torque = 0;
		for (int i=0; i<np_mobile; i++) {
			if (max_total_force < forceResultant[i].sq_norm()){
				max_total_force = forceResultant[i].sq_norm();
			}
			if (max_total_torque < torqueResultant[i].sq_norm()){
				max_total_torque = torqueResultant[i].sq_norm();
			}
		}
		cerr << "force balance: " << sqrt(max_total_force) << endl;
		cerr << "torque balance: " << sqrt(max_total_torque) << endl;
		if (p.simulation_mode >= 10 && p.simulation_mode <= 20) {
			int i_np_1 = np_mobile+np_wall1;
			// inner wheel
			// Positions of wall particles are at r =
			force_normal_wall1 = 0;
			double torque_wall1 = 0;
			for (int i=np_mobile; i<i_np_1; i++) {
				vec3d unitvec_out = origin_of_rotation-position[i];
				unitvec_out.y = 0;
				unitvec_out.unitvector();
				force_normal_wall1 += dot(forceResultant[i], unitvec_out);
				vec3d torque_tmp = cross(position[i]-origin_of_rotation, forceResultant[i]);
				torque_wall1 += torque_tmp.y+torqueResultant[i].y;
			}
			force_tang_wall1 = torque_wall1/(radius_in-radius_wall_particle);
			// outer wheel
			force_normal_wall2 = 0;
			double torque_wall2 = 0;
			for (int i=i_np_1; i<np; i++) {
				vec3d unitvec_out = position[i]-origin_of_rotation;
				unitvec_out.y = 0;
				unitvec_out.unitvector();
				force_normal_wall2 += dot(forceResultant[i], unitvec_out);
				vec3d torque_tmp = cross(position[i]-origin_of_rotation, forceResultant[i]);
				torque_wall2 += torque_tmp.y+torqueResultant[i].y;
			}
			force_tang_wall2 = torque_wall2/(radius_out+radius_wall_particle);
			cerr << " normal:" << force_normal_wall1 << ' ' << force_normal_wall2 << endl;
			cerr << " tangential:" << force_tang_wall1 << ' ' << force_tang_wall2 << ' ' << torque_wall1 << ' ' << torque_wall2 << endl;
		} else if (p.simulation_mode > 40) {
			int i_np_1 = np_mobile+np_wall1;
			// bottom wall
			force_tang_wall1 = 0;
			force_normal_wall1 = 0;
			for (int i=np_mobile; i<i_np_1; i++) { // bottom
					force_tang_wall1   += forceResultant[i].x;
					force_normal_wall1 += forceResultant[i].z;
			}
			// top wall
			force_tang_wall2 = 0;
			force_normal_wall2 = 0;
			for (int i=i_np_1; i<np; i++) {
				force_tang_wall2   += forceResultant[i].x;
				force_normal_wall2 += forceResultant[i].z;
			}
			cerr << "Ft " <<   force_tang_wall1 << ' ' <<   force_tang_wall2 << endl;
			cerr << "Fn " << force_normal_wall1 << ' ' << force_normal_wall2 << endl;
		}
	}
}

void System::forceResultantReset()
{
	for (int i=0; i<np; i++) {
		forceResultant[i].reset();
		torqueResultant[i].reset();
	}
}

void System::checkForceBalance()
{
	// 1st way: does not work: forceResultand != 0
	forceResultantReset();
	forceResultantInterpaticleForces();
	unsigned int i, j;
	for (const auto &inter: interaction) {
		std::tie(i, j) = inter.get_par_num();
		vec3d lubforce_i = inter.lubrication.getTotalForce();
		forceResultant[i] += lubforce_i;
		forceResultant[j] -= lubforce_i;
	}
	// 2nd way: works
	forceResultantReset();
	forceResultantInterpaticleForces();
	forceResultantLubricationForce();
}

void System::timeEvolutionEulersMethod(bool calc_stress,
									   double time_end,
									   double strain_end)
{
	/**
	 \brief One full time step, Euler's method.

	 This method is never used when running a Brownian simulation.
	 */
	in_predictor = true;
	in_corrector = true;
	if (!pairwise_resistance) {
		computeVelocitiesStokesDrag();
	} else {
		computeVelocities(calc_stress);
	}
	if (wall_rheology && calc_stress) {
		forceResultantReset();
		forceResultantInterpaticleForces();
	}
	if (calc_stress) {
		if (wall_rheology) {
			wallForces();
		}
		calcStressPerParticle();
		if (wall_rheology) {
			calcStress();
		}
		if (!p.output.out_particle_stress.empty() || couette_stress) {
			calcTotalStressPerParticle();
		}
	}
	timeStepMove(time_end, strain_end);
	for (int i=0; i<np; i++) {
		na_disp[i] += na_velocity[i]*dt;
	}
	if (eventLookUp != NULL) {
		(this->*eventLookUp)();
	}
}

/****************************************************************************************************
 ******************************************** Mid-Point Scheme ***************************************
 ****************************************************************************************************/

void System::timeEvolutionPredictorCorrectorMethod(bool calc_stress,
												   double time_end,
												   double strain_end)
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
	 - \f$ \bm{U}^{+} = \bm{A}^{-1}( \bm{X}' ) \bm{F} ( \bm{X}' )  \f$
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

	if (pairwise_resistance) {
		computeVelocities(calc_stress); // divided velocities for stress calculation
	} else {
		computeVelocitiesStokesDrag();
	}

	if (wall_rheology && calc_stress) {
		forceResultantReset();
		forceResultantInterpaticleForces();
	}

	if (calc_stress) {
		if (wall_rheology) {
			wallForces();
		}
		calcStressPerParticle(); // stress compornents
	}
	timeStepMovePredictor(time_end, strain_end);

	/* corrector */
	in_predictor = false;
	in_corrector = true;
	if (pairwise_resistance) {
		computeVelocities(calc_stress);
	} else {
		computeVelocitiesStokesDrag();
	}
	if (calc_stress) {
		calcStressPerParticle(); // stress compornents
		if (wall_rheology || brownian_dominated) {
			calcStress();
		}
		if (!p.output.out_particle_stress.empty() || couette_stress) {
			calcTotalStressPerParticle();
		}
	}
	timeStepMoveCorrector();
	for (int i=0; i<np; i++) {
		na_disp[i] += na_velocity[i]*dt;
	}
}

void System::calcOrderParameter()
{
	if (twodimension == false) {
		cerr << "The bond-orientation order parameter can be calculated only 2D simulation." << endl;
		throw runtime_error("  The bond-orientation order parameter is only for 2D simulation. out_bond_order_parameter6 must be false for 3D simulation\n");
	}
	for (auto phi6_ : phi6) {
		phi6_.real(0);
		phi6_.imag(0);
	}
	vector<int> n_neighobr(np, 0);
	for (auto &inter: interaction) {
		double phase = inter.nvec.angle_elevation_xz()+M_PI; // nvec = p1 - p0
		complex<double> bond_phase1(cos(6*phase), sin(6*phase));
		/*
		 complex<double> bond_phase2(cos(6*(phase+M_PI)), sin(6*(phase+M_PI)));
		 bond_phase2 = bond_phase1
		 */
		phi6[inter.get_p0()] += bond_phase1;
		phi6[inter.get_p1()] += bond_phase1;
		n_neighobr[inter.get_p0()] ++;
		n_neighobr[inter.get_p1()] ++;
	}
	for (int i=0; i<np; i++) {
		phi6[i] *= 1.0/n_neighobr[i];
	}
}

void System::adaptTimeStep()
{
	/**
	 \brief Adapt the time step so that the maximum relative displacement is p.disp_max .
	 */
	if (max_velocity > 0 || max_sliding_velocity > 0) { // small density system can have na_velocity=0
		if (max_velocity > max_sliding_velocity) {
			dt = p.disp_max/max_velocity;
		} else {
			dt = p.disp_max/max_sliding_velocity;
		}
	} else {
		dt = p.disp_max/shear_rate;
	}
	if (dt*shear_rate > p.disp_max) { // cases where na_velocity < \dotgamma*radius
		dt = p.disp_max/shear_rate;
	}
	if (p.dt_max > 0) {
		if (dt > p.dt_max) {
			dt = p.dt_max;
		}
	}
	if (p.dt_min > 0) {
		if (dt < p.dt_min) {
			dt = p.dt_min;
		}
	}
}

void System::adaptTimeStep(double time_end, double strain_end)
{
	/**
	 \brief Adapt the time step so that (a) the maximum relative displacement is p.disp_max, and (b) time or strain does not get passed the end value.
	 */
	adaptTimeStep();
	// To stop exactly at t == time_end or strain == strain_end,
	// whatever comes first
	if (strain_end >= 0) {
		if (fabs(dt*shear_rate) > strain_end-clk.cumulated_strain) {
			dt = fabs((strain_end-clk.cumulated_strain)/shear_rate);
		}
	}
	if (time_end >= 0) {
		if (get_time()+dt > time_end) {
			dt = time_end-get_time();
		}
	}
}

void System::retrimProcess()
{
	cerr << "retrim" << endl;
	strain_retrim += strain_retrim_interval;
	updateH();
	for (int i=0; i<np; i++) {
		retrim(position[i]);
		boxset.box(i);
		velocity[i] -= u_inf[i];
		u_inf[i] = grad_u*position[i];
		velocity[i] +=  u_inf[i];
	}
	boxset.updateExtFlow();
}

void System::timeStepMove(double time_end, double strain_end)
{
	/**
	 \brief Moves particle positions according to previously computed velocities, Euler method step.
	 */
	/* @@@@ NOTE
	 * shear_rate = 1 for both simple shear and extensional flow
	 * dot_epsion = shear_rate / 2 is always true.
	 * clk.cumulated_strain = shear_rate * t for both simple shear and extensional flow.
	 */
	/* Adapt dt to get desired p.disp_max	 */
	if (!p.fixed_dt) {
		adaptTimeStep(time_end, strain_end);
	}
	clk.time_ += dt;
	total_num_timesteps ++;
	/* evolve PBC */
	timeStepBoxing();

	/* move particles */
	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (angle_output) {
		for (int i=0; i<np; i++) {
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	if (retrim_ext_flow) {
		cerr << "clk.cumulated_strain = " << clk.cumulated_strain << endl;
		retrimProcess();
	}
	checkNewInteraction();
	updateInteractions();
}

void System::timeStepMovePredictor(double time_end, double strain_end)
{
	/**
	 \brief Moves particle positions according to previously computed velocities, predictor step.
	 */
	if (!brownian) { // adaptative time-step for non-Brownian cases
 		if (!p.fixed_dt) {
 			adaptTimeStep(time_end, strain_end);
 		}
 	}
	clk.time_ += dt;
	total_num_timesteps ++;
	/* evolve PBC
	 * The periodic boundary condition is updated in predictor.
	 * It must not be updated in corrector.
	 */
	timeStepBoxing();
	for (int i=0; i<np; i++) {
		displacement(i, velocity[i]*dt);
	}
	if (angle_output) {
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
		velocity[i] = 0.5*(velocity[i]+velocity_predictor[i]); // real velocity, in predictor and in corrector
		ang_velocity[i] = 0.5*(ang_velocity[i]+ang_velocity_predictor[i]);
	}
	for (int i=0; i<np; i++) {
		displacement(i, (velocity[i]-velocity_predictor[i])*dt);
	}
	if (angle_output) {
		for (int i=0; i<np; i++) {
			angle[i] += (ang_velocity[i].y-ang_velocity_predictor[i].y)*dt; // no cross_shear in 2d
		}
	}
	if (retrim_ext_flow) {
		cerr << "clk.cumulated_strain = " << clk.cumulated_strain << endl;
		retrimProcess();
	}
	checkNewInteraction();
	updateInteractions();
}

bool System::keepRunning(double time_end, double strain_end)
{
	if (clk.cumulated_strain > strain_end-1e-8 && strain_end >= 0) {
		return false;
	}
	if (get_time() > time_end-1e-8 && time_end>=0) {
		return false;
	}
	if (!events.empty()) {
		return false;
	}
	return true;
}

void System::timeEvolution(double time_end, double strain_end)
{
	/**
	 \brief Main time evolution routine: evolves the system until time_end

	 This method essentially loops the appropriate one time step
	 method method, according to the Euler vs predictor-corrector or
	 strain rate vs stress controlled choices. On the last time step,
	 the stress is computed.
	 (In the case of low Peclet simulations, the stress is computed at every time step.)
	 r
	 \param time_end Time to reach.
	 */
	static bool firsttime = true;
	in_predictor = false;
	in_corrector = false;
	if (firsttime) {
		double dt_bak = dt; // to avoid stretching contact spring
		dt = 0;
		checkNewInteraction();
		in_predictor = true;
		updateInteractions();
		in_predictor = false;
		dt = dt_bak;
		firsttime = false;
	}
	bool calc_stress = false;
	if (brownian_dominated) {
		calc_stress = true;
	}
	retrim_ext_flow = false;
	avg_dt = 0;
	avg_dt_nb = 0;
	double dt_bak = -1;
	double loop_time_adjust = 0;
	if (p.fixed_dt == true) {
		loop_time_adjust = dt;
	}
	while (keepRunning(time_end - loop_time_adjust, strain_end - loop_time_adjust)) {
		retrim_ext_flow = false; // used in ext_flow simulation
		if (!brownian && !p.fixed_dt) { // adaptative time-step for non-Brownian cases
			adaptTimeStep(time_end, strain_end);
		}
		if (ext_flow) {
			if (clk.cumulated_strain+shear_rate*dt > strain_retrim) {
				cerr << "strain_retrim = " << strain_retrim << endl;
				dt_bak = dt;
				dt = (strain_retrim-clk.cumulated_strain)/shear_rate;
				retrim_ext_flow = true;
				strain_end = strain_retrim;
				calc_stress = true;
			}
		}
		(this->*timeEvolutionDt)(calc_stress, time_end, strain_end); // no stress computation except at low Peclet
		avg_dt += dt;
		avg_dt_nb++;
		if (dt_bak != -1){
			dt = dt_bak;
		}
		if (sig_caught == SIGINT) { // return to Simulation immediatly for checkpointing
			return;
		}
	};
	if (avg_dt_nb > 0) {
		avg_dt /= avg_dt_nb;
	} else {
		avg_dt = dt;
	}
	if (events.empty() && retrim_ext_flow == false) {
		calc_stress = true;
		(this->*timeEvolutionDt)(calc_stress, time_end, strain_end); // last time step, compute the stress
	}
	if (p.auto_determine_knkt
		&& clk.cumulated_strain > p.start_adjust) {
		adjustContactModelParameters();
	}
}

void System::createNewInteraction(int i, int j, double scaled_interaction_range)
{
	// new interaction
	Interaction inter(this, i, j, scaled_interaction_range);
	interaction.push_back(inter); // could emplace_back if Interaction gets a move ctor
	// tell i and j their new partner
	interaction_partners[i].push_back(j);
	interaction_partners[j].push_back(i);
}

bool System::hasNeighbor(int i, int j)
{
	for (int k : interaction_partners[i]) {
		if (j == k) {
			return true;
		}
	}
	return false;
}

void System::removeNeighbors(int i, int j)
{
	// return interaction_partners[i].find(j) != interaction_partners[i].end();
	vector<int> &neighi = interaction_partners[i];
	vector<int> &neighj = interaction_partners[j];
	int l = neighi.back();
	if (l != j) {
		for (unsigned int k=0; k<neighi.size(); k++) {
			if (neighi[k] == j) {
				neighi[k] = l;
				break;
			}
		}
	}
	neighi.pop_back();
	l = neighj.back();
	if (l != i) {
		for (unsigned int k=0; k<neighj.size(); k++) {
			if (neighj[k] == i) {
				neighj[k] = l;
				break;
			}
		}
	}
	neighj.pop_back();
}

void System::checkNewInteraction()
{
	/**
	 \brief Checks if there are new pairs of interacting particles. If so, creates and sets up the corresponding Interaction objects.

	 To be called after particle moved.
	 */
	vec3d pos_diff;
	double sq_dist;
	if (!ext_flow) {
		for (int i=0; i<np-1; i++) {
			for (auto j : boxset.neighborhood(i)) {
				if (j > i) {
					if (!hasNeighbor(i, j)) {
						pos_diff = position[j]-position[i];
						periodizeDiff(pos_diff);
						sq_dist = pos_diff.sq_norm();
						double scaled_interaction_range = (this->*calcInteractionRange)(i, j);
						double sq_dist_lim = scaled_interaction_range*scaled_interaction_range;
						if (sq_dist < sq_dist_lim) {
							createNewInteraction(i, j, scaled_interaction_range);
						}
					}
				}
			}
		}
	} else {
		vec3d pd_shift;
		for (int i=0; i<np-1; i++) {
			for (auto j : boxset.neighborhood(i)) {
				if (j > i) {
					if (!hasNeighbor(i, j)) {
						pos_diff = position[j]-position[i];
						periodizeDiffExtFlow(pos_diff, pd_shift, i, j);
						sq_dist = pos_diff.sq_norm();
						double scaled_interaction_range = (this->*calcInteractionRange)(i, j);
						double sq_dist_lim = scaled_interaction_range*scaled_interaction_range;
						if (sq_dist < sq_dist_lim) {
							createNewInteraction(i, j, scaled_interaction_range);
						}
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
	for (unsigned int k=0; k<interaction.size(); k++) {
		bool deactivated = false;
		interaction[k].updateState(deactivated);
		if (friction && interaction[k].contact.is_active()) {
			double sq_sliding_velocity = interaction[k].contact.relative_surface_velocity_sqnorm;
			if (sq_sliding_velocity > sq_max_sliding_velocity) {
				sq_max_sliding_velocity = sq_sliding_velocity;
			}
		}
		if (deactivated) {
			auto ij = interaction[k].get_par_num();
			removeNeighbors(ij.first, ij.second);
			interaction[k] = interaction[interaction.size()-1];
			interaction.pop_back();
		}
	}
	max_sliding_velocity = sqrt(sq_max_sliding_velocity);
}

void System::declareResistance(int p0, int p1)
{
	if (p1 < np_mobile) { // i and j mobile
		nb_blocks_mm[p0]++;
	} else if (p0 >= np_mobile) { // i and j fixed
		nb_blocks_ff[p0-np_mobile]++;
	} else {
		nb_blocks_mf[p0]++;
	}
	pairwise_resistance_changed = true;
}

void System::eraseResistance(int p0, int p1)
{
	if (p1 < np_mobile) { // i and j mobile
		nb_blocks_mm[p0]--;
	} else if (p0 >= np_mobile) { // i and j fixed
		nb_blocks_ff[p0-np_mobile]--;
	} else {
		nb_blocks_mf[p0]--;
	}
	pairwise_resistance_changed = true;
}

void System::buildResistanceMatrix()
{
	/**
	 \brief Builds the resistance matrix

	 The built matrix \f$R_{\mathrm{FU}}\f$ (in Bossis and Brady \cite
	 brady_stokesian_1988 notations) contains pairwise resistances,
	 coming from lubrication or contact dashpots.
	 */
	unsigned int size_mm = 0;
	for (auto bnb: nb_blocks_mm) {
		size_mm += bnb;
	}
	unsigned int size_mf = 0;
	for (auto bnb: nb_blocks_mf) {
		size_mf += bnb;
	}
	unsigned int size_ff = 0;
	for (auto bnb: nb_blocks_ff) {
		size_ff += bnb;
	}

	// create a new resistance matrix in stokes_solver
	/* [note]
	 * The resistance matrix is reset with resistance_matrix_dblock,
	 * which is calculated at the beginning.
	 */
	stokes_solver.resetResistanceMatrix(size_mm, size_mf, size_ff,
										resistance_matrix_dblock, pairwise_resistance_changed);
	pairwise_resistance_changed = false;
	for (int i=0; i<np-1; i ++) {
		stokes_solver.startNewColumn();
		for (auto& inter : interaction_list[i]) {
			int j = inter->partner(i);
			if (j > i) {
				if (inter->hasPairwiseResistance()) { // Range of interaction can be larger than range of lubrication
					stokes_solver.addResistanceBlocks(i, j,
													  inter->RFU_DBlocks(),
													  inter->RFU_ODBlock());
				}
			}
		}
	}
	stokes_solver.matrixFillingDone();
}

void System::computeForcesOnWallParticles()
{
	/**
		\brief This method computes the force (and torque, for now, but it might be dropped)
		on the fixed particles.

		It is designed with simple shear with walls under stress controlled conditions in mind,
		so it decomposes the force in a rate-proportional part and a rate-independent part.

		*/
	throw runtime_error(" Control stress with walls disabled for now .\n");
	if (!zero_shear) {
		throw runtime_error(" Stress-control with walls requires zero_shear==true .\n");
	}
	vector<vec3d> force (p.np_fixed);
	vector<vec3d> torque (p.np_fixed);

	// Compute the part of the velocity of mobile particles
	// that is not coming from the wall velocities
	vector<vec3d> na_velocity_mobile (np_mobile);
	vector<vec3d> na_ang_velocity_mobile (np_mobile);
	const auto &vel_contact = na_velo_components["contact"].vel;
	const auto &ang_vel_contact = na_velo_components["contact"].ang_vel;

	for (int i=0; i<np_mobile; i++) {
		na_velocity_mobile[i] = vel_contact[i];
		na_ang_velocity_mobile[i] = ang_vel_contact[i];
	}
	if (repulsiveforce) {
		const auto &vel_repulsive = na_velo_components["repulsion"].vel;
		const auto &ang_vel_repulsive = na_velo_components["repulsion"].ang_vel;
		for (int i=0; i<np_mobile; i++) {
			na_velocity_mobile[i] += vel_repulsive[i];
			na_ang_velocity_mobile[i] += ang_vel_repulsive[i];
		}
	}
	if (delayed_adhesion) {
		const auto &vel_adhesion = na_velo_components["delayed_adhesion"].vel;
		const auto &ang_vel_adhesion = na_velo_components["delayed_adhesion"].ang_vel;
		for (int i=0; i<np_mobile; i++) {
			na_velocity_mobile[i] += vel_adhesion[i];
			na_ang_velocity_mobile[i] += ang_vel_adhesion[i];
		}
	}
	// from this, we can compute the hydro force on the wall that does *not* depend on the wall velocity
	stokes_solver.multiply_by_RFU_fm(na_velocity_mobile, na_ang_velocity_mobile, force, torque);

	// Now we sum up this hydro part with the other non-rate dependent forces (contact, etc)
	const auto &contact_force = force_components["contact"].force;
	const auto &contact_torque = force_components["contact"].torque;
	const auto &repulsive_force = force_components["repulsion"].force;
	const auto &adhesion_force = force_components["delayed_adhesion"].force;
	throw std::runtime_error("computeForcesOnWallParticles : ongoing work.... sorry");
	for (int i=0; i<p.np_fixed; i++) {
		non_rate_proportional_wall_force[i] = -force[i];
		non_rate_proportional_wall_torque[i] = -torque[i];
		non_rate_proportional_wall_force[i] += contact_force[i+np_mobile];
		non_rate_proportional_wall_torque[i] += contact_torque[i+np_mobile];
		if (repulsiveforce) {
			non_rate_proportional_wall_force[i] += repulsive_force[i+np_mobile];
		}
	}

	// Now the part proportional to the wall speed

	// From the mobile particles
	const auto &vel_hydro_from_fixed = na_velo_components["from_fixed"].vel;
	const auto &ang_vel_hydro_from_fixed = na_velo_components["from_fixed"].ang_vel;
	for (int i=0; i<np_mobile; i++) {
		na_velocity_mobile[i] = vel_hydro_from_fixed[i];
		na_ang_velocity_mobile[i] = ang_vel_hydro_from_fixed[i];
	}
	stokes_solver.multiply_by_RFU_fm(na_velocity_mobile, na_ang_velocity_mobile, force, torque);
	for (int i=0; i<p.np_fixed; i++) {
		rate_proportional_wall_force[i] = -force[i];
		rate_proportional_wall_torque[i] = -torque[i];
	}

	// From the fixed particles themselves. This should be zero if these particles form a wall
	// (i.e. they move with zero relative velocity) and if the Stokes drag is zero (which is controlled by sd_coeff)
	// As we do not want to make too many assumptions here (especially regarding the Stokes drag)
	// we compute it. [Probably a p.no_stokes_drag should be introduced at some point.]
	vector<vec3d> na_velocity_fixed (p.np_fixed);
	vector<vec3d> na_ang_velocity_fixed (p.np_fixed);
	for (int i=0; i<p.np_fixed; i++) {
		na_velocity_fixed[i] = na_velocity[i+np_mobile];
		na_ang_velocity_fixed[i] = na_ang_velocity[i+np_mobile];
	}
	stokes_solver.multiply_by_RFU_ff(na_velocity_fixed, na_ang_velocity_fixed, force, torque);
	for (int i=0; i<p.np_fixed; i++) {
		rate_proportional_wall_force[i] -= force[i];
		rate_proportional_wall_torque[i] -= torque[i];
	}
}

void System::forceResultantLubricationForce()
{
	/* Only F = R_FU U is calculated, but R_FE E is not implemented yet.
	 * So, we cannot check the force balance with E^{inf} yet.
	 */
	/*
	 *  F^{M} = R_FU^{MM} U^{M}
	 */
	vector<double> force_m_to_m (6*np_mobile);
	vector<double> minus_mobile_velocities (6*np_mobile);
	for (int i=0; i<np_mobile; i++) {
		int i6 = 6*i;
		minus_mobile_velocities[i6  ] = -na_velocity[i].x;
		minus_mobile_velocities[i6+1] = -na_velocity[i].y;
		minus_mobile_velocities[i6+2] = -na_velocity[i].z;
		minus_mobile_velocities[i6+3] = -na_ang_velocity[i].x;
		minus_mobile_velocities[i6+4] = -na_ang_velocity[i].y;
		minus_mobile_velocities[i6+5] = -na_ang_velocity[i].z;
	}
	stokes_solver.multiply_by_RFU_mm(minus_mobile_velocities, force_m_to_m);
	for (int i=0; i<np_mobile; i++) {
		int i6 = 6*i;
		forceResultant[i].x += force_m_to_m[i6];
		forceResultant[i].y += force_m_to_m[i6+1];
		forceResultant[i].z += force_m_to_m[i6+2];
		torqueResultant[i].x += force_m_to_m[i6+3];
		torqueResultant[i].y += force_m_to_m[i6+4];
		torqueResultant[i].z += force_m_to_m[i6+5];
	}
	if (mobile_fixed) {
		/*
		 *  F^{M} += R_FU^{MF} U^{F}
		 */
		vector<double> force_f_to_m (6*np_mobile);
		vector<double> minus_fixed_velocities (6*p.np_fixed);
		for (int i=0; i<p.np_fixed; i++) {
			int i6 = 6*i;
			int i_fixed = i+np_mobile;
			minus_fixed_velocities[i6  ] = -na_velocity[i_fixed].x;
			minus_fixed_velocities[i6+1] = -na_velocity[i_fixed].y;
			minus_fixed_velocities[i6+2] = -na_velocity[i_fixed].z;
			minus_fixed_velocities[i6+3] = -na_ang_velocity[i_fixed].x;
			minus_fixed_velocities[i6+4] = -na_ang_velocity[i_fixed].y;
			minus_fixed_velocities[i6+5] = -na_ang_velocity[i_fixed].z;
		}
		stokes_solver.multiply_by_RFU_mf(minus_fixed_velocities, force_f_to_m);
		for (int i=0; i<np_mobile; i++) {
			int i6 = 6*i;
			forceResultant[i].x += force_f_to_m[i6];
			forceResultant[i].y += force_f_to_m[i6+1];
			forceResultant[i].z += force_f_to_m[i6+2];
			torqueResultant[i].x += force_f_to_m[i6+3];
			torqueResultant[i].y += force_f_to_m[i6+4];
			torqueResultant[i].z += force_f_to_m[i6+5];
		}
		/*
		 *  F^{F} += R_FU^{FM} U^{M}
		 */
		vector<double> force_m_to_f (6*p.np_fixed);
		for (int i=0; i<np_mobile; i++) {
			int i6 = 6*i;
			minus_mobile_velocities[i6  ] = -na_velocity[i].x;
			minus_mobile_velocities[i6+1] = -na_velocity[i].y;
			minus_mobile_velocities[i6+2] = -na_velocity[i].z;
			minus_mobile_velocities[i6+3] = -na_ang_velocity[i].x;
			minus_mobile_velocities[i6+4] = -na_ang_velocity[i].y;
			minus_mobile_velocities[i6+5] = -na_ang_velocity[i].z;
		}
		stokes_solver.multiply_by_RFU_fm(minus_mobile_velocities, force_m_to_f);
		for (int i=np_mobile; i<np; i++) {
			int i6 = 6*(i-np_mobile);
			forceResultant[i].x += force_m_to_f[i6];
			forceResultant[i].y += force_m_to_f[i6+1];
			forceResultant[i].z += force_m_to_f[i6+2];
			torqueResultant[i].x += force_m_to_f[i6+3];
			torqueResultant[i].y += force_m_to_f[i6+4];
			torqueResultant[i].z += force_m_to_f[i6+5];
		}
		/*
		 *  F^{F} += R_FU^{FF} U^{F}
		 */
		vector<double> force_f_to_f (6*p.np_fixed);
		for (int i=0; i<p.np_fixed; i++) {
			int i6 = 6*i;
			int i_fixed = i+np_mobile;
			minus_fixed_velocities[i6  ] = -na_velocity[i_fixed].x;
			minus_fixed_velocities[i6+1] = -na_velocity[i_fixed].y;
			minus_fixed_velocities[i6+2] = -na_velocity[i_fixed].z;
			minus_fixed_velocities[i6+3] = -na_ang_velocity[i_fixed].x;
			minus_fixed_velocities[i6+4] = -na_ang_velocity[i_fixed].y;
			minus_fixed_velocities[i6+5] = -na_ang_velocity[i_fixed].z;
		}
		stokes_solver.multiply_by_RFU_ff(minus_fixed_velocities, force_f_to_f);
		for (int i=np_mobile; i<np; i++) {
			int i6 = 6*(i-np_mobile);
			forceResultant[i].x += force_f_to_f[i6];
			forceResultant[i].y += force_f_to_f[i6+1];
			forceResultant[i].z += force_f_to_f[i6+2];
			torqueResultant[i].x += force_f_to_f[i6+3];
			torqueResultant[i].y += force_f_to_f[i6+4];
			torqueResultant[i].z += force_f_to_f[i6+5];
		}
	}
}

void System::setBrownianForceToParticle(vector<vec3d> &force,
										vector<vec3d> &torque)
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

	 */
	if(!in_predictor) { // The Brownian force must be the same in the predictor and the corrector
		return;
	}
	if (mobile_fixed) {
		throw runtime_error("Brownian algorithm with fixed particles not implemented yet.\n");
	}
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	double sqrt_2_dt_amp = sqrt(2*p.brownian/dt);
	for (unsigned int i=0; i<force.size(); i++) {
		force[i].x = sqrt_2_dt_amp*GRANDOM; // \sqrt(2kT/dt) * random vector A (force and torque)
		force[i].y = sqrt_2_dt_amp*GRANDOM;
		force[i].z = sqrt_2_dt_amp*GRANDOM;
		torque[i].x = sqrt_2_dt_amp*GRANDOM; // \sqrt(2kT/dt) * random vector A (force and torque)
		torque[i].y = sqrt_2_dt_amp*GRANDOM;
		torque[i].z = sqrt_2_dt_amp*GRANDOM;
	}

	if (pairwise_resistance) {
		/* L*L^T = RFU
		 */
		stokes_solver.setRHS(force, torque);
		stokes_solver.compute_LTRHS(force, torque); // F_B = \sqrt(2kT/dt) * L^T * A
		stokes_solver.resetRHS();
	} else {
		/*
		 *  F_B = \sqrt(2kT/dt) * L^T * A
		 *  U_B = (RFU)^{-1} F_B
		 *  In Stokes drag simulation
		 *  F_B = \sqrt(2kT/dt) * sqrt(RFU) * A
		 *  U_B = (RFU)^{-1} F_B
		 *	    = (RFU)^{-1} \sqrt(2kT/dt) * sqrt(RFU) * A
		 *      = \sqrt(2kT/dt) * A / sqrt(RFU)
		 *  In order to reduce trivial calculations,
		 *  Here, sqrt(RFU) is not included in F_B
		 *  F_B = \sqrt(2kT/dt) * A
		 *  In the function computeVelocitiesStokesDrag(),
		 *  U_B = F_B / sqrt(RFU)
		 */
	}
}

void System::setDashpotForceToParticle(vector<vec3d> &force,
									   vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	if (!zero_shear) {
		vec3d GEi, GEj, HEi, HEj;
		unsigned int i, j;
		for (const auto &inter: interaction) {
			if (inter.contact.is_active() && inter.contact.dashpot.is_active()) {
				std::tie(i, j) = inter.get_par_num();
				std::tie(GEi, GEj, HEi, HEj) = inter.contact.dashpot.getRFU_Uinf(u_inf[i], u_inf[j], omega_inf);
				force[i] += GEi;
				force[j] += GEj;
				torque[i] += HEi;
				torque[j] += HEj;
			}
		}
	}
}

void System::setHydroForceToParticle_squeeze(vector<vec3d> &force,
											 vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	if (!zero_shear) {
		vec3d GEi, GEj;
		unsigned int i, j;
		for (const auto &inter: interaction) {
			if (inter.lubrication.is_active()) {
				std::tie(i, j) = inter.get_par_num();
				std::tie(GEi, GEj) = inter.lubrication.calcGE_squeeze(E_infinity); // G*E_\infty term
				force[i] += GEi;
				force[j] += GEj;
			}
		}
	}
}

void System::setHydroForceToParticle_squeeze_tangential(vector<vec3d> &force,
														vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	if (!zero_shear) {
		vec3d GEi, GEj, HEi, HEj;
		unsigned int i, j;
		for (const auto &inter: interaction) {
			if (inter.lubrication.is_active()) {
				std::tie(i, j) = inter.get_par_num();
				std::tie(GEi, GEj, HEi, HEj) = inter.lubrication.calcGEHE_squeeze_tangential(E_infinity); // G*E_\infty term, no gamma dot
				force[i] += GEi;
				force[j] += GEj;
				torque[i] += HEi;
				torque[j] += HEj;
			}
		}
	}
}

void System::setContactForceToParticle(vector<vec3d> &force,
									   vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	for (const auto &inter: interaction) {
		if (inter.contact.is_active()) {
			inter.contact.addUpForceTorque(force, torque);
		}
	}
}

void System::setRepulsiveForceToParticle(vector<vec3d> &force,
										 vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	for (const auto &inter: interaction) {
		inter.repulsion.addUpForce(force);
	}
}

void System::setTActAdhesionForceToParticle(vector<vec3d> &force,
								            vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	unsigned int i, j;
	for (const auto &inter: interaction) {
		std::tie(i, j) = inter.get_par_num();
		inter.delayed_adhesion->addUpForce(force[i], force[j]);
	}
}

void System::setFixedParticleForceToParticle(vector<vec3d> &force,
											 vector<vec3d> &torque)
{
	vector<double> force_torque_from_fixed (6*np_mobile);
	// @@ TODO: avoid copy of the velocities and forces
	vector<double> minus_fixed_velocities (6*p.np_fixed);
	for (int i=0; i<p.np_fixed; i++) {
		int i6 = 6*i;
		int i_fixed = i+np_mobile;
		minus_fixed_velocities[i6  ] = -na_velocity[i_fixed].x;
		minus_fixed_velocities[i6+1] = -na_velocity[i_fixed].y;
		minus_fixed_velocities[i6+2] = -na_velocity[i_fixed].z;
		minus_fixed_velocities[i6+3] = -na_ang_velocity[i_fixed].x;
		minus_fixed_velocities[i6+4] = -na_ang_velocity[i_fixed].y;
		minus_fixed_velocities[i6+5] = -na_ang_velocity[i_fixed].z;
	}
	stokes_solver.multiply_by_RFU_mf(minus_fixed_velocities, force_torque_from_fixed); // -R_FU^mf*fixed_velocities
	for (unsigned int i=0; i<force.size(); i++) {
		auto i6 = 6*i;
		force[i].x = force_torque_from_fixed[i6  ];
		force[i].y = force_torque_from_fixed[i6+1];
		force[i].z = force_torque_from_fixed[i6+2];
		torque[i].x = force_torque_from_fixed[i6+3];
		torque[i].y = force_torque_from_fixed[i6+4];
		torque[i].z = force_torque_from_fixed[i6+5];
	}
}

void System::setSolverRHS(const ForceComponent &fc)
{
	if (fc.has_torque) {
		for (int i=0; i<np; i++) {
			stokes_solver.setRHSForce(i, fc.force[i]);
			stokes_solver.setRHSTorque(i, fc.torque[i]);
		}
	} else {
		for (int i=0; i<np; i++) {
			stokes_solver.setRHSForce(i, fc.force[i]);
		}
		stokes_solver.resetRHStorque();
	}
}

void System::addToSolverRHS(const ForceComponent &fc)
{
	if (fc.has_torque) {
		for (int i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, fc.force[i]);
			stokes_solver.addToRHSTorque(i, fc.torque[i]);
		}
	} else {
		for (int i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, fc.force[i]);
		}
	}
}

void System::computeMaxNAVelocity()
{
	/**
	 \brief Compute the maximum non-affine velocity

	 Note: it does \b not compute the velocities, just takes the maximum.
	 */
	double sq_max_na_velocity = 0;
	double sq_na_velocity;
	for (int i=0; i<np; i++) {
		sq_na_velocity = na_velocity[i].sq_norm();
		if (sq_na_velocity > sq_max_na_velocity) {
			sq_max_na_velocity = sq_na_velocity;
		}
	}
	max_velocity = sqrt(sq_max_na_velocity);
	if (brownian) {
		double sq_max_na_velocity_br = 0;
		double sq_na_velocity_br;
		for (int i=0; i<np; i++) {
			sq_na_velocity_br = na_velo_components["brownian"].vel[i].sq_norm();
			if (sq_na_velocity_br > sq_max_na_velocity_br) {
				sq_max_na_velocity_br = sq_na_velocity_br;
			}
		}
		max_velocity_brownian = sqrt(sq_max_na_velocity_br);
		double sq_max_na_velocity_con = 0;
		double sq_na_velocity_con;
		for (int i=0; i<np; i++) {
			sq_na_velocity_con = na_velo_components["contact"].vel[i].sq_norm();
			if (sq_na_velocity_con > sq_max_na_velocity_con) {
				sq_max_na_velocity_con = sq_na_velocity_con;
			}
		}
		max_velocity_contact = sqrt(sq_max_na_velocity_con);
	}
	//	double sq_na_ang_velocity;
	//	for (int i=0; i<np; i++) {
	//		sq_na_ang_velocity = na_ang_velocity[i].sq_norm()*radius_squared[i];
	//		if (sq_max_na_velocity < sq_na_ang_velocity) {
	//			sq_max_na_velocity = sq_na_ang_velocity;
	//		}
	//	}
}

void System::computeVelocityWithoutComponents()
{
	buildResistanceMatrix();
	for (auto &fc: force_components) {
		if (fc.first != "brownian" || in_predictor) {
			CALL_MEMBER_FN(*this, fc.second.getForceTorque)(fc.second.force, fc.second.torque);
		}
		addToSolverRHS(fc.second);
	}
	stokes_solver.solve(na_velocity, na_ang_velocity); // get V
	if (brownian && twodimension) {
		rushWorkFor2DBrownian(na_velocity, na_ang_velocity);
	}
}

void System::computeVelocityByComponents()
{
	/**
	 \brief Compute velocities component by component.
	 */
	buildResistanceMatrix();
	for (auto &fc: force_components) {
		if (fc.first != "brownian" || in_predictor) {
			CALL_MEMBER_FN(*this, fc.second.getForceTorque)(fc.second.force, fc.second.torque);
		}
		setSolverRHS(fc.second);
		stokes_solver.solve(na_velo_components[fc.first].vel,
							na_velo_components[fc.first].ang_vel);
	}
	if (brownian && twodimension) {
		rushWorkFor2DBrownian(na_velo_components["brownian"].vel,
							  na_velo_components["brownian"].ang_vel);
	}
}

void System::setVelocityDifference()
{
	vel_difference = 2*dot(E_infinity, {0, 0, lz});
}

void System::set_shear_rate(double shear_rate_)
{
	shear_rate = shear_rate_;
	omega_inf = omegahat_inf*shear_rate;
	E_infinity = Ehat_infinity*shear_rate;
	if (!ext_flow) {
		setVelocityDifference();
	} else {
		grad_u = shear_rate*grad_u_hat;
	}
}

void System::setImposedFlow(Sym2Tensor EhatInfty, vec3d OhatInfty)
{
	Ehat_infinity = EhatInfty;
	omegahat_inf = OhatInfty;
	if (twodimension) {
		if (fabs(Ehat_infinity.elm[3]) > 1e-15 || fabs(Ehat_infinity.elm[4]) > 1e-15) {
			throw runtime_error(" System:: Error: 2d simulation with Einf_{y?} != 0");
		} else {
			Ehat_infinity.elm[3] = 0;
			Ehat_infinity.elm[4] = 0;
		}
		if (fabs(omegahat_inf.x) > 1e-15 || fabs(omegahat_inf.z) > 1e-15) {
			throw runtime_error(" System:: Error: 2d simulation with Omega_inf not along y");
		} else {
			omegahat_inf.x = 0;
			omegahat_inf.z = 0;
		}
	}
	//	omega_inf = omegahat_inf*shear_rate;
	//E_infinity = Ehat_infinity*shear_rate;
}

void System::setShearDirection(double theta_shear) // will probably be deprecated soon
{
	if (!ext_flow) {
		p.theta_shear = theta_shear;
		double costheta_shear = cos(theta_shear);
		double sintheta_shear = sin(theta_shear);
		setImposedFlow({0, 0, costheta_shear/2, sintheta_shear/2, 0, 0},
					   {-0.5*sintheta_shear, 0.5*costheta_shear, 0});
	}
}

void System::computeShearRate()
{
	/**
	 \brief Compute the shear rate under stress control conditions.
	 */
	assert(abs(shear_rate-1) < 1e-15);
	calcStressPerParticle();
	Sym2Tensor rate_prop_stress;
	Sym2Tensor rate_indep_stress;
	gatherStressesByRateDependencies(rate_prop_stress, rate_indep_stress);
	/*
	 *  target_stress = rate_indep_stress_ + rate * rate_prop_stress_at_1
	 *  rate = (target_stress - rate_indep_stress)/rate_prop_stress_at_1
	 *
	 *  kappa = (Sigma.E)/(E.E)
	 *  sigama = -p*I + 2 eta*E + ... = -p*I + kappa*E + ...
	 *  E.E = 1/2 (if shear rate = 1)
	 *  eta = kappa / 2 = (Sigma.E)/(2*E.E) = Sigma.E
	 *  sigma_target = rate_prop_shearstress(rate=1)*(rate/1) + rate_indep_shearstress
	 *  rate = (sigma_target - rate_indep_shearstress) / rate_prop_shearstress(rate=1)
	 */
	double rate_prop_shearstress_rate1 = doubledot(rate_prop_stress, getEinfty()); // computed with rate=1, o here it is also the viscosity.
	double rate_indep_shearstress = doubledot(rate_indep_stress, getEinfty());
	if (brownian) {
		rate_prop_shearstress_rate1_ave.update(rate_prop_shearstress_rate1, get_time());
		rate_indep_shearstress_ave.update(rate_indep_shearstress, get_time());
		rate_prop_shearstress_rate1 = rate_prop_shearstress_rate1_ave.get();
		rate_indep_shearstress = rate_indep_shearstress_ave.get();
	}
	if (rate_prop_shearstress_rate1 != 0) {
		double rate = (target_stress-rate_indep_shearstress)/rate_prop_shearstress_rate1;
		set_shear_rate(rate);
	}
	if (clk.cumulated_strain < init_strain_shear_rate_limit) {
		if (shear_rate > init_shear_rate_limit) {
			set_shear_rate(init_shear_rate_limit);
		}
	}
}

void System::computeShearRateWalls()
{
	/**
	 \brief Compute the coefficient to give to the velocity of the fixed particles under stress control conditions.
	 */

	computeForcesOnWallParticles();

	double total_rate_dep_wall_shear_stress = 0;
	double total_rate_indep_wall_shear_stress = 0;

	for (int i=0; i<p.np_fixed; i++) {
		total_rate_dep_wall_shear_stress += dot(fixed_velocities[i], rate_proportional_wall_force[i]);
		total_rate_indep_wall_shear_stress += dot(fixed_velocities[i], non_rate_proportional_wall_force[i]);
	}
	double wall_surface;
	if (twodimension) {
		wall_surface = lx;
	} else {
		wall_surface = lx*ly;
	}

	total_rate_dep_wall_shear_stress /= wall_surface;
	total_rate_indep_wall_shear_stress /= wall_surface;

	// // the total_rate_dep_wall_shear_stress is computed above with shear_rate=1, so here it is also a viscosity.
	set_shear_rate((-target_stress-total_rate_indep_wall_shear_stress)/total_rate_dep_wall_shear_stress);

	if (clk.cumulated_strain < init_strain_shear_rate_limit) {
		if (shear_rate > init_shear_rate_limit) {
			set_shear_rate(init_shear_rate_limit);
		}
	}
	if (p.simulation_mode == 31) {
		force_upwall.reset();
		force_downwall.reset();
		for (int i=0; i<p.np_fixed; i++) {
			if (fixed_velocities[i].x > 0) {
				force_upwall += shear_rate*rate_proportional_wall_force[i]+non_rate_proportional_wall_force[i];
			}
			if (fixed_velocities[i].x < 0) {
				force_downwall += shear_rate*rate_proportional_wall_force[i]+non_rate_proportional_wall_force[i];
			}
		}
	}
}

void System::tmpMixedProblemSetVelocities()
{
	if (p.simulation_mode == 1) {
		/* Shear reversal simulation
		 */
		static double time_next = 16;
		static double direction = 1;
		if (get_time() > time_next) {
			direction *= -1;
			time_next += 16;
			cerr << direction << endl;
		}
		for (int i=np_mobile; i<np; i++) { // temporary: particles perfectly advected
			na_velocity[i].reset();
			na_ang_velocity[i].reset();
		}
		na_velocity[np_mobile].x = direction;
	} else if (p.simulation_mode == 2) {
		/* ????, yeah what happened here and below?? Haha. I guess whoever created simulation_mode 2/3 does not use it often :)
		 */
		for (int i=np_mobile; i<np; i++) { // temporary: particles perfectly advected
			na_velocity[i].reset();
			na_ang_velocity[i].reset();
		}
		na_ang_velocity[np_mobile].y = 2*shear_rate;
	} else if (p.simulation_mode == 3) {
		/* ????
		 */

		for (int i=np_mobile; i<np; i++) { // temporary: particles perfectly advected
			na_velocity[i].reset();
			na_ang_velocity[i].reset();
		}
		na_velocity[np_mobile].x = 1;
	} else if (p.simulation_mode == 4) {
		static double time_next = 10;
		if (get_time() > time_next) {
			if (zero_shear == true) {
				zero_shear = false;
			} else {
				zero_shear = true;
			}
			time_next += 10;
		}
		if (zero_shear) {
			for (int i=np_mobile; i<np; i++) { // temporary: particles perfectly advected
				na_velocity[i].reset();
				na_ang_velocity[i].reset();
			}
		}
	} else if (p.simulation_mode >= 10 && p.simulation_mode < 20) {
		int i_np_in = np_mobile+np_wall1;
		// inner wheel
		for (int i=np_mobile; i<i_np_in; i++) { // temporary: particles perfectly advected
			na_velocity[i] = {-omega_wheel_in*(position[i].z-origin_of_rotation.z),
				0,
				omega_wheel_in*(position[i].x-origin_of_rotation.x)};
			na_ang_velocity[i] = {0, -omega_wheel_in, 0};
		}
		// outer wheel
		for (int i=i_np_in; i<np; i++) { // temporary: particles perfectly advected
			na_velocity[i] = {-omega_wheel_out*(position[i].z-origin_of_rotation.z),
				0,
				omega_wheel_out*(position[i].x-origin_of_rotation.x)};
			na_ang_velocity[i] = {0, -omega_wheel_out, 0};
		}
	} else if (p.simulation_mode == 21) {
		static double time_next = p.strain_reversal;
		if (get_time() > time_next) {
			p.theta_shear += M_PI;
			setShearDirection(p.theta_shear);
			time_next += p.strain_reversal;
		}
	} else if (p.simulation_mode == 31) {
		auto &vel_from_fixed = na_velo_components["from_fixed"];
		for (int i=np_mobile; i<np; i++) {
			vel_from_fixed.vel[i] = shear_rate*fixed_velocities[i-np_mobile];
			vel_from_fixed.ang_vel[i].reset();
			na_velocity[i] = vel_from_fixed.vel[i];
			na_ang_velocity[i] = vel_from_fixed.ang_vel[i];
		}
	} else if (p.simulation_mode == 41) {
		int i_np_wall1 = np_mobile+np_wall1;
		double wall_velocity = shear_rate*system_height;
		for (int i=np_mobile; i<i_np_wall1; i++) {
			na_velocity[i] = {-wall_velocity/2, 0, 0};
			na_ang_velocity[i].reset();
		}
		for (int i=i_np_wall1; i<np; i++) {
			na_velocity[i] = {wall_velocity/2, 0, 0};
			na_ang_velocity[i].reset();
		}
	} else if (p.simulation_mode == 42) {
		int i_np_wall1 = np_mobile+np_wall1;
		double wall_velocity = shear_rate*system_height;
		for (int i=np_mobile; i<i_np_wall1; i++) {
			na_velocity[i].reset();
			na_ang_velocity[i].reset();
		}
		for (int i=i_np_wall1; i<np; i++) {
			na_velocity[i] = {wall_velocity, 0, 0};
			na_ang_velocity[i].reset();
		}
	} else if (p.simulation_mode == 51) {
		int i_np_in = np_mobile+np_wall1;
		// inner wheel
		double l = lx/2;
		vec3d origin_of_rotation2(lx/2, 0, l);
		vec3d origin_of_rotation3(  lx, 0, 0);
		double x1 = l/sqrt(2);
		double x2 = x1+radius_in*sqrt(2);
		for (int i=i_np_in; i<np; i++) {
			if (position[i].x < x1) {
				na_velocity[i] = {-omega_wheel_out*(position[i].z),
					0,
					omega_wheel_out*(position[i].x)};
				na_ang_velocity[i] = {0, -omega_wheel_out, 0};
			} else if (position[i].x < x2) {
				na_velocity[i] = {-omega_wheel_in*(position[i].z-origin_of_rotation2.z),
					0,
					omega_wheel_in*(position[i].x-origin_of_rotation2.x)};
				na_ang_velocity[i] = {0, -omega_wheel_in, 0};
			} else {
				na_velocity[i] = {-omega_wheel_out*(position[i].z-origin_of_rotation3.z),
					0,
					omega_wheel_out*(position[i].x-origin_of_rotation3.x)};
				na_ang_velocity[i] = {0, -omega_wheel_out, 0};
			}
		}
	}
}

void System::sumUpVelocityComponents()
{
	for (int i=0; i<np_mobile; i++) {
		na_velocity[i].reset();
		na_ang_velocity[i].reset();
	}
	for (const auto &vc: na_velo_components) {
		const auto &vel = vc.second.vel;
		const auto &ang_vel = vc.second.ang_vel;
		for (int i=0; i<np_mobile; i++) {
			na_velocity[i] += vel[i];
			na_ang_velocity[i] += ang_vel[i];
		}
	}
}

void System::setFixedParticleVelocities()
{
	if (p.simulation_mode == 0) {
		for (int i=np_mobile; i<np; i++) { // temporary: particles perfectly advected
			na_velocity[i].reset();
			na_ang_velocity[i].reset();
		}
	} else if (p.simulation_mode > 0) {
		tmpMixedProblemSetVelocities();
	}
}

void System::rescaleVelHydroStressControlled()
{
	for (auto &vc: na_velo_components) {
		if (vc.second.rate_dependence == RATE_PROPORTIONAL) {
			vc.second *= shear_rate;
		}
	}
}

void System::computeVelocities(bool divided_velocities)
{
	/**
	 \brief Compute velocities in the current configuration.

	 \param divided_velocities Divide the velocities in components
	 (hydro, contacts, Brownian, ...). (Note that in Brownian
	 simulations the Brownian component is always computed explicitely, independently of the values of divided_velocities.)
	 */
	stokes_solver.resetRHS();
	if (divided_velocities || control == Parameters::ControlVariable::stress) {
		if (control == Parameters::ControlVariable::stress) {
			set_shear_rate(1);
		}
		computeUInf();
		setFixedParticleVelocities();
		computeVelocityByComponents();
		if (control == Parameters::ControlVariable::stress) {
			if (p.simulation_mode != 31) {
				computeShearRate();
			} else {
				computeShearRateWalls();
			}
			rescaleVelHydroStressControlled();
		}
		sumUpVelocityComponents();
		// checkForceBalance();
	} else {
		computeUInf();
		setFixedParticleVelocities();
		computeVelocityWithoutComponents();
	}
	/*
	 * The max velocity is used to find dt from max displacement
	 * at each time step.
	 */
	if (!p.fixed_dt && in_predictor) {
		computeMaxNAVelocity();
	}
	adjustVelocityPeriodicBoundary();
	if (divided_velocities && wall_rheology) {
		if (in_predictor) {
			forceResultantLubricationForce();
		}
	}
	stokes_solver.solvingIsDone();
}

void System::computeVelocitiesStokesDrag()
{
	/**
	 \brief Compute velocities in Stokes-drag simulation.

	 Note: Velocities of particles are simply proportional to the total forces acting on respective particles.
	 When the contact model includes dashpots, Stokes-drag simulation cannot be used.
	 */
	for (int i=0; i<np; i++) {
		na_velocity[i].reset();
		na_ang_velocity[i].reset();
	}
	for (auto &fc: force_components) {
		if (fc.first != "brownian" || in_predictor) {
			CALL_MEMBER_FN(*this, fc.second.getForceTorque)(fc.second.force, fc.second.torque);
		}
		if (fc.first != "brownian") {
			for (unsigned int i=0; i<na_velocity.size(); i++) {
				na_velocity[i] += fc.second.force[i]/stokesdrag_coeff_f[i];
				na_ang_velocity[i] += fc.second.torque[i]/stokesdrag_coeff_t[i];
			}
		} else {
			for (unsigned int i=0; i<na_velocity.size(); i++) {
				/* See the comment given in generateBrownianForces()
				 */
				na_velocity[i] += fc.second.force[i]/stokesdrag_coeff_f_sqrt[i];
				na_ang_velocity[i] += fc.second.torque[i]/stokesdrag_coeff_t_sqrt[i];
			}
		}
	}
	if (brownian && twodimension) {
		rushWorkFor2DBrownian(na_velocity, na_ang_velocity);
	}
	adjustVelocityPeriodicBoundary();
}

void System::computeUInf()
{
	for (int i=0; i<np; i++) {
		u_inf[i].reset();
	}
	if (!zero_shear) {
		for (int i=0; i<np; i++) {
			u_inf[i] = dot(E_infinity, position[i]) + cross(omega_inf, position[i]);
		}
	}
}

void System::adjustVelocityPeriodicBoundary()
{
	if (control == Parameters::ControlVariable::stress) { // in rate control it is already done in computeVelocities()
		computeUInf();
	}
	for (int i=0; i<np; i++) {
		velocity[i] = na_velocity[i];
		ang_velocity[i] = na_ang_velocity[i];
	}
	if (!zero_shear) {
		if (!ext_flow) {
			/**** simple shear flow ****/
			for (int i=0; i<np; i++) {
				velocity[i] += u_inf[i];
				ang_velocity[i] += omega_inf;
			}
		} else {
			/**** extensional flow ****/
			for (int i=0; i<np; i++) {
				velocity[i] += u_inf[i];
			}
		}
	}
}

void System::rushWorkFor2DBrownian(vector<vec3d> &vel, vector<vec3d> &ang_vel)
{
	/* [note]
	 * Native 2D simulation is not implemented yet.
	 * As a quick implementation, the velocity elements for extra dimension are set to zero
	 */
	if (p.monolayer) {
		/* Particle (3D sphere) cannot move along y-direction.
		 * All other degrees of freedom exist.
		 */
		for (int i=0; i<np; i++) {
			vel[i].y = 0; // @@ To be checked
		}
	} else {
		/* Particle (2D disk) can rotate only along y-axis.
		 */
		for (int i=0; i<np; i++) {
			vel[i].y = 0; // @@ To be checked
			ang_vel[i].x = 0;
			ang_vel[i].z = 0;
		}
	}
}

void System::displacement(int i, const vec3d& dr)
{
	position[i] += dr;
	/* Note:
	 * When the position of the particle is periodized,
	 * we need to modify the velocity, which was already evaluated.
	 * The position and velocity will be used to calculate the contact forces.
	 */
	if (!ext_flow) {
		/**** simple shear flow ****/
		int z_shift = periodize(position[i]);
		if (z_shift != 0) {
			velocity[i] += z_shift*vel_difference;
		}
	} else {
		/**** extensional flow ****/
		bool pd_transport = false;
		periodizeExtFlow(i, pd_transport);
		if (!zero_shear && pd_transport) {
			velocity[i] -= u_inf[i];
			u_inf[i] = grad_u*position[i];
			velocity[i] +=  u_inf[i];
		}
	}
	boxset.box(i);
}

void System::periodizeExtFlow(const int &i, bool &pd_transport)
{
	if (boxset.boxType(i) != 1) {
		vec3d s = deform_backward*position[i];
		pd_transport = false;
		if (s.z >= lz) {
			s.z -= lz;
			pd_transport = true;
		} else if (s.z < 0) {
			s.z += lz;
			pd_transport = true;
		}
		if (s.x >= lx) {
			s.x -= lx;
			pd_transport = true;
		} else if (s.x < 0) {
			s.x += lx;
			pd_transport = true;
		}
		if (pd_transport) {
			position[i] = deform_forward*s;
		}
	}
	if (!twodimension) {
		if (position[i].y >= ly) {
			position[i].y -= ly;
		} else if (position[i].y < 0) {
			position[i].y += ly;
		}
	}
}

// [0,l]
int System::periodize(vec3d& pos)
{
	/* Lees-Edwards boundary condition
	 *
	 */
	int z_shift = 0;
	if (pos.z >= lz) {
		pos.z -= lz;
		pos -= shear_disp;
		z_shift = -1;
	} else if (pos.z < 0) {
		pos.z += lz;
		pos += shear_disp;
		z_shift = 1;
	}
	while (pos.x >= lx) {
		pos.x -= lx;
	}
	while (pos.x < 0) {
		pos.x += lx;
	}
	if (!twodimension) {
		if (pos.y >= ly) {
			pos.y -= ly;
		} else if (pos.y < 0) {
			pos.y += ly;
		}
	}
	return z_shift;
}

// [0,l]
vec3d System::periodized(const vec3d &pos)
{
	vec3d periodized_pos = pos;
	periodize(periodized_pos);
	return periodized_pos;
}

int System::periodizeDiff(vec3d& pos_diff)
{
	/** Periodize a separation vector with Lees-Edwards boundary condition

		On return pos_diff is the separation vector corresponding to the closest copies,
		and velocity_offset contains the velocity difference produced by Lees-Edwards between the 2 points.
	 */
	int zshift = 0;
	if (pos_diff.z > lz_half) {
		pos_diff.z -= lz;
		pos_diff -= shear_disp;
		zshift = -1;
	} else if (pos_diff.z < -lz_half) {
		pos_diff.z += lz;
		pos_diff += shear_disp;
		zshift = 1;
	}
	while (pos_diff.x > lx_half) {
		pos_diff.x -= lx;
	}
	while (pos_diff.x < -lx_half) {
		pos_diff.x += lx;
	}
	if (!twodimension) {
		if (pos_diff.y > ly_half) {
			pos_diff.y -= ly;
		} else if (pos_diff.y < -ly_half) {
			pos_diff.y += ly;
		}
	}
	return zshift;
}

void System::periodizeDiffExtFlow(vec3d& pos_diff, vec3d& pd_shift, const int &i, const int &j)
{
	pd_shift = boxset.periodicDiffShift(i, j);
	if (twodimension == false) {
		if (pos_diff.y > ly_half) {
			pd_shift.y -= ly;
		} else if (pos_diff.y < -ly_half) {
			pd_shift.y += ly;
		}
	}
	if (pd_shift.is_nonzero()) {
		pos_diff += pd_shift;
	}
}

void System::setSystemVolume()
{
	string indent = "  System::\t";
	if (z_top == -1) {
		system_height = lz;
	} else {
		/* wall particles are at z = z_bot - a and z_top + a
		 */
		system_height = z_top-z_bot;
	}
	if (twodimension) {
		system_volume = lx*system_height;
		cout << indent << "lx = " << lx << " lz = " << lz << " system_height = " << system_height << endl;
	} else {
		system_volume = lx*ly*system_height;
		cout << indent << "lx = " << lx << " lz = " << lz << " ly = " << ly << endl;
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

	double overlap = -evaluateMinGap(*this);
	overlap_avg.update(overlap, clk.cumulated_strain);
	double max_disp_tan = evaluateMaxDispTan(*this);
	max_disp_tan_avg.update(max_disp_tan, clk.cumulated_strain);
	kn_avg.update(p.kn, clk.cumulated_strain);
	kt_avg.update(p.kt, clk.cumulated_strain);

	static double previous_cumulated_strain = 0;
	double deltagamma = (clk.cumulated_strain-previous_cumulated_strain);
	double kn_target = kn_avg.get()*overlap_avg.get()/p.overlap_target;
	double dkn = (kn_target-p.kn)*deltagamma/p.memory_strain_k;

	p.kn += dkn;
	if (p.kn < p.min_kn_auto_det) {
		p.kn = p.min_kn_auto_det;
	}
	if (p.kn > p.max_kn_auto_det) {
		p.kn = p.max_kn_auto_det;
	}
	if (p.disp_tan_target != -1) {
		double kt_target = kt_avg.get()*max_disp_tan_avg.get()/p.disp_tan_target;
		double dkt = (kt_target-p.kt)*deltagamma/p.memory_strain_k;
		p.kt += dkt;
		if (p.kt < p.min_kt_auto_det) {
			p.kt = p.min_kt_auto_det;
		}
		if (p.kt > p.max_kt_auto_det) {
			p.kt = p.max_kt_auto_det;
		}
	} else {
		p.kt = p.kn;
	}
	adaptTimeStep();
	if (dt < p.min_dt_auto_det) {
		dt = p.min_dt_auto_det;
	}
	if (dt > p.max_dt_auto_det) {
		dt = p.max_dt_auto_det;
	}
	previous_cumulated_strain = clk.cumulated_strain;
	resetContactModelParameer();
}

void System::resetContactModelParameer()
{
	for (auto &inter: interaction) {
		inter.contact.setSpringConstants();
	}
}

double System::calcInteractionRangeDefault(int i, int j)
{
	return p.interaction_range*0.5*(radius[i]+radius[j]);
}

double System::calcLubricationRange(int i, int j)
{
	double rad_ratio = radius[i]/radius[j];
	if (rad_ratio < 2 && rad_ratio > 0.5) {
		return (2+p.lub_max_gap)*0.5*(radius[i]+radius[j]);
	} else {
		double minradius = (radius[i]<radius[j] ? radius[i] : radius[j]);
		return radius[i]+radius[j]+p.lub_max_gap*minradius;
	}
}

void System::retrim(vec3d& pos)
{
	while (pos.z >= lz) {
		pos.z -= lz;
	}
	while (pos.z < 0) {
		pos.z += lz;
	}
	while (pos.x >= lx) {
		pos.x -= lx;
	}
	while (pos.x < 0) {
		pos.x += lx;
	}
}

void System::updateH()
{
	/*
	 * strainH is twice of extensional strain.
	 *
	 */
	/**** extensional flow ****/
	/* H_orig = (exp(epsilon_dot t) 0 0 / 0 1 0 / 0 0 exp(-epsilon_dot t))
	 * H = R(-magicangle) H_orig R(magicangle)
	 */
	double diff_ext_strain = 0.5*(clk.cumulated_strain-(strain_retrim-strain_retrim_interval));
	double exp_strain_x = exp(diff_ext_strain); // extensional strain = 0.5*(shear strain)
	double exp_strain_z = 1.0/exp_strain_x;
	/******** Set H     *****************************************/
	// 00, 01, 02, 12, 11, 22
	deform_forward.setSymmetric(exp_strain_x*sq_cos_ma+exp_strain_z*sq_sin_ma, // 00
								0, // 01
								-(exp_strain_x-exp_strain_z)*cos_ma_sin_ma, // 02
								0, // 12
								1, // 11
								exp_strain_x*sq_sin_ma+exp_strain_z*sq_cos_ma); // 22
	deform_backward = deform_forward.inverse();
	/************** for boxing **************************************/
	box_axis1 = lx*deform_forward.getLine(0); // 6--7 = 10--11
	box_axis2 = lz*deform_forward.getLine(2); // 6--10 = 7--11
	box_diagonal_6_11 = box_axis1+box_axis2; // 6--11
	box_diagonal_7_10 = box_axis2-box_axis1; // 7--10
}

void System::yaplotBoxing(std::ofstream &fout_boxing)
{
	boxset.yaplotBox(fout_boxing);
	vec3d e1;
	vec3d e2;
	vec3d dy(0, 0, 0);
	if (twodimension) {
		dy.y = 0.01;
	}
	e1 = lx*deform_forward.getLine(0);
	e2 = lz*deform_forward.getLine(2);
	fout_boxing << "@ 0" << endl;
	fout_boxing << "r 0.2\n";
	fout_boxing << "s 0 -0.01 0 " << e1 -dy<< endl;
	fout_boxing << "s 0 -0.01 0 " << e2 -dy << endl;
	fout_boxing << "s " << e1 -dy << ' ' << e1+e2 -dy << endl;
	fout_boxing << "s " << e2 -dy << ' ' << e1+e2 -dy << endl;
	fout_boxing << "@ 0" << endl;

	fout_boxing << "r 0.5\n";
	for (int i=0; i<np; i++) {
		//fout_boxing << "@ " << boxset.boxType(i) << endl;
		fout_boxing << "c " << position[i]-2*dy << endl;
	}
	if (1) {
		fout_boxing << "@ 2" << endl;
		int i = 0;
		fout_boxing << "r 1\n";
		fout_boxing << "c " << position[i] << endl;
		fout_boxing << "@ 4" << endl;
		for (auto j : boxset.neighborhood(i)) {
			if (i != j) {
				fout_boxing << "c " << position[j] << endl;
			}
		}
	}
	fout_boxing << "y 5" << endl;
	fout_boxing << "@ 5" << endl;
	fout_boxing << "r 1\n";
	for (auto op: overlap_particles) {
		fout_boxing << "c " << position[op]-3*dy << endl;
	}
	overlap_particles.clear(); // @@@ for debuging

	if (/* DISABLES CODE */ (0)) {
		fout_boxing << "@ 8" << endl;
		for (unsigned int k=0; k<interaction.size(); k++) {
			int p0 = interaction[k].get_p0();
			int p1 = interaction[k].get_p1();
			vec3d nvec = interaction[k].nvec;
			fout_boxing << "l " << position[p0]-3*dy << ' ' << position[p0]-3*dy + nvec  << endl;
			fout_boxing << "l " << position[p1]-3*dy << ' ' << position[p1]-3*dy - nvec  << endl;
		}
	}

	fout_boxing << endl;
}
