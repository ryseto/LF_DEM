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
#include "SolventFlow.h"
#include "AlmostEqual.h"
#include "KraynikReinelt.h"
#include "LeesEdwards.h"
#include "Random.h"
#include "global.h"

#ifdef SIGINT_CATCH
extern volatile sig_atomic_t sig_caught;
#endif

System::System(State::BasicCheckpoint chkp):
np(0),
clk(chkp.clock),
shear_type(ShearType::simple_shear),
twodimension(false),
control(Parameters::ControlVariable::rate),
wall_rheology(false),
mobile_fixed(false),
couette_stress(false),
dt(0),
avg_dt(0),
target_stress(0),
init_strain_shear_rate_limit(0),
init_shear_rate_limit(999),
z_bot(-1),
z_top(-1),
eventLookUp(NULL)
//forbid_displacement(false)
{
}

void System::allocateRessources()
{
	if (np <= 0) {
		throw runtime_error("System::allocateRessources() :  np is 0 or negative, cannot allocate this.");
	}
	radius_cubed.resize(np);
	if (!res_solver) {
		// Stokes-drag simulation
		stokesdrag_coeff_f.resize(np);
		stokesdrag_coeff_f_sqrt.resize(np);
		stokesdrag_coeff_t.resize(np);
		stokesdrag_coeff_t_sqrt.resize(np);
	}
	// Velocity
	velocity = ParticleVelocity(np, VelocityType::total, RateDependence::dependent);
	na_velocity = ParticleVelocity(np, VelocityType::nonaffine, RateDependence::dependent);
	if (p.integration_method == 1) {
		velocity_predictor = ParticleVelocity(np, VelocityType::total, RateDependence::dependent);
	}
	vel_bg = ParticleVelocity(np, VelocityType::total, RateDependence::proportional);
	velgrad_bg = ParticleVelocityGrad(np);

	na_disp.resize(np);
	if (mobile_fixed) {
		non_rate_proportional_wall_force.resize(p.np_fixed);
		non_rate_proportional_wall_torque.resize(p.np_fixed);
		rate_proportional_wall_force.resize(p.np_fixed);
		rate_proportional_wall_torque.resize(p.np_fixed);
	}
	
	// Forces and Stress
	forceResultant.resize(np);
	torqueResultant.resize(np);
	total_stress_pp.resize(np);
	phi6.resize(np);
	n_contact.resize(np);
	if (p.solvent_flow) {
		sflow = new SolventFlow;
		phi_local.resize(np);
		//forbid_displacement = true;
	}
}

void System::declareForceComponents()
{
	// Only declare in force components the forces on the rhs of
	// R_FU*(U-U_inf) = RHS
	// These forces will be used to compute the na_velo_components,
	// a set of components that must add up exactly to the total non-affine velocity

	bool torque = true;

	if (is_brownian()) {
		force_components["brownian"] = ForceComponent(np, RateDependence::independent, torque);
		declared_forces.push_back("brownian");
	}
	if (p.bodyforce > 0) {
		force_components["body_force"] = ForceComponent(np, RateDependence::independent, torque);
		declared_forces.push_back("body_force");
	}
	/********** Force R_FU^{mf}*(U^f-U^f_inf)  *************/
	if (mobile_fixed) {
		throw std::runtime_error("Mobile-fixed to be fixed");
		// rate proportional with walls, but this can change
		if (p.lub.model == "normal") {
			force_components["from_fixed"] = ForceComponent(np, RateDependence::proportional, !torque);
			declared_forces.push_back("from_fixed");
		}
		if (p.lub.model == "tangential") {
			force_components["from_fixed"] = ForceComponent(np, RateDependence::proportional, torque);
			declared_forces.push_back("from_fixed");
		}
	}
	if (p.confinement.on) {
		force_components["confinement"] = ForceComponent(np, RateDependence::independent, !torque);
		declared_forces.push_back("confinement");
	}
    if (p.magnetic_field_type != 0) {
        force_components["magnetic_force"] = ForceComponent(np, RateDependence::independent, torque);
        declared_forces.push_back("magnetic_force");
    }
}

void System::setForceToParticle(const string &component, vector<vec3d> &force, vector<vec3d> &torque)
{
	if (component == "confinement") {
		setConfinementForce(force, torque);
	} else if (component == "body_force") {
		setBodyForce(force, torque);
	} else if (component == "brownian") {
		setBrownianForceToParticle(force, torque);
	} else if (component == "magnetic_force") {
        setMagneticForce(force, torque);
    }
}


void System::declareVelocityComponents()
{
	// Only declare in na_velo_components a set of components that add up
	// exactly to the non-affine velocity
	for (const auto &fc: force_components) {
		na_velo_components[fc.first] = ParticleVelocity(fc.second.force.size(), VelocityType::nonaffine, fc.second.rate_dependence);
	}
}

void System::setConfiguration(const vector <vec3d>& initial_positions,
							  const vector <double>& radii,
							  const vector <double>& angles)
{
	/**
		\brief Set positions of the particles for initialization.
	 */
	string indent = "  System::\t";
	set_np(initial_positions.size());

	conf = std::make_shared<ParticleConfig> (np);
	conf->position = initial_positions;
	conf->radius = radii;
	if (twodimension) {
		conf->angle = angles;
	}
	np_mobile = np-p.np_fixed;
	if (np_mobile <= 0) {
		throw runtime_error("np_fixed>=np");
	}
	if (p.np_fixed > 0) {
		mobile_fixed = true;
		cerr << "np fixed " << p.np_fixed << endl;
	}
	radius_wall_particle = conf->radius[np-1];

	particle_volume = 0;
	if (twodimension) {
		for (auto r: conf->radius) {
			particle_volume += r*r;
		}
		particle_volume *= M_PI;
	} else {
		for (auto r: conf->radius) {
			particle_volume += r*r*r;
		}
		particle_volume *= 4*M_PI/3.;
	}
}

void System::setupParameters()
{
	string indent = "  System::\t";

	if (p.integration_method > 1) {
		ostringstream error_str;
		error_str << indent << "integration_method = " << p.integration_method << endl << indent << "The integration method is not impremented yet." << endl;
		throw runtime_error(error_str.str());
	}
	
	if (p.auto_determine_knkt) {
		kn_avg.setRelaxationTime(p.memory_strain_avg);
		kt_avg.setRelaxationTime(p.memory_strain_avg);
		overlap_avg.setRelaxationTime(p.memory_strain_avg);
		max_disp_tan_avg.setRelaxationTime(p.memory_strain_avg);
	}

	if (is_brownian()) {
		double stress_avg_relaxation_parameter = 10*p.output.time_interval_output_data.value; // 0 --> no average
		stress_avg.setRelaxationTime(stress_avg_relaxation_parameter);
		rate_prop_shearstress_rate1_ave.setRelaxationTime(p.brownian_relaxation_time);
		rate_indep_shearstress_ave.setRelaxationTime(p.brownian_relaxation_time);
		// if ((Interactions::hasPairwiseResistanceStdInteractionStdInteraction(p) || Interactions::hasDimer(p))
		if (Interactions::hasPairwiseResistanceStdInteraction(p)
			&& p.integration_method != 1) {
			ostringstream error_str;
			error_str << "Brownian simulation with multiplicative noise needs to use the Predictor-Corrector method." << endl;
			error_str << "Modify the parameter file." << endl;
			throw runtime_error(error_str.str());
		}
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
	if (wall_rheology) {
		throw std::runtime_error("wall_rheology commented out in System.");
	}
}

template<typename T>
void System::setupGenericConfiguration(T config, Parameters::ControlVariable control_)
{
	string indent = "  System::\t";
	cout << indent << "Setting up System... " << endl;
	
	np = config.position.size();
	np_mobile = np - p.np_fixed;

	control = control_;

	container = {config.lx, config.ly, config.lz};
	twodimension = config.ly == 0;

	setupParameters();
	setConfiguration(config.position, config.radius, config.angle);

	auto max_range = Interactions::maxRangeStdInteraction(p, conf->radius);
	if (shear_type == ShearType::extensional_flow) {
		kr = std::make_shared<BC::KraynikReineltBC>(container, p.magic_angle, imposed_flow);
		pairconf = std::make_shared<Geometry::KraynikReineltPairwiseConfig>(conf, kr, max_range);

	} else {
		lees = std::make_shared<BC::LeesEdwardsBC>(container, config.lees_edwards_disp, !p.keep_input_strain, imposed_flow);
		pairconf = std::make_shared<Geometry::LeesEdwardsPairwiseConfig>(conf, lees, max_range);
	}

	if (Interactions::hasPairwiseResistanceStdInteraction(p)) {
		if (!mobile_fixed) {
			res_solver = std::make_shared<Dynamics::PairwiseResistanceVelocitySolver>(p.sd_coeff, conf->radius);
		} else {
			throw std::runtime_error("Mobile-fixed to be fixed");
		}
	}
	
	// Memory
	allocateRessources();

	total_num_timesteps = 0;

	auto pair_manager = std::make_shared<Interactions::PairManager> (np);
	interaction = std::make_shared<Interactions::StdInteractionManager>(np, 
																		pair_manager,
																		conf.get(), 
																		&vel_bg, 
																		&velgrad_bg, 
																		pairconf, 
																		&p,
																		res_solver,
																		config.contact_states);
	interaction->declareForceComponents(force_components);
	declareForceComponents();
	declareVelocityComponents();
	declareStressComponents();

	setupSystemPostConfiguration();

	cout << indent << "Setting up System... [ok]" << endl;
}

void System::setupConfiguration(struct base_shear_configuration config, Parameters::ControlVariable control_)
{
	setupGenericConfiguration(config, control_);
}

void System::setupConfiguration(struct fixed_velo_configuration config, Parameters::ControlVariable control_)
{
	np_wall1 = config.np_wall1;
	np_wall2 = config.np_wall2;
	//p.np_fixed = conf.fixed_velocities.size();
	p.np_fixed = np_wall1+np_wall2;
	z_bot = config.z_bot;
	z_top = config.z_top;
	setupGenericConfiguration(config, control_);
	throw std::runtime_error("Mobile-fixed to be fixed");	
	// setFixedVelocities(conf.fixed_velocities);
}

void System::setupConfiguration(struct circular_couette_configuration config, Parameters::ControlVariable control_)
{
	p.np_fixed = config.np_wall1 + config.np_wall2;
	np_wall1 = config.np_wall1;
	np_wall2 = config.np_wall2;
	radius_in = config.radius_in;
	radius_out = config.radius_out;
	setupGenericConfiguration(config, control_);
}

void System::addDimers(const std::vector<Interactions::Dimer::DimerState> &dimers,
					   const std::vector <struct Interactions::contact_state>& cs)
{
	// reinitialize interaction machinery, start by dimers, to not create StdInteractions for dimer pairs
	auto pair_manager = std::make_shared<Interactions::PairManager> (np);
	dimer_manager = std::make_shared<Interactions::Dimer::DimerManager>(np, 
																		pair_manager,
																		conf.get(),
																		&vel_bg,
																		pairconf,
																		&p,
																		res_solver,
																		dimers);
	interaction = std::make_shared<Interactions::StdInteractionManager>(np, 
																		pair_manager,
																		conf.get(), 
																		&vel_bg, 
																		&velgrad_bg, 
																		pairconf, 
																		&p,
																		res_solver,
																		cs);
	
	dimer_manager->declareForceComponents(force_components);
	interaction->declareForceComponents(force_components);
	declareVelocityComponents();
	declareStressComponents();
}

void System::addDimers(const std::vector<Interactions::Dimer::UnloadedDimerState> &dimers,
					   const std::vector <struct Interactions::contact_state>& cs)
{
	auto pair_manager = std::make_shared<Interactions::PairManager> (np);
	dimer_manager = std::make_shared<Interactions::Dimer::DimerManager>(np, 
																		pair_manager,
																		conf.get(),
																		&vel_bg,
																		pairconf,
																		&p,
																		res_solver,
																		dimers);
	interaction = std::make_shared<Interactions::StdInteractionManager>(np, 
																		pair_manager,
																		conf.get(), 
																		&vel_bg, 
																		&velgrad_bg, 
																		pairconf, 
																		&p,
																		res_solver,
																		cs);
	
	dimer_manager->declareForceComponents(force_components);
	interaction->declareForceComponents(force_components);
	declareVelocityComponents();
	declareStressComponents();
}

struct base_shear_configuration confConvertBase2Shear(const struct base_configuration &config,
													  vec3d lees_edwards_disp)
{
	struct base_shear_configuration base_shear;
	base_shear.lx = config.lx;
	base_shear.ly = config.ly;
	base_shear.lz = config.lz;
	base_shear.volume_or_area_fraction = config.volume_or_area_fraction;
	base_shear.position = config.position;
	base_shear.radius = config.radius;
	base_shear.angle = config.angle;
	base_shear.contact_states = config.contact_states;

	base_shear.lees_edwards_disp = lees_edwards_disp;
	return base_shear;
}

// void System::setupConfiguration(const struct delayed_adhesion_configuration &config,
// 								Parameters::ControlVariable control_)
// {
// 	setupGenericConfiguration(confConvertBase2Shear(config.base, config.lees_edwards_disp), control_);
// 	Interactions::TActAdhesion::setupInteractions(*interaction, config.adhesion_states, get_time());
// }

void System::setupConfiguration(const struct activated_adhesion_configuration &config,
								Parameters::ControlVariable control_)
{
	setupGenericConfiguration(confConvertBase2Shear(config.base, config.lees_edwards_disp), control_);
	Interactions::ActAdhesion::setupInteractions(*interaction, config.adhesion_states);
}

void System::setupSystemPostConfiguration()
{
	for (int i=0; i<np; i++) {
		radius_cubed[i] = pow(conf->radius[i], 3);
	}
	omega_wheel_in  = 0;
	omega_wheel_out = 0;
	if (p.simulation_mode >= 10 && p.simulation_mode <= 20) {
		origin_of_rotation = {0.5*container.lx, 0, 0.5*container.lz};
		for (int i=np_mobile; i<np; i++) {
			conf->angle[i] = -atan2(conf->position[i].z-origin_of_rotation.z,
							 		 conf->position[i].x-origin_of_rotation.x);
		}
		double omega_wheel = (radius_out-radius_in)*imposed_flow->getRate()/radius_in;
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
		double omega_wheel = (radius_out-radius_in)*imposed_flow->getRate()/radius_out;
		omega_wheel_out = -omega_wheel;
		omega_wheel_in  = omega_wheel*radius_out/radius_in;
	}
	dt = p.dt;
	if (p.fixed_dt) {
		avg_dt = dt;
	}
}

struct base_configuration System::getBaseConfiguration() const
{
	base_configuration config;
	config.lx = container.lx;
	config.ly = container.ly;
	config.lz = container.lz;
	config.volume_or_area_fraction = particle_volume/getSystemVolume();

	config.position = conf->position;
	config.radius = conf->radius;
	if (twodimension) {
		config.angle = conf->angle;
	}
	config.contact_states = interaction->getContacts();
	return config;
}

void System::timeStepBoxing()
{
	/**
		\brief Apply a strain step to the boxing system.
	 */
	if (!imposed_flow->zero_shear()) {
		double strain_increment = imposed_flow->getRate()*dt;
		clk.cumulated_strain += strain_increment;
		if (shear_type == ShearType::simple_shear) {
			lees->incrementStrain(dt);
		}
		if (shear_type == ShearType::extensional_flow) {
			throw std::runtime_error("timeStepBoxing:: Extensional flow to fix");
		}
	} else {
		if (wall_rheology || p.simulation_mode == 31) {
			vec3d strain_increment = 2*dot(imposed_flow->getSymGradU(), {0, 0, 1})*dt;
			clk.cumulated_strain += strain_increment.norm();
			// shear_strain += strain_increment;
			angle_wheel += dt*(omega_wheel_in-omega_wheel_out);
		}
	}
	pairconf->updateAfterDeformation();
}

void System::eventShearJamming()
{
	/**
	 \brief Create an event when the shear rate is less than p.sj_shear_rate
	*/
	static int cnt_jamming = 0;
	if (abs(imposed_flow->getRate()) < p.sj_shear_rate && computeMaxNAVelocity() < p.sj_velocity) {
		cnt_jamming ++;
		if (cnt_jamming > p.sj_check_count) {
			Event ev;
			ev.type = "jammed_shear_rate";
			events.push_back(Event(ev));
		}
	} else {
		cnt_jamming = 0;
	}
}

void System::forceResultantInterpaticleForces()
{
	throw std::runtime_error(" forceResultantInterpaticleForces deactivated ");
// 	auto &contact_force = force_components["contact"].force;
// 	int np_tmp = np;
// 	for (int i=0; i<np_tmp; i++) {
// 		forceResultant[i] += contact_force[i];
// 	}
// 	if (Interactions::has_friction(p.contact.friction_model)) {
// 		auto &contact_torque = force_components["contact"].torque;
// 		for (int i=0; i<np_tmp; i++) {
// 			torqueResultant[i] += contact_torque[i];
// 		}
// 	}
// 	if (repulsiveforce) {
// 		auto &repulsive_force = force_components["repulsion"].force;
// 		for (int i=0; i<np_tmp; i++) {
// 			forceResultant[i] += repulsive_force[i];
// 		}
// 	}
// 	if (delayed_adhesion) {
// 		auto &adhesion_force = force_components["repulsion"].force;
// 		for (int i=0; i<np_tmp; i++) {
// 			forceResultant[i] += adhesion_force[i];
// 		}
// 	}
//    if (magnetic_interaction) {
//        auto &adhesion_force = force_components["repulsion"].force;
//        for (int i=0; i<np_tmp; i++) {
//            forceResultant[i] += adhesion_force[i];
//        }
//    }
}

void System::wallForces()
{
	throw std::runtime_error(" wallForces deactivated ");

// 	if (wall_rheology) {
// 		double max_total_force = 0;
// 		double max_total_torque = 0;
// 		for (int i=0; i<np_mobile; i++) {
// 			if (max_total_force < forceResultant[i].sq_norm()){
// 				max_total_force = forceResultant[i].sq_norm();
// 			}
// 			if (max_total_torque < torqueResultant[i].sq_norm()){
// 				max_total_torque = torqueResultant[i].sq_norm();
// 			}
// 		}
// 		cerr << "force balance: " << sqrt(max_total_force) << endl;
// 		cerr << "torque balance: " << sqrt(max_total_torque) << endl;
// 		if (p.simulation_mode >= 10 && p.simulation_mode <= 20) {
// 			int i_np_1 = np_mobile+np_wall1;
// 			// inner wheel
// 			// Positions of wall particles are at r =
// 			force_normal_wall1 = 0;
// 			double torque_wall1 = 0;
// 			for (int i=np_mobile; i<i_np_1; i++) {
// 				vec3d unitvec_out = origin_of_rotation-position[i];
// 				unitvec_out.y = 0;
// 				unitvec_out.unitvector();
// 				force_normal_wall1 += dot(forceResultant[i], unitvec_out);
// 				vec3d torque_tmp = cross(position[i]-origin_of_rotation, forceResultant[i]);
// 				torque_wall1 += torque_tmp.y+torqueResultant[i].y;
// 			}
// 			force_tang_wall1 = torque_wall1/(radius_in-radius_wall_particle);
// 			// outer wheel
// 			force_normal_wall2 = 0;
// 			double torque_wall2 = 0;
// 			for (int i=i_np_1; i<np; i++) {
// 				vec3d unitvec_out = position[i]-origin_of_rotation;
// 				unitvec_out.y = 0;
// 				unitvec_out.unitvector();
// 				force_normal_wall2 += dot(forceResultant[i], unitvec_out);
// 				vec3d torque_tmp = cross(position[i]-origin_of_rotation, forceResultant[i]);
// 				torque_wall2 += torque_tmp.y+torqueResultant[i].y;
// 			}
// 			force_tang_wall2 = torque_wall2/(radius_out+radius_wall_particle);
// 			cerr << " normal:" << force_normal_wall1 << ' ' << force_normal_wall2 << endl;
// 			cerr << " tangential:" << force_tang_wall1 << ' ' << force_tang_wall2 << ' ' << torque_wall1 << ' ' << torque_wall2 << endl;
// 		} else if (p.simulation_mode > 40) {
// 			int i_np_1 = np_mobile+np_wall1;
// 			// bottom wall
// 			force_tang_wall1 = 0;
// 			force_normal_wall1 = 0;
// 			for (int i=np_mobile; i<i_np_1; i++) { // bottom
// 				force_tang_wall1   += forceResultant[i].x;
// 				force_normal_wall1 += forceResultant[i].z;
// 			}
// 			// top wall
// 			force_tang_wall2 = 0;
// 			force_normal_wall2 = 0;
// 			for (int i=i_np_1; i<np; i++) {
// 				force_tang_wall2   += forceResultant[i].x;
// 				force_normal_wall2 += forceResultant[i].z;
// 			}
// 			cerr << "Ft " <<   force_tang_wall1 << ' ' <<   force_tang_wall2 << endl;
// 			cerr << "Fn " << force_normal_wall1 << ' ' << force_normal_wall2 << endl;
// 		}
// 	}
}

void System::forceResultantReset()
{
	throw std::runtime_error(" forceResultantReset deactivated ");
// 	for (int i=0; i<np; i++) {
// 		forceResultant[i].reset();
// 		torqueResultant[i].reset();
// 	}
}

// void System::checkForceBalance()
// {
// 	// 1st way: does not work: forceResultand != 0
// 	forceResultantReset();
// 	forceResultantInterpaticleForces();
// 	unsigned int i, j;
// 	for (const auto &inter: *interaction) {
// 		std::tie(i, j) = inter.get_par_num();
// 		vec3d lubforce_i = inter.lubrication.getTotalForce();
// 		forceResultant[i] += lubforce_i;
// 		forceResultant[j] -= lubforce_i;
// 	}
// 	// 2nd way: works
// 	forceResultantReset();
// 	forceResultantInterpaticleForces();
// 	forceResultantLubricationForce();
// }

// void System::checkStaticForceBalance()
// {
// 	vector< vector<vec3d> > contact_force;
// 	contact_force.resize(np);
// 	for (auto &inter: *interaction) {
// 		if (inter.contact.is_active()) {
// 			unsigned int i = inter.get_p0();
// 			unsigned int j = inter.get_p1();
// 			if (n_contact[i] >= 2 && n_contact[j] >= 2) {
// 				contact_force[i].push_back(inter.contact.getSpringForce());
// 				contact_force[j].push_back(-inter.contact.getSpringForce());
// 			}
// 		}
// 	}
// 	double max_imbalance = 0;
// 	vec3d total_force;
// 	double total_imbalance = 0;
// 	int cnt = 0;
// 	for (auto &cf : contact_force) {
// 		if (cf.size() >= 2) {
// 			total_force.reset();
// 			double abs_total_force = 0;
// 			for (auto &ff : cf) {
// 				total_force += ff;
// 				abs_total_force += ff.norm();
// 			}
// 			if (abs_total_force > 0) {
// 				double imbalance = total_force.norm()/abs_total_force;
// 				if (imbalance > max_imbalance) {
// 					max_imbalance = imbalance;
// 				}
// 				total_imbalance += imbalance;
// 				cnt++;
// 			}
// 		}
// 	}
// 	if (cnt > 0) {
// 		cerr << "imbalance = " << max_imbalance << ' ' << total_imbalance/cnt << endl;
// 	}
// 	max_force_imbalance = max_imbalance;
// }

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
	if (shear_type == ShearType::simple_shear || shear_type == ShearType::extensional_flow) {
		if (!Interactions::hasPairwiseResistanceStdInteraction(p)) {
			computeNonAffineVelocitiesStokesDrag();
		} else {
			computeNonAffineVelocities(calc_stress, true);
		}
		computeTotalVelocity();
	} else {
		if (!Interactions::hasPairwiseResistanceStdInteraction(p)) {
			computeNonAffineVelocitiesStokesDrag();
		} else {
			computeNonAffineVelocities(calc_stress, true);
		}
		sflow->update();
		adjustVelocitySolventFlow();
	}
	if (wall_rheology && calc_stress) { // @@@@ calc_stress remove????
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
		if (!p.output.out_particle_stress.empty() || couette_stress || p.output.out_gsd) {
			calcTotalStressPerParticle();
		}
	}
	// if (p.output.recording_interaction_history) {
	// 	recordHistory();
	// }
	timeStepMove(time_end, strain_end);
	for (int i=0; i<np; i++) {
		na_disp[i] += na_velocity.vel[i]*dt;
	}
	if (eventLookUp != NULL) {
		(this->*eventLookUp)();
	}
}

//void System::sflowIteration(bool calc_stress)
//{
//	double diff_u = 0;
//	int cnt = 0;
//	do {
//		if (!pairwise_resistance) {
//			computeVelocitiesStokesDrag();
//		} else {
//			computeVelocities(calc_stress, true);
//		}
//		diff_u = sflow->update();
//		if (cnt == 1 || cnt % 100 == 0) {
//			cerr << cnt << ' '  << diff_u << ' '  << sflow->get_pressure_grad_x() << endl;
//		}
//		cnt++;
//		if (cnt > 1000) {
//			break;
//		}
//	} while (diff_u > 1e-3);
//	sflow->pressureController();
//	if (cnt > 1) {
//		cerr << cnt << endl;
//	}
//	adjustVelocitySolventFlow();
//	return ;
//}

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
	 + \f$ \bm{U}^{+} = \bm{A}^{-1}( \bm{X}' ) ( \bm{F}_\mathrm{B} + \bm{F} ( \bm{X}' ) )  \f$ (\b same \f$\bm{F}_\mathrm{B}\f$ as in the first step)
	 + \f$ \bm{X}(t + dt) = \bm{X}(t) + \frac{1}{2}(\bm{U}^{+}+\bm{U}^{-})dt =  \bm{X}' + \frac{1}{2}(\bm{U}^{+}-\bm{U}^{-})dt \f$

	 */
	/* predictor */
	in_predictor = true;
	in_corrector = false;

	if (Interactions::hasPairwiseResistanceStdInteraction(p)) {
		computeNonAffineVelocities(calc_stress, true); // divided velocities for stress calculation
	} else {
		computeNonAffineVelocitiesStokesDrag();
	}
	if (!p.solvent_flow) {
		computeTotalVelocity();
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
	if (Interactions::hasPairwiseResistanceStdInteraction(p)) {
		computeNonAffineVelocities(calc_stress, true);
	} else {
		computeNonAffineVelocitiesStokesDrag();
	}
	if (!p.solvent_flow) {
		computeTotalVelocity();
	}
	if (calc_stress) {
		calcStressPerParticle(); // stress compornents
		if (wall_rheology) {
			calcStress();
		}
		if (!p.output.out_particle_stress.empty() || couette_stress || p.output.out_gsd) {
			calcTotalStressPerParticle();
		}
	}
	timeStepMoveCorrector();
	for (int i=0; i<np; i++) {
		na_disp[i] += na_velocity.vel[i]*dt;
	}
	if (eventLookUp != NULL) {
		(this->*eventLookUp)();
	}
}

void System::adaptTimeStepWithVelocities()
{
	/**
	 \brief Adapt the time step so that the maximum relative displacement is p.disp_max .
	 */

	/*
	 * The max velocity is used to find dt from max displacement
	 * at each time step.
	 */
	double max_na_velocity = computeMaxNAVelocity();
	double max_interaction_vel = interaction->getMaxRelativeVelocity(&velocity);
	if (dimer_manager) {
		double max_dimer_vel = dimer_manager->getMaxRelativeVelocity(&velocity);
		if (max_dimer_vel > max_interaction_vel) {
			max_interaction_vel = max_dimer_vel;
		}
	}
	if (max_na_velocity > 0 || max_interaction_vel > 0) { // small density system can have na_velocity=0
		if (max_na_velocity > max_interaction_vel) {
			dt = p.disp_max/max_na_velocity;
		} else {
			dt = p.disp_max/max_interaction_vel;
		}
	} else {
		dt = p.disp_max/imposed_flow->getRate();
	}
	if (dt*imposed_flow->getRate() > p.disp_max) { // cases where na_velocity < \dotgamma*radius
		dt = p.disp_max/imposed_flow->getRate();
	}
	if (p.dt_max > 0) {
		if (dt > p.dt_max) {
			dt = p.dt_max;
		}
	}
	if (p.solvent_flow) {
		double dt_sflow = p.disp_max/max_na_velocity;
		if (dt_sflow < dt) {
			dt = dt_sflow;
		}
	}
}

void System::adaptTimeStepWithBounds(double time_end, double strain_end)
{
	/**
	 \brief Adapt the time step so that (a) the maximum relative displacement is p.disp_max, and (b) time or strain does not get passed the end value.
	 */
	// To stop exactly at t == time_end or strain == strain_end,
	// whatever comes first
	if (shear_type != ShearType::solvent_flow) {
		if (strain_end >= 0) {
			if (fabs(dt*imposed_flow->getRate()) > strain_end-clk.cumulated_strain) {
				dt = fabs((strain_end-clk.cumulated_strain)/imposed_flow->getRate());
			}
		}
		if (time_end >= 0) {
			if (get_time()+dt > time_end) {
				dt = time_end-get_time();
			}
		}
	}
}

void System::timeStepMove(double time_end, double strain_end)
{
	/**
	 \brief Moves particle positions according to previously computed velocities, Euler method step.
	 */
	/* @@@@ NOTE
	 * shear_rate = 1 for both simple shear and extensional flow (???)
	 * dot_epsion = shear_rate / 2 is always true.
	 * clk.cumulated_strain = shear_rate * t for both simple shear and extensional flow.
	 */
	if (!p.fixed_dt) {
		adaptTimeStepWithVelocities();
	}
	adaptTimeStepWithBounds(time_end, strain_end);
	clk.time_ += dt;
	total_num_timesteps ++;
	/* evolve PBC */
	timeStepBoxing();
	
	/* move particles */
	//	if (!forbid_displacement) {
	for (int i=0; i<np; i++) {
		displacement(i, velocity.vel[i]*dt);
	}
	if (twodimension) {
		for (int i=0; i<np; i++) {
			conf->angle[i] += velocity.ang_vel[i].y*dt;
		}
	}
	if (shear_type == ShearType::extensional_flow) {
		if (almost_equal(clk.cumulated_strain, kr->getStrainRetrim(), 2)) {
			cerr << "clk.cumulated_strain = " << clk.cumulated_strain << endl;
			kr->retrimProcess(conf->position, clk.cumulated_strain);
			for (int i=0; i<np; i++) {
				velocity.vel[i] -= vel_bg.vel[i];
				vel_bg.vel[i] = imposed_flow->getGradU()*conf->position[i];
				velocity.vel[i] += vel_bg.vel[i];
			}
		}
	}
	interaction->checkNewInteractions();
	interaction->updateInteractions(dt, &velocity);
	if (dimer_manager) {
		dimer_manager->updateInteractions(dt, &velocity);
	}
}

void System::timeStepMovePredictor(double time_end, double strain_end)
{
	/**
	 \brief Moves particle positions according to previously computed velocities, predictor step.
	 */
	if (!p.fixed_dt) {
		adaptTimeStepWithVelocities();
	}
	adaptTimeStepWithBounds(time_end, strain_end);

	clk.time_ += dt;
	total_num_timesteps ++;
	/* evolve PBC
	 * The periodic boundary condition is updated in predictor.
	 * It must not be updated in corrector.
	 */
	timeStepBoxing();
	for (int i=0; i<np; i++) {
		displacement(i, velocity.vel[i]*dt);
	}
	
	if (twodimension) {
		for (int i=0; i<np; i++) {
			conf->angle[i] += velocity.ang_vel[i].y*dt;
		}
	}
	interaction->saveState();
	interaction->updateInteractions(dt, &velocity);
	if (dimer_manager) {
		dimer_manager->saveState();
		dimer_manager->updateInteractions(dt, &velocity);
	}
	/*
	 * Keep V^{-} to use them in the corrector.
	 */
	for (int i=0; i<np; i++) {
		velocity_predictor.vel[i] = velocity.vel[i];
		velocity_predictor.ang_vel[i] = velocity.ang_vel[i];
	}
}

void System::timeStepMoveCorrector()
{
	/**
	 \brief Moves particle positions according to previously computed velocities, corrector step.
	 */
	for (int i=0; i<np; i++) {
		velocity.vel[i] = 0.5*(velocity.vel[i]+velocity_predictor.vel[i]); // real velocity, in predictor and in corrector
		velocity.ang_vel[i] = 0.5*(velocity.ang_vel[i]+velocity_predictor.ang_vel[i]);
	}
	for (int i=0; i<np; i++) {
		displacement(i, (velocity.vel[i]-velocity_predictor.vel[i])*dt);
	}
	if (twodimension) {
		for (int i=0; i<np; i++) {
			conf->angle[i] += (velocity.ang_vel[i].y-velocity_predictor.ang_vel[i].y)*dt; // no cross_shear in 2d
		}
	}
	if (shear_type == ShearType::extensional_flow) {
		if (almost_equal(clk.cumulated_strain, kr->getStrainRetrim(), 2)) {
			cerr << "clk.cumulated_strain = " << clk.cumulated_strain << endl;
			kr->retrimProcess(conf->position, clk.cumulated_strain);
			for (int i=0; i<np; i++) {
				velocity.vel[i] -= vel_bg.vel[i];
				vel_bg.vel[i] = imposed_flow->getGradU()*conf->position[i];
				velocity.vel[i] += vel_bg.vel[i];
			}
		}
	}
	interaction->restoreState();
	interaction->checkNewInteractions();
	interaction->updateInteractions(dt, &velocity);
	if (dimer_manager) {
		dimer_manager->restoreState();
		dimer_manager->updateInteractions(dt, &velocity);
	}
}

bool System::keepRunning(double time_end, double strain_end)
{
	if (clk.cumulated_strain > strain_end-1e-8 && strain_end >= 0) {
		return false;
	}
	if (get_time() > time_end-1e-8 && time_end >= 0) {
		return false;
	}
	if (!events.empty()) {
		return false;
	}
	return true;
}

void System::calculateForces()
{
	interaction->checkNewInteractions();
	in_predictor = true;
	interaction->updateInteractions(0, &velocity);
	if (dimer_manager) {
		dimer_manager->updateInteractions(0, &velocity);
	}
	in_predictor = false;
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

	double loop_time_adjust = 0;

	if (p.fixed_dt == true && !firsttime) {
		loop_time_adjust = dt;
	}

	if (firsttime) {
		calculateForces();
		firsttime = false;
	}

	bool calc_stress = false;
	if (is_brownian()) {
		calc_stress = true;
	}
	avg_dt = 0;
	avg_dt_nb = 0;
	while (keepRunning(time_end-loop_time_adjust, strain_end-loop_time_adjust)) {
		if (p.fixed_dt) {
			dt = p.dt;
		}
		double time_bound, strain_bound;
		if (is_brownian()) {
			time_bound = -1;
			strain_bound = -1;
		} else if (shear_type == ShearType::extensional_flow) {
			time_bound = -1;
			strain_bound = kr->getStrainRetrim();
		} else {
			time_bound = time_end;
			strain_bound = strain_end;
		}
		if (p.integration_method == 0) {
			timeEvolutionEulersMethod(calc_stress, time_bound, strain_bound);
		} else if (p.integration_method == 1) {
			timeEvolutionPredictorCorrectorMethod(calc_stress, time_bound, strain_bound);
		}
		avg_dt += dt;
		avg_dt_nb++;
#ifdef SIGINT_CATCH
		if (sig_caught == SIGINT) { // return to Simulation immediatly for checkpointing
			return;
		}
#endif
	};
	if (avg_dt_nb > 0) {
		avg_dt /= avg_dt_nb;
	} else {
		avg_dt = dt;
	}
	if (!(!events.empty() || 
		  (shear_type == ShearType::extensional_flow && almost_equal(clk.cumulated_strain, kr->getStrainRetrim(), 2) ))) {
		if (shear_type != ShearType::solvent_flow) {
			calc_stress = true;
		}
		if (p.integration_method == 0) {
			timeEvolutionEulersMethod(calc_stress, time_end, strain_end);
		} else if (p.integration_method == 1) {
			timeEvolutionPredictorCorrectorMethod(calc_stress, time_end, strain_end);
		}
	}
	if (p.auto_determine_knkt
		&& clk.cumulated_strain > p.start_adjust) {
		adjustContactModelParameters();
	}
}

// void System::computeForcesOnWallParticles()
// {
// 	/**
// 		\brief This method computes the force (and torque, for now, but it might be dropped)
// 		on the fixed particles.

// 		It is designed with simple shear with walls under stress controlled conditions in mind,
// 		so it decomposes the force in a rate-proportional part and a rate-independent part.

// 		*/
// 	throw runtime_error(" Control stress with walls disabled for now .\n");
// 	if (!zero_shear) {
// 		throw runtime_error(" Stress-control with walls requires zero_shear==true .\n");
// 	}
// 	vector<vec3d> force (p.np_fixed);
// 	vector<vec3d> torque (p.np_fixed);

// 	// Compute the part of the velocity of mobile particles
// 	// that is not coming from the wall velocities
// 	vector<vec3d> na_velocity_mobile (np_mobile);
// 	vector<vec3d> na_ang_velocity_mobile (np_mobile);
// 	const auto &vel_contact = na_velo_components["contact"].vel;
// 	const auto &ang_vel_contact = na_velo_components["contact"].ang_vel;

// 	for (int i=0; i<np_mobile; i++) {
// 		na_velocity_mobile[i] = vel_contact[i];
// 		na_ang_velocity_mobile[i] = ang_vel_contact[i];
// 	}
// 	if (repulsiveforce) {
// 		const auto &vel_repulsive = na_velo_components["repulsion"].vel;
// 		const auto &ang_vel_repulsive = na_velo_components["repulsion"].ang_vel;
// 		for (int i=0; i<np_mobile; i++) {
// 			na_velocity_mobile[i] += vel_repulsive[i];
// 			na_ang_velocity_mobile[i] += ang_vel_repulsive[i];
// 		}
// 	}
// 	if (delayed_adhesion) {
// 		const auto &vel_adhesion = na_velo_components["delayed_adhesion"].vel;
// 		const auto &ang_vel_adhesion = na_velo_components["delayed_adhesion"].ang_vel;
// 		for (int i=0; i<np_mobile; i++) {
// 			na_velocity_mobile[i] += vel_adhesion[i];
// 			na_ang_velocity_mobile[i] += ang_vel_adhesion[i];
// 		}
// 	}
// 	// from this, we can compute the hydro force on the wall that does *not* depend on the wall velocity
// 	stokes_solver.multiply_by_RFU_fm(na_velocity_mobile, na_ang_velocity_mobile, force, torque);

// 	// Now we sum up this hydro part with the other non-rate dependent forces (contact, etc)
// 	const auto &contact_force = force_components["contact"].force;
// 	const auto &contact_torque = force_components["contact"].torque;
// 	const auto &repulsive_force = force_components["repulsion"].force;
// 	const auto &adhesion_force = force_components["delayed_adhesion"].force;
// 	throw std::runtime_error("computeForcesOnWallParticles : ongoing work.... sorry");
// 	for (int i=0; i<p.np_fixed; i++) {
// 		non_rate_proportional_wall_force[i] = -force[i];
// 		non_rate_proportional_wall_torque[i] = -torque[i];
// 		non_rate_proportional_wall_force[i] += contact_force[i+np_mobile];
// 		non_rate_proportional_wall_torque[i] += contact_torque[i+np_mobile];
// 		if (repulsiveforce) {
// 			non_rate_proportional_wall_force[i] += repulsive_force[i+np_mobile];
// 		}
// 	}

// 	// Now the part proportional to the wall speed

// 	// From the mobile particles
// 	const auto &vel_hydro_from_fixed = na_velo_components["from_fixed"].vel;
// 	const auto &ang_vel_hydro_from_fixed = na_velo_components["from_fixed"].ang_vel;
// 	for (int i=0; i<np_mobile; i++) {
// 		na_velocity_mobile[i] = vel_hydro_from_fixed[i];
// 		na_ang_velocity_mobile[i] = ang_vel_hydro_from_fixed[i];
// 	}
// 	stokes_solver.multiply_by_RFU_fm(na_velocity_mobile, na_ang_velocity_mobile, force, torque);
// 	for (int i=0; i<p.np_fixed; i++) {
// 		rate_proportional_wall_force[i] = -force[i];
// 		rate_proportional_wall_torque[i] = -torque[i];
// 	}

// 	// From the fixed particles themselves. This should be zero if these particles form a wall
// 	// (i.e. they move with zero relative velocity) and if the Stokes drag is zero (which is controlled by sd_coeff)
// 	// As we do not want to make too many assumptions here (especially regarding the Stokes drag)
// 	// we compute it. [Probably a p.no_stokes_drag should be introduced at some point.]
// 	vector<vec3d> na_velocity_fixed (p.np_fixed);
// 	vector<vec3d> na_ang_velocity_fixed (p.np_fixed);
// 	for (int i=0; i<p.np_fixed; i++) {
// 		na_velocity_fixed[i] = na_velocity[i+np_mobile];
// 		na_ang_velocity_fixed[i] = na_ang_velocity[i+np_mobile];
// 	}
// 	stokes_solver.multiply_by_RFU_ff(na_velocity_fixed, na_ang_velocity_fixed, force, torque);
// 	for (int i=0; i<p.np_fixed; i++) {
// 		rate_proportional_wall_force[i] -= force[i];
// 		rate_proportional_wall_torque[i] -= torque[i];
// 	}
// }

void System::forceResultantLubricationForce()
{
	throw std::runtime_error(" forceResultantLubricationForce deactivated ");
// 	/* Only F = R_FU U is calculated, but R_FE E is not implemented yet.
// 	 * So, we cannot check the force balance with E^{inf} yet.
// 	 */
// 	/*
// 	 *  F^{M} = R_FU^{MM} U^{M}
// 	 */
// 	vector<double> force_m_to_m (6*np_mobile);
// 	vector<double> minus_mobile_velocities (6*np_mobile);
// 	for (int i=0; i<np_mobile; i++) {
// 		int i6 = 6*i;
// 		minus_mobile_velocities[i6  ] = -na_velocity[i].x;
// 		minus_mobile_velocities[i6+1] = -na_velocity[i].y;
// 		minus_mobile_velocities[i6+2] = -na_velocity[i].z;
// 		minus_mobile_velocities[i6+3] = -na_ang_velocity[i].x;
// 		minus_mobile_velocities[i6+4] = -na_ang_velocity[i].y;
// 		minus_mobile_velocities[i6+5] = -na_ang_velocity[i].z;
// 	}
// 	stokes_solver.multiply_by_RFU_mm(minus_mobile_velocities, force_m_to_m);
// 	for (int i=0; i<np_mobile; i++) {
// 		int i6 = 6*i;
// 		forceResultant[i].x += force_m_to_m[i6];
// 		forceResultant[i].y += force_m_to_m[i6+1];
// 		forceResultant[i].z += force_m_to_m[i6+2];
// 		torqueResultant[i].x += force_m_to_m[i6+3];
// 		torqueResultant[i].y += force_m_to_m[i6+4];
// 		torqueResultant[i].z += force_m_to_m[i6+5];
// 	}
// 	if (mobile_fixed) {
// 		/*
// 		 *  F^{M} += R_FU^{MF} U^{F}
// 		 */
// 		vector<double> force_f_to_m (6*np_mobile);
// 		vector<double> minus_fixed_velocities (6*p.np_fixed);
// 		for (int i=0; i<p.np_fixed; i++) {
// 			int i6 = 6*i;
// 			int i_fixed = i+np_mobile;
// 			minus_fixed_velocities[i6  ] = -na_velocity[i_fixed].x;
// 			minus_fixed_velocities[i6+1] = -na_velocity[i_fixed].y;
// 			minus_fixed_velocities[i6+2] = -na_velocity[i_fixed].z;
// 			minus_fixed_velocities[i6+3] = -na_ang_velocity[i_fixed].x;
// 			minus_fixed_velocities[i6+4] = -na_ang_velocity[i_fixed].y;
// 			minus_fixed_velocities[i6+5] = -na_ang_velocity[i_fixed].z;
// 		}
// 		stokes_solver.multiply_by_RFU_mf(minus_fixed_velocities, force_f_to_m);
// 		for (int i=0; i<np_mobile; i++) {
// 			int i6 = 6*i;
// 			forceResultant[i].x += force_f_to_m[i6];
// 			forceResultant[i].y += force_f_to_m[i6+1];
// 			forceResultant[i].z += force_f_to_m[i6+2];
// 			torqueResultant[i].x += force_f_to_m[i6+3];
// 			torqueResultant[i].y += force_f_to_m[i6+4];
// 			torqueResultant[i].z += force_f_to_m[i6+5];
// 		}
// 		/*
// 		 *  F^{F} += R_FU^{FM} U^{M}
// 		 */
// 		vector<double> force_m_to_f (6*p.np_fixed);
// 		for (int i=0; i<np_mobile; i++) {
// 			int i6 = 6*i;
// 			minus_mobile_velocities[i6  ] = -na_velocity[i].x;
// 			minus_mobile_velocities[i6+1] = -na_velocity[i].y;
// 			minus_mobile_velocities[i6+2] = -na_velocity[i].z;
// 			minus_mobile_velocities[i6+3] = -na_ang_velocity[i].x;
// 			minus_mobile_velocities[i6+4] = -na_ang_velocity[i].y;
// 			minus_mobile_velocities[i6+5] = -na_ang_velocity[i].z;
// 		}
// 		stokes_solver.multiply_by_RFU_fm(minus_mobile_velocities, force_m_to_f);
// 		for (int i=np_mobile; i<np; i++) {
// 			int i6 = 6*(i-np_mobile);
// 			forceResultant[i].x += force_m_to_f[i6];
// 			forceResultant[i].y += force_m_to_f[i6+1];
// 			forceResultant[i].z += force_m_to_f[i6+2];
// 			torqueResultant[i].x += force_m_to_f[i6+3];
// 			torqueResultant[i].y += force_m_to_f[i6+4];
// 			torqueResultant[i].z += force_m_to_f[i6+5];
// 		}
// 		/*
// 		 *  F^{F} += R_FU^{FF} U^{F}
// 		 */
// 		vector<double> force_f_to_f (6*p.np_fixed);
// 		for (int i=0; i<p.np_fixed; i++) {
// 			int i6 = 6*i;
// 			int i_fixed = i+np_mobile;
// 			minus_fixed_velocities[i6  ] = -na_velocity[i_fixed].x;
// 			minus_fixed_velocities[i6+1] = -na_velocity[i_fixed].y;
// 			minus_fixed_velocities[i6+2] = -na_velocity[i_fixed].z;
// 			minus_fixed_velocities[i6+3] = -na_ang_velocity[i_fixed].x;
// 			minus_fixed_velocities[i6+4] = -na_ang_velocity[i_fixed].y;
// 			minus_fixed_velocities[i6+5] = -na_ang_velocity[i_fixed].z;
// 		}
// 		stokes_solver.multiply_by_RFU_ff(minus_fixed_velocities, force_f_to_f);
// 		for (int i=np_mobile; i<np; i++) {
// 			int i6 = 6*(i-np_mobile);
// 			forceResultant[i].x += force_f_to_f[i6];
// 			forceResultant[i].y += force_f_to_f[i6+1];
// 			forceResultant[i].z += force_f_to_f[i6+2];
// 			torqueResultant[i].x += force_f_to_f[i6+3];
// 			torqueResultant[i].y += force_f_to_f[i6+4];
// 			torqueResultant[i].z += force_f_to_f[i6+5];
// 		}
// 	}
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
	if (!in_predictor) { // The Brownian force must be the same in the predictor and the corrector
		return;
	}
	if (mobile_fixed) {
		cerr << "Brownian algorithm with fixed particles not implemented yet.\n";
//		throw runtime_error("Brownian algorithm with fixed particles not implemented yet.\n");
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

	if (Interactions::hasPairwiseResistanceStdInteraction(p)) {
		/* L*L^T = RFU
		 */
		res_solver->setSolverRHS(force, torque);
		res_solver->compute_LTRHS(force, torque); // F_B = \sqrt(2kT/dt) * L^T * A
		res_solver->resetRHS();
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
		 *  In the function computeNonAffineVelocitiesStokesDrag(),
		 *  U_B = F_B / sqrt(RFU)
		 */
	}
}

void System::setBodyForce(vector<vec3d> &force,
						  vector<vec3d> &torque)
{
	for (auto &t: torque) {
		t.reset();
	}
	double angle = M_PI*p.body_force_angle/180;
	double bf_x = p.bodyforce*cos(angle); // cos(angle);
	double bf_z = -p.bodyforce*sin(angle);
	for (int i=0; i<np_mobile; i++) {
		force[i].set(radius_cubed[i]*bf_x, 0 , radius_cubed[i]*bf_z);
	}
	for (int i=np_mobile; i<np; i++) {
		force[i].reset();
	}
}

void System::setConfinementForce(vector<vec3d> &force,
								  vector<vec3d> &torque)
{
	for (auto &t: torque) {
		t.reset();
	}
	for (int i=0; i<np_mobile; i++) {
		if (conf->position[i].y > p.confinement.y_max) {
			force[i].set(0, -p.confinement.k.value*(conf->position[i].y-p.confinement.y_max), 0);
		} else if (conf->position[i].y < p.confinement.y_min) {
			force[i].set(0, -p.confinement.k.value*(conf->position[i].y-p.confinement.y_min), 0);
		} else {
			force[i].reset();
		}
	}
}

// void System::setFixedParticleForceToParticle(vector<vec3d> &force,
// 											 vector<vec3d> &torque)
// {
// 	vector<double> force_torque_from_fixed (6*np_mobile);
// 	// @@ TODO: avoid copy of the velocities and forces
// 	vector<double> minus_fixed_velocities (6*p.np_fixed);
// 	for (int i=0; i<p.np_fixed; i++) {
// 		int i6 = 6*i;
// 		int i_fixed = i+np_mobile;
// 		minus_fixed_velocities[i6  ] = -na_velocity[i_fixed].x;
// 		minus_fixed_velocities[i6+1] = -na_velocity[i_fixed].y;
// 		minus_fixed_velocities[i6+2] = -na_velocity[i_fixed].z;
// 		minus_fixed_velocities[i6+3] = -na_ang_velocity[i_fixed].x;
// 		minus_fixed_velocities[i6+4] = -na_ang_velocity[i_fixed].y;
// 		minus_fixed_velocities[i6+5] = -na_ang_velocity[i_fixed].z;
// 	}
// 	stokes_solver.multiply_by_RFU_mf(minus_fixed_velocities, force_torque_from_fixed); // -R_FU^mf*fixed_velocities
// 	for (unsigned int i=0; i<np_mobile; i++) {
// 		auto i6 = 6*i;
// 		force[i].x = force_torque_from_fixed[i6  ];
// 		force[i].y = force_torque_from_fixed[i6+1];
// 		force[i].z = force_torque_from_fixed[i6+2];
// 		torque[i].x = force_torque_from_fixed[i6+3];
// 		torque[i].y = force_torque_from_fixed[i6+4];
// 		torque[i].z = force_torque_from_fixed[i6+5];
// 	}
// }

void System::setMagneticForce(vector<vec3d> &force, vector<vec3d> &torque)
{
    magnetic_field_gradient.resize(3,3);
    if (p.magnetic_field_type == 1) {
        magnetic_field.set(0,0,1);
    } else if (p.magnetic_field_type == 2) {
        magnetic_field.set(0,0,1);                      // Need to modify to satisfy Maxwell's equation
        magnetic_field_gradient[0].set(0,0,1);          // Better to read external file including magnetic field and gradient
        magnetic_field_gradient[1].set(0,0,0);
        magnetic_field_gradient[2].set(0,0,0);
    } else if (p.magnetic_field_type == 3) {
        magnetic_field.set(0,0,cos(p.magnetic_field_freq*get_time()));
    } else if (p.magnetic_field_type == 4) {
        magnetic_field.set(sin(p.magnetic_field_freq*get_time()),0,cos(p.magnetic_field_freq*get_time()));
    }
    
    for (int i=0; i<np; i++)
    {
        if (twodimension) {
            magnetic_field_gradient.resize(np,3);
            magnetic_dipole_moment[i].set(cos(conf->angle[i]),0,sin(conf->angle[i]));
        }
        if (p.magnetic_field_type == 2) {
            double force_x = p.langevin_parameter*dot(magnetic_dipole_moment[i],magnetic_field_gradient[0]);
            double force_y = p.langevin_parameter*dot(magnetic_dipole_moment[i],magnetic_field_gradient[1]);
            double force_z = p.langevin_parameter*dot(magnetic_dipole_moment[i],magnetic_field_gradient[2]);
            force[i].set(force_x,force_y,force_z);
        } else {
            force[i].set(0,0,0);
        }
        torque[i] = p.langevin_parameter*cross(magnetic_dipole_moment[i],magnetic_field);
    }
}

vec3d System::get_shear_strain() const
{
	if (shear_type == ShearType::simple_shear) {
		return lees->getShearStrain();
	}
	if (shear_type == ShearType::extensional_flow) {
		throw std::runtime_error("getShearStrain for Kraynik-Reinelt?");
	}
	return 0;
}


double System::computeMaxNAVelocity()
{
	/**
	 \brief Compute the maximum non-affine velocity

	 Note: it does \b not compute the velocities, just takes the maximum.
	 */

	double sq_max_na_velocity = 0;
	for (int i=0; i<np; i++) {
		auto sq_na_velocity = na_velocity.vel[i].sq_norm();
		if (sq_na_velocity > sq_max_na_velocity) {
			sq_max_na_velocity = sq_na_velocity;
		}
	}
	return sqrt(sq_max_na_velocity);
}

double System::computeMaxVelocity()
{
	/**
	 \brief Compute the maximum non-affine velocity
	 
	 Note: it does \b not compute the velocities, just takes the maximum.
	 */
	
	double sq_max_velocity = 0;
	for (int i=0; i<np; i++) {
		auto sq_velocity = velocity.vel[i].sq_norm();
		if (sq_velocity > sq_max_velocity) {
			sq_max_velocity = sq_velocity;
		}
	}
	return sqrt(sq_max_velocity);
}

void System::computeVelocityWithoutComponents(bool rebuild)
{
	if (rebuild) {
		if (!dimer_manager) {
			res_solver->buildResistanceMatrix(*interaction);
		} else {
			res_solver->buildResistanceMatrix(*interaction, *dimer_manager);
		}
	}
	for (const auto &component: interaction->declared_forces) {
		interaction->setForceToParticle(component, force_components[component].force, force_components[component].torque);
		res_solver->addToSolverRHS(force_components[component]);
	}
	if (dimer_manager) {
		for (const auto &component: dimer_manager->declared_forces) {
			dimer_manager->setForceToParticle(component, force_components[component].force, force_components[component].torque);
			res_solver->addToSolverRHS(force_components[component]);
		}
	}
	for (const auto &component: declared_forces) {
		if (component != "brownian" || in_predictor) {
			setForceToParticle(component, force_components[component].force, force_components[component].torque);
		}
		res_solver->addToSolverRHS(force_components[component]);
	}
	res_solver->solve(na_velocity.vel, na_velocity.ang_vel); // get V
	if (is_brownian() && twodimension) {
		rushWorkFor2DBrownian(na_velocity.vel, na_velocity.ang_vel);
	}
}

void System::computeVelocityByComponents()
{
	/**
	 \brief Compute velocities component by component.
	 */
	if (!dimer_manager) {
		res_solver->buildResistanceMatrix(*interaction);
	} else {
		res_solver->buildResistanceMatrix(*interaction, *dimer_manager);
	}
	for (const auto &component: interaction->declared_forces) {
		interaction->setForceToParticle(component, force_components[component].force, force_components[component].torque);
		res_solver->setSolverRHS(force_components[component]);
		res_solver->solve(na_velo_components[component].vel,
							na_velo_components[component].ang_vel);
	}
	if (dimer_manager) {
		for (const auto &component: dimer_manager->declared_forces) {
			dimer_manager->setForceToParticle(component, force_components[component].force, force_components[component].torque);
			res_solver->setSolverRHS(force_components[component]);
			res_solver->solve(na_velo_components[component].vel,
								na_velo_components[component].ang_vel);
		}
	}
	for (const auto &component: declared_forces) {
		if (component != "brownian" || in_predictor) {
			setForceToParticle(component, force_components[component].force, force_components[component].torque);
		}
		res_solver->setSolverRHS(force_components[component]);
		res_solver->solve(na_velo_components[component].vel,
							na_velo_components[component].ang_vel);
	}
	if (is_brownian() && twodimension) {
		rushWorkFor2DBrownian(na_velo_components["brownian"].vel,
							  na_velo_components["brownian"].ang_vel);
	}
}

void System::computeShearRate()
{
	/**
	 \brief Compute the shear rate under stress control conditions.
	 */
	assert(abs(imposed_flow->getRate()-1) < 1e-15);
	calcStressPerParticle();
	Sym2Tensor rate_prop_stress;
	Sym2Tensor rate_indep_stress;
	gatherStressesByRateDependencies(rate_prop_stress, rate_indep_stress);
	/*
	 *  target_stress = rate_indep_stress + rate * rate_prop_stress_at_1
	 *  rate = (target_stress - rate_indep_stress)/rate_prop_stress_at_1
	 *
	 *  kappa = (Sigma.E)/(E.E)
	 *  sigama = -p*I + 2 eta*E + ... = -p*I + kappa*E + ...
	 *  E.E = 1/2 (if shear rate = 1)
	 *  eta = kappa / 2 = (Sigma.E)/(2*E.E) = Sigma.E
	 *  sigma_target = rate_prop_shearstress(rate=1)*(rate/1) + rate_indep_shearstress
	 *  rate = (sigma_target - rate_indep_shearstress) / rate_prop_shearstress(rate=1)
	 */
	double rate_prop_shearstress_rate1 = doubledot(rate_prop_stress, imposed_flow->getSymGradU()); // computed with rate=1, o here it is also the viscosity.
	double rate_indep_shearstress = doubledot(rate_indep_stress, imposed_flow->getSymGradU());
	if (is_brownian()) {
		rate_prop_shearstress_rate1_ave.update(rate_prop_shearstress_rate1, get_time());
		rate_indep_shearstress_ave.update(rate_indep_shearstress, get_time());
		rate_prop_shearstress_rate1 = rate_prop_shearstress_rate1_ave.get();
		rate_indep_shearstress = rate_indep_shearstress_ave.get();
	}
	if (rate_prop_shearstress_rate1 != 0) {
		double rate = (target_stress-rate_indep_shearstress)/rate_prop_shearstress_rate1;
		imposed_flow->setRate(rate);
	} else {
		cerr << "rate_prop_shearstress_rate1 = 0 ???" << endl;
		exit (1);
	}
	if (clk.cumulated_strain < init_strain_shear_rate_limit) {
		if (imposed_flow->getRate() > init_shear_rate_limit) {
			imposed_flow->setRate(init_shear_rate_limit);
		}
	}
}

void System::computeShearRateWalls()
{
	throw std::runtime_error(" computeShearRateWalls deactivated ");
// 	/**
// 	 \brief Compute the coefficient to give to the velocity of the fixed particles under stress control conditions.
// 	 */

// 	computeForcesOnWallParticles();

// 	double total_rate_dep_wall_shear_stress = 0;
// 	double total_rate_indep_wall_shear_stress = 0;
// 	cerr << "np_fixed =" << p.np_fixed << endl;
// 	for (int i=0; i<p.np_fixed; i++) {
// 		total_rate_dep_wall_shear_stress += dot(fixed_velocities[i], rate_proportional_wall_force[i]);
// 		total_rate_indep_wall_shear_stress += dot(fixed_velocities[i], non_rate_proportional_wall_force[i]);
// 	}
// 	double wall_surface;
// 	if (twodimension) {
// 		wall_surface = lx;
// 	} else {
// 		wall_surface = lx*ly;
// 	}

// 	total_rate_dep_wall_shear_stress /= wall_surface;
// 	total_rate_indep_wall_shear_stress /= wall_surface;

// 	// // the total_rate_dep_wall_shear_stress is computed above with shear_rate=1, so here it is also a viscosity.
// 	imposed_flow->setRate((-target_stress-total_rate_indep_wall_shear_stress)/total_rate_dep_wall_shear_stress);

// 	if (clk.cumulated_strain < init_strain_shear_rate_limit) {
// 		if (shear_rate > init_shear_rate_limit) {
// 			imposed_flow->setRate(init_shear_rate_limit);
// 		}
// 	}
// 	if (p.simulation_mode == 31) {
// 		force_upwall.reset();
// 		force_downwall.reset();
// 		for (int i=0; i<p.np_fixed; i++) {
// 			if (fixed_velocities[i].x > 0) {
// 				force_upwall += shear_rate*rate_proportional_wall_force[i]+non_rate_proportional_wall_force[i];
// 			}
// 			if (fixed_velocities[i].x < 0) {
// 				force_downwall += shear_rate*rate_proportional_wall_force[i]+non_rate_proportional_wall_force[i];
// 			}
// 		}
// 	}
// }

// void System::tmpMixedProblemSetVelocities()
// {
// 	if (p.simulation_mode == 1) {
// 		/* Shear reversal simulation
// 		 */
// 		static double time_next = 16;
// 		static double direction = 1;
// 		if (get_time() > time_next) {
// 			direction *= -1;
// 			time_next += 16;
// 			cerr << direction << endl;
// 		}
// 		for (int i=np_mobile; i<np; i++) { // temporary: particles perfectly advected
// 			na_velocity[i].reset();
// 			na_ang_velocity[i].reset();
// 		}
// 		na_velocity[np_mobile].x = direction;
// 	} else if (p.simulation_mode == 4) {
// 		static double time_next = 10;
// 		if (get_time() > time_next) {
// 			if (zero_shear == true) {
// 				zero_shear = false;
// 			} else {
// 				zero_shear = true;
// 			}
// 			time_next += 10;
// 		}
// 		if (zero_shear) {
// 			for (int i=np_mobile; i<np; i++) { // temporary: particles perfectly advected
// 				na_velocity[i].reset();
// 				na_ang_velocity[i].reset();
// 			}
// 		}
// 	} else if (p.simulation_mode >= 10 && p.simulation_mode < 20) {
// 		int i_np_in = np_mobile+np_wall1;
// 		// inner wheel
// 		for (int i=np_mobile; i<i_np_in; i++) { // temporary: particles perfectly advected
// 			na_velocity[i] = {-omega_wheel_in*(position[i].z-origin_of_rotation.z),
// 				0,
// 				omega_wheel_in*(position[i].x-origin_of_rotation.x)};
// 			na_ang_velocity[i] = {0, -omega_wheel_in, 0};
// 		}
// 		// outer wheel
// 		for (int i=i_np_in; i<np; i++) { // temporary: particles perfectly advected
// 			na_velocity[i] = {-omega_wheel_out*(position[i].z-origin_of_rotation.z),
// 				0,
// 				omega_wheel_out*(position[i].x-origin_of_rotation.x)};
// 			na_ang_velocity[i] = {0, -omega_wheel_out, 0};
// 		}
// 	} else if (p.simulation_mode == 21) {
// 		static double time_next = p.strain_reversal;
// 		if (get_time() > time_next) {
// 			p.theta_shear += M_PI;
// 			setShearDirection(p.theta_shear);
// 			time_next += p.strain_reversal;
// 		}
// 	} else if (p.simulation_mode == 31) {
// 		auto &vel_from_fixed = na_velo_components["from_fixed"];
// 		for (int i=np_mobile; i<np; i++) {
// 			vel_from_fixed.vel[i] = shear_rate*fixed_velocities[i-np_mobile];
// 			vel_from_fixed.ang_vel[i].reset();
// 			na_velocity[i] = vel_from_fixed.vel[i];
// 			na_ang_velocity[i] = vel_from_fixed.ang_vel[i];
// 		}
// 	} else if (p.simulation_mode == 41) {
// 		int i_np_wall1 = np_mobile+np_wall1;
// 		double wall_velocity = shear_rate*system_height;
// 		for (int i=np_mobile; i<i_np_wall1; i++) {
// 			na_velocity[i] = {-wall_velocity/2, 0, 0};
// 			na_ang_velocity[i].reset();
// 		}
// 		for (int i=i_np_wall1; i<np; i++) {
// 			na_velocity[i] = {wall_velocity/2, 0, 0};
// 			na_ang_velocity[i].reset();
// 		}
// 	} else if (p.simulation_mode == 42) {
// 		int i_np_wall1 = np_mobile+np_wall1;
// 		double wall_velocity = shear_rate*system_height;
// 		for (int i=np_mobile; i<i_np_wall1; i++) {
// 			na_velocity[i].reset();
// 			na_ang_velocity[i].reset();
// 		}
// 		for (int i=i_np_wall1; i<np; i++) {
// 			na_velocity[i] = {wall_velocity, 0, 0};
// 			na_ang_velocity[i].reset();
// 		}
// 	} else if (p.simulation_mode == 51) {
// 		int i_np_in = np_mobile+np_wall1;
// 		// inner wheel
// 		double l = lx/2;
// 		vec3d origin_of_rotation2(lx/2, 0, l);
// 		vec3d origin_of_rotation3(  lx, 0, 0);
// 		double x1 = l/sqrt(2);
// 		double x2 = x1+radius_in*sqrt(2);
// 		for (int i=i_np_in; i<np; i++) {
// 			if (position[i].x < x1) {
// 				na_velocity[i] = {-omega_wheel_out*(position[i].z),
// 					0,
// 					omega_wheel_out*(position[i].x)};
// 				na_ang_velocity[i] = {0, -omega_wheel_out, 0};
// 			} else if (position[i].x < x2) {
// 				na_velocity[i] = {-omega_wheel_in*(position[i].z-origin_of_rotation2.z),
// 					0,
// 					omega_wheel_in*(position[i].x-origin_of_rotation2.x)};
// 				na_ang_velocity[i] = {0, -omega_wheel_in, 0};
// 			} else {
// 				na_velocity[i] = {-omega_wheel_out*(position[i].z-origin_of_rotation3.z),
// 					0,
// 					omega_wheel_out*(position[i].x-origin_of_rotation3.x)};
// 				na_ang_velocity[i] = {0, -omega_wheel_out, 0};
// 			}
// 		}
// 	} else if (p.simulation_mode == 60) {
// 		if (mobile_fixed) {
// 			int i_np_wall1 = np_mobile+np_wall1;
// 			for (int i=np_mobile; i<i_np_wall1; i++) {
// 				na_velocity[i].reset();
// 				na_ang_velocity[i].reset();
// 			}
// 			for (int i=i_np_wall1; i<np; i++) {
// 				na_velocity[i].reset();
// 				na_ang_velocity[i].reset();
// 			}
// 		}
// 	}
}

void System::sumUpVelocityComponents()
{
	for (int i=0; i<np_mobile; i++) {
		na_velocity.vel[i].reset();
		na_velocity.ang_vel[i].reset();
	}
	for (const auto &vc: na_velo_components) {
		const auto &vel = vc.second.vel;
		const auto &ang_vel = vc.second.ang_vel;
		for (int i=0; i<np_mobile; i++) {
			na_velocity.vel[i] += vel[i];
			na_velocity.ang_vel[i] += ang_vel[i];
		}
	}
}

void System::setFixedParticleVelocities()
{
	throw std::runtime_error("Mobile-fixed to be fixed");
// 	if (p.simulation_mode == 0) {
// 		for (int i=np_mobile; i<np; i++) { // temporary: particles perfectly advected
// 			na_velocity[i].reset();
// 			na_ang_velocity[i].reset();
// 		}
// 	} else if (p.simulation_mode > 0) {
// 		tmpMixedProblemSetVelocities();
// 	}
}

void System::rescaleVelHydroStressControlled()
{
	for (auto &vc: na_velo_components) {
		if (vc.second.rate_dependence == RateDependence::proportional) {
			vc.second *= imposed_flow->getRate();
		}
	}
}

void System::computeNonAffineVelocities(bool divided_velocities, bool mat_rebuild)
{
	/**
	 \brief Compute velocities in the current configuration.

	 \param divided_velocities Divide the velocities in components
	 (hydro, contacts, Brownian, ...). (Note that in Brownian
	 simulations the Brownian component is always computed explicitely, independently of the values of divided_velocities.)
	 */
	res_solver->resetRHS();
	if (divided_velocities || control == Parameters::ControlVariable::stress) {
		if (control == Parameters::ControlVariable::stress) {
			imposed_flow->setRate(1);
		}
		computeUInf(); // note: after imposed_flow->setRate(1);
		if (mobile_fixed) {
			setFixedParticleVelocities();
		}
		computeVelocityByComponents();
		if (control == Parameters::ControlVariable::stress) {
			// Stress-controlled simulation
			if (p.simulation_mode != 31) {
				computeShearRate();
			} else {
				computeShearRateWalls();
			}
			rescaleVelHydroStressControlled();
		}
		sumUpVelocityComponents();
	} else {
		computeUInf();
		if (mobile_fixed) {
			setFixedParticleVelocities();
		}
		computeVelocityWithoutComponents(mat_rebuild);
	}
	if (in_predictor) {
		if (eventLookUp != NULL) {
			computeMaxNAVelocity();
		}
	}
	if (divided_velocities && wall_rheology) {
		if (in_predictor) {
			forceResultantLubricationForce();
		}
	}
	res_solver->solvingIsDone();
}

void System::computeNonAffineVelocitiesStokesDrag()
{
	/**
	 \brief Compute velocities in Stokes-drag simulation.

	 Note: Velocities of particles are simply proportional to the total forces acting on respective particles.
	 When the contact model includes dashpots, Stokes-drag simulation cannot be used.
	 */
	for (int i=0; i<np; i++) {
		na_velocity.vel[i].reset();
		na_velocity.ang_vel[i].reset();
	}
	for (const auto &component: interaction->declared_forces) {
		interaction->setForceToParticle(component, force_components[component].force, force_components[component].torque);
		for (unsigned int i=0; i<na_velocity.vel.size(); i++) {
			na_velocity.vel[i]     += force_components[component].force[i] /stokesdrag_coeff_f[i];
			na_velocity.ang_vel[i] += force_components[component].torque[i]/stokesdrag_coeff_t[i];
		}
	}
	if (dimer_manager) {
		for (const auto &component: dimer_manager->declared_forces) {
			dimer_manager->setForceToParticle(component, force_components[component].force, force_components[component].torque);
			for (unsigned int i=0; i<na_velocity.vel.size(); i++) {
				na_velocity.vel[i]     += force_components[component].force[i] /stokesdrag_coeff_f[i];
				na_velocity.ang_vel[i] += force_components[component].torque[i]/stokesdrag_coeff_t[i];
			}
		}
	}
	for (const auto &component: declared_forces) {
		if (component != "brownian" || in_predictor) {
			setForceToParticle(component, force_components[component].force, force_components[component].torque);
		}
		for (unsigned i=0; i<na_velocity.vel.size(); i++) {
			/* See the comment given in generateBrownianForces()
			 */
			na_velocity.vel[i]     += force_components[component].force[i] /stokesdrag_coeff_f_sqrt[i];
			na_velocity.ang_vel[i] += force_components[component].torque[i]/stokesdrag_coeff_t_sqrt[i];
		}
	}
	if (is_brownian() && twodimension) {
		rushWorkFor2DBrownian(na_velocity.vel, na_velocity.ang_vel);
	}
}

void System::computeUInf()
{
	for (int i=0; i<np; i++) {
		vel_bg.vel[i].reset();
		vel_bg.ang_vel[i].reset();
	}
	if (!imposed_flow->zero_shear()) {
		for (int i=0; i<np; i++) {
			vel_bg.vel[i] = dot(imposed_flow->getSymGradU(), conf->position[i]) + cross(imposed_flow->getAntiSymGradU(), conf->position[i]);
			vel_bg.ang_vel[i] = imposed_flow->getAntiSymGradU();
			velgrad_bg.E[i] = imposed_flow->getSymGradU();
		}
	}
}

void System::computeTotalVelocity()
{
	if (control == Parameters::ControlVariable::stress) { // in rate control it is already done in computeNonAffineVelocities()
		computeUInf();
	}
	for (int i=0; i<np; i++) {
		velocity.vel[i] = na_velocity.vel[i] + vel_bg.vel[i];
		velocity.ang_vel[i] = na_velocity.ang_vel[i] + vel_bg.ang_vel[i];
	}
}

void System::adjustVelocitySolventFlow()
{
	vector<double> st_tens(3);
	for (int i=0; i<np_mobile; i++) {
		sflow->localFlow(conf->position[i], vel_bg.vel[i], vel_bg.ang_vel[i], st_tens);
		velgrad_bg.E[i].set(st_tens[0], 0, st_tens[1], 0, 0, st_tens[2]);// xx xy xz yz yy zz
	}
	for (int i=0; i<np_mobile; i++) {
		velocity.vel[i] = na_velocity.vel[i] + vel_bg.vel[i];
		velocity.ang_vel[i] = na_velocity.ang_vel[i] + vel_bg.ang_vel[i];
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
	conf->position[i] += dr;
	/* Note:
	 * When the position of the particle is periodized,
	 * we need to modify the velocity, which was already evaluated.
	 * The position and velocity will be used to calculate the contact forces.
	 */
	if (shear_type == ShearType::simple_shear) {
		/**** simple shear flow ****/
		lees->periodize(conf->position[i]);
	} else if (shear_type == ShearType::extensional_flow) {
		/**** extensional flow ****/
		kr->periodize(i, conf->position[i]);
	} else {
		throw std::runtime_error(" What to do here?");
		// periodize(position[i]);
	}
	pairconf->updateAfterParticleMove(i);

	if (!imposed_flow->zero_shear()) {
		velocity.vel[i] -= vel_bg.vel[i];
		vel_bg.vel[i] = imposed_flow->getGradU()*conf->position[i];
		velocity.vel[i] +=  vel_bg.vel[i];
	}
}


double System::getSystemVolume() const
{
	string indent = "  System::\t";
	double system_height, system_volume;
	if (z_top == -1) {
		system_height = container.lz;
	} else {
		/* wall particles are at z = z_bot - a and z_top + a
		 */
		system_height = z_top-z_bot;
	}
	if (twodimension) {
		system_volume = container.lx*system_height;
	} else {
		system_volume = container.lx*container.ly*system_height;
	}
	return system_volume;
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
	kn_avg.update(p.contact.kn, clk.cumulated_strain);
	kt_avg.update(p.contact.kt, clk.cumulated_strain);

	static double previous_cumulated_strain = 0;
	double deltagamma = (clk.cumulated_strain-previous_cumulated_strain);
	double kn_target = kn_avg.get()*overlap_avg.get()/p.overlap_target;
	double dkn = (kn_target-p.contact.kn)*deltagamma/p.memory_strain_k;

	p.contact.kn += dkn;
	if (p.contact.kn < p.min_kn_auto_det) {
		p.contact.kn = p.min_kn_auto_det;
	}
	if (p.contact.kn > p.max_kn_auto_det) {
		p.contact.kn = p.max_kn_auto_det;
	}
	if (p.disp_tan_target != -1) {
		double kt_target = kt_avg.get()*max_disp_tan_avg.get()/p.disp_tan_target;
		double dkt = (kt_target-p.contact.kt)*deltagamma/p.memory_strain_k;
		p.contact.kt += dkt;
		if (p.contact.kt < p.min_kt_auto_det) {
			p.contact.kt = p.min_kt_auto_det;
		}
		if (p.contact.kt > p.max_kt_auto_det) {
			p.contact.kt = p.max_kt_auto_det;
		}
	} else {
		p.contact.kt = p.contact.kn;
	}
	adaptTimeStepWithVelocities();
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
	throw std::runtime_error("System:: resetContactModelParameer() broken.");
	// for (auto &inter: *interaction) {
	// 	inter->contact->setSpringConstants(p.contact);
	// 	inter->contact->setDashpotConstants();
	// }
}

// void System::yaplotBoxing(std::ofstream &fout_boxing)
// {
// 	boxset.yaplotBox(fout_boxing);
// 	vec3d e1;
// 	vec3d e2;
// 	vec3d dy(0, 0, 0);
// 	if (twodimension) {
// 		dy.y = 0.01;
// 	}
// 	e1 = lx*deform_forward.getLine(0);
// 	e2 = lz*deform_forward.getLine(2);
// 	fout_boxing << "@ 0" << endl;
// 	fout_boxing << "r 0.2\n";
// 	fout_boxing << "s 0 -0.01 0 " << e1 -dy<< endl;
// 	fout_boxing << "s 0 -0.01 0 " << e2 -dy << endl;
// 	fout_boxing << "s " << e1 -dy << ' ' << e1+e2 -dy << endl;
// 	fout_boxing << "s " << e2 -dy << ' ' << e1+e2 -dy << endl;
// 	fout_boxing << "@ 0" << endl;

// 	fout_boxing << "r 0.5\n";
// 	for (int i=0; i<np; i++) {
// 		//fout_boxing << "@ " << boxset.boxType(i) << endl;
// 		fout_boxing << "c " << position[i]-2*dy << endl;
// 	}
// 	if (1) {
// 		fout_boxing << "@ 2" << endl;
// 		int i = 0;
// 		fout_boxing << "r 1\n";
// 		fout_boxing << "c " << position[i] << endl;
// 		fout_boxing << "@ 4" << endl;
// 		for (auto j : boxset.neighborhood(i)) {
// 			if (i != j) {
// 				fout_boxing << "c " << position[j] << endl;
// 			}
// 		}
// 	}
// 	fout_boxing << "y 5" << endl;
// 	fout_boxing << "@ 5" << endl;
// 	fout_boxing << "r 1\n";
// 	for (auto op: overlap_particles) {
// 		fout_boxing << "c " << position[op]-3*dy << endl;
// 	}
// 	overlap_particles.clear(); // @@@ for debuging

// 	if (/* DISABLES CODE */ (0)) {
// 		fout_boxing << "@ 8" << endl;
// 		for (unsigned int k=0; k<interaction.size(); k++) {
// 			unsigned int p0 = interaction[k].get_p0();
// 			unsigned int p1 = interaction[k].get_p1();
// 			vec3d nvec = interaction[k].nvec;
// 			fout_boxing << "l " << position[p0]-3*dy << ' ' << position[p0]-3*dy + nvec  << endl;
// 			fout_boxing << "l " << position[p1]-3*dy << ' ' << position[p1]-3*dy - nvec  << endl;
// 		}
// 	}

// 	fout_boxing << endl;
// }

// void System::recordHistory()
// {
// 	for (unsigned int k=0; k<interaction.size(); k++) {
// //		if (interaction[k].lubrication.is_active()) {
// //			interaction[k].lubrication.calcLubricationForce();
// //		} else {
// //			interaction[k].lubrication.force = 0;
// //		}
// 		interaction[k].recordHistory();
// 	}
// }

void System::countContactNumber()
{
	/*
	 * This counts effective contact number per particle.
	 * It considers only contacts to effective particles in the force balance.
	 * Coordination number after excluding such particles is relevent for isostatic conditions.
	 */
	n_contact.clear();
	n_contact.resize(np, 0);
	for (auto &inter: *interaction) {
		if (inter->contact) {
			n_contact[inter->get_p0()] ++;
			n_contact[inter->get_p1()] ++;
		}
	}
	bool cn_change;
	do {
		cn_change = false;
		for (auto &inter: *interaction) {
			if (inter->contact) {
				if (n_contact[inter->get_p0()] == 1
					&& n_contact[inter->get_p1()] == 1) {
					/* If two particles are in contact with one bond,
					 * these particles are isolated.
					 * The contact will disappear.
					 * Residual overlaps can be considered as numerical artifact.
					 */
					n_contact[inter->get_p0()] = 0;
					n_contact[inter->get_p1()] = 0;
					cn_change = true;
				} else if (n_contact[inter->get_p0()] == 1
						   || n_contact[inter->get_p1()] == 1) {
					/*
					 * If one particle has only one contacting neighbor,
					 * the particle is a dead-end branch.
					 * Removing them.
					 * This iteration algortim removes them from the tip one by one.
					 */
					n_contact[inter->get_p0()] --;
					n_contact[inter->get_p1()] --;
					cn_change = true;
				}
			}
		}
	} while (cn_change);
	int num_node_contact_network = 0;
	int num_bond_contact_network = 0;
	for (int i = 0; i < np; i++) {
		if (n_contact[i] >= 2) {
			num_node_contact_network ++;
			num_bond_contact_network += n_contact[i];
		}
	}
	effective_coordination_number = num_bond_contact_network*(1.0/num_node_contact_network);
}

void System::initSolventFlow(string simulation_type)
{
	sflow->init(this, simulation_type);
}

vec3d System::meanParticleVelocity()
{
	vec3d mean_velocity(0);
	for (int i=0; i<np_mobile; i++) {
		mean_velocity += velocity.vel[i];
	}
	return mean_velocity/np_mobile;
}

vec3d System::meanParticleAngVelocity()
{
	vec3d mean_ang_velocity(0);
	for (int i=0; i<np_mobile; i++) {
		mean_ang_velocity += velocity.ang_vel[i];
	}
	return mean_ang_velocity/np_mobile;
	
}

