//
//  System.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class System
 \brief Central class holding the suspension's configuration and the methods to evolve the dynamics
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__System__
#define __LF_DEM__System__
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <fstream>
#include <queue>
#include <list>
#include <string>
#include <tuple>
#include <map>
#include <complex>
#include "global.h"
#include "Configuration.h"
#include "States.h"
#include "Sym2Tensor.h"
#include "vec3d.h"
#include "Matrix.h"
#include "ParameterSet.h"
#include "Averager.h"
#include "Events.h"
#include "StressComponent.h"
#include "SolventFlow.h"
#include "Box3d.h"
#include "KraynikReinelt.h"
#include "LeesEdwards.h"
#include "PairwiseConfig.h"
#include "ParticleConfig.h"
#include "ShearType.h"
#include "StdInteractionManager.h"
#include "PairwiseResistanceVelocitySolver.h"


class MTRand;

#ifdef USE_DSFMT
#include "dSFMT-src-2.2.3/dSFMT.h"
#endif

class Simulation;
// class Interaction;
class BoxSet;
class SolventFlow;

class System {
private:

	int np; ///< number of particles
	std::vector<unsigned int> nb_blocks_mm;
	std::vector<unsigned int> nb_blocks_mf;
	std::vector<unsigned int> nb_blocks_ff;
	//bool forbid_displacement;
	int total_num_timesteps;
	Geometry::box3d container;
	struct State::Clock clk;
	double angle_wheel; // rotational angle of rotary couette geometory
	double particle_volume;
	std::vector <vec3d> na_disp;
	std::vector< std::string > declared_forces;

	/* data */
	bool keepRunning(double time_end, double strain_end);
	bool keepRunning(const std::string& time_or_strain, const double& value_end);
	void timeEvolutionEulersMethod(bool calc_stress,
								   double time_end,
								   double strain_end);
	void timeEvolutionPredictorCorrectorMethod(bool calc_stress,
											   double time_end,
											   double strain_end);
	void timeStepMove(double time_end, double strain_end);
	void timeStepMoveCorrector();
	void timeStepMovePredictor(double time_end, double strain_end);
	void timeStepBoxing();
	void adaptTimeStepWithVelocities();
	void adaptTimeStepWithBounds(double time_end, double strain_end);
	void setBodyForce(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setConfinementForce(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setBrownianForceToParticle(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setForceToParticle(const std::string &component, std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void computeNonAffineVelocities(bool divided_velocities, bool mat_rebuild);
	void computeNonAffineVelocitiesStokesDrag();
	void computeVelocityWithoutComponents(bool rebuild);
	void computeVelocityByComponents();
	void computeVelocityByComponentsFixedParticles();
	void sumUpVelocityComponents();
	void setFixedParticleVelocities();
	void computeBrownianVelocities();
	void tmpMixedProblemSetVelocities();
	void computeTotalVelocity();
	void adjustVelocitySolventFlow();
	void rushWorkFor2DBrownian(std::vector<vec3d> &vel, std::vector<vec3d> &ang_vel); // We need to implement real 2D simulation.
	void computeUInf();
	void computeShearRate();
	void computeShearRateWalls();
	void computeForcesOnWallParticles();
	void computeVelocityCoeffFixedParticles();
	void rescaleVelHydroStressControlled();
	void sflowFiniteRe(bool calc_stress);
	void sflowIteration(bool calc_stress);
	void computeMaxNAVelocity();
	void computeMaxVelocity();
	void forceResultantReset();
	void forceResultantLubricationForce();
	void forceResultantInterpaticleForces();
	void checkForceBalance();
	void wallForces();
	void smoothStressTransition();
#ifndef USE_DSFMT
	MTRand *r_gen;
#endif
#ifdef USE_DSFMT
	dsfmt_t r_gen;
#endif
	std::vector<double> radius_cubed;
	std::vector<double> stokesdrag_coeff_f;
	std::vector<double> stokesdrag_coeff_t;
	std::vector<double> stokesdrag_coeff_f_sqrt;
	std::vector<double> stokesdrag_coeff_t_sqrt;
	
	void adjustContactModelParameters();
	Averager<double> kn_avg;
	Averager<double> kt_avg;
	Averager<double> overlap_avg;
	Averager<double> max_disp_tan_avg;

	void declareStressComponents();
	void declareVelocityComponents();
	void declareForceComponents();

	template<typename T> void setupGenericConfiguration(T conf, Parameters::ControlVariable control_);
	void setupBrownian();
	void setupParameters();
	void setupSystemPostConfiguration();
	void setConfiguration(const std::vector <vec3d>& initial_positions,
						  const std::vector <double>& radii,
						  const std::vector <double>& angles);
	bool hasPairwiseResistance();

 protected:
 public:
	System(struct State::BasicCheckpoint = State::zero_time_basicchkp);

	std::shared_ptr<Parameters::ParameterSet> p;
	int np_mobile; ///< number of mobile particles
	ShearType shear_type;
	std::shared_ptr<BC::LeesEdwardsBC> lees;
	std::shared_ptr<BC::KraynikReineltBC> kr;
	std::shared_ptr<Geometry::ImposedDeformation> imposed_flow;
	std::shared_ptr<Geometry::PairwiseConfig> pairconf;
	std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> res_solver;

	// Simulation parameters
	bool twodimension;
	Parameters::ControlVariable control;
	bool wall_rheology;
	bool mobile_fixed;
	bool couette_stress;
	double system_height;
	bool in_predictor;
	bool in_corrector;

	std::shared_ptr<ParticleConfig> conf;
	ParticleVelocity vel_bg; // u_inf or u_local -> solvent velocity at particle center if particle was absent
	ParticleVelocityGrad velgrad_bg; // solvent velocity gradient at particle center if particle was absent

	std::vector<double> phi_local;
	std::vector<vec3d> forceResultant;
	std::vector<vec3d> torqueResultant;
	std::vector<vec3d> non_rate_proportional_wall_force;
	std::vector<vec3d> non_rate_proportional_wall_torque;
	std::vector<vec3d> rate_proportional_wall_force;
	std::vector<vec3d> rate_proportional_wall_torque;

	SolventFlow *sflow;

	// std::vector<double> mu; // friction coeffient
	ParticleVelocity velocity;
	ParticleVelocity velocity_predictor;
	ParticleVelocity na_velocity;
	std::vector<Sym2Tensor> total_stress_pp; // per particle
	std::vector<std::complex<double>> phi6;
	Sym2Tensor total_stress;
	std::vector<int> n_contact;
	std::ofstream fout_history;

	std::shared_ptr<Interactions::StdInteractionManager> interaction;
	void gatherStressesByRateDependencies(Sym2Tensor &rate_prop_stress,
										  Sym2Tensor &rate_indep_stress);
	std::map<std::string, ForceComponent> force_components;
	std::map<std::string, Sym2Tensor> total_stress_groups;
	std::map<std::string, StressComponent> stress_components;
	std::map<std::string, ParticleVelocity> na_velo_components;
	Averager<Sym2Tensor> stress_avg;
	Averager<double> rate_prop_shearstress_rate1_ave;
	Averager<double> rate_indep_shearstress_ave;
	double dt;
	double avg_dt;
	int avg_dt_nb;
	double system_volume;
	vec3d shear_disp; // lees-edwards shift between top and bottom. only shear_disp.x, shear_disp.y is used
	double max_na_velocity;
	double max_velocity;
	double max_force_imbalance;
	double target_stress;
	double init_strain_shear_rate_limit;
	double init_shear_rate_limit;

	/**** temporal circular gap setup ***********/
	vec3d origin_of_rotation;
	double omega_wheel_in;
	double omega_wheel_out;
	int np_wall1;
	int np_wall2;
	double radius_wall_particle;
	double radius_in;  // wall particles are at r = radius_in - radius_wall_particle;
	double radius_out; // wall particles are at r = radius_out + radius_wall_particle;
	double z_bot;
	double z_top;
	double force_tang_wall1;
	double force_tang_wall2;
	double force_normal_wall1;
	double force_normal_wall2;
	double shearstress_wall1;
	double shearstress_wall2;
	double normalstress_wall1;
	double normalstress_wall2;
	vec3d force_upwall;
	vec3d force_downwall;
	double effective_coordination_number;
	double stress_transition_target;
	/**** pipe flow setup ***********/
	double pressure_difference;
	/****************************************************************************************************
	 * Extensional flow using Kraynik-Reinelt Method was originally implemented                         *
	 * by Antonio Martiniello and Giulio Giuseppe Giusteri from Auguest to November 2016 at OIST.       *
	 ****************************************************************************************************/

	std::vector <int> overlap_particles;

	std::list <Event> events;
	
	/****************************************/
	double getSystemVolume();
	void setContacts(const std::vector <struct contact_state>& cs);
	std::vector <struct contact_state> getContacts() const;
	struct base_configuration getBaseConfiguration() const;

	void setupConfiguration(struct base_shear_configuration c, Parameters::ControlVariable control_);
	void setupConfiguration(struct fixed_velo_configuration c, Parameters::ControlVariable control_);
	void setupConfiguration(struct circular_couette_configuration c, Parameters::ControlVariable control_);
	void setupConfiguration(const struct delayed_adhesion_configuration &conf, Parameters::ControlVariable control_);
	void resetContactModelParameer();
	void allocateRessources();
	void timeEvolution(double time_end, double strain_end);
	void displacement(int i, const vec3d& dr);
	void calculateForces(); //
	void calcStress();
	void calcStressPerParticle();
	void calcContactXFPerParticleStressControlled();
	void gatherVelocitiesByRateDependencies(ParticleVelocity &rateprop_vel,
												ParticleVelocity &rateindep_vel) const;
	void calcTotalStressPerParticle();
	void getStressCouette(int i,
						  double &stress_rr,
						  double &stress_thetatheta,
						  double &stress_rtheta);
	void initializeBoxing();

	/*************************************************************/
	void (System::*eventLookUp)();
	void eventShearJamming();
	// void yaplotBoxing(std::ofstream &fout_boxing); // Extensional flow Periodic Boundary condition
	// void recordHistory();
	// void openHistoryFile(std::string rec_filename);

	void countContactNumber();
	void checkStaticForceBalance();
	void initSolventFlow(std::string simulation_type);
	vec3d meanParticleVelocity();
	vec3d meanParticleAngVelocity();

	bool is_brownian() const
	{
		return p->brownian > 0;
	}
	
	bool has_body_force() const
	{
		return p->body_force > 0;
	}

	double get_lx() const
	{
		return container.lx;
	}

	double get_ly() const
	{
		return container.ly;
	}

	double get_lz() const
	{
		return container.lz;
	}

	double get_time() const
	{
		return clk.time_;
	}

	void set_np(int val)
	{
		np = val;
	}

	void set_np_mobile(int val)
	{
		np_mobile = val;
	}

	unsigned get_np() const
	{
		return np;
	}

	vec3d get_shear_strain()
	{
		if (shear_type == ShearType::simple_shear) {
			return lees->getShearStrain();
		}
		if (shear_type == ShearType::extensional_flow) {
			throw std::runtime_error("getShearStrain for Kraynik-Reinelt?");
		}
		return 0;
	}

	double get_cumulated_strain() const
	{
		return clk.cumulated_strain;
	}
	
	void reset_cumulated_strain()
	{
		clk.cumulated_strain = 0;
	}

	double get_angle_wheel()
	{
		return angle_wheel;
	}

	double get_omega_wheel()
	{
		return omega_wheel_in-omega_wheel_out;
	}

	std::size_t get_nb_interactions() const
	{
		return interaction->size();
	}

	int get_total_num_timesteps()
	{
		return total_num_timesteps;
	}

	const std::vector <vec3d> & getNonAffineDisp()
	{
		return na_disp;
	}
	
	void resetNonAffineDispData()
	{
		for (auto &elm: na_disp) {
			elm.reset();
		}
	}
};
#endif /* defined(__LF_DEM__System__) */
