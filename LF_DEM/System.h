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
#include <vector>
#include <fstream>
#include <queue>
#include <list>
#include <string>
#include <tuple>
#include <map>
#include "StressTensor.h"
#include "Interaction.h"
#include "vec3d.h"
#include "BoxSet.h"
#include "StokesSolver.h"
#include "ParameterSet.h"
#include "Averager.h"
#include "Events.h"
#include "StressComponent.h"
#include "VelocityComponent.h"
#include "ForceComponent.h"
#include "cholmod.h"
#ifndef USE_DSFMT
#include "MersenneTwister.h"
#endif
#ifdef USE_DSFMT
#include "dSFMT-src-2.2.3/dSFMT.h"
#endif

class Simulation;
// class Interaction;
class BoxSet;

struct ForceAmplitudes
{
	double repulsion;
	double temperature;
	double sqrt_temperature;
	double contact;
	double cohesion;
	double critical_normal_force;
	double ft_max;
};

class System{
private:
	int np; ///< number of particles
	int maxnb_interactionpair;
	int maxnb_interactionpair_per_particle;
	int nb_of_active_interactions_mm;
	int nb_of_active_interactions_mf;
	int nb_of_active_interactions_ff;
	int nb_of_pairwise_resistances_mm;
	int nb_of_pairwise_resistances_mf;
	int nb_of_pairwise_resistances_ff;
	int nb_of_contacts_mm;
	int nb_of_contacts_mf;
	int nb_of_contacts_ff;
	bool pairwise_resistance_changed;
	int total_num_timesteps;
	double time_; ///< time elapsed since beginning of the time evolution.
	double time_in_simulation_units; ///< time elapsed since beginning of the time evolution. \b note: this is measured in Simulation (output) units, not in internal System units.
	double lx;
	double ly;
	double lz;
	double lx_half; // =lx/2
	double ly_half; // =ly/2
	double lz_half; // =lz/2
	double shear_strain;
	double angle_wheel; // rotational angle of rotary couette geometory
	int linalg_size;
	double costheta_shear;
	double sintheta_shear;
	double shear_rate;

	vec3d omega_inf;
	std::vector <vec3d> u_inf;
	/* data */
	bool keepRunning(double time_end, double strain_end);
	bool keepRunning(const std::string& time_or_strain, const double& value_end);
	void (System::*timeEvolutionDt)(bool, double, double);
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
	void adaptTimeStep();
	void adaptTimeStep(double time_end, double strain_end);
	void setContactForceToParticle(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setRepulsiveForceToParticle(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setFixedParticleForceToParticle(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setDashpotForceToParticle(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setHydroForceToParticle_squeeze(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setHydroForceToParticle_squeeze_tangential(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void buildHydroTerms();
	void (System::*buildLubricationTerms)();
	void buildLubricationTerms_squeeze(); // lubrication_model = 1
	void buildLubricationTerms_squeeze_tangential(); // lubrication_model = 2
	void setBrownianForceToParticle(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setSolverRHS(const ForceComponent &fc);
	void addToSolverRHS(const ForceComponent &fc);
	void resetForceComponents();
	void computeVelocities(bool divided_velocities);
	void computeVelocitiesStokesDrag();
	void computeVelocityWithoutComponents();
	void computeVelocityByComponents();
	void computeVelocityByComponentsFixedParticles();
	void sumUpVelocityComponents();
	void setFixedParticleVelocities();
	void computeBrownianVelocities();
	void tmpMixedProblemSetVelocities();
	void adjustVelocitiesLeesEdwardsPeriodicBoundary();
	void rushWorkFor2DBrownian(std::vector<vec3d> &vel, std::vector<vec3d> &ang_vel); // We need to implement real 2D simulation.
	void computeUInf();
	void computeShearRate();
	void computeShearRateWalls();
	void computeForcesOnWallParticles();
	void computeVelocityCoeffFixedParticles();
	void rescaleVelHydroStressControlled();
	void addUpInteractionStressGU(std::vector<StressTensor> &stress_comp,
	                              const std::vector<vec3d> &non_affine_vel,
	                              const std::vector<vec3d> &non_affine_ang_vel);
	void addUpInteractionStressME(std::vector<StressTensor> &stress_comp);

	void computeMaxNAVelocity();
	double (System::*calcInteractionRange)(int, int);
	void forceResultantReset();
	void forceResultantLubricationForce();
	void forceResultantInterpaticleForces();
	void checkForceBalance();
	void wallForces();
	bool hasNeighbor(int i, int j);
	void setVelocityDifference();
#ifndef USE_DSFMT
	MTRand *r_gen;
#endif
#ifdef USE_DSFMT
	dsfmt_t r_gen;
#endif
	bool angle_output;
	std::vector<double> radius_cubed;
	std::vector<double> radius_squared;
	std::vector<double> stokesdrag_coeff_f;
	std::vector<double> stokesdrag_coeff_t;
	std::vector<double> stokesdrag_coeff_f_sqrt;
	std::vector<double> stokesdrag_coeff_t_sqrt;
	std::vector <struct DBlock> resistance_matrix_dblock;

	void adjustContactModelParameters();
	Averager<double> kn_avg;
	Averager<double> kt_avg;
	Averager<double> overlap_avg;
	Averager<double> max_disp_tan_avg;
	std::list <Event>& events;

protected:
 public:
	System(ParameterSet& ps, std::list <Event>& ev);
	~System(){};
	ParameterSet& p;
	int np_mobile; ///< number of mobile particles
	// Interaction types
	bool brownian;
	bool friction;
	bool rolling_friction;
	bool repulsiveforce;
	bool cohesion;
	bool critical_load;
	bool lowPeclet;
	bool lubrication;
	bool pairwise_resistance;
	// Simulation parameters
	bool twodimension;
	bool rate_controlled;
	bool stress_controlled;
	bool zero_shear;
	bool wall_rheology;
	bool mobile_fixed;
	bool couette_stress;
	double system_height;
	bool in_predictor;
	bool in_corrector;
	std::vector<vec3d> position;
	std::vector<vec3d> forceResultant;
	std::vector<vec3d> torqueResultant;
	std::vector<vec3d> non_rate_proportional_wall_force;
	std::vector<vec3d> non_rate_proportional_wall_torque;
	std::vector<vec3d> rate_proportional_wall_force;
	std::vector<vec3d> rate_proportional_wall_torque;

	std::vector<Interaction> interaction;
	BoxSet boxset;
	std::vector<double> radius;
	std::vector<double> angle; // for 2D visualization

	std::vector<vec3d> velocity;
	std::vector<vec3d> velocity_predictor;
	std::vector<vec3d> na_velocity;
	std::vector<vec3d> ang_velocity;
	std::vector<vec3d> ang_velocity_predictor;
	std::vector<vec3d> na_ang_velocity;
	std::vector<vec3d> fixed_velocities;
	std::vector<StressTensor> total_stress_pp; // per particle
	StressTensor total_stress;

	void gatherStressesByRateDependencies(StressTensor &rate_prop_stress,
	                                      StressTensor &rate_indep_stress);

	std::map<std::string, ForceComponent> force_components;
	std::map<std::string, StressTensor> total_stress_groups;
	std::map<std::string, StressComponent> stress_components;
	void declareStressComponents();
	std::map<std::string, VelocityComponent> velocity_components;
	void declareVelocityComponents();
	void declareForceComponents();


	Averager<StressTensor> stress_avg;
	double dt;
	double avg_dt;
	int avg_dt_nb;
	double system_volume;
	std::vector < std::set <Interaction*, compare_interaction> > interaction_list;
	std::vector < std::vector<int> > interaction_partners;
	int nb_interaction;
	vec3d shear_disp; // lees-edwards shift between top and bottom. only shear_disp.x, shear_disp.y is used
	double max_velocity;
	double max_sliding_velocity;
	std::queue<int> deactivated_interaction;
	double target_stress;
	double init_strain_shear_rate_limit;
	double init_shear_rate_limit;
	/* Velocity difference between top and bottom
	 * in Lees-Edwards boundary condition
	 * vel_difference = shear_rate * lz
	 */
	vec3d vel_difference;
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


	double *ratio_unit_time; // to convert System time in Simulation time

	/****************************************/
	void setSystemVolume();
	void setConfiguration(const std::vector <vec3d>& initial_positions,
						  const std::vector <double>& radii);
	void setFixedVelocities(const std::vector <vec3d>& vel);
	void setContacts(const std::vector <struct contact_state>& cs);
	void getContacts(std::vector <struct contact_state>& cs);
	void setInteractions_GenerateInitConfig();
	void setupSystemPreConfiguration(std::string control, bool is2d);
	void setupSystemPostConfiguration();
	void allocateRessourcesPreConfiguration();
	void allocateRessourcesPostConfiguration();
	void timeEvolution(double time_end, double strain_end);
	void displacement(int i, const vec3d& dr);
	void checkNewInteraction();
	void createNewInteraction(int i, int j, double scaled_interaction_range);
	void removeNeighbors(int i, int j);
	void updateNumberOfInteraction(int p0, int p1, int val);
	void updateNumberOfPairwiseResistances(int p0, int p1, int val);
	void updateNumberOfContacts(int p0, int p1, int val);
	void updateInteractions();

	int periodize(vec3d&);
	int periodize_diff(vec3d&);
	vec3d periodized(const vec3d&);
	void calcStress();
	void calcStressPerParticle();
	void calcContactXFPerParticleStressControlled();
	void gatherVelocitiesByRateDependencies(std::vector<vec3d> rateprop_vel,
                                          std::vector<vec3d> rateprop_ang_vel,
                                          std::vector<vec3d> rateindep_vel,
                                          std::vector<vec3d> rateindep_ang_vel) const;
	void calcTotalStressPerParticle();
	void getStressCouette(int i,
						  double &stress_rr,
						  double &stress_thetatheta,
						  double &stress_rtheta);
	StokesSolver stokes_solver;
	void initializeBoxing();
	/*************************************************************/
	double calcInteractionRangeDefault(int, int);
	double calcLubricationRange(int, int);
	void (System::*eventLookUp)();
	void eventShearJamming();

	void setBoxSize(double lx_, double ly_, double lz_)
	{
		lx = lx_;
		lx_half = 0.5*lx;
		ly = ly_;
		ly_half = 0.5*ly;
		lz = lz_;
		lz_half = 0.5*lz;
	}

	double get_lx()
	{
		return lx;
	}

	double get_ly()
	{
		return ly;
	}

	double get_lz()
	{
		return lz;
	}

	double Lx_half()
	{
		return lx_half;
	}

	double Ly_half()
	{
		return ly_half;
	}

	double Lz_half()
	{
		return lz_half;
	}

	double get_time_in_simulation_units()
	{
		return time_in_simulation_units;
	}

	double get_time()
	{
		return time_;
	}

	double get_shear_rate()
	{
		return shear_rate;
	}

	void set_shear_rate(double sr);

	vec3d get_vel_difference()
	{
		return vel_difference;
	}

	void set_np(int val)
	{
		np = val;
	}

	void set_np_mobile(int val)
	{
		np_mobile = val;
	}

	int get_np() const
	{
		return np;
	}

	double get_shear_strain()
	{
		return shear_strain;
	}

	double get_angle_wheel()
	{
		return angle_wheel;
	}

	double get_omega_wheel()
	{
		return omega_wheel_in-omega_wheel_out;
	}

	void set_dt(double val)
	{
		dt = val;
	}

	double get_nb_of_active_interactions()
	{
		return nb_of_active_interactions_mm + nb_of_active_interactions_mf + nb_of_active_interactions_ff;
	}

	int get_total_num_timesteps()
	{
		return total_num_timesteps;
	}

	std::tuple<double,double> getCosSinShearAngle()
	{
		return std::make_tuple(costheta_shear, sintheta_shear);
	}

	void setShearDirection(double theta_shear){
		costheta_shear = cos(theta_shear);
		sintheta_shear = sin(theta_shear);
		setVelocityDifference();
	}
	struct ForceAmplitudes amplitudes;
};
#endif /* defined(__LF_DEM__System__) */
