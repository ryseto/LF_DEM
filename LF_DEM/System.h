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
#include "StressTensor.h"
#include "Interaction.h"
#include "vec3d.h"
#include "BoxSet.h"
#include "StokesSolver.h"
#include "ParameterSet.h"
#include "Averager.h"
#include "Events.h"
#include "cholmod.h"
#ifndef USE_DSFMT
#include "MersenneTwister.h"
#endif
#ifdef USE_DSFMT
#include "dSFMT-src-2.2.3/dSFMT.h"
#endif

class Simulation;
class Interaction;
class BoxSet;

struct ForceAmplitudes
{
	double repulsion;
	double sqrt_temperature;
	double contact;
	double cohesion;
	double magnetic;
	double critical_normal_force;
	double ft_max;
};

class System{
private:
	int np; ///< nb of particles
	int nmobile;
	int maxnb_interactionpair;
	int maxnb_interactionpair_per_particle;
	int nb_of_active_interactions_mm;
	int nb_of_active_interactions_mf;
	int nb_of_active_interactions_ff;

	int total_num_timesteps;
	double time; ///< time elapsed since beginning of the time evolution.
	double time_in_simulation_units; ///< time elapsed since beginning of the time evolution. \b note: this is measured in Simulation (output) units, not in internal System units.
	double shear_rate;
	double lx;
	double ly;
	double lz;
	double lx_half; // =lx/2
	double ly_half; // =ly/2
	double lz_half; // =lz/2
	double system_volume;
	double shear_strain;
	double total_energy;
	int linalg_size;
	double costheta_shear;
	double sintheta_shear;
	/* data */
	bool keepRunning(const std::string& time_or_strain, const double& value_end);
	void (System::*timeEvolutionDt)(bool, const std::string&, const double&);
	void timeEvolutionEulersMethod(bool calc_stress, const std::string& time_or_strain, const double& value_end);
	void timeEvolutionPredictorCorrectorMethod(bool calc_stress,
											   const std::string& time_or_strain,
											   const double& value_end);
	void timeStepMove(const std::string& time_or_strain,
					  const double& value_end);
	void timeStepMoveCorrector();
	void timeStepMovePredictor(const std::string& time_or_strain,
							   const double& value_end);
	void timeStepBoxing(const double strain_increment);
	void adaptTimeStep();
	void adaptTimeStep(const std::string& time_or_strain,
					   const double& value_end);
	void setContactForceToParticle();
	void setRepulsiveForceToParticle();
	void setMagneticForceToParticle();
	void buildHydroTerms(bool, bool);
	void (System::*buildLubricationTerms)(bool, bool);
	void buildLubricationTerms_squeeze(bool mat, bool rhs); // lubrication_model = 1
	void buildLubricationTerms_squeeze_tangential(bool mat, bool rhs); // lubrication_model = 2
	void generateBrownianForces();
	void buildContactTerms(bool);
	void buildRepulsiveForceTerms(bool);
	void buildMagneticForceTerms(bool);
	void computeVelocities(bool divided_velocities);
	void computeVelocitiesStokesDrag();
	void computeVelocityComponents();
	void computeBrownianVelocities();
	void adjustVelocitiesLeesEdwardsPeriodicBoundary();
	void rushWorkFor2DBrownian(); // We need to implement real 2D simulation.
	void computeShearRate();
	void stressReset();
	void computeMaxNAVelocity();
	double (System::*calcInteractionRange)(const int&, const int&);
	double evaluateMinGap();
	double evaluateMaxContactGap();
	double evaluateMaxDispTan();
	double evaluateMaxDispRolling();
	double evaluateMaxFcNormal();
	double evaluateMaxFcTangential();
	void evaluateMaxContactVelocity();
	double evaluateMaxVelocity();
	double evaluateMaxAngVelocity();
	void countNumberOfContact();
	vec3d randUniformSphere(double);
#ifndef USE_DSFMT
	MTRand *r_gen;
	MTRand rand_gen;
#endif
#ifdef USE_DSFMT
	dsfmt_t r_gen;
	dsfmt_t rand_gen;
	unsigned long hash(time_t, clock_t);
#endif
	bool angle_output;
	double *radius_cubed;
	double *radius_squared;
	double *stokesdrag_coeff_f;
	double *stokesdrag_coeff_t;
	double *stokesdrag_coeff_f_sqrt;
	double *stokesdrag_coeff_t_sqrt;
	void adjustContactModelParameters();
	Averager<double> *kn_avg;
	Averager<double> *kt_avg;
	Averager<double> *overlap_avg;
	Averager<double> *max_disp_tan_avg;
	std::list <Event>& events;
	/*
	 * Simulation for magnetic particles
	 */
	double sq_magnetic_interaction_range;
	std::vector<std::pair<vec3d, std::pair<int,int> > > magnetic_force_stored;
	std::vector<std::vector<int> > magnetic_pair;
	void updateMagneticPair();
	double time_update_magnetic_pair;

 protected:
 public:
	System(ParameterSet& ps, std::list <Event>& ev);
	~System();
	ParameterSet& p;
	int test_simulation; //@@@ This test simulation may be temporal to debug the mix problem.
	// Interaction types
	bool brownian;
	bool friction;
	bool rolling_friction;
	bool repulsiveforce;
	bool magnetic;
	bool cohesion;
	bool critical_load;
	bool lowPeclet;
	// Simulation parameters
	bool twodimension;
	bool rate_controlled;
	bool stress_controlled;
	bool zero_shear; ///< To be used for relaxation to generate initial configuration.
	double volume_fraction;
	bool in_predictor;
	bool in_corrector;
	vec3d *position;
	Interaction *interaction;
	BoxSet boxset;
	double *radius;
	double *angle; // for 2D visualization
	std::vector <struct DBlock> resistance_matrix_dblock;

	vec3d *velocity;
	vec3d *velocity_predictor;
	vec3d *na_velocity;
	vec3d *ang_velocity;
	vec3d *ang_velocity_predictor;
	vec3d *na_ang_velocity;
	vec3d *vel_repulsive;
	vec3d *ang_vel_repulsive;
	vec3d *vel_contact;
	vec3d *ang_vel_contact;
	vec3d *vel_hydro;
	vec3d *ang_vel_hydro;
	vec3d *vel_brownian;
	vec3d *ang_vel_brownian;
	vec3d *vel_magnetic;
	vec3d *ang_vel_magnetic;
	vec3d *contact_force;
	vec3d *contact_torque;
	vec3d *repulsive_force;
	vec3d *magnetic_moment;
	vec3d *magnetic_force;
	vec3d *magnetic_torque;
	std::vector<double> magnetic_susceptibility;
	std::vector<vec3d> brownian_force;
	StressTensor* lubstress; // G U + M E
	StressTensor* contactstressGU; // per particle
	StressTensor* contactstressXF; // per particle
	StressTensor* repulsivestressGU; // per particle
	StressTensor* repulsivestressXF; // per particle
	StressTensor* brownianstressGU; // per particle
	StressTensor* brownianstressGU_predictor; // per particle
	StressTensor* magneticstressGU; // per particle
	std::vector<StressTensor> magneticstressXF;
	StressTensor total_stress;
	StressTensor total_hydro_stress;
	StressTensor total_contact_stressXF;
	StressTensor total_contact_stressGU;
	StressTensor total_repulsive_stressXF;
	StressTensor total_repulsive_stressGU;
	StressTensor total_brownian_stressGU;
	StressTensor total_magnetic_stressXF;
	StressTensor total_magnetic_stressGU;
	Averager<StressTensor> *stress_avg;
	double dt;
	/* double kn; */
	/* double kt; */
	/* double kr; */
	double kn_master;
	double kt_master;
	double kr_master;
	double lub_coeff_contact;
	double einstein_stress;
	double einstein_viscosity;
	// resistance coeffient for normal mode
	double log_lub_coeff_contact_tan_dashpot;
	double log_lub_coeff_contact_tan_lubrication;
	double log_lub_coeff_contact_tan_total;
	std::set <Interaction*> *interaction_list;
	std::unordered_set <int> *interaction_partners;
	int nb_interaction;
	double *ratio_unit_time; // to convert System time in Simulation time
	vec3d shear_disp; // lees-edwards shift between top and bottom. only shear_disp.x, shear_disp.y is used
	/* For non-Brownian suspension:
	 * dimensionless_number = 6*pi*mu*a^2*shear_rate/F_repulsive(0)
	 * For Brownian suspension, it should be Peclet number
	 */
	//	double dimensionless_number;
	/* Velocity difference between top and bottom
	 * in Lees-Edwards boundary condition
	 * vel_difference = shear_rate * lz
	 */
	vec3d vel_difference;
	double max_velocity;
	double max_relative_velocity;
	double max_sliding_velocity;
	double max_ang_velocity;
	double min_reduced_gap;
	double max_contact_gap;
	double max_disp_tan;
	double max_disp_rolling;
	std::queue<int> deactivated_interaction;
	double max_contact_velo_tan;
	double max_contact_velo_normal;
	double ave_contact_velo_tan;
	double ave_contact_velo_normal;
	double ave_sliding_velocity;
	int contact_nb; // gap < 0
	int fric_contact_nb; // fn > f* in the critical load model
	double max_fc_normal;
	double max_fc_tan;
	std::string simu_name;
	double target_stress;
	double init_strain_shear_rate_limit;
	double init_shear_rate_limit;
	double new_contact_gap; // When gel structure is imported it needs to be larger than 0 at the begining.
	void setSystemVolume(double depth = 0);
	void setConfiguration(const std::vector <vec3d>& initial_positions,
						  const std::vector <double>& radii,
						  double lx_, double ly_, double lz_);
	void setContacts(const std::vector <struct contact_state>& cs);
	void getContacts(std::vector <struct contact_state>& cs);
	void setInteractions_GenerateInitConfig();
	void setupSystem(std::string control);
	void allocatePositionRadius();
	void allocateRessources();
	void timeEvolution(const std::string& time_or_strain,
					   const double& value_end);
//	void timeEvolution(double value_end); // @@@ DEPRECATED
	void displacement(int i, const vec3d& dr);
	void checkNewInteraction();
	void createNewInteraction(int i, int j, double scaled_interaction_range);
	void destroyInteraction(int k);
	void updateInteractions();
	void updateMagneticInteractions();
	void updateUnscaledContactmodel();
	int periodize(vec3d&);
	void periodize_diff(vec3d&, int&);
	void periodize_diff(vec3d&);
	void calcStress();
	void calcStressPerParticle();
	void analyzeState();
	double evaluateAvgContactGap();
	StokesSolver stokes_solver;
	void initializeBoxing();
	void calcLubricationForce(); // for visualization of force chains
	void calcPotentialEnergy();
	/*
	 * Simulation for magnetic particles
	 */
	bool magnetic_rotation_active;
	double magnetic_dd_energy; // Magnetic dipole-dipole energy per particle
	double angle_external_magnetic_field;
	vec3d external_magnetic_field;
	void setMagneticConfiguration(const std::vector <vec3d>& magnetic_moment,
								  const std::vector <double>& magnetic_susceptibility);
	void setInducedMagneticMoment();
	void setMagneticMomentZero();
	void calcMagneticEnergy();
	/*************************************************************/
	double calcInteractionRangeDefault(const int&, const int&);
	double calcLubricationRange(const int& i, const int& j);
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

	unsigned int getTotalNumberOfContacts()
	{
		return contact_nb;
	}

	double getContactNumber()
	{
		return (double)2*contact_nb/np;
	}

	double getFrictionalContactNumber()
	{
		return (double)2*fric_contact_nb/np;
	}

	double get_lx()
	{
		return lx;
	}

	double get_ly()
	{
		return ly;
	}

	inline double get_lz()
	{
		return lz;
	}

	inline double Lx_half()
	{
		return lx_half;
	}

	inline double Ly_half()
	{
		return ly_half;
	}

	inline double Lz_half()
	{
		return lz_half;
	}

	double get_time_in_simulation_units()
	{
		return time_in_simulation_units;
	}

	double get_time()
	{
		return time;
	}

	inline double get_shear_rate()
	{
		return shear_rate;
	}

	inline void set_shear_rate(double sr)
	{
		shear_rate = sr;
	}

	inline void set_np(int val)
	{
		np = val;
	}

	inline int get_np()
	{
		return np;
	}

	inline double get_shear_strain()
	{
		return shear_strain;
	}

	inline void set_dt(double val)
	{
		dt = val;
	}

	inline double get_nb_of_active_interactions()
	{
		return nb_of_active_interactions_mm + nb_of_active_interactions_mf + nb_of_active_interactions_ff;
	}

	int get_total_num_timesteps()
	{
		return total_num_timesteps;
	}

	double get_total_energy()
	{
		return total_energy;
	}

	std::tuple<double,double> getCosSinShearAngle()
	{
		return std::make_tuple(costheta_shear, sintheta_shear);
	}
	struct ForceAmplitudes amplitudes;
};
#endif /* defined(__LF_DEM__System__) */
