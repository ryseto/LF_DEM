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
#include "StressTensor.h"
#include "Interaction.h"
#include "vec3d.h"
#include "BoxSet.h"
#include "StokesSolver.h"
#include "ParameterSet.h"
#include "Averager.h"
#include "cholmod.h"
#ifndef USE_DSFMT
#include "MersenneTwister.h"
#endif
#ifdef USE_DSFMT
#include "dSFMT-src-2.2.3/dSFMT.h"
#endif

using namespace std;

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
	int maxnb_interactionpair;
	int maxnb_interactionpair_per_particle;
	int nb_of_active_interactions;
	int total_num_timesteps;
	double time;
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
	int dof;
	//double repulsiveforce_length; // repulsive force length
	//int integration_method; // 0: Euler's method 1: PredictorCorrectorMethod
	/* data */
	void (System::*timeEvolutionDt)(bool);
	void timeEvolutionEulersMethod(bool calc_stress);
	void timeEvolutionPredictorCorrectorMethod(bool calc_stress);
	void timeStepMove();
	void timeStepMoveCorrector();
	void timeStepMovePredictor();
	void timeStepBoxing(const double strain_increment);
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
	void computeVelocityComponents();
	void computeShearRate();
	void forceReset();
	void torqueReset();
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
	double *radius_cubed;
	double *radius_squared;
	void adjustContactModelParameters();
	Averager<double> *kn_avg;
	Averager<double> *kt_avg;
	Averager<double> *overlap_avg;
	Averager<double> *max_disp_tan_avg;
	bool fixed_dt;
	/*
	 * Simulation for magnetic particles
	 */
	double num_magnetic_particles;
	double sq_magnetic_interaction_range;
	vector<pair<vec3d, pair<int,int>>> magnetic_force_stored;
	vector<vector<int>> magnetic_pair;
	void updateMagneticPair();
	double time_update_magnetic_pair;
	
 protected:
 public:
	System();
	~System();
	void importParameterSet(ParameterSet &ps);

	ParameterSet p;
	// Interaction types
	bool brownian;
	bool friction;
	bool rolling_friction;
	bool repulsiveforce;
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
	double *resistance_matrix_dblock;
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
	vector <double> magnetic_moment_norm;
	vector <double> magnetic_susceptibility;
	double *brownian_force;
	StressTensor* lubstress; // G U + M E
	StressTensor* contactstressGU; // by particle
	StressTensor* repulsivestressGU; // by particle
	StressTensor* brownianstressGU; // by particle
	StressTensor* brownianstressGU_predictor; // by particle
	StressTensor total_stress;
	StressTensor total_hydro_stress;
	StressTensor total_contact_stressXF_normal;
	StressTensor total_contact_stressXF_tan;
	StressTensor total_contact_stressXF;
	StressTensor total_contact_stressGU;
	StressTensor total_repulsive_stressXF;
	StressTensor total_repulsive_stressGU;
	StressTensor total_repulsive_stress;
	StressTensor total_brownian_stressGU;
	Averager<StressTensor> *stress_avg;
	double dt;
	double kn;
	double kt;
	double kr;
	double kn_master;
	double kt_master;
	double kr_master;
	double ft_max;
	double mu_static; // static friction coefficient
	double mu_dynamic; // dynamic friction coefficient
	double mu_rolling; // rolling friction coeffient
	double lub_coeff_contact;
	double magnetic_coeffient; // (3*mu0)/(4*M_PI)
	double einstein_stress;
	double einstein_viscosity;
	// resistance coeffient for normal mode
	double log_lub_coeff_contact_tan_dashpot;
	double log_lub_coeff_contact_tan_lubrication;
	double log_lub_coeff_contact_tan_total;
	set <Interaction*> *interaction_list;
	set <int> *interaction_partners;
	int nb_interaction;
	double shear_disp;
	/* For non-Brownian suspension:
	 * dimensionless_number = 6*pi*mu*a^2*shear_rate/F_repulsive(0)
	 * For Brownian suspension, it should be Peclet number
	 */
	//	double dimensionless_number;
	/* Velocity difference between top and bottom
	 * in Lees-Edwards boundary condition
	 * vel_difference = shear_rate * lz
	 */
	double vel_difference;
	double max_velocity;
	double max_relative_velocity;
	double max_sliding_velocity;
	double max_ang_velocity;
	double min_reduced_gap;
	double max_contact_gap;
	double max_disp_tan;
	double max_disp_rolling;
	queue<int> deactivated_interaction;
	double max_contact_velo_tan;
	double max_contact_velo_normal;
	double ave_contact_velo_tan;
	double ave_contact_velo_normal;
	double ave_sliding_velocity;
	int contact_nb; // gap < 0
	int fric_contact_nb; // fn > f* in the critical load model
	double max_fc_normal;
	double max_fc_tan;
	string simu_name;
	double target_stress;
	double init_strain_shear_rate_limit;
	double init_shear_rate_limit;
	double new_contact_gap; // When gel structure is imported it needs to be larger than 0 at the begining.
	/*
	 * Simulation for magnetic particles
	 */
	bool magnetic_rotation_active;
	double magnetic_energy;
	double magnetic_field_square;
	double angle_external_magnetic_field;
	vec3d external_magnetic_field;

	/////////////////////////////////
	void setSystemVolume(double depth = 0);
	void setConfiguration(const vector <vec3d> &initial_positions,
						  const vector <double> &radii,
						  double lx_, double ly_, double lz_);
	void setMagneticConfiguration(const vector <vec3d> &magnetic_moment,
								  const vector <double> &magnetic_susceptibility);
	void setMagneticMomentExternalField();

	void setInteractions_GenerateInitConfig();
	void setupSystem(string control);
	void setupBrownian();
	void allocatePositionRadius();
	void allocateRessources();
	void timeEvolution(double time_output_data);
	void displacement(int i, const vec3d &dr);
	void checkNewInteraction();
	void checkInteractionEnd();
	void updateInteractions();
	void updateMagneticInteractions();
	void updateUnscaledContactmodel();
	double lubricationForceFactor(int i, int j);
	int periodize(vec3d &);
	void periodize_diff(vec3d &, int &);
	void periodize_diff_unsheared(vec3d &);
	void stressBrownianReset();
	void calcStress();
	void calcStressPerParticle();
	void analyzeState();
	StokesSolver stokes_solver;
	void lubricationStress(int i, int j);
	void initializeBoxing();
	void calcLubricationForce(); // for visualization of force chains
	void calcPotentialEnergy();
	void calcMagneticEnergy();
	/*************************************************************/
	double calcInteractionRangeDefault(const int&, const int&);
	double calcLubricationRange(const int& i, const int& j);

	
	void setBoxSize(double lx_, double ly_, double lz_)
	{
		lx = lx_;
		lx_half = 0.5*lx;
		ly = ly_;
		ly_half = 0.5*ly;
		lz = lz_;
		lz_half = 0.5*lz;
		cerr << "box: " << lx << ' ' << ly << ' ' << lz << endl;
	}

	double getContactNumber()
	{
		return (double)2*contact_nb/np;
	}

	double getFrictionalContactNumber()
	{
		return (double)2*contact_nb/np;
	}
	
	double get_lx()
	{
		return lx;
	}

	double get_ly()
	{
		return ly;
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
		return nb_of_active_interactions;
	}

	int get_total_num_timesteps()
	{
		return total_num_timesteps;
	}

	double get_total_energy()
	{
		return total_energy;
	}
	
	struct ForceAmplitudes amplitudes;
};
#endif /* defined(__LF_DEM__System__) */
