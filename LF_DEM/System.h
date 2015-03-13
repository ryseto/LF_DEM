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
#include "MersenneTwister.h"

using namespace std;

class Simulation;
class Interaction;
class BoxSet;

struct ForceAmplitudes
{
	double repulsion;
	double sqrt_temperature;
	double contact;
};

class System{
private:
	int np;
	int maxnb_interactionpair;
	int maxnb_interactionpair_per_particle;
	int nb_of_active_interactions;
	double time;
	double shear_rate;
	double disp_max;
	double lx;
	double ly;
	double lz;
	double lx_half; // =lx/2
	double ly_half; // =ly/2
	double lz_half; // =lz/2
	double system_volume;
	double shear_strain;
	double lub_max_gap;
	int linalg_size;
	int linalg_size_per_particle;
	int dof;
	int max_lub_int;
	double repulsiveforce_length; // repulsive force length (dimensionless)
	int integration_method; // 0: Euler's method 1: PredictorCorrectorMethod
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
	void buildHydroTerms(bool, bool);
	void (System::*buildLubricationTerms)(bool, bool);
	void buildLubricationTerms_squeeze(bool mat, bool rhs); // lubrication_model = 1
	void buildLubricationTerms_squeeze_tangential(bool mat, bool rhs); // lubrication_model = 2
	void generateBrownianForces();
	void buildContactTerms(bool);
	void buildRepulsiveForceTerms(bool);
	double (System::*calcInteractionRange)(const int&, const int&);
	double calcLubricationRange(const int& i, const int& j);
	double evaluateMinGap();
	double evaluateMaxDispTan();
	double evaluateMaxDispRolling();
	double evaluateMaxFcNormal();
	double evaluateMaxFcTangential();
	void evaluateMaxContactVelocity();
	double evaluateMaxVelocity();
	double evaluateMaxAngVelocity();
	void countNumberOfContact();
	MTRand *r_gen;
	double *radius_cubed;
	double *radius_squared;
	ParameterSet p;
	void adjustContactModelParameters();
	Averager<double> *kn_avg;
	Averager<double> *kt_avg;
	Averager<double> *overlap_avg;
	Averager<double> *max_disp_tan_avg;
	bool lowPeclet;
	
 protected:
 public:
	System();
	~System();
	void importParameterSet(ParameterSet &ps);
	bool brownian;	
	bool in_predictor;
	bool in_corrector;
	bool twodimension;
	bool rate_controlled;
	bool stress_controlled;
	/* zero_shear:
	 * To be used for relaxation to generate initial configuration.
	 */
	bool zero_shear;
	int shear_direction; // 1: positive shear rate, -1: negative shear rate
	/* When Pe is small, the simulation parameters need to be adjusted differently.
	 * For high Pe, spring constants are scaled with shear rate being proportional to Pe.
	 * In dimensionless simulation, kn and kt appear the same.
	 * For low Pe, they need to be scaled with Pe.
	 *
	 */
	double shear_strain_end;
	double critical_normal_force;
	double scale_factor_SmallPe;
	double volume_fraction;
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
	vec3d *contact_force;
	vec3d *contact_torque;
	vec3d *repulsive_force;
	double *brownian_force;
	StressTensor* lubstress; // G U + M E
	StressTensor* contactstressGU; // by particle
	StressTensor* repulsivestressGU; // by particle
	StressTensor* brownianstressGU; // by particle
	StressTensor* brownianstressGU_predictor; // by particle
	StressTensor contactstressXF_normal;
	StressTensor contactstressXF_tan;
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
	double dt; // <=== It should be called d_strain.
	double kn;
	double kt;
	double kr;
	double kn_master;
	double kt_master;
	double kr_master;
	int friction_model;
	bool friction;
	bool rolling_friction;
	bool repulsiveforce;
	bool cohesion;
	double mu_static; // static friction coefficient
	double mu_dynamic; // dynamic friction coefficient
	double mu_rolling; // rolling friction coeffient
	double dimensionless_cohesive_force;
	double lub_coeff_contact;
	/* sd_coeff:
	 * Full Stokes drag is given by sd_coeff = 1.
	 * sd_coeff = 0 makes the resistance matrix singular.
	 * sd_coeff = 1e-3 may be reasonable way to remove effect of
	 *
	 */
	double sd_coeff;
	double einstein_stress;
	double cohesive_force;
	// resistance coeffient for normal mode
	double log_lub_coeff_contact_tan_dashpot;
	double log_lub_coeff_contact_tan_lubrication;
	double log_lub_coeff_contact_tan_total;
	set <Interaction*> *interaction_list;
	set <int> *interaction_partners;
	/*
	 * Lubrication model
	 * 0 no lubrication
	 * 1 1/xi lubrication (only squeeze mode)
	 * 2 log(1/xi) lubrication
	 * 3 1/xi lubrication for h>0 and tangential dashpot.
	 */
	int lubrication_model;
	int nb_interaction;
	/*
	 * Leading term of lubrication force is 1/reduced_gap,
	 * with reduced_gap the gap
	 * reduced_gap = 2r/(a0+a1) - 2.
	 * we set a cutoff for the lubrication interaction,
	 * such that the lub term is proportional to:
	 *
	 * 1/(reduced_gap+lub_reduce_parameter)
	 * when reduced_gap > 0.
	 */
	double lub_reduce_parameter;
	double shear_disp;
	/* For non-Brownian suspension:
	 * dimensionless_number = 6*pi*mu*a^2*shear_rate/F_repulsive(0)
	 * For Brownian suspension, it should be Peclet number
	 */
	double dimensionless_number;
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
	double average_fc_normal;
	double max_fc_normal;
	double max_fc_tan;
	string simu_name;
	double target_stress_input;
	double target_stress;
	double init_strain_shear_rate_limit;
	double init_shear_rate_limit;
	double new_contact_gap; // When gel structure is imported it needs to be larger than 0 at the begining.
	void setSystemVolume(double depth = 0);
	void setConfiguration(const vector <vec3d> &initial_positions,
						  const vector <double> &radii,
						  double lx_, double ly_, double lz_);
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
	void updateUnscaledContactmodel();
	double lubricationForceFactor(int i, int j);
	int periodize(vec3d &);
	void periodize_diff(vec3d &, int &);
	void computeVelocities(bool divided_velocities);
	void computeVelocitiesStressControlled();
	void forceReset();
	void torqueReset();
	void stressReset();
	void stressBrownianReset();
	void calcStress();
	void calcStressPerParticle();
	void analyzeState();
	StokesSolver stokes_solver;
	void lubricationStress(int i, int j);
	void initializeBoxing();
	void calcLubricationForce(); // for visualization of force chains
	void setupShearFlow(bool activate)
	{
		if (activate) {
			/* In dimensionless simulations for non Browninan simulation,
			 * shear rate is always 1.
			 */
			vel_difference = lz;
		} else {
			vel_difference = 0;
		}
	}
	/*************************************************************/
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
	double getParticleContactNumber()
	{
		return (double)2*contact_nb/np;
	}
	void set_integration_method(int val)
	{
		integration_method = val;
	}
	void set_lubrication_model(int val)
	{
		lubrication_model = val;
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
	inline void set_lub_max_gap(double val)
	{
		lub_max_gap = val;
		if (lub_max_gap >= 1) {
			cerr << "lub_max_gap must be smaller than 1\n";
			exit(1);
		}
	}
	inline double get_lub_max_gap()
	{
		return lub_max_gap;
	}
	inline void set_dt(double val)
	{
		dt = val;
	}
	void set_disp_max(double val)
	{
		disp_max = val;
	}
	void set_repulsiveforce_length(double val)
	{
		repulsiveforce_length = val;
	}
	void set_sd_coeff(double val)
	{
		sd_coeff = val;
	}
	inline double get_repulsiveforce_length()
	{
		return repulsiveforce_length;
	}
//	void set_mu_static(double val)
//	{
//		mu_static = val;
//	}
//	inline double get_mu_static()
//	{
//		return mu_static;
//	}
	//inline double get_lub_coeff_contact(){return lub_coeff_contact;}
	//	inline double get_log_lub_coeff_dynamicfriction(){
	/* In a sliding state, the resistance coeffient is the sum of
	 * lubrication and dashpot.
	 * In our standard model, we do not set the dashpot.
	 * Only tangential lubrications are considered.
	 */
	//		return log_lub_coeff_contact_tan_total;
	//	}
	inline double get_nb_of_active_interactions()
	{
		return nb_of_active_interactions;
	}
	struct ForceAmplitudes amplitudes;
};
#endif /* defined(__LF_DEM__System__) */
