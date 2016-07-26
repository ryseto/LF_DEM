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
// class Interaction;
class BoxSet;

struct ForceAmplitudes
{
	double repulsion;
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
	double shear_strain;
	double angle_wheel; // rotational angle of rotary couette geometory
	double total_energy;
	int linalg_size;
	double costheta_shear;
	double sintheta_shear;
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
	void setContactForceToParticle();
	void setRepulsiveForceToParticle();
	void buildHydroTerms(bool, bool);
	void (System::*buildLubricationTerms)(bool, bool);
	void buildLubricationTerms_squeeze(bool mat, bool rhs); // lubrication_model = 1
	void buildLubricationTerms_squeeze_tangential(bool mat, bool rhs); // lubrication_model = 2
	void buildHydroTermsFromFixedParticles();
	std::vector<double> computeForceFromFixedParticles();
	void generateBrownianForces();
	void buildContactTerms(bool);
	void buildRepulsiveForceTerms(bool);
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
	void rushWorkFor2DBrownian(); // We need to implement real 2D simulation.
	void computeShearRate();
	void computeShearRateWalls();
	void computeForcesOnWallParticles();
	void computeVelocityCoeffFixedParticles();
	void rescaleVelHydroStressControlled();
	void rescaleVelHydroStressControlledFixed();
	void stressReset();
	void addLubricationStress(Interaction &);
	void computeMaxNAVelocity();
	double (System::*calcInteractionRange)(int, int);
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
	void forceResultantReset();
	void forceResultantLubricationForce();
	void forceResultantInterpaticleForces();
	void checkForceBalance();
	void wallForces();
	bool hasNeighbor(int i, int j);

#ifndef USE_DSFMT
	MTRand *r_gen;
#endif
#ifdef USE_DSFMT
	dsfmt_t r_gen;
	unsigned long hash(time_t, clock_t);
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
	~System();
	ParameterSet& p;
	int np_mobile; ///< number of mobile particles
	int test_simulation; //@@@ This test simulation may be temporal to debug the mix problem.
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

	std::vector<vec3d>velocity;
	std::vector<vec3d>velocity_predictor;
	std::vector<vec3d> na_velocity;
	std::vector<vec3d>ang_velocity;
	std::vector<vec3d>ang_velocity_predictor;
	std::vector<vec3d> na_ang_velocity;
	std::vector<vec3d> vel_repulsive;
	std::vector<vec3d> ang_vel_repulsive;
	std::vector<vec3d> vel_contact;
	std::vector<vec3d> ang_vel_contact;
	std::vector<vec3d> vel_hydro;
	std::vector<vec3d> ang_vel_hydro;
	std::vector<vec3d> vel_brownian;
	std::vector<vec3d> ang_vel_brownian;
	std::vector<vec3d> vel_hydro_from_fixed;
	std::vector<vec3d> ang_vel_hydro_from_fixed;
	std::vector<vec3d> fixed_velocities;
	std::vector<vec3d>contact_force;
	std::vector<vec3d>contact_torque;
	std::vector<vec3d>repulsive_force;
	std::vector<vec3d> brownian_force_torque;
	std::vector<StressTensor> lubstress; // G U + M E
	std::vector<StressTensor> contactstressGU; // per particle
	std::vector<StressTensor> contactstressXF; // per particle
	std::vector<StressTensor> repulsivestressGU; // per particle
	std::vector<StressTensor> repulsivestressXF; // per particle
	std::vector<StressTensor> brownianstressGU; // per particle
	std::vector<StressTensor> brownianstressGU_predictor; // per particle
	std::vector<StressTensor> total_stress_pp; // per particle
	std::vector<StressTensor> hydrofromfixedstressGU; // per particle
	StressTensor total_stress;
	StressTensor total_hydro_stress;
	StressTensor total_contact_stressXF;
	StressTensor total_contact_stressGU;
	StressTensor total_repulsive_stressXF;
	StressTensor total_repulsive_stressGU;
	StressTensor total_brownian_stressGU;
	StressTensor total_hydrofromfixed_stressGU;
	Averager<StressTensor> stress_avg;
	double dt;
	double avg_dt;
	int avg_dt_nb;
	/* double kn; */
	/* double kt; */
	/* double kr; */
	double kn_master;
	double kt_master;
	double kr_master;
	double lub_coeff_contact;
	double system_volume;
	// resistance coeffient for normal mode
	double log_lub_coeff_contact_tan_dashpot;
	double log_lub_coeff_contact_tan_lubrication;
	double log_lub_coeff_contact_tan_total;
	std::vector < std::set <Interaction*> > interaction_list;
	std::vector < std::vector<int> > interaction_partners;
	// std::unordered_set <int> *interaction_partners;
	int nb_interaction;
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
	void allocatePositionRadius();
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

	void updateUnscaledContactmodel();
	int periodize(vec3d&);
	void periodize_diff(vec3d&, int&);
	void periodize_diff(vec3d&);
	void calcStress();
	void calcStressPerParticle();
	void calcTotalStressPerParticle();
	void getStressCouette(int i,
						  double &stress_rr,
						  double &stress_thetatheta,
						  double &stress_rtheta);
	void analyzeState();
	double evaluateAvgContactGap();
	StokesSolver stokes_solver;
	void initializeBoxing();
    //void calcLubricationForce(); // for visualization of force chains
	void calcPotentialEnergy();
	/*************************************************************/
	double calcInteractionRangeDefault(int, int);
	double calcLubricationRange(int, int);
	void (System::*eventLookUp)();
	void eventShearJamming();
	// std::pair<vec3d,vec3d> checkForceOnWalls();

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

    inline void set_np_mobile(int val)
    {
        np_mobile = val;
    }

	inline int get_np()
	{
		return np;
	}

	inline double get_shear_strain()
	{
		return shear_strain;
	}

	inline double get_angle_wheel()
	{
		return angle_wheel;
	}

	inline double get_omega_wheel()
	{
		return omega_wheel_in-omega_wheel_out;
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
