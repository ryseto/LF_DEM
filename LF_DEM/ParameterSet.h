//
//  ParameterSet.h
//  LF_DEM
//
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__ParameterSet__
#define __LF_DEM__ParameterSet__
#include <string>

struct ParameterSet
{
	/*******************************************************
	 SIMULATION
	 *******************************************************/
	int simulation_mode;
	/*******************************************************
	 INTERACTIONS
	********************************************************/
	double repulsion;				///< Amplitude of the repulsive force [0]
	double critical_load;			///< Amplitude of the critical load [0]
	double cohesion;				///< Amplitude of the cohesion [0]
	double brownian;				///< Amplitude of the Brownian force [0]
	double repulsive_length;				///< "Debye" screering length for the repulsive force [0.05]
	double repulsive_max_length;            ///< Maximum length until which the repulsive force can reach. If -1, no limit. (e.g. length of polymer brush) [-1]
	double interaction_range;		///< maximum range (center-to-center) for interactions (repulsive force, etc.). If -1, lub_max_gap is used as cutoff [-1]
	int np_fixed;
	/*******************************************************
	 HYDRODYNAMICS
	********************************************************/
	/*
	* Stokes drag coeffient
	*/
	/* sd_coeff:
	* Full Stokes drag is given by sd_coeff = 1.
	* sd_coeff = 0 makes the resistance matrix singular.
	* sd_coeff = 1e-3 may be reasonable to remove the effect of the drag
	*/
	double sd_coeff;                         ///< Stokes drag coeffient. Full drag is 1. [1]
	/*
	* Lubrication model
	* 0 no lubrication (no dashpot as well)
	* 1 1/xi lubrication (only squeeze mode)
	* 2 log(1/xi) lubrication
	* 3 dashpot in the contact model (not implemented yet)
	*   In this case lub_max_gap should be automatically zero. (not yet implemented)
	*/
	std::string lubrication_model;                   ///< Lubrication type. "none": no lubrication, "normal": 1/xi lubrication (only squeeze mode), "tangential": "normal" plus log(1/xi) terms (shear and pump modes) ["tangential"]
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
	double lub_reduce_parameter;             ///< Lubrication regularization length ("roughness length") [1e-3]
	double lub_max_gap;                          ///< Lubrication range (in interparticle gap distance) [0.5]
	/*******************************************************
	 CONTACTS
	********************************************************/
	/*
	 * 0 No friction
	 * 1 Linear friction law Ft < mu Fn
	 * 2 Threshold friction without repulsive force
	 */
	double friction_model;                   ///< Friction model. 0: No friction. 1: Coulomb. 2: Coulomb with threshold. [1]
	double mu_static;                        ///< friction coefficient (static) [1]
	double mu_dynamic;                        ///< friction coefficient (dynamic). If -1, mu_dynamic = mu_static [-1]
	double mu_rolling;                        ///< friction coefficient (rolling) [0]
	double ft_max;							///< max tangential force in friction_model = 5 [1]
	/*
	 * Contact force parameters
	 * kn: normal spring constant
	 * kt: tangential spring constant
	 */
	double kn;                               ///< Particle stiffness: normal spring constant [2000input_units]
	double kt;                               ///< Particle stiffness: tangential spring constant [0.5kn]
	double kr;                               ///< Particle stiffness: rolling spring constant [0.5kn]
		/*
		 * contact_relaxation_factor:
		 *
		 * This gives the coeffient of the resistance term for h < 0.
		 * - If the value is negative, the value of 1/lub_reduce_parameter is used.
		 *
		 */
	double contact_relaxation_time;          ///< Relaxation time (normal) of the contact model. Sets the normal dashpot. If <0, use normal lubrication at contact as normal dashpot. [1e-3input_unit]
	double contact_relaxation_time_tan;      ///< Relaxation time (tangential) of the contact model. Sets the tangential dashpot. If <0, use tangential lubrication at contact as tangential dashpot. [-1input_unit]

	/*******************************************************
	 INTEGRATOR
	********************************************************/
	double time_end;                 ///< Time of the simulation. [10 time unit]

	int integration_method;                  ///< Integrator. 0: Euler's Method, 1: predictor-corrector. [1]
	bool fixed_dt;							///< Use constant dt [false]
	double dt;                           ///< When fixed_dt == false: initial time step value. When fixed_dt == true: time step value. [1e-4 time unit]
	double disp_max;                         ///< When fixed_dt == false only: maximum displacement at each time step, the time step size dt is determined from disp_max at every step. [2e-3 length unit]



	/*******************************************************
	 OUTPUT
	********************************************************/

	double time_interval_output_data;      ///< Output interval for outputing data_* file [0.01 time unit]
	double time_interval_output_config;    ///< Output interval for outputing int_* and par_* files [0.1 time unit]
	bool log_time_interval;                ///< Output in logarithmic time [false]
	double initial_log_time;               ///< Initial output time in log time mode [1e-4]
	int nb_output_data_log_time;           ///< Nb of data output in log time mode [100]
	int nb_output_config_log_time;           ///< Nb of config output in log time mode (must be <= nb_output_data_log_time) [100]
	bool origin_zero_flow;                   ///< Output: the middle height of the simulation box is set to the flow zero level. [true]

	bool out_data_particle;                  ///< Output par_* file [true]
	bool out_data_interaction;               ///< Output int_* file [true]
	bool out_binary_conf;					///< Output binary configurations conf_*.bin files [false]
	std::string out_particle_stress;				///< Output stress per particle in pst_* file, indicating which component ("c" for contact, "r" for repulsion, "b" for Brownian, "t" for total, "l" for lubrication) by a string, e.g "tc" for total stress and contact stress [""]
	bool out_data_vel_components;						///< Output velocity components in the par* file [false]

	/*******************************************************
	 CONTACT PARAMETERS AUTO-DETERMINATION
	********************************************************/
	bool auto_determine_knkt;                ///< auto-determine stiffnesses knowing overlap and tangential displacement targets, see System::adjustContactModelParameters for details [false]
	double overlap_target;                   ///< max overlap to reach when auto-determining stiffness  [0.05]
	double disp_tan_target;                  ///< max tangential displacement to reach when auto-determining stiffness [0.05]
	double memory_strain_k;                  ///< relaxation strain for the stiffness determination [0.02]
	double memory_strain_avg;                ///< averaging strain for the stiffness determination [0.01]
	double min_kn;                           ///< min normal spring constant when auto-determining stiffness [1000]
	double max_kn;                           ///< max normal spring constant when auto-determining stiffness [1000000]
	double min_kt;                           ///< min tangential spring constant when auto-determining stiffness [1000]
	double max_kt;                           ///< max tangential spring constant when auto-determining stiffness [1000000]
	double min_dt;                           ///< min time step when auto-determining stiffness [1e-7]
	double max_dt;                           ///< max time step when auto-determining stiffness [1e-3]
	double start_adjust;                     ///< strain after which aut-determination of stiffnesses starts [0.2]


	/*******************************************************
	OTHER PARAMETERS
	********************************************************/
	double Pe_switch;                        ///< Value of Peclet below which Brownian units are used [5]
	bool monolayer;							///< Particle movements are confined in monolayer. 3D rotations are allowed. [false]
	double rest_threshold; ///< criteria to judge saturation of deformation, i.e. jammed state etc. [1e-4]
	std::string event_handler;  ///< Select event handler [""]
 	bool cross_shear;        ///< Allow arbitrary shear in xy plane (actual direction set by theta_shear) [false]
	double theta_shear;  ///< Shear direction, in degress, 0 is shear along x, 90 is shear along y [0]
	double strain_reversal;  ///< for test_simulation = 21 (rtest1)
};

#endif/* defined(__LF_DEM__ParameterSet__) */
