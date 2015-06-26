//
//  ParameterSet.h
//  LF_DEM
//
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__ParameterSet__
#define __LF_DEM__ParameterSet__

struct ParameterSet
{
	
	/*******************************************************
	 INTERACTIONS 
	********************************************************/
	double repulsion_amplitude;				///< Amplitude of the repulsive force [0]
	double critical_load_amplitude;			///< Amplitude of the critical load [0]
	double cohesion_amplitude;				///< Amplitude of the cohesion [0]
	double brownian_amplitude;				///< Amplitude of the Brownian force [0]
	double magnetic_amplitude;				///< Amplitude of the magnetic force [0]
	double repulsive_length;				///< "Debye" length for the repulsive force [0.05]
	int magnetic_type;						///< Magnetic, 1: fixed magnetic dipole, 2: magnetic susceptible particles
	double interaction_range;		///< maximum range (center-to-center) for interactions (repulsive force, magnetic force, etc.). If -1, lub_max_gap is used as cutoff [-1]
	double magnetic_interaction_range; ///< [20]
	double init_angle_external_magnetic_field;  ///< Initial angle of external magnetic field
	double rot_step_external_magnetic_field;  ///<
	double step_interval_external_magnetic_field;  ///<

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
	* 0 no lubrication
	* 1 1/xi lubrication (only squeeze mode)
	* 2 log(1/xi) lubrication
	* 3 ???
	*/
	int lubrication_model;                   ///< Lubrication type. 0: no lubrication, 1: 1/xi lubrication (only squeeze mode), 2: log(1/xi) lubrication. [2]	
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
	 * 3 Threshold friction without repulsion + mu inf
	 */
	double friction_model;                   ///< Friction model. 0: No friction. 1: Coulomb. 2: Coulomb with threshold. 3 infinite mu Coulomb with threshold [1]
	double mu_static;                        ///< friction coefficient (static) [1]
	double mu_dynamic;                        ///< friction coefficient (dynamic). If -1, mu_dynamic = mu_static [-1]
	double mu_rolling;                        ///< friction coefficient (rolling) [0]
	double ft_max;							///< max tangential force in friction_model = 5 [1]	
	/*
	 * Contact force parameters
	 * kn: normal spring constant
	 * kt: tangential spring constant
	 */
	bool unscaled_contactmodel;              ///< Scale the particles' stiffness with force scale [true under stress control, false under rate control]
	double kn;                               ///< Particle stiffness: normal spring constant [2000 under stress control, 10000 under rate control]
	double kt;                               ///< Particle stiffness: tangential spring constant [1000 under stress control, 6000 under rate control]
	double kr;                               ///< Particle stiffness: rolling spring constant [1000 under stress control, 6000 under rate control]	
		/*
		 * contact_relaxation_factor:
		 *
		 * This gives the coeffient of the resistance term for h < 0.
		 * - If the value is negative, the value of 1/lub_reduce_parameter is used.
		 *
		 */
	double contact_relaxation_time;          ///< Relaxation time (normal) of the contact model (in units of p.dt) [1e-3]
	double contact_relaxation_time_tan;      ///< Relaxation time (tangential) of the contact model (in units of p.dt) [0]
		
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
	double timeinterval_update_magnetic_pair; ///< [0.02]

	double time_interval_output_data;      ///< Output interval for outputing data_* file [0.01 time unit]
	double time_interval_output_config;    ///< Output interval for outputing int_* and par_* files [0.1 time unit]
	bool origin_zero_flow;                   ///< Output: the middle height of the simulation box is set to the flow zero level. [true]

	bool out_data_particle;                  ///< Output par_* file [true]
	bool out_data_interaction;               ///< Output int_* file [true]
	
	
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
	double start_adjust;                     ///< strain after which aut-determination of stiffnesses starts [0.2]
	

	/*******************************************************
	OTHER PARAMETERS
	********************************************************/
	double Pe_switch;                        ///< Value of Peclet below which Brownian units are used [5]	
	bool monolayer;							///< Particle movements are confined in monolayer. 3D rotations are allowed. [false]
	double rest_threshold; ///< criteria to judge saturation of deformation, i.e. jammed state etc. [1e-4]
	
};

#endif/* defined(__LF_DEM__ParameterSet__) */
