//
//  ParameterSet.h
//  LF_DEM
//
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__ParameterSet__
#define __LF_DEM__ParameterSet__

struct ParameterSet{
	double Pe_switch = 5; ///< Value of Peclet below which low Peclet mode is enabled
	double dt_max = 1e-4; ///< [Euler]: initial time/strain step value. [Pedictor/Corrector or Brownian]: time/strain step value
	double dt_lowPeclet = 1e-4; ///< [Brownian]: time step (not strain) in low Peclet mode
	double disp_max = 2e-3;///< [Euler]: maximum displacement at each time step. Time step size dt is determined from disp_max at every step.
	
	int integration_method = 1; ///< Integrator. 0: Euler's Method, 1: predictor-corrector

	/*
	 * Stokes drag coeffient
	 */
	double sd_coeff = 1;  ///< Stokes drag coeffient. Full drag is 1
	/*
	 * Lubrication model
	 * 0 no lubrication
	 * 1 1/xi lubrication (only squeeze mode)
	 * 2 log(1/xi) lubrication
	 * 3 ???
	 */
	int lubrication_model = 2; ///< Lubrication type. 0: no lubrication, 1: 1/xi lubrication (only squeeze mode), 2: log(1/xi) lubrication
	/*
	 * 0 No friction
	 * 1 Linear friction law Ft < mu Fn
	 * 2 Threshold friction without repulsive force
	 * 3 Threshold friction without repulsion + mu inf
	 */
	double friction_model = 1; ///< Friction model. 0: No friction. 1: Coulomb. 2: Coulomb with threshold. 3 infinite mu Coulomb with threshold

	bool rolling_friction = false;///< Activate rolling friction.
	double shear_strain_end = 10;///< Length of the simulation, in strain units.
	double lub_max = 2.5;///< Lubrication range (center-to-center distance)
	/*
	 * gap_nondim_min: gives reduced lubrication (maximum coeeffient).
	 *
	 */
	double lub_reduce_parameter = 1e-3;///< Lubrication regularization length ("roughness length")
	/*
	 * contact_relaxation_factor:
	 *
	 * This gives the coeffient of the resistance term for h < 0.
	 * - If the value is negative, the value of 1/lub_reduce_parameter is used.
	 *
	 */
	double contact_relaxation_time = 1e-3;///< Relaxation time (normal) of the contact model
	double contact_relaxation_time_tan = 0;///< Relaxation time (tangential) of the contact model
	/*
	 * Contact force parameters
	 * kn: normal spring constant
	 * kt: tangential spring constant
	 */
	bool unscaled_contactmodel;///< Scale the particles' stiffness with shear rate
	double kn;///< Particle stiffness: normal spring constant
	double kt;///< Particle stiffness: tangential spring constant
	double kr;///< Particle stiffness: rolling spring constant
	double kn_lowPeclet;///< Particle stiffness: normal spring constant in low Peclet mode
	double kt_lowPeclet;///< Particle stiffness: tangential spring constant in low Peclet mode
	double kr_lowPeclet;///< Particle stiffness: rolling spring constant in low Peclet mode

	double overlap_target = 0.05;///< max overlap to reach when auto-determining stiffness
	double disp_tan_target = 0.05;///< max tangential displacement to reach when auto-determining stiffness
	double max_kn = 1000000;///< max normal spring constant when auto-determining stiffness (auto-determination exits with failure return if kn>max_kn)

	double repulsive_length = 0.05;///< max normal spring constant when auto-determining stiffness (auto-determination exits with failure return if kn>max_kn)

	double mu_static = 1;///< friction coefficient (static)

	double strain_interval_output_data = 0.01;///< Output interval for outputing rheo_* file
	double strain_interval_output_config = 0.1;///< Output interval for outputing int_* and par_* files
	bool origin_zero_flow = true;///< Output: the middle height of the simulation box is set to the flow zero level.

	bool out_data_particle = true;///< Output par_* file
	bool out_data_interaction = true;///< Output int_* file
};

#endif/* defined(__LF_DEM__ParameterSet__) */
