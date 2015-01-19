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
	
	int _integration_method = 1; ///< Integrator. 0: Euler's Method, 1: predictor-corrector

	/*
	 * Stokes drag coeffient
	 */
	double _sd_coeff = 1;  ///< Stokes drag coeffient. Full drag is 1
	/*
	 * Lubrication model
	 * 0 no lubrication
	 * 1 1/xi lubrication (only squeeze mode)
	 * 2 log(1/xi) lubrication
	 * 3 ???
	 */
	int _lubrication_model = 2; ///< Lubrication type. 0: no lubrication, 1: 1/xi lubrication (only squeeze mode), 2: log(1/xi) lubrication
	/*
	 * 0 No friction
	 * 1 Linear friction law Ft < mu Fn
	 * 2 Threshold friction without repulsive force
	 * 3 Threshold friction without repulsion + mu inf
	 */
	double friction_model = 1; ///< Friction model. 0: No friction. 1: Coulomb. 2: Coulomb with threshold. 3 infinite mu Coulomb with threshold

	bool rolling_friction = false;
	/*
	 * Shear flow
	 *  shear_rate: shear rate
	 *  strain(): total strain (length of simulation)
	 *
	 */
	double shear_strain_end = 10;
	/*
	 * Lubrication force
	 * lub_max: reduced large cutoff distance for lubrication
	 * I think lub_max = 2.5 and 3 generate different results.
	 * We should give suffiently larger value.
	 * The value 3 or 3.5 should be better (To be checked.)
	 */
	double _lub_max = 2.5;
	/*
	 * gap_nondim_min: gives reduced lubrication (maximum coeeffient).
	 *
	 */
	double lub_reduce_parameter = 1e-3;
	/*
	 * contact_relaxation_factor:
	 *
	 * This gives the coeffient of the resistance term for h < 0.
	 * - If the value is negative, the value of 1/lub_reduce_parameter is used.
	 *
	 */
	double contact_relaxation_time = 1e-3;
	double contact_relaxation_time_tan = 0;
	/*
	 * Contact force parameters
	 * kn: normal spring constant
	 * kt: tangential spring constant
	 */
	bool unscaled_contactmodel;
	double kn;
	double kt;
	double kr;
	double kn_lowPeclet;
	double kt_lowPeclet;
	double kr_lowPeclet;

	double overlap_target = 0.05;
	double disp_tan_target = 0.05;
	double max_kn = 1000000;
	/*
	 * repulsive force parameter
	 * Short range repulsion is assumed.
	 * cf_amp_dl0: cf_amp_dl at shearrate = 1
	 */
	double repulsive_length;
	/*
	 * mu_static: static friction coeffient
	 * mu_dynamic: dynamic friction coeffient
	 */
	double _mu_static = 1;
	/*
	 * Output interval:
	 * strain_interval_output_data is for outputing rheo_...
	 * strain_interval_output is for outputing int_... and par_...
	 */
	double strain_interval_output_data = 0.01;
	double strain_interval_output_config = 0.1;
	/*
	 *  Data output
	 */
	/*
	 * The middle height of the simulation box is set to the flow zero level.
	 */
	bool origin_zero_flow = true;
	/*
	 * position and interaction data
	 */
	bool out_data_particle = true;
	bool out_data_interaction = true;
};

#endif/* defined(__LF_DEM__ParameterSet__) */
