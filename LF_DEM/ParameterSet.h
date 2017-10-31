//
//  ParameterSet.h
//  LF_DEM
//
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__ParameterSet__
#define __LF_DEM__ParameterSet__
#include <string>
#include <functional>
#include <vector>
#include "DimensionalQty.h"



struct ParameterSet
{
	/*******************************************************
	 SIMULATION
	 *******************************************************/
	int simulation_mode;
	std::string flow_type;
	double magic_angle; // magic angle for extensinoal flow
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
	int friction_model;                   ///< Friction model. 0: No friction. 1: Coulomb. 2: Coulomb with threshold. [1]
	double mu_static;                        ///< friction coefficient (static) [1]
	double mu_dynamic;                        ///< friction coefficient (dynamic). If -1, mu_dynamic = mu_static [-1]
	double mu_rolling;                        ///< friction coefficient (rolling) [0]
	double ft_max;                            ///< max tangential force in friction_model = 5 [1]
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
	double dt_max;                           ///< max time step for adaptive dt [-1]
	double dt_min;                           ///< min time step for adaptive dt [-1]



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
	bool out_bond_order_parameter6;        ///< Output amplitudes and arguments of 6-fold bond orientation order parameters in the par* file [false]

	/*******************************************************
	 CONTACT PARAMETERS AUTO-DETERMINATION
	********************************************************/
	bool auto_determine_knkt;                ///< auto-determine stiffnesses knowing overlap and tangential displacement targets, see System::adjustContactModelParameters for details [false]
	double overlap_target;                   ///< max overlap to reach when auto-determining stiffness  [0.05]
	double disp_tan_target;                  ///< max tangential displacement to reach when auto-determining stiffness [0.05]
	double memory_strain_k;                  ///< relaxation strain for the stiffness determination [0.02]
	double memory_strain_avg;                ///< averaging strain for the stiffness determination [0.01]
	double min_kn_auto_det;                           ///< min normal spring constant when auto-determining stiffness [1000]
	double max_kn_auto_det;                           ///< max normal spring constant when auto-determining stiffness [1000000]
	double min_kt_auto_det;                           ///< min tangential spring constant when auto-determining stiffness [1000]
	double max_kt_auto_det;                           ///< max tangential spring constant when auto-determining stiffness [1000000]
	double min_dt_auto_det;                           ///< min time step when auto-determining stiffness [1e-7]
	double max_dt_auto_det;                           ///< max time step when auto-determining stiffness [1e-3]
	double start_adjust;                     ///< strain after which aut-determination of stiffnesses starts [0.2]


	/*******************************************************
	OTHER PARAMETERS
	********************************************************/
	double Pe_switch;                        ///< Value of Peclet below which Brownian units are used [5]
	bool monolayer;							///< Particle movements are confined in monolayer. 3D rotations are allowed. [false]
	double rest_threshold; ///< criteria to judge saturation of deformation, i.e. jammed state etc. [1e-4]
	std::string event_handler;  ///< Select event handler [""]
	double theta_shear;  ///< Shear direction, in degress, 0 is shear along x, 90 is shear along y [0]
	double strain_reversal;  ///< for test_simulation = 21 (rtest1)
	bool keep_input_strain;  ///< Use as initial strain value the strain from initial Lees-Edwards displacement [false]
	double brownian_relaxation_time; ///< Averaging time scale in the stress controlled simulation for Brownian [1]
};


inline void setFromKeyValue(ParameterSet &p, const std::string &key, const Dimensional::DimensionalQty<double> &dim_value)
{
	double value = dim_value.value;
	if (key == "contact_relaxation_time") {
		p.contact_relaxation_time = value;
	} else if (key == "contact_relaxation_time_tan") {
		p.contact_relaxation_time_tan = value;
	} else if (key == "time_end") {
		p.time_end = value;
	} else if (key == "time_interval_output_config") {
		p.time_interval_output_config = value;
	} else if (key == "time_interval_output_data") {
		p.time_interval_output_data = value;
	} else if (key == "initial_log_time") {
		p.initial_log_time = value;
	} else if (key == "min_kn_auto_det") {
		p.min_kn_auto_det = value;
	} else if (key == "max_kn_auto_det") {
		p.max_kn_auto_det = value;
	} else if (key == "min_kt_auto_det") {
		p.min_kt_auto_det = value;
	} else if (key == "max_kt_auto_det") {
		p.max_kt_auto_det = value;
	} else {
		throw std::runtime_error("Unknown dimensional parameter "+key);
	}
}

template <typename T>
inline void setFromMap(ParameterSet &p,
                       const std::map <std::string, T> &param_map)
{
	for (const auto & param: param_map) {
		setFromKeyValue(p, param.first, param.second);
	}
}



template<typename T>
struct InputParameter
{
	std::string name_str;
	std::function<void(ParameterSet &, InputParameter<T>)> exportToParameterSet;
	T value;
};

#define PARAM_INIT(name, default_value) {#name,  [](ParameterSet &p, InputParameter<decltype(ParameterSet::name)> in) {p.name = in.value;}, default_value}
#define PARAM_INIT_DIMVAL(name, type, default_value) {#name,  [](ParameterSet &p, InputParameter<Dimensional::DimensionalValue<decltype(ParameterSet::name)>> in) {p.name = in.value.value;}, default_value}

static std::vector<InputParameter<bool>> InputBoolParameterSet = \
{
	PARAM_INIT(fixed_dt, false),
	PARAM_INIT(keep_input_strain, false),
	PARAM_INIT(monolayer, false),
	PARAM_INIT(auto_determine_knkt, false),
	PARAM_INIT(out_bond_order_parameter6, false),
	PARAM_INIT(out_data_vel_components, false),
	PARAM_INIT(out_binary_conf, false),
	PARAM_INIT(out_data_interaction, true),
	PARAM_INIT(out_data_particle, true),
	PARAM_INIT(origin_zero_flow, true),
	PARAM_INIT(log_time_interval, false),
};

static std::vector<InputParameter<double>> InputDoubleParameterSet = \
{
	PARAM_INIT(brownian_relaxation_time, 1),
	PARAM_INIT(strain_reversal, 1),
	PARAM_INIT(theta_shear, 0),
	PARAM_INIT(rest_threshold, 1e-4),
	PARAM_INIT(Pe_switch, 5),
	PARAM_INIT(start_adjust, 0.2),
	PARAM_INIT(max_dt_auto_det, 1e-3),
	PARAM_INIT(min_dt_auto_det, 1e-7),
	PARAM_INIT(memory_strain_avg, 0.01),
	PARAM_INIT(memory_strain_k, 0.02),
	PARAM_INIT(disp_tan_target, 0.05),
	PARAM_INIT(overlap_target, 0.05),
	PARAM_INIT(dt_min, -1),
	PARAM_INIT(dt_max, -1),
	PARAM_INIT(dt, 1e-4),
	PARAM_INIT(disp_max, 2e-3),
	PARAM_INIT(lub_max_gap, 0.5),
	PARAM_INIT(lub_reduce_parameter, 1e-3),
	PARAM_INIT(sd_coeff, 1),
	PARAM_INIT(interaction_range, -1),
	PARAM_INIT(repulsive_length, 0.05),
	PARAM_INIT(repulsive_max_length, -1),
	PARAM_INIT(magic_angle, 0),
	PARAM_INIT(mu_static, 1),
	PARAM_INIT(mu_dynamic, -1),
	PARAM_INIT(mu_rolling, 0),
};


static std::vector<InputParameter<int>> InputIntParameterSet = \
{
	PARAM_INIT(nb_output_data_log_time, 100),
	PARAM_INIT(nb_output_config_log_time, 100),
	PARAM_INIT(integration_method, 1),
	PARAM_INIT(friction_model, 1),
	PARAM_INIT(np_fixed, 0),
	PARAM_INIT(simulation_mode, 0)
};

static std::vector<InputParameter<str::string>> InputStrParameterSet = \
{
	PARAM_INIT(flow_type, ""),
	PARAM_INIT(event_handler, ""),
	PARAM_INIT(lubrication_model, "tangential"),
};

static std::vector<InputParameter<str::string>> InputDimValParameterSet = \
{
	PARAM_INIT(flow_type, ""),
	PARAM_INIT(lubrication_model, "tangential"),
};




void Simulation::setDefaultParameters(Dimensional::DimensionalQty<double> control_value)
{

	/**
	 \brief Set default values for ParameterSet parameters.
	 */
	auto input_scale = control_value.unit;
	auto input_scale_str = Dimensional::unit2suffix(input_scale);

	autoSetParameters("time_end", "10h");
	autoSetParameters("contact_relaxation_time", "1e-3"+input_scale_str);
	autoSetParameters("contact_relaxation_time_tan", "-1"+input_scale_str);
	if (input_scale != Dimensional::Unit::kn) {
		autoSetParameters("kn", "2000"+input_scale_str);
		autoSetParameters("min_kn_auto_det", "1000"+input_scale_str);
		autoSetParameters("max_kn_auto_det", "1000000"+input_scale_str);
	}
	if (input_scale != Dimensional::Unit::kt) {
		autoSetParameters("kt", "0.5kn");
		autoSetParameters("min_kt_auto_det", "1000"+input_scale_str);
		autoSetParameters("max_kt_auto_det", "1000000"+input_scale_str);
	}
	if (input_scale != Dimensional::Unit::kr) {
		autoSetParameters("kr", "0kn");
	}
	autoSetParameters("time_interval_output_data", "1e-2h");
	autoSetParameters("time_interval_output_config", "1e-1h");
}


#endif/* defined(__LF_DEM__ParameterSet__) */
