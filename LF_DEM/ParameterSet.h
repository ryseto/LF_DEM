//
//  ParameterSet.h
//  LF_DEM
//
//  Copyright (c) 2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__ParameterSet__
#define __LF_DEM__ParameterSet__
#include <string>
#include "ActivatedAdhesion_Params.h"
#include "ConfinementParams.h"
#include "LubricationParams.h"
#include "RepulsiveForceParams.h"
#include "ContactParams.h"
#include "AgeingContactParams.h"
#include "VanDerWaalsParams.h"
#include "DimerParams.h"
#include "DimensionalQty.h"
#include "MagneticInteractionParams.h"

/*===================================================
 =            Adding a parameter to LF_DEM           =
 =====================================================
 
 
 
 **** Non force-defining parameter *******
 
 Adding a parameter which is *not* a force scale is done in 2 steps:
 
 1. Add the parameter in the ParameterSet structure below,
 documentation usually provides the parameter role and its default value within brackets.
 If you plan to add several related parameters for e.g. a new interaction, bundling them
 in a dedicated struct will help organizing ParameterSet.
 
 2. Parameter assignement on LF_DEM start up is provided by the ParameterSetFactory class.
 In ParameterSetFactory.cpp, modify setDefaultValues() function to declare
 the new parameter in one of the std::vectors BoolParams, IntParams, etc,
 according to the parameter type. These vectors are storing the information necessary
 to match a string in a parameter file to a member of ParameterSet.
 If both string and member have the same name
 (i.e. if you want ParameterSet.AParam to be set in the parameter file with "AParam = value;"),
 you can use the macros PARAM_INIT* to save a few typos.
 
 
 **** Force-defining parameter *******
 
 If the parameter is a defining a new force scale, follow step 1 and 2 as above.
 Then:
 
 3. Add the new force scale name in the Dimensional::Unit enum (in DimensionalQty.h)
 
 
 =====  End of Adding a parameter  ======*/

namespace Parameters {
	
	struct outputParams
	{
		Dimensional::DimensionalQty<double> time_interval_output_data;		///< Output interval for outputing data_* file [0.01 time unit]
		Dimensional::DimensionalQty<double> time_interval_output_config;	///< Output interval for outputing int_* and par_* files [0.1 time unit]
		bool log_time_interval;				///< Output in logarithmic time [false]
		bool new_material_functions;		///< Output new material functions (eta, lambda0, lambda3, insted of standard viscometric functions (eta, N1, N2) [false]
		Dimensional::DimensionalQty<double> initial_log_time;		///< Initial output time in log time mode [1e-4]
		int nb_output_data_log_time;		///< Nb of data output in log time mode [100]
		int nb_output_config_log_time;		///< Nb of config output in log time mode (must be <= nb_output_data_log_time) [100]
		bool origin_zero_flow;				///< Output: the middle height of the simulation box is set to the flow zero level. [true]
		bool out_data_particle;				///< Output par_* file [true]
		bool out_data_interaction;			///< Output int_* file [true]
		bool out_binary_conf;				///< Output binary configurations conf_*.bin files [false]
		bool out_gsd;						///< Output gsd data file [true]
		bool out_gsd_na_velocity;    		///< Output gsd data file [true]
        int gsd_nb_mark;                    ///< Number of markers in each particle to visualize rotation [0]
        double gsd_size_mark;               ///< Radius of marks to visualize rotation [0.2]
		std::string out_particle_stress;	///< Output stress per particle in pst_* file, indicating which component ("c" for contact, "r" for repulsion, "b" for Brownian, "t" for total, "l" for lubrication) by a string, e.g "tc" for total stress and contact stress [""]
		bool out_data_vel_components;		///< Output velocity components in the par* file [false]
		bool out_na_vel;					///< Output non-affine velocity components in the par* file [false]
		bool out_na_disp;					///< Output non-affine displacements since last time step in the par* file [false]
		double recording_start;				///< Hisotry is recorded after this strain [1]
		bool effective_coordination_number; ///< Count and output effective coordination number [false]
	};
	
	struct ParameterSet
	{
		/*******************************************************
		 SIMULATION
		 *******************************************************/
		int simulation_mode;
		std::string flow_type;
		double magic_angle; // magic angle for extensinoal flow
		int sflow_boundary_conditions; // boundary conditions for flow fields 0: all PD 1: x-PD and z-Wall  [0]
		bool solvent_flow; ///< [false]
		/*******************************************************
		 INTERACTIONS
		 ********************************************************/
		double brownian;				///< Amplitude of the Brownian force [0]
		struct Interactions::RepulsiveForceParams repulsion;
		struct Interactions::vanDerWaalsForceParams vdw;
		double bodyforce;              ///< Amplitude of the body force [0]
		int np_fixed;
		// struct Interactions::TActAdhesion::Parameters TA_adhesion; ///< Time delayed adhesion params (.adhesion_range [1e-2], .adhesion_max_force [0h], .activation_time [0h])
		struct Interactions::ActAdhesion::Params activated_adhesion; ///< Time delayed adhesion params (.adhesion_range [1e-2], .adhesion_max_force [0h], .activation_time [0h])
		struct Interactions::Dimer::DimerParams dimer;  ///< see DimerParams.h
        struct Interactions::MagneticInteractionParams params_magnetic_int;
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
		double sd_coeff;			///< Stokes drag coeffient. Full drag is 1. [1]
		/*
		 * Lubrication model
		 * 0 no lubrication (no dashpot as well)
		 * 1 1/xi lubrication (only squeeze mode)
		 * 2 log(1/xi) lubrication
		 * 3 dashpot in the contact model (not implemented yet)
		 *   In this case lub_max_gap should be automatically zero. (not yet implemented)
		 */
		Interactions::Lub::LubParams lub;         ///< see LubricationParams.h
		
		/*******************************************************
		 CONTACTS
		 ********************************************************/
		/*
		 * 0 No friction
		 * 1 Coulomb's friction law Ft < mu Fn
		 * 2 Critical load model (Threshold friction without repulsive force)
		 * 3 Critical load model with mu = infinity
		 * 4 Coulomb's friction law with mu = infinity
		 * 5 Constant maximum tangenetial force (no normal force dependence)
		 * 6 Coulomb's friction law with a maximum tangential force.
		 */
		struct Interactions::ContactParams contact;

		/*******************************************************
		 INTEGRATOR
		 ********************************************************/
		Dimensional::DimensionalQty<double> time_end;		///< Time of the simulation. [10 time unit]
		int integration_method;			///< Integrator. 0: Euler's Method, 1: predictor-corrector. [1]
		bool fixed_dt;					///< Use constant dt [false]
		double dt;				///< When fixed_dt == false: initial time step value. When fixed_dt == true: time step value. [1e-4 time unit]
		double disp_max;		///< When fixed_dt == false only: maximum displacement at each time step, the time step size dt is determined from disp_max at every step. [2e-3 length unit]
		double dt_max;          ///< When fixed_dt == false: forbid too large time step (if -1, no limitation) [-1]
		double dt_jamming;      ///< fixed_dt == true  [0.001s]
		/*******************************************************
		 OUTPUT
		 ********************************************************/
		outputParams output;
		
		/*******************************************************
		 CONTACT PARAMETERS AUTO-DETERMINATION
		 ********************************************************/
		bool auto_determine_knkt;		///< auto-determine stiffnesses knowing overlap and tangential displacement targets, see System::adjustContactModelParameters for details [false]
		double overlap_target;			///< max overlap to reach when auto-determining stiffness  [0.05]
		double disp_tan_target;			///< max tangential displacement to reach when auto-determining stiffness [0.05]
		double memory_strain_k;			///< relaxation strain for the stiffness determination [0.02]
		double memory_strain_avg;		///< averaging strain for the stiffness determination [0.01]
		double min_kn_auto_det;			///< min normal spring constant when auto-determining stiffness [0]
		double max_kn_auto_det;			///< max normal spring constant when auto-determining stiffness [0]
		double min_kt_auto_det;			///< min tangential spring constant when auto-determining stiffness [0]
		double max_kt_auto_det;			///< max tangential spring constant when auto-determining stiffness [0]
		double min_dt_auto_det;			///< min time step when auto-determining stiffness [1e-7]
		double max_dt_auto_det;			///< max time step when auto-determining stiffness [1e-3]
		double start_adjust;			///< strain after which aut-determination of stiffnesses starts [0.2]

		/*******************************************************
		 OTHER PARAMETERS
		 ********************************************************/
		double Pe_switch;						///< Value of Peclet below which Brownian units are used [5]
		bool monolayer;							///< Particle movements are confined in monolayer. 3D rotations are allowed. [false]
		double rest_threshold; ///< criteria to judge saturation of deformation, i.e. jammed state etc. [1e-4]
		std::string event_handler;  ///< Select event handler [""]
		double sj_disp_max_goal;	///< Maximum displacment for shear jamming [1e-6]
		double sj_disp_max_shrink_factor;	///< rescaling factor for negative shear rate for shear jamming [1.1]
		double sj_shear_rate;		///< Shear rate to judge shear jamming [0]
		double sj_velocity;         ///< Velocity to judge shear jamming [1e-3]
		int sj_check_count;			///< Jamming is judeged after counting this number [500] 
		int sj_reversal_repetition;			///< Repetition number for shear reversal [2]
		std::string sj_program_file;        ///< [""]
		double theta_shear;  ///< Shear direction, in degrees, 0 is shear along x, 90 is shear along y [0]
		double strain_reversal;  ///< for test_simulation = 21 (rtest1)
		bool keep_input_strain;  ///< Use as initial strain value the strain from initial Lees-Edwards displacement [false]
		double brownian_relaxation_time; ///< Averaging time scale in the stress controlled simulation for Brownian [1]
		bool check_static_force_balance;
		double body_force_angle;  ///< parallel to wall 0 and vertical to wall 90 [0]
		double sflow_dx; // mesh size in unit of particle radius [5]
		double sflow_smooth_length; // mesh size in unit of particle radius [3]
		int sflow_Darcy_power; // lambda for phi/(1-phi)^{lambda} (Modified Darcy's law) [0]
		double sflow_Darcy_coeff; // [1]
//		double sflow_ReNum; // [0.1]
		double sflow_ReNum_p; // [0.001]
		double sflow_pcontrol_increment; // To adjust pressure diffrence to fix flux value along x-direction [1e-4]
		double sflow_pcontrol_rtime; // [0.1]
		double sflow_pcontrol_damper; // [0.1]
		double sflow_target_flux; // Target flux along x-direction controlled by pressure difference [0]

		Confinement::Parameters confinement;
		struct Interactions::AgeingContactParams ageing_contact;
        
        int magnetic_field_type;                        // Type of exeternally applied magnetic field
        double langevin_parameter;                      // Langevin parameter (alpha) = \miu_0*m*H_0/k_B*T
        double magnetic_field_freq;                     // Dimensionless angular frequency of magnetic field = 2*\pi*f/Dr
	};
	
} // namespace Parameters

#endif/* defined(__LF_DEM__ParameterSet__) */
