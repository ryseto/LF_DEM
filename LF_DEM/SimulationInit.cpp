//
//  SimulationInit.cpp
//  LF_DEM
//
//  Created by Romain Mari on 08/10/15.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "Simulation.h"
#include <string>
#include <sstream>
#include <stdexcept>
#include "Configuration.h"
using namespace std;


void Simulation::importPreSimulationData(string filename)
{
	// @@@ DEPRECATED?
	ifstream fin_PreSimulationData;
	fin_PreSimulationData.open(filename.c_str());
	if (!fin_PreSimulationData) {
		ostringstream error_str;
		error_str  << " Pre-simulation data file '" << filename << "' not found." << endl;
		throw runtime_error(error_str.str());
	}
	double stress_, shear_rate_;
	while (fin_PreSimulationData >> stress_ >> shear_rate_) {
		if (stress_ == target_stress_input) {
			break;
		}
	}
}

void Simulation::echoInputFiles(string in_args,
								vector<string>& input_files)
{
	/**
	 \brief Print the entire information needed to reproduce the simulation in Simulation::fout_input
	 */
	fout_input << "# LF_DEM version " << GIT_VERSION << ", called with:" << endl;
	fout_input << in_args << endl << endl;
	for (const string& in_file : input_files) {
		ifstream in_f;
		string line;
		in_f.open(in_file.c_str());
		if (in_f.is_open()) {
			fout_input << "********** File " << in_file << " ************" << endl << endl;
			while (in_f.good()) {
				getline(in_f, line);
				fout_input << line << endl;
			}
			fout_input << endl << endl;
		}
		in_f.close();
	}
	fout_input.close();
}


void Simulation::exportForceAmplitudes()
{
	/**
	 \brief Copy the input force alues in the ForceAmplitude struct of the System class
	 */
	string indent = "  Simulation::\t";
	cout << indent+"Forces used:" << endl;
	indent += "\t";

	auto forces = system_of_units.getForceTree();
	sys.repulsiveforce = forces.count(Dimensional::Unit::repulsion) > 0;
	if (sys.repulsiveforce) {
		sys.p.repulsion = forces.at(Dimensional::Unit::repulsion).value;
		cout << indent+"Repulsive force (in \"" << Dimensional::unit2suffix(forces.at(Dimensional::Unit::repulsion).unit) << "\" units): " << sys.p.repulsion << endl;
	}

	sys.critical_load = forces.count(Dimensional::Unit::critical_load) > 0;
	if (sys.critical_load) {
		sys.p.critical_load = forces.at(Dimensional::Unit::critical_load).value;
		cout << indent+"Critical Load (in \"" << Dimensional::unit2suffix(forces.at(Dimensional::Unit::critical_load).unit) << "\" units): " << sys.p.critical_load << endl;
	}

	sys.cohesion = forces.count(Dimensional::Unit::cohesion) > 0;
	if (sys.cohesion) {
		sys.p.cohesion = forces.at(Dimensional::Unit::cohesion).value;
		cout << indent+"Cohesion (in \"" << Dimensional::unit2suffix(forces.at(Dimensional::Unit::cohesion).unit) << "\" units): " << sys.p.cohesion << endl;
	}

	bool is_ft_max = forces.count(Dimensional::Unit::ft_max) > 0;
	if (is_ft_max) {
		sys.p.ft_max = forces.at(Dimensional::Unit::ft_max).value;
		cout << indent+"Max Tangential Load (in \"" << Dimensional::unit2suffix(forces.at(Dimensional::Unit::ft_max).unit) << "\" units): " << sys.p.ft_max << endl;
	}

	sys.brownian = forces.count(Dimensional::Unit::brownian) > 0;
	if (sys.brownian) {
		sys.p.brownian = forces.at(Dimensional::Unit::brownian).value;
		cout << indent+"Brownian force (in \"" << Dimensional::unit2suffix(forces.at(Dimensional::Unit::brownian).unit) << "\" units): " << sys.p.brownian << endl;
	}
	sys.p.kn = forces.at(Dimensional::Unit::kn).value;
	cout << indent+"Normal contact stiffness (in \"" << Dimensional::unit2suffix(forces.at(Dimensional::Unit::kn).unit) << "\" units): " << sys.p.kn << endl;
	sys.p.kt = forces.at(Dimensional::Unit::kt).value;
	cout << indent+"Sliding contact stiffness (in \"" << Dimensional::unit2suffix(forces.at(Dimensional::Unit::kt).unit) << "\" units): " << sys.p.kt << endl;
	sys.p.kr = forces.at(Dimensional::Unit::kr).value;
	cout << indent+"Rolling contact stiffness (in \"" << Dimensional::unit2suffix(forces.at(Dimensional::Unit::kr).unit) << "\" units): " << sys.p.kr << endl;

	if (forces.count(Dimensional::Unit::hydro) > 0) { // == if rate controlled
		sys.set_shear_rate(forces.at(Dimensional::Unit::hydro).value);
	}
}

void Simulation::setupNonDimensionalization(Dimensional::DimensionalQty<double> control_value){
	/**
	 \brief Non-dimensionalize the simulation.

		This function determines the most appropriate unit scales to use in the System class depending on the input parameters (Brownian/non-Brownian, shear rate, stress/rate controlled), and converts all the input values in these units.
	 */
	string indent = "  Simulation::\t";
	if (control_var == Parameters::ControlVariable::rate) {
		input_rate = control_value.value; // @@@ Renaming is required?
	}
	if (control_var == Parameters::ControlVariable::rate) {// || control_var == Parameters::ControlVariable::viscnb) {
		if (input_rate != 0) {
			system_of_units.add(Dimensional::Unit::hydro, control_value);
			if (internal_unit == Dimensional::Unit::none){
				internal_unit = system_of_units.getLargestUnit();
			}
			if (internal_unit == Dimensional::Unit::brownian) {
				sys.brownian_dominated = true;
			}
		} else {
			if (control_value.unit == Dimensional::Unit::brownian) {
				cout << indent << "Brownian at Pe = 0 " << endl;
				internal_unit = Dimensional::Unit::brownian;
				sys.brownian_dominated = true;
				sys.zero_shear = true;
			} else if (control_value.unit == Dimensional::Unit::repulsion) {
				cout << indent << "non-Brownain at rate = 0 " << endl;
				internal_unit = Dimensional::Unit::repulsion;
				sys.zero_shear = true;
			}
		}
	} else if (control_var == Parameters::ControlVariable::stress) {
		system_of_units.add(Dimensional::Unit::stress, control_value);
		internal_unit = control_value.unit;
	}
	output_unit = control_value.unit;
	cout << indent << "internal units = " << Dimensional::unit2suffix(internal_unit) << endl;
	cout << indent << "output units = " << Dimensional::unit2suffix(output_unit) << endl;
	system_of_units.setInternalUnit(internal_unit);
	exportForceAmplitudes();
	for (auto &dimval: dimensional_input_params) {
		system_of_units.convertToInternalUnit(dimval.second);
	}
}


void Simulation::assertParameterCompatibility()
{
	// test for incompatibilities
	if (sys.brownian == true) {
		if (sys.pairwise_resistance && p.integration_method != 1) {
			 // @@@@ This test is broken as System has not yet set pairwise resistance. For now test is duplicated later on in System
			ostringstream error_str;
			error_str << "Brownian simulation needs to use the Predictor-Corrector method." << endl;
			error_str << "Modify the parameter file." << endl;
			throw runtime_error(error_str.str());
		}
	}
	if (control_var == Parameters::ControlVariable::stress) {
		if (p.integration_method != 0) {
			cerr << "Warning : use of the Predictor-Corrector method for the stress controlled simulation is experimental." << endl;
		}
		//p.integration_method = 0;
	}
	if (sys.critical_load) {
		p.friction_model = 2;
		cerr << "Warning : critical load simulation -> switched to friction_model=2" << endl;
	}
}

void Simulation::resolveTimeOrStrainParameters(const map<string, Dimensional::DimensionalQty<double>> &dim_params)
{
	/**
		\brief Interpret time units.

		We have to treat times as a special case, because we sometimes want to
		use times expressed in units of inverse shear rate (i.e. time=strain).
		Because the shear rate might change in the simulation, it is not possible
		to perform a simple change of unit at any time in the simulation.

		Ex 1: I am running a stress controlled simulation in repulsive time units.
		I set kn = 100r and time_end = 1kn;
		Then I know that I need to stop the simulation when sys.time()==1/100
		(i.e. when the time expressed in repulsive units reaches 1/100)

		Ex 2: I am running a stress controlled simulation in repulsive time units.
		I set kn = 100r but this time time_end = 1h;
		There is no way to know what 1h corresponds to in repulsive units,
		as the shear rate is evolving with time (stopping at sys.time()==???).
		The only way is to stop when sys.shear_strain()==1.

		Because these 2 cases lead to 2 different tests, we have to inform
		the System class which test to perform. If time_end > 0, System tests with
		System::time()==time_end. If time_end==-1, the test is done with
		System::shear_strain()==strain_end.

		We have to do this not only for time_end, but also for every time defined
		in the parameters.
	 */
	if (dim_params.at("time_end").dimension == Dimensional::Dimension::Strain) {
		time_end = -1;
		strain_end = dim_params.at("time_end").value;
	} else {
		time_end = dim_params.at("time_end").value;
	}
	if (p.log_time_interval) {
		if (dim_params.at("time_end").dimension != dim_params.at("initial_log_time").dimension &&
			(dim_params.at("time_end").dimension == Dimensional::Dimension::Strain || dim_params.at("initial_log_time").dimension == Dimensional::Dimension::Strain)) {
			throw runtime_error(" If one of time_end or initial_log_time is a strain (\"h\" unit), than both must be.\n");
		}
	}
}

void Simulation::setConfigToSystem(bool binary_conf, const std::string &filename)
{
	if (binary_conf) {
		auto format = getBinaryConfigurationFileFormat(filename);
		switch(format) {
			case ConfFileFormat::bin_format_base_new:
				{
					auto conf = readBinaryBaseConfiguration(filename);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case ConfFileFormat::bin_format_fixed_vel:
				{
					auto conf = readBinaryFixedVeloConfiguration(filename);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			default:
				throw std::runtime_error("Unable to read config binary format "+to_string(static_cast<int>(format)));
		}
	} else {
		auto format = getTxtConfigurationFileFormat(filename);
		switch(format) {
			case ConfFileFormat::txt_format_base_old:
				{
					auto conf = readTxtBaseConfiguration(filename);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case ConfFileFormat::txt_format_base_new:
				{
					auto conf = readTxtBaseConfiguration(filename);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case ConfFileFormat::txt_format_fixed_vel:
				{
					auto conf = readTxtFixedVeloConfiguration(filename);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case ConfFileFormat::txt_format_circular_couette:
				{
					auto conf = readTxtCircularCouetteConfiguration(filename);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			default:
				throw std::runtime_error("Unable to read config text format "+to_string(static_cast<int>(format)));
		}
	}
}


void Simulation::setupFlow(Dimensional::DimensionalQty<double> control_value)
{
	// @@@ This function is quite messy, should be fixed 
	// when we rewrite shear and extensional in a consistent manner


	/* dot_gamma = 1 --> dot_epsilon = 0;
	 *
	 */
	if (control_value.value != 0) {
		double dimensionless_deformation_rate = 0.5;
		if (!sys.ext_flow) {
			/* simple shear flow
			 * shear_rate = 2*dot_epsilon
			 */
			Einf_base.set(0, 0, dimensionless_deformation_rate, 0, 0, 0);
			Omegainf_base.set(0, dimensionless_deformation_rate, 0);
			sys.setImposedFlow(Einf_base, Omegainf_base);
			stress_basis_0 = {-dimensionless_deformation_rate/2, 0, 0, 0,
				dimensionless_deformation_rate, -dimensionless_deformation_rate/2};
			stress_basis_3 = {-dimensionless_deformation_rate, 0, 0, 0, 0, dimensionless_deformation_rate};
		} else {
			/* extensional flow
			 *
			 */
			p.magic_angle = atan(0.5*(sqrt(5)-1)); // simulation box needs to be tilted in this angle.
			matrix grad_u_orig(dimensionless_deformation_rate, 0, 0,
							   0, 0, 0,
							   0, 0, -dimensionless_deformation_rate);
			matrix rotation, rotation_inv;
			rotation.set_rotation(-p.magic_angle, 'y');
			rotation_inv.set_rotation(p.magic_angle, 'y');
			sys.grad_u_hat = rotation_inv*grad_u_orig*rotation;
			sys.grad_u = rotation_inv*grad_u_orig*rotation;
			Einf_base.setSymmetrize(sys.grad_u);
			Omegainf_base.set(0, 0, 0);
			sys.setImposedFlow(Einf_base, Omegainf_base);
			matrix mat_stress_basis_0(-dimensionless_deformation_rate/2, 0, 0,
									  0, dimensionless_deformation_rate, 0,
									  0, 0, -dimensionless_deformation_rate/2);
			matrix mat_stress_basis_3(0, 0, dimensionless_deformation_rate,
									  0, 0, 0,
									  dimensionless_deformation_rate, 0, 0);
			mat_stress_basis_0 = rotation_inv*mat_stress_basis_0*rotation;
			mat_stress_basis_3 = rotation_inv*mat_stress_basis_3*rotation;
			stress_basis_0.setSymmetrize(mat_stress_basis_0);
			stress_basis_3.setSymmetrize(mat_stress_basis_3);
		}
	} else {
		cerr << " dimensionlessnumber = " << control_value.value << endl;
		Sym2Tensor Einf_common = {0, 0, 0, 0, 0, 0};
		vec3d Omegainf(0, 0, 0);
		sys.setImposedFlow(Einf_common, Omegainf);
		sys.zero_shear = true;
	}
}

void Simulation::setupSimulation(string in_args,
                                 vector<string>& input_files,
                                 bool binary_conf,
                                 Dimensional::DimensionalQty<double> control_value,
                                 string flow_type,
                                 string simu_identifier)
{
	/**
	 \brief Set up the simulation.

		This function is intended to be generically used to set up the simulation. It processes the input parameters, non-dimensionalizes the system and starts up a System class with the relevant parameters.
	 */
	string indent = "  Simulation::\t";
	cout << indent << "Simulation setup starting... " << endl;
	string filename_import_positions = input_files[0];
	string filename_parameters = input_files[1];
	sys.p.flow_type = flow_type; // shear or extension or mix (not implemented yet)
	/*
	 * @@@@ This way to prepare relaxed initial configuration should be changed.
	 */
	if (filename_parameters.find("init_relax", 0) != string::npos) {
		cout << "init_relax" << endl;
		sys.zero_shear = true;
	} else {
		sys.zero_shear = false;
	}
	
	setupFlow(control_value);
	setDefaultParameters(control_value);
	readParameterFile(filename_parameters);
	if (sys.ext_flow) {
		p.origin_zero_flow = false;
	}
	setupOptionalSimulation(indent);
	setupNonDimensionalization(control_value);

	assertParameterCompatibility();

	if (input_files[3] != "not_given") {
		throw runtime_error("pre-simulation data deprecated?");
	}
	resolveTimeOrStrainParameters(dimensional_input_params);
	setFromMap(p, dimensional_input_params);

	setConfigToSystem(binary_conf, filename_import_positions);
	 //@@@@ temporary repair
	if (input_files[2] != "not_given") {
		if (sys.brownian && !p.auto_determine_knkt) {
			contactForceParameterBrownian(input_files[2]);
		} else {
			contactForceParameter(input_files[2]);
		}
	}

	p_initial = p;
	sys.resetContactModelParameer(); //@@@@ temporary repair


	if (!sys.ext_flow) {
		// simple shear
		sys.setVelocityDifference();
	} else {
		// extensional flow
		sys.vel_difference.reset();
	}
	if (!restart_from_chkp) {
		simu_name = prepareSimulationName(binary_conf, filename_import_positions, filename_parameters,
		                                  simu_identifier, control_value);
	}
	openOutputFiles();
	echoInputFiles(in_args, input_files);
	cout << indent << "Simulation setup [ok]" << endl;
}

void Simulation::autoSetParameters(const string &keyword, const string &value)
{
	/**
	 \brief Parse an input parameter
	 */
	string numeral, suffix;
	if (keyword == "simulation_mode") {
		p.simulation_mode = atoi(value.c_str());
	} else if (keyword == "lubrication_model") {
		p.lubrication_model = value;
	} else if (keyword == "friction_model") {
		p.friction_model = atoi(value.c_str());
	} else if (keyword == "repulsion") {
		system_of_units.add(Dimensional::Unit::repulsion, 
							Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword));
	} else if (keyword == "cohesion") {
		system_of_units.add(Dimensional::Unit::cohesion, 
							Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword));
	} else if (keyword == "brownian") {
		system_of_units.add(Dimensional::Unit::brownian, 
							Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword));
	} else if (keyword == "critical_load") {
		system_of_units.add(Dimensional::Unit::critical_load, 
							Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword));
	} else if (keyword == "monolayer") {
		p.monolayer = str2bool(value);
	} else if (keyword == "repulsive_length") {
		p.repulsive_length = atof(value.c_str());
	} else if (keyword == "repulsive_max_length") {
		p.repulsive_max_length = atof(value.c_str());
	} else if (keyword == "lub_reduce_parameter") {
		p.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::Time, value, keyword);
	} else if (keyword == "contact_relaxation_time_tan"){
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::Time, value, keyword);
	} else if (keyword == "disp_max") {
		p.disp_max = atof(value.c_str());
	} else if (keyword == "dt_max") {
		p.dt_max = atof(value.c_str());
	} else if (keyword == "dt_min") {
		p.dt_min = atof(value.c_str());
	} else if (keyword == "time_end") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::Time, value, keyword);
	} else if (keyword == "integration_method") {
		p.integration_method = atoi(value.c_str());
	} else if (keyword == "lub_max_gap") {
		p.lub_max_gap = atof(value.c_str());
	} else if (keyword == "interaction_range") {
		p.interaction_range = atof(value.c_str());
	} else if (keyword == "sd_coeff") {
		p.sd_coeff = atof(value.c_str());
	} else if (keyword == "kn") {
		system_of_units.add(Dimensional::Unit::kn, 
							Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword));
	} else if (keyword == "kt") {
		system_of_units.add(Dimensional::Unit::kt, 
							Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword));
	} else if (keyword == "kr") {
		system_of_units.add(Dimensional::Unit::kr,
							Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword));
	} else if (keyword == "dt") {
		p.dt = atof(value.c_str());
	} else if (keyword == "Pe_switch") {
		p.Pe_switch = atof(value.c_str());
	} else if (keyword == "mu_static") {
		p.mu_static = atof(value.c_str());
	} else if (keyword == "mu_dynamic") {
		p.mu_dynamic = atof(value.c_str());
	} else if (keyword == "mu_rolling") {
		p.mu_rolling = atof(value.c_str());
	} else if (keyword == "time_interval_output_config") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::TimeOrStrain, value, keyword);
	} else if (keyword == "time_interval_output_data") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::TimeOrStrain, value, keyword);
	} else if (keyword == "log_time_interval") {
		p.log_time_interval = str2bool(value);
	} else if (keyword == "initial_log_time") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::TimeOrStrain, value, keyword);
	} else if (keyword == "nb_output_data_log_time") {
		p.nb_output_data_log_time = atoi(value.c_str());
	} else if (keyword == "nb_output_config_log_time") {
		p.nb_output_config_log_time = atoi(value.c_str());
	} else if (keyword == "out_data_particle") {
		p.out_data_particle = str2bool(value);
	} else if (keyword == "out_data_interaction") {
		p.out_data_interaction = str2bool(value);
	} else if (keyword == "out_data_vel_components") {
		p.out_data_vel_components = str2bool(value);
	} else if (keyword == "origin_zero_flow") {
		p.origin_zero_flow = str2bool(value);
	} else if (keyword == "auto_determine_knkt") {
		p.auto_determine_knkt = str2bool(value.c_str());
	} else if (keyword == "overlap_target") {
		p.overlap_target = atof(value.c_str());
	} else if (keyword == "disp_tan_target") {
		p.disp_tan_target = atof(value.c_str());
	} else if (keyword == "memory_strain_avg") {
		p.memory_strain_avg = atof(value.c_str());
	} else if (keyword == "memory_strain_k") {
		p.memory_strain_k = atof(value.c_str());
	} else if (keyword == "start_adjust") {
		p.start_adjust = atof(value.c_str());
	} else if (keyword == "min_kn_auto_det") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword);
	} else if (keyword == "max_kn_auto_det") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword);
	} else if (keyword == "min_kt_auto_det") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword);
	} else if (keyword == "max_kt_auto_det") {
		dimensional_input_params[keyword] = Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword);
	} else if (keyword == "min_dt_auto_det") {
		p.min_dt_auto_det = atof(value.c_str());
	} else if (keyword == "max_dt_auto_det") {
		p.max_dt_auto_det = atof(value.c_str());
	} else if (keyword == "rest_threshold") {
		p.rest_threshold = atof(value.c_str());
	} else if (keyword == "ft_max") {
		system_of_units.add(Dimensional::Unit::ft_max, Dimensional::str2DimensionalQty(Dimensional::Dimension::Force, value, keyword));
	} else if (keyword == "fixed_dt") {
		p.fixed_dt = str2bool(value);
	} else if (keyword == "theta_shear") {
		p.theta_shear = atof(value.c_str());
		p.theta_shear *= M_PI/180;  // convert in radians
	} else if (keyword == "strain_reversal") {
		p.strain_reversal = atof(value.c_str());
	} else if (keyword == "event_handler") {
		p.event_handler = value;
		p.event_handler.erase(remove(p.event_handler.begin(), p.event_handler.end(), '\"' ), p.event_handler.end());
	} else if (keyword == "out_particle_stress") {
		p.out_particle_stress = value;
		p.out_particle_stress.erase(remove(p.out_particle_stress.begin(), p.out_particle_stress.end(), '\"' ), p.out_particle_stress.end());
	} else if (keyword == "out_bond_order_parameter6") {
		p.out_bond_order_parameter6 = str2bool(value);
	} else if (keyword == "out_binary_conf") {
		p.out_binary_conf = str2bool(value);
	} else if (keyword == "np_fixed") {
		p.np_fixed = atoi(value.c_str());
	} else if (keyword == "keep_input_strain") {
		p.keep_input_strain = str2bool(value);
	} else if (keyword == "brownian_relaxation_time") {
		p.brownian_relaxation_time = atof(value.c_str());
	} else {
		ostringstream error_str;
		error_str  << "keyword " << keyword << " is not associated with an parameter" << endl;
		throw runtime_error(error_str.str());
	}
}

void Simulation::readParameterFile(const string& filename_parameters)
{
	/**
	 \brief Read and parse the parameter file
	 */
	ifstream fin;
	fin.open(filename_parameters.c_str());
	if (!fin) {
		ostringstream error_str;
		error_str  << " Parameter file '" << filename_parameters << "' not found." <<endl;
		throw runtime_error(error_str.str());
	}
	string keyword, value;
	while (!fin.eof()) {
		string line;
		if (!getline(fin, line, ';')) {
			break;
		}
		if (fin.eof()) {
			break;
		}
		string str_parameter;
		removeBlank(line);
		str_parameter = line;
		string::size_type begin_comment;
		string::size_type end_comment;
		do {
			begin_comment = str_parameter.find("/*");
			end_comment = str_parameter.find("*/");
			if (begin_comment > 10000) {
				break;
			}
			str_parameter = str_parameter.substr(end_comment+2);
		} while (true);
		if (begin_comment > end_comment) {
			cerr << str_parameter.find("/*") << endl;
			cerr << str_parameter.find("*/") << endl;
			throw runtime_error("syntax error in the parameter file.");
		}
		string::size_type pos_slashslash = str_parameter.find("//");
		if (pos_slashslash != string::npos) {
			throw runtime_error(" // is not the syntax to comment out. Use /* comment */");
		}
		Str2KeyValue(str_parameter, keyword, value);
		autoSetParameters(keyword, value);
	}
	fin.close();
	return;
}

void Simulation::setDefaultParameters(Dimensional::DimensionalQty<double> control_value)
{

	/**
	 \brief Set default values for ParameterSet parameters.
	 */
	auto input_scale = control_value.unit;
	auto input_scale_str = Dimensional::unit2suffix(input_scale);

	autoSetParameters("Pe_switch", "5");
	autoSetParameters("dt", "1e-4");
	autoSetParameters("disp_max", "1e-3");
	autoSetParameters("dt_max", "-1");
	autoSetParameters("dt_max", "-1");
	autoSetParameters("monolayer", "false");
	autoSetParameters("rest_threshold", "1e-4");
	autoSetParameters("integration_method", "1");
	autoSetParameters("np_fixed", "0");
	autoSetParameters("sd_coeff", "1");
	autoSetParameters("lubrication_model", "tangential");
	autoSetParameters("friction_model", "1");
	autoSetParameters("time_end", "10h");
	autoSetParameters("lub_max_gap", "0.5");
	autoSetParameters("interaction_range", "-1");
	autoSetParameters("lub_reduce_parameter", "1e-3");
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
	autoSetParameters("min_dt_auto_det", "1e-7");
	autoSetParameters("max_dt_auto_det", "1e-3");
	if (input_scale != Dimensional::Unit::kr) {
		autoSetParameters("kr", "0kn");
	}
	autoSetParameters("auto_determine_knkt", "false");
	autoSetParameters("overlap_target", "0.05");
	autoSetParameters("disp_tan_target", "0.05");
	autoSetParameters("memory_strain_avg", "0.01");
	autoSetParameters("memory_strain_k", "0.02");
	autoSetParameters("start_adjust", "0.2");
	autoSetParameters("repulsive_length", "0.05");
	autoSetParameters("repulsive_max_length", "-1");
	autoSetParameters("mu_static", "1");
	autoSetParameters("mu_dynamic", "-1");
	autoSetParameters("mu_rolling", "0");
	autoSetParameters("time_interval_output_data", "1e-2h");
	autoSetParameters("time_interval_output_config", "1e-1h");
	autoSetParameters("log_time_interval", "false");
	autoSetParameters("initial_log_time", "1e-4h");
	autoSetParameters("nb_output_data_log_time", "100");
	autoSetParameters("nb_output_config_log_time", "100");
	autoSetParameters("origin_zero_flow", "true");
	autoSetParameters("out_data_particle", "true");
	autoSetParameters("out_data_interaction", "true");
	autoSetParameters("out_particle_stress", "");
	autoSetParameters("out_binary_conf", "false");
	autoSetParameters("out_data_vel_components", "false");
	autoSetParameters("fixed_dt", "false");
	autoSetParameters("dt_min", "-1");
	autoSetParameters("dt_max", "-1");
	autoSetParameters("theta_shear", "0");
	autoSetParameters("event_handler", "");
	autoSetParameters("simulation_mode", "0");
	autoSetParameters("keep_input_strain", "false");
	autoSetParameters("brownian_relaxation_time", "1");
	autoSetParameters("out_bond_order_parameter6", "false");
}

//inline string columnDefinition(int &cnb, const string &type, const string &name)
//{
//	stringstream defs;
//	if (type == "vec3d") {
//		array<string, 3> xyz = {"x", "y", "z"};
//		for (auto &u : xyz) {
//			stringstream col_def_complement;
//			defs << "#" << cnb << ": "<< name << " " << u << "\n";
//			cnb ++;
//		}
//	} else if (type=="scalar") {
//		defs << "#" << cnb << ": "<< name << "\n";
//	} else {
//		throw runtime_error(" unknown type for column def\n");
//	}
//	return defs.str();
//}

void Simulation::openOutputFiles()
{
	/**
	 \brief Set up the output files

		This function determines a simulation name from the parameters, opens the output files with the corresponding name and prints their header.
	 */

	stringstream data_header;
	createDataHeader(data_header);
	outdata.setFile("data_"+simu_name+".dat",
	                data_header.str(), force_to_run, restart_from_chkp);
	outdata_st.setFile("st_"+simu_name+".dat",
	                   data_header.str(), force_to_run, restart_from_chkp);

	if (!p.out_particle_stress.empty()) {
		outdata_pst.setFile("pst_"+simu_name+".dat",
		                    data_header.str(), force_to_run, restart_from_chkp);

	}
	string time_filename = "t_"+simu_name+".dat";
	fout_time.open(time_filename.c_str());
	string input_filename = "input_"+simu_name+".dat";
	if (!restart_from_chkp) {
		fout_input.open(input_filename.c_str());
	} else {
		fout_input.open(input_filename.c_str(), fstream::out | fstream::app);
	}
	if (p.out_data_particle) {
		outdata_par.setFile("par_"+simu_name+".dat",
		                    data_header.str(), force_to_run, restart_from_chkp);

	}
	if (p.out_data_interaction) {
		outdata_int.setFile("int_"+simu_name+".dat",
		                    data_header.str(), force_to_run, restart_from_chkp);
	}
	//string box_name = "box_"+simu_name+".dat";
	//fout_boxing.open(box_name);
	return;
}

string Simulation::prepareSimulationName(bool binary_conf,
			                             const std::string& filename_import_positions,
			                             const std::string& filename_parameters,
			                             const std::string& simu_identifier,
			                             Dimensional::DimensionalQty<double> control_value)
{
	/**
	 \brief Determine simulation name.
	 */
	ostringstream ss_simu_name;
	string::size_type pos_name_end = filename_import_positions.find_last_of(".");
	string::size_type param_name_end = filename_parameters.find_last_of(".");
	string::size_type pos_name_start = filename_import_positions.find_last_of("/");
	string::size_type param_name_start = filename_parameters.find_last_of("/");
	if (pos_name_start == std::string::npos) {
		pos_name_start = -1;
	}
	if (param_name_start == std::string::npos) {
		param_name_start = -1;
	}
	pos_name_start += 1;
	param_name_start += 1;
	ss_simu_name << filename_import_positions.substr(pos_name_start, pos_name_end-pos_name_start);
	ss_simu_name << "_";
	ss_simu_name << filename_parameters.substr(param_name_start, param_name_end-param_name_start);
	ostringstream string_control_parameters;
	if (control_var == Parameters::ControlVariable::rate) {
		string_control_parameters << "_" << "rate";
	}
	if (control_var == Parameters::ControlVariable::stress) {
		string_control_parameters << "_" << "stress";
	}
	// if (control_var == Parameters::ControlVariable::viscnb) {
	// 	string_control_parameters << "_" << "viscnb";
	// }
	string_control_parameters << control_value.value << Dimensional::unit2suffix(control_value.unit);
	ss_simu_name << string_control_parameters.str();
	ss_simu_name << "_" << sys.p.flow_type;
	if (simu_identifier != "") {
		ss_simu_name << "_";
		ss_simu_name << simu_identifier;
	}
	string indent = "  Simulation::\t";
	cout << indent << "filename: " << ss_simu_name.str() << endl;

	return ss_simu_name.str();
}

TimeKeeper Simulation::initTimeKeeper()
{
	TimeKeeper tk;
	if (p.log_time_interval) {
		tk.addClock("data", LogClock(p.initial_log_time,
									 p.time_end,
									 p.nb_output_data_log_time,
									 dimensional_input_params["time_end"].dimension == Dimensional::Dimension::Strain));
	} else {
		tk.addClock("data", LinearClock(p.time_interval_output_data,
										dimensional_input_params["time_interval_output_data"].dimension == Dimensional::Dimension::Strain));
	}
	if (p.log_time_interval) {
		tk.addClock("config", LogClock(p.initial_log_time,
									   p.time_end,
									   p.nb_output_config_log_time,
									   dimensional_input_params["time_end"].dimension == Dimensional::Dimension::Strain));
	} else {
		tk.addClock("config", LinearClock(p.time_interval_output_config,
			 							  dimensional_input_params["time_interval_output_data"].dimension == Dimensional::Dimension::Strain));
	}
	return tk;
}
