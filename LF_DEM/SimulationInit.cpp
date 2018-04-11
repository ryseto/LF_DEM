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
#include "ParameterSetFactory.h"

using namespace std;

void Simulation::contactForceParameter(string filename)
{
	/**
	 \brief Load a file containing spring constants and time steps as functions of the volume fraction.

	 Input file must be formatted as:
	 phi kn kt dt
	 */
	auto conf = sys.getBaseConfiguration();
	ifstream fin_knktdt;
	fin_knktdt.open(filename.c_str());
	if (!fin_knktdt) {
		ostringstream error_str;
		error_str  << " Contact parameter file '" << filename << "' not found." << endl;
		throw runtime_error(error_str.str());
	}
	// temporal variables to keep imported values.
	double phi_, kn_, kt_, dt_;
	// To find parameters for considered volume fraction phi.
	bool found = false;
	while (fin_knktdt >> phi_ >> kn_ >> kt_ >> dt_) {
		if (abs(phi_-conf.volume_or_area_fraction) < 1e-10) {
			found = true;
			break;
		}
	}
	fin_knktdt.close();
	if (found) {
		// Set the parameter object
		p.kn = kn_, p.kt = kt_, p.dt = dt_;
		string indent = "  Simulation::\t";
		cout << indent << "Input for kn, kt, dt = " << phi_ << ' ' << kn_ << ' ' << kt_ << ' ' << dt_ << endl;
	} else {
		ostringstream error_str;
		error_str  << " Error: file " << filename.c_str() << " contains no data for vf = " << conf.volume_or_area_fraction << endl;
		throw runtime_error(error_str.str());
	}
}

void Simulation::contactForceParameterBrownian(string filename)
{
	/**
	 \brief Load a file containing spring constants and time steps as functions of the volume fraction and Peclet number.

	 Input file must be formatted as:
	 phi peclet kn kt dt
	 */
	auto conf = sys.getBaseConfiguration();
	ifstream fin_knktdt;
	fin_knktdt.open(filename.c_str());
	if (!fin_knktdt) {
		ostringstream error_str;
		error_str  << " Contact parameter file '" << filename << "' not found." << endl;
		throw runtime_error(error_str.str());
	}
	// temporal variables to keep imported values.
	double phi_, peclet_, kn_, kt_, dt_;
	bool found = false;
	auto forces = system_of_units.getForceScales();
	auto peclet = 1./forces.at(Dimensional::Unit::brownian).value;
	while (fin_knktdt >> phi_ >> peclet_ >> kn_ >> kt_ >> dt_) {
		if (abs(phi_-conf.volume_or_area_fraction) < 1e-10 && peclet_ == peclet) {
			found = true;
			break;
		}
	}
	fin_knktdt.close();

	if (found) {
		p.kn = kn_, p.kt = kt_, p.dt = dt_;
		string indent = "  Simulation::\t";
		cout << indent << "Input for vf = " << phi_ << " and Pe = " << peclet_ << " : kn = " << kn_ << ", kt = " << kt_ << " and dt = " << dt_ << endl;
	} else {
		ostringstream error_str;
		error_str  << " Error: file " << filename.c_str() << " contains no data for vf = " << conf.volume_or_area_fraction << " and Pe = " << peclet << endl;
		throw runtime_error(error_str.str());
	}
}


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


void Simulation::setupNonDimensionalization(Dimensional::DimensionalQty<double> control_value, 
											Parameters::ParameterSetFactory &PFact)
{
	/**
	 \brief Non-dimensionalize the simulation.

		This function determines the most appropriate unit scales to use in the System class depending on the input parameters (Brownian/non-Brownian, shear rate, stress/rate controlled), and converts all the input values in these units.
	 */
	string indent = "  Simulation::\t";

	// feed in the force scales to the UnitSystem solver
	for (auto &fs: PFact.getForceScales()) {
		system_of_units.add(fs.type, fs.dim_qty);
	}

	// determine the internal unit to be used
	Dimensional::Unit internal_unit = Dimensional::Unit::hydro;
	if (control_var == Parameters::ControlVariable::rate) {// || control_var == Parameters::ControlVariable::viscnb) {
		if (control_value.value != 0) {
			system_of_units.add(Dimensional::Unit::hydro, control_value);
			system_of_units.setInternalUnit(Dimensional::Unit::hydro);
			internal_unit = Dimensional::Unit::hydro;
			double largest_force_val = 1;
			for (auto &fs: PFact.getForceScales()) {
				if (fs.dim_qty.value > largest_force_val && 
					fs.type != Dimensional::Unit::kn &&
					fs.type != Dimensional::Unit::kt &&
					fs.type != Dimensional::Unit::kr) {
					largest_force_val = fs.dim_qty.value;
					internal_unit = fs.type;
				}
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

	// set the internal unit to actually determine force and parameter non-dimensionalized values 
	system_of_units.setInternalUnit(internal_unit);
	PFact.setSystemOfUnits(system_of_units);
	cout << indent << "internal units = " << Dimensional::unit2suffix(internal_unit) << endl;

	// set the output unit
	output_unit = control_value.unit;
	cout << indent << "output units = " << Dimensional::unit2suffix(output_unit) << endl;
	
	// when there is a hydro force, its value is the non-dimensionalized shear rate.
	auto forces = system_of_units.getForceScales();
	if (forces.count(Dimensional::Unit::hydro) > 0) { // == if rate controlled
		sys.set_shear_rate(forces.at(Dimensional::Unit::hydro).value);
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

void Simulation::setConfigToSystem(bool binary_conf, const std::string &filename)
{
	if (binary_conf) {
		auto format = getBinaryConfigurationFileFormat(filename);
		switch(format) {
			case ConfFileFormat::bin_format_base_shear:
			case ConfFileFormat::bin_format_base_new:
				{
					auto conf = readBinaryBaseShearConfiguration(filename);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case ConfFileFormat::bin_format_fixed_vel_shear:
				{
					auto conf = readBinaryFixedVeloConfiguration(filename);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case ConfFileFormat::bin_delayed_adhesion:
				{
					auto conf = readBinaryDelayedAdhesionConfiguration(filename);
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
	/*
	 * @@@@ This way to prepare relaxed initial configuration should be changed.
	 */
	if (filename_parameters.find("init_relax", 0) != string::npos) {
		cout << "init_relax" << endl;
		sys.zero_shear = true;
	} else {
		sys.zero_shear = false;
	}
	
	Parameters::ParameterSetFactory PFactory;
	PFactory.setFromFile(filename_parameters);
	setupNonDimensionalization(control_value, PFactory);
	
	if (control_var == Parameters::ControlVariable::stress) {
		target_stress_input = control_value.value; //@@@ Where should we set the target stress???
		sys.target_stress = target_stress_input/6/M_PI; //@@@
	}
		
	p = PFactory.getParameterSet();

    setupFlow(control_value); // Including parameter p setting.
    
	p.flow_type = flow_type; // shear or extension or mix (not implemented yet)

	if (sys.ext_flow) {
		p.output.origin_zero_flow = false;
	}
	setupOptionalSimulation(indent);

	assertParameterCompatibility();

	if (input_files[3] != "not_given") {
		throw runtime_error("pre-simulation data deprecated?");
	}

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
	if (simu_name.empty()) {
		simu_name = prepareSimulationName(binary_conf, filename_import_positions, filename_parameters,
		                                  simu_identifier, control_value);
	}
	openOutputFiles();
	echoInputFiles(in_args, input_files);
	cout << indent << "Simulation setup [ok]" << endl;
}


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

	if (!p.output.out_particle_stress.empty()) {
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
	if (p.output.out_data_particle) {
		outdata_par.setFile("par_"+simu_name+".dat",
		                    data_header.str(), force_to_run, restart_from_chkp);

	}
	if (p.output.out_data_interaction) {
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
	ss_simu_name << "_" << p.flow_type;
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
	if (p.output.log_time_interval) {
		tk.addClock("data", LogClock(p.output.initial_log_time.value,
									 p.time_end.value,
									 p.output.nb_output_data_log_time,
									 p.time_end.dimension == Dimensional::Dimension::Strain));
	} else {
		tk.addClock("data", LinearClock(p.output.time_interval_output_data.value,
										p.output.time_interval_output_data.dimension == Dimensional::Dimension::Strain));
	}
	if (p.output.log_time_interval) {
		tk.addClock("config", LogClock(p.output.initial_log_time.value,
									   p.time_end.value,
									   p.output.nb_output_config_log_time,
									   p.time_end.dimension == Dimensional::Dimension::Strain));
	} else {
		tk.addClock("config", LinearClock(p.output.time_interval_output_config.value,
			 							  p.output.time_interval_output_config.dimension == Dimensional::Dimension::Strain));
	}
	return tk;
}
