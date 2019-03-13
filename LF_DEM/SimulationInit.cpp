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

void Simulation::setupNonDimensionalization(Parameters::ParameterSetFactory &PFact)
{
	/**
	 \brief Non-dimensionalize the simulation.
	 
	 This function
	 (1) determines the most appropriate unit scales to use in the System class depending on the input parameters (Brownian/non-Brownian, shear rate, stress/rate controlled)
	 and
	 (2) converts all the input values in these units.
	 */
	indent = "  Simulation::\t";
	Dimensional::Unit internal_unit = Dimensional::Unit::hydro;
	internal_unit = determineUnit(PFact);
	convertForces(internal_unit, PFact);
}

Dimensional::Unit Simulation::determineUnit(Parameters::ParameterSetFactory &PFact)
{
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
					fs.type != Dimensional::Unit::kr &&
					fs.type != Dimensional::Unit::adhesion) {
					largest_force_val = fs.dim_qty.value;
					internal_unit = fs.type;
				}
			}
			if (internal_unit == Dimensional::Unit::brownian) {
				sys.brownian = true; // @@@@ to be checked
				sys.brownian_dominated = true;
			}
		} else {
			if (control_value.unit == Dimensional::Unit::brownian) {
				cout << indent << "Brownian at Pe = 0 " << endl;
				internal_unit = Dimensional::Unit::brownian;
				sys.brownian = true; // @@@@ to be checked
				sys.brownian_dominated = true;
				sys.zero_shear = true;
			} else if (control_value.unit == Dimensional::Unit::repulsion) {
				cout << indent << "non-Brownain at rate = 0 " << endl;
				internal_unit = Dimensional::Unit::repulsion;
				sys.zero_shear = true;
			} else {
				cout << indent << "non-Brownain at rate = 0 " << endl;
				internal_unit = Dimensional::Unit::kn;
				sys.zero_shear = true;
			}
		}
	} else if (control_var == Parameters::ControlVariable::stress) {
		system_of_units.add(Dimensional::Unit::stress, control_value);
		internal_unit = control_value.unit;
		//		internal_unit = Dimensional::Unit::stress;
	}
	return internal_unit;
}

void Simulation::convertForces(Dimensional::Unit &internal_unit,
							   Parameters::ParameterSetFactory &PFact)
{
	// set the internal unit to actually determine force and parameter non-dimensionalized values
	system_of_units.setInternalUnit(internal_unit);
	PFact.setSystemOfUnits(system_of_units);
	cout << indent << "internal units = " << Dimensional::unit2suffix(internal_unit) << endl;

	// set the output unit
	output_unit = control_value.unit;
	cout << indent << "output units = " << Dimensional::unit2suffix(output_unit) << endl;

	// when there is a hydro force, its value is the non-dimensionalized shear rate.
	auto forces = system_of_units.getForceScales();
	if (control_var == Parameters::ControlVariable::rate) {
		if (!sys.zero_shear) {
			sys.set_shear_rate(forces.at(Dimensional::Unit::hydro).value);
		}
	}
	if (control_var == Parameters::ControlVariable::stress) {
		sys.target_stress = forces.at(Dimensional::Unit::stress).value;
	}
}

void Simulation::assertParameterCompatibility()
{
	// test for incompatibilities
	if (sys.brownian == true) {
		if (sys.pairwise_resistance && sys.p.integration_method != 1) {
			 // @@@@ This test is broken as System has not yet set pairwise resistance. For now test is duplicated later on in System
			ostringstream error_str;
			error_str << "Brownian simulation needs to use the Predictor-Corrector method." << endl;
			error_str << "Modify the parameter file." << endl;
			throw runtime_error(error_str.str());
		}
	}
	if (control_var == Parameters::ControlVariable::stress) {
		if (sys.p.integration_method != 0) {
			cerr << "Warning : use of the Predictor-Corrector method for the stress controlled simulation is experimental." << endl;
		}
		//p.integration_method = 0;
	}
	if (sys.critical_load_model) {
		sys.p.friction_model = 2;
		cerr << "Warning : critical load simulation -> switched to friction_model=2" << endl;
	}
	if (sys.p.output.recording_interaction_history) {
		cerr << "Interaction history recording needs to use the Euler's Method." << endl;
		sys.p.integration_method = 0;
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

void Simulation::setupFlow()
{
	// @@@ This function is quite messy, should be fixed 
	// when we rewrite shear and extensional in a consistent manner
	/* dot_gamma = 1 --> dot_epsilon = 0;
	 *
	 */
	if (sys.p.flow_type == "extension") {
		sys.ext_flow = true;
	} else {
		sys.ext_flow = false;
	}
	if (control_value.value != 0) {
		double dimensionless_deformation_rate = 0.5;
		if (!sys.ext_flow) {
			/* simple shear flow
			 * shear_rate = 2*dot_epsilon
			 *
			 * The basis tensorial bases, D, E, G are given in the simple shear coordinate.
			 * sigma = -p I + 2*eta*D + 2*lambda0*E + 2*lambda3*G
			 * D = ((0, 0, 1/2), (0, 0, 0), (1/2, 0, 0))
			 * E = ((-1/4, 0, 0), (0, 1/2, 0), (0, 0, -1/4))
			 * G = ((-1/2, 0, 0), (0, 0, 0), (0, 0, 1/2)
			 */
			Einf_base.set(0, 0, 1, 0, 0, 0); // = D
			Omegainf_base.set(0, 1, 0);
			sys.setImposedFlow(dimensionless_deformation_rate*Einf_base, dimensionless_deformation_rate*Omegainf_base);
			stress_basis_0 = {-dimensionless_deformation_rate/2, 0, 0, 0,
				dimensionless_deformation_rate, -dimensionless_deformation_rate/2}; // = E
			stress_basis_3 = {-dimensionless_deformation_rate, 0, 0, 0, 0, dimensionless_deformation_rate}; // = G
		} else {
			/* extensional flow
			 *
			 */
			sys.p.magic_angle = atan(0.5*(sqrt(5)-1)); // simulation box needs to be tilted in this angle.
			matrix grad_u_orig(dimensionless_deformation_rate, 0, 0,
							   0, 0, 0,
							   0, 0, -dimensionless_deformation_rate);
			matrix rotation, rotation_inv;
			rotation.set_rotation(-sys.p.magic_angle, 'y');
			rotation_inv.set_rotation(sys.p.magic_angle, 'y');
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
		Einf_base.set(0, 0, 1, 0, 0, 0);
		Omegainf_base.set(0, 1, 0);
		Sym2Tensor Einf_zero = {0, 0, 0, 0, 0, 0};
		vec3d Omegainf_zero(0, 0, 0);
		sys.setImposedFlow(Einf_zero, Omegainf_zero);
		sys.zero_shear = true;
	}
}

void Simulation::setupSimulation(string in_args,
								 vector<string>& input_files,
								 bool binary_conf,
								 string simu_identifier)
{
	/**
	 \brief Set up the simulation.

		This function is intended to be generically used to set up the simulation. It processes the input parameters, non-dimensionalizes the system and starts up a System class with the relevant parameters.
	 */
	cout << indent << "Simulation setup starting... " << endl;
	string filename_import_positions = input_files[0];
	string filename_parameters = input_files[1];
	Dimensional::Unit guarranted_unit; // a unit we're sure will mean something, for ParameterSetFactory to set default dimensional qties.
	if (control_var == Parameters::ControlVariable::rate) {
		if (control_value.value != 0) {
			guarranted_unit = Dimensional::Unit::hydro;
		} else {
			guarranted_unit = control_value.unit;
		}
	} else if (control_var == Parameters::ControlVariable::stress) {
		guarranted_unit = control_value.unit;
	} else {
		ostringstream error_str;
		error_str  << "control_var is not set properly." << endl;
		throw runtime_error(error_str.str());
	}

	Parameters::ParameterSetFactory PFactory(guarranted_unit);
	PFactory.setFromFile(filename_parameters);
	setupNonDimensionalization(PFactory);
	if (control_var == Parameters::ControlVariable::stress) {
		target_stress_input = control_value.value; //@@@ Where should we set the target stress???
		sys.target_stress = target_stress_input/6/M_PI; //@@@
	}
	sys.p = PFactory.getParameterSet();
	setupFlow(); // Including parameter p setting.
	if (sys.ext_flow) {
		sys.p.output.origin_zero_flow = false;
	}
	if (sys.p.output.relative_position_view) {
		sys.p.output.origin_zero_flow = false;
	}
	setupOptionalSimulation(); // @@@ To be removed

	assertParameterCompatibility();

	setConfigToSystem(binary_conf, filename_import_positions);

	p_initial = sys.p;
	
	sys.resetContactModelParameer(); //@@@@ temporary repair

	if (!sys.ext_flow) {
		// simple shear
		cerr << "simple shear " << endl;
		sys.setVelocityDifference();
	} else {
		// extensional flow
		cerr << "extensional flow " << endl;
		sys.vel_difference.reset();
	}
	if (simu_name.empty()) {
		simu_name = prepareSimulationName(binary_conf, filename_import_positions, filename_parameters,
										  simu_identifier);
	}
	openOutputFiles();
	echoInputFiles(in_args, input_files);
	checkDispersionType();
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

	if (!sys.p.output.out_particle_stress.empty()) {
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
	if (sys.p.output.out_data_particle) {
		outdata_par.setFile("par_"+simu_name+".dat",
							data_header.str(), force_to_run, restart_from_chkp);
	}
	if (sys.p.output.out_data_interaction) {
		outdata_int.setFile("int_"+simu_name+".dat",
							data_header.str(), force_to_run, restart_from_chkp);
	}
	if (sys.p.output.out_gsd) {
		string gsd_filename = simu_name+".gsd";
		gsd_create(gsd_filename.c_str(), "CIL", "hoomd", gsd_make_version(1, 1));
		gsd_open(&gsdOut, gsd_filename.c_str() , GSD_OPEN_APPEND);
	}
	//string box_name = "box_"+simu_name+".dat";
	//fout_boxing.open(box_name);
	return;
}

string Simulation::prepareSimulationName(bool binary_conf,
										 const std::string& filename_import_positions,
										 const std::string& filename_parameters,
										 const std::string& simu_identifier)
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
	if (sys.p.output.log_time_interval) {
		tk.addClock("data", LogClock(sys.p.output.initial_log_time.value,
									 sys.p.time_end.value,
									 sys.p.output.nb_output_data_log_time,
									 sys.p.time_end.dimension == Dimensional::Dimension::Strain));
	} else {
		tk.addClock("data", LinearClock(sys.p.output.time_interval_output_data.value,
										sys.p.output.time_interval_output_data.dimension == Dimensional::Dimension::Strain));
	}
	if (sys.p.output.log_time_interval) {
		tk.addClock("config", LogClock(sys.p.output.initial_log_time.value,
									   sys.p.time_end.value,
									   sys.p.output.nb_output_config_log_time,
									   sys.p.time_end.dimension == Dimensional::Dimension::Strain));
	} else {
		tk.addClock("config", LinearClock(sys.p.output.time_interval_output_config.value,
										  sys.p.output.time_interval_output_config.dimension == Dimensional::Dimension::Strain));
	}
	return tk;
}

void Simulation::checkDispersionType()
{
	int cnt_type = 0;
	np1 = sys.get_np();
	for (int i=0; i<sys.get_np()-1; i++) {
		if (sys.radius[i+1] != sys.radius[i]) {
			cnt_type ++;
			np1 = i+1;
		}
	}
	if (cnt_type == 0) {
		dispersion_type = DispersionType::mono;
	} else if (cnt_type == 1) {
		dispersion_type = DispersionType::bi;
	} else {
		dispersion_type = DispersionType::bi;
	}
}

