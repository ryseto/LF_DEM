//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Romain Mari on 08/10/15.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "Simulation.h"
#include <string>
#include <sstream>

using namespace std;

void Simulation::contactForceParameter(string filename)
{
	ifstream fin_knktdt;
	fin_knktdt.open(filename.c_str());
	if (!fin_knktdt) {
		cerr << " Contact parameter file '" << filename << "' not found." << endl;
		exit(1);
	}
	
	// temporal variables to keep imported values.
	double phi_, kn_, kt_, dt_;
	// To find parameters for considered volume fraction phi.
	bool found = false;
	while (fin_knktdt >> phi_ >> kn_ >> kt_ >> dt_) {
		if (phi_ == volume_or_area_fraction) {
			found = true;
			break;
		}
	}
	fin_knktdt.close();
	
	if (found) {
		// Set the parameter object
		p.kn = kn_, p.kt = kt_, p.dt = dt_;
		cout << " Input for kn, kt, dt = " << phi_ << ' ' << kn_ << ' ' << kt_ << ' ' << dt_ << endl;
	} else {
		cerr << " Error: file " << filename.c_str() << " contains no data for vf = " << phi_ << endl;
		exit(1);
	}
}

void Simulation::contactForceParameterBrownian(string filename)
{
	ifstream fin_knktdt;
	fin_knktdt.open(filename.c_str());
	if (!fin_knktdt) {
		cerr << " Contact parameter file '" << filename << "' not found." <<endl;
		exit(1);
	}
	
	// temporal variables to keep imported values.
	double phi_, peclet_, kn_, kt_, dt_;
	bool found = false;
	while (fin_knktdt >> phi_ >> peclet_ >> kn_ >> kt_ >> dt_) {
		if (phi_ == volume_or_area_fraction && peclet_ == dimensionless_numbers["hydro/thermal"]) {
			found = true;
			break;
		}
	}
	fin_knktdt.close();
	
	if (found) {
		p.kn = kn_, p.kt = kt_, p.dt = dt_;
		cout << "Input for vf = " << phi_ << " and Pe = " << peclet_ << " : kn = " << kn_ << ", kt = " << kt_ << " and dt = " << dt_ << endl;
	} else {
		cerr << " Error: file " << filename.c_str() << " contains no data for vf = " << volume_or_area_fraction << " and Pe = " << dimensionless_numbers["hydro/thermal"] << endl;
		exit(1);
	}
}

void Simulation::importPreSimulationData(string filename)
{
	ifstream fin_PreSimulationData;
	fin_PreSimulationData.open(filename.c_str());
	if (!fin_PreSimulationData) {
		cerr << " Pre-simulation data file '" << filename << "' not found." << endl;
		exit(1);
	}
	
	double stress_, shear_rate_;
	while (fin_PreSimulationData >> stress_ >> shear_rate_) {
		if (stress_ == target_stress_input) {
			break;
		}
	}
	shear_rate_expectation = shear_rate_;
}

void Simulation::echoInputFiles(string in_args, vector<string> &input_files)
{
	fout_input << "# LF_DEM version " << GIT_VERSION << ", called with:" << endl;
	fout_input << in_args << endl << endl;
	for (const string &in_file : input_files) {
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

void Simulation::resolveUnitSystem(string unit) // can we express all forces in unit "unit"?
{
	values[unit] = 1;
	suffixes[unit] = unit;
	
	// now resolve the other force units
	set <string> resolved_units;
	resolved_units.clear();
	resolved_units.insert(unit);
	unsigned int resolved = resolved_units.size();
	unsigned int previous_resolved;
	
	do {
		previous_resolved = resolved;
		for (const auto& f: suffixes) {
			string force_type = f.first;
			string suffix = f.second;
			if (resolved_units.find(suffix) != resolved_units.end()) {  // then we know how to convert to unit
				values[force_type] *= values[suffix];
				suffixes[force_type] = unit;
				resolved_units.insert(force_type);
			}
		}
		resolved = resolved_units.size();
	} while (previous_resolved < resolved);
	
	// check we found everyone
	if (resolved < suffixes.size()) {
		for (const auto& f: suffixes) {
			string force_type = f.first;
			string suffix = f.second;
			if (resolved_units.find(suffix) == resolved_units.end()) {
				cerr << "Error: force type \"" << force_type << "\" has an unknown scale \"" << suffix << "\"" << endl;
			}
		}
		exit(1);
	}
	
	// determine the dimensionless_numbers
	for (const auto& f1: suffixes) {
		string force1_type = f1.first;
		for (const auto& f2: suffixes) {
			string force2_type = f2.first;
			dimensionless_numbers[force2_type+'/'+force1_type] = values[force2_type]/values[force1_type];
			dimensionless_numbers[force1_type+'/'+force2_type] = 1/dimensionless_numbers[force2_type+'/'+force1_type];
		}
	}
}

void Simulation::convertInputForcesStressControlled(double dimensionlessnumber, string rate_unit){
	/**
	 \brief Chooses units for the simulation and convert the forces to this unit (stress controlled case).
	 
	 The strategy is the following:
	 1. Convert all the forces in the force unit "w" given by the input stress (LF_DEM -s a_numberw), by solving recursively, ie:
		a. w = 1w
		b. y = mw
		c. z = oy = omw
		d. etc...
		
	 In the future, we may allow other unit scale than the one given by the input stress.
	 */
	
	string force_type = rate_unit; // our force defining the shear rate
	
	if (force_type == "hydro") {
		cerr << " Error: please give a stress in non-hydro units." << endl;
		exit(1);
		/*
		 Note:
		 Although it is in some cases possible to run under stress control without any non-hydro force scale,
		 it is not always possible and as a consequence it is a bit dangerous to let the user do so.
		 
		 With only hydro, the problem is that the target stress \tilde{S} cannot take any possible value, as
		 \tilde{S} = S/(\eta_0 \dot \gamma) = \eta/\eta_0
		 --> It is limited to the available range of viscosity.
		 If you give a \tilde{S} outside this range (for example \tilde{S}=0.5), you run into troubles.
		 */
	}
	if (force_type == "thermal") {
		cerr << " Error: stress controlled Brownian simulations are not yet implemented." << endl;
		exit(1);
	}
	
	sys.set_shear_rate(1);
	// we take as a unit scale the one given by the user with the stress
	// TODO: other choices may be better when several forces are used.
	internal_unit_scales = force_type;
	target_stress_input = dimensionlessnumber;
	sys.target_stress = target_stress_input/6/M_PI;
	
	// convert all other forces to internal_unit_scales
	resolveUnitSystem(internal_unit_scales);
}

void Simulation::convertInputForcesRateControlled(double dimensionlessnumber, string rate_unit)
{
	/**
	 \brief Chooses units for the simulation and convert the forces to this unit (rate controlled case).
	 
	 The strategy is the following:
	 1. Convert all the forces in hydro force unit, by solving recursively, ie:
		a. h = 1h
		b. "LF_DEM -r nx" --> x = (1/n)h (the rate given in input is also the ratio of forces h/x)
		c. y = mx = (m/n)h
		d. z = oy = (om/n)h
		e. etc...
	 2. Decide the unit "w" for the simulation
	 3. Convert all the forces to this unit (by multiplying every force by h/w)
	 */
	string force_type = rate_unit; // our force defining the shear rate
	if (force_type == "hydro") {
		cerr << "Error: cannot define the shear rate in hydro unit." << endl;
		exit(1);
	}
	
	if (values[force_type] > 0) {
		cerr << "Error: redefinition of the rate (given both in the command line and in the parameter file with \"" << force_type << "\" force)" << endl;
		exit(1);
	}
	// switch this force in hydro units
	values[force_type] = 1/dimensionlessnumber;
	suffixes[force_type] = "hydro";
	
	// convert all other forces to hydro
	resolveUnitSystem("hydro");
	
	//	chose simulation unit
	setUnitScaleRateControlled();
	
	// convert from hydro scale to chosen scale
	convertForceValues(internal_unit_scales);
	
	bool is_brownian = dimensionless_numbers.find("hydro/thermal") != dimensionless_numbers.end();
	if (is_brownian) {
		sys.brownian = true;
		p.brownian_amplitude = values["thermal"];
		cout << "Brownian, Peclet number " << dimensionless_numbers["hydro/thermal"] << endl;
	} else {
		cout << "non-Brownian" << endl;
	}
}

void Simulation::convertInputForcesMagnetic(double dimensionlessnumber, string rate_unit)
{
	string force_type = rate_unit;
	if (force_type != "thermal") {
		cerr << "unit needs to be thermal" << endl;
		exit(1);
	}
	if (values[force_type] > 0) {
		cerr << "Error: redefinition of the rate (given both in the command line and in the parameter file with \"" << force_type << "\" force)" << endl;
		exit(1);
	}
	// switch this force in hydro units
	values[force_type] = 1/dimensionlessnumber;
	suffixes[force_type] = "magnetic";
	// convert all other forces to hydro
	resolveUnitSystem("magnetic");
	//	chose simulation unit
	cerr << "setUnitScaleRateControlled" << endl;
	setUnitScaleMagnetic();
	// convert from hydro scale to chosen scale
	convertForceValues(internal_unit_scales);
	sys.brownian = true;
	p.brownian_amplitude = values["thermal"];
	cerr << "Brownian, Peclet number " << dimensionless_numbers["magnetic/thermal"] << endl;
}

void Simulation::setLowPeclet()
{
	sys.lowPeclet = true;
	double scale_factor_SmallPe = p.Pe_switch/dimensionless_numbers["hydro/thermal"];
	p.memory_strain_k /= scale_factor_SmallPe;
	p.memory_strain_avg /= scale_factor_SmallPe;
	p.start_adjust /= scale_factor_SmallPe;
	p.dt *= p.Pe_switch; // to make things continuous at Pe_switch
}

void Simulation::convertForceValues(string new_unit)
{
	for (auto& f: suffixes) {
		string force_type = f.first;
		string old_unit = f.second;
		if (old_unit != new_unit) {
			values[force_type] *= dimensionless_numbers[old_unit+'/'+new_unit];
			f.second = new_unit;
		}
	}
}

void Simulation::setUnitScaleRateControlled()
{
	bool is_brownian;
	
	if (dimensionless_numbers.find("hydro/thermal") != dimensionless_numbers.end()
		|| dimensionless_numbers.find("magnetic/thermal") != dimensionless_numbers.end()) {
		is_brownian = true;
	} else {
		is_brownian = false;
	}

	if (is_brownian) {
		if (dimensionless_numbers["hydro/thermal"] > p.Pe_switch && !sys.zero_shear) { // hydro units
			internal_unit_scales = "hydro";
			sys.amplitudes.sqrt_temperature = 1/sqrt(dimensionless_numbers["hydro/thermal"]);
			sys.set_shear_rate(1);
		} else { // low Peclet mode
			internal_unit_scales = "thermal";
			sys.amplitudes.sqrt_temperature = 1;
			sys.set_shear_rate(dimensionless_numbers["hydro/thermal"]);
			setLowPeclet();
		}
	} else {
		internal_unit_scales = "hydro";
		sys.set_shear_rate(1);
	}
	
	if (control_var == "magnetic") {
		sys.amplitudes.sqrt_temperature = 0.3;
	}
}

void Simulation::setUnitScaleMagnetic()
{
	internal_unit_scales = "thermal";
	sys.amplitudes.sqrt_temperature = 1;
}

void Simulation::exportForceAmplitudes()
{
	bool is_repulsive = values.find("repulsive") != values.end();
	if (is_repulsive) {
		sys.repulsiveforce = true;
		sys.amplitudes.repulsion = values["repulsive"];
		cout << " Repulsive force (in \"" << suffixes["repulsive"] << "\" units): " << sys.amplitudes.repulsion << endl;
	}
	bool is_critical_load = values.find("critical_load") != values.end();
	if (is_critical_load) {
		sys.critical_load = true;
		sys.amplitudes.critical_normal_force = values["critical_load"];
		cout << " Critical Load (in \"" << suffixes["critical_load"] << "\" units): " << sys.amplitudes.critical_normal_force << endl;
	}
	bool is_cohesive = values.find("cohesive") != values.end();
	if (is_cohesive) {
		sys.cohesion = true;
		sys.amplitudes.cohesion = values["cohesive"];
		cout << " Cohesion (in \"" << suffixes["cohesive"] << "\" units): " << sys.amplitudes.cohesion << endl;
	}
	
	//	bool is_magnetic = values.find("m") != values.end();
	//	if (is_magnetic) {
	//		sys.amplitudes.magnetic = values["m"];
	//		cout << " Magnetic force (in \"" << suffixes["m"] << "\" units): " << p.magnetic_amplitude << endl; // unused now, should map to a quantity in sys.amplitudes
	//		cout << " values[m] = "  << values["m"] << endl;
	//	}
	bool is_ft_max = values.find("ft") != values.end();
	if (is_ft_max) {
		sys.amplitudes.ft_max = values["ft"];
		cout << " Max tangential load (in \"" << suffixes["ft"] << "\" units): " << sys.amplitudes.ft_max << endl;
	}
}

void Simulation::convertInputValues(string new_unit)
{
	for (auto& inv: input_values) {
		string old_unit = inv.unit;
		if (old_unit != new_unit) {
			if (old_unit != "hydro" && values.find(old_unit) == values.end()) {
				cerr << " Error: trying to convert " << inv.name << " from an unknown unit \"" << inv.unit 	<< "\"" << endl;
				exit(1);
			}
			if (inv.type == "stiffness") {
				*(inv.value) *= dimensionless_numbers[old_unit+'/'+new_unit];
				inv.unit = new_unit;
			}
			if (inv.type == "time") {
				if (old_unit == "hydro") { // = it is a strain, better keeping it that way
					inv.unit = "strain";
				} else {
					*(inv.value) *= dimensionless_numbers[new_unit+'/'+old_unit];
					inv.unit = new_unit;
				}
			}
		}
		cout << inv.name << " (in \"" << inv.unit << "\" units): " << *(inv.value) << endl;
	}
}

void Simulation::setupNonDimensionalization(double dimensionlessnumber, string input_scale){
	/*
	 * Force unit
	 */
	input_scale = unit_longname[input_scale];
	if (control_var == "rate") {
		input_rate = dimensionlessnumber;
		input_rate_unit = input_scale;
	}
	if (control_var == "rate") {
		convertInputForcesRateControlled(dimensionlessnumber, input_scale);
	} else if (control_var == "stress") {
		convertInputForcesStressControlled(dimensionlessnumber, input_scale);
		p.unscaled_contactmodel = true;
	} else if (control_var == "magnetic") {
		convertInputForcesMagnetic(dimensionlessnumber, input_scale);
	} else {
		exit(1);
	}
	exportForceAmplitudes();
	cerr << "internal_unit_scales = " << internal_unit_scales << endl;

	convertInputValues(internal_unit_scales);

	output_unit_scales = input_scale;

}

void Simulation::setupSimulation(string in_args,
								 vector<string> &input_files,
								 bool binary_conf,
								 double dimensionlessnumber,
								 string input_scale)
{
	filename_import_positions = input_files[0];
	filename_parameters = input_files[1];

	if (filename_parameters.find("init_relax", 0) != string::npos) {
		cerr << "init_relax" << endl;
		sys.zero_shear = true;
	} else if (control_var == "magnetic") {
		sys.zero_shear = true;
	} else {
		sys.zero_shear = false;
	}
	setDefaultParameters();
	readParameterFile();
	for (const auto& f: suffixes) {
		string_control_parameters << "_" << f.first << values[f.first] << f.second;
	}
	if (control_var == "rate") {
		string_control_parameters << "_rate";
	} else if (control_var == "stress") {
		string_control_parameters << "_stress";
	} else if (control_var == "magnetic") {
		string_control_parameters << "_magnetic";
	}
	string_control_parameters << dimensionlessnumber << input_scale;

	setupNonDimensionalization(dimensionlessnumber, input_scale);

	/*
	 * Simulate parameters
	 */
	// test for incompatibilities
	if (sys.brownian == true) {
		if (p.integration_method != 1) {
			cerr << "Brownian simulation needs to use the Predictor-Corrector method." << endl;
			cerr << "Modify the parameter file: " << filename_parameters << endl;
			exit(1);
		}
	}
	if (control_var == "stress") {
		if (p.integration_method != 0) {
			cerr << "Need to be so far the Euler method for the stress controlled simulation." << endl;
		}
		p.integration_method = 0;
	}
	if (sys.critical_load) {
		p.friction_model = 2;
	}

	if (binary_conf) {
		importConfigurationBinary();
	} else {
		importConfiguration();
	}
	sys.shear_disp = initial_lees_edwards_disp;
	
	if (input_files[2] != "not_given") {
		if (sys.brownian && !p.auto_determine_knkt) {
			contactForceParameterBrownian(input_files[2]);
		} else {
			contactForceParameter(input_files[2]);
		}
	}

	for (const auto& inv: input_values) {
		if (inv.name == "time_end") {
			if (inv.unit == "strain") {
				time_end = -1;
				strain_end = p.time_end;
				cout << "  strain_end = " << strain_end << endl;
			} else {
				time_end = p.time_end;
				cout << "  time_end = " << time_end << endl;
			}
		}
	}
	if (control_var == "magnetic") {
		time_end = p.time_end;
	}

	if (input_files[3] != "not_given") {
		importPreSimulationData(input_files[3]);
		time_interval_output_data = p.time_interval_output_data/shear_rate_expectation;
		time_interval_output_config = p.time_interval_output_config/shear_rate_expectation;
	} else {
		for (const auto& inv: input_values) {
			if (inv.name == "time_interval_output_data") {
				if (inv.unit == "strain") {
					time_interval_output_data = -1;
					strain_interval_output_data = p.time_interval_output_data;
					cout << "  strain_interval_output_data = " << strain_interval_output_data << endl;
				} else {
					time_interval_output_data = p.time_interval_output_data;
					cout << "  time_interval_output_data = " << time_interval_output_data << endl;
				}
			}
			if (inv.name == "time_interval_output_config") {
				if (inv.unit == "strain") {
					time_interval_output_config = -1;
					strain_interval_output_config = p.time_interval_output_config;
					cout << "  strain_interval_output_config = " << strain_interval_output_config << endl;
				} else {
					time_interval_output_config = p.time_interval_output_config;
					cout << "  time_interval_output_config = " << time_interval_output_config << endl;
				}
			}
		}
	}
	
	sys.setupSystem(control_var);

	openOutputFiles(binary_conf);
	echoInputFiles(in_args, input_files);
}

void Simulation::autoSetParameters(const string &keyword, const string &value)
{
	string numeral, suffix;
	bool caught_suffix = true;
	if (keyword == "lubrication_model") {
		p.lubrication_model = atoi(value.c_str());
	} else if (keyword == "friction_model") {
		if (p.friction_model == 2) {
			cerr << "!!Neglected friction_model in parameter file!!" << endl;
		} else {
			p.friction_model = atoi(value.c_str());
		}
	} else if (keyword == "repulsion_amplitude") {
		caught_suffix = getSuffix(value, numeral, suffix);
		suffix = unit_longname[suffix];
		p.repulsion_amplitude = atof(numeral.c_str());
		suffixes["repulsive"] = suffix;
		values["repulsive"] = atof(numeral.c_str());
	} else if (keyword == "cohesion_amplitude") {
		caught_suffix = getSuffix(value, numeral, suffix);
		suffix = unit_longname[suffix];
		suffixes["cohesive"] = suffix;
		values["cohesive"] = atof(numeral.c_str());
	} else if (keyword == "brownian_amplitude") {
		caught_suffix = getSuffix(value, numeral, suffix);
		suffix = unit_longname[suffix];
		suffixes["thermal"] = suffix;
		values["thermal"] = atof(numeral.c_str());
	} else if (keyword == "critical_load_amplitude") {
		caught_suffix = getSuffix(value, numeral, suffix);
		suffix = unit_longname[suffix];
		suffixes["critical_load"] = suffix;
		values["critical_load"] = atof(numeral.c_str());
	} else if (keyword == "magnetic_amplitude") {
		caught_suffix = getSuffix(value, numeral, suffix);
		suffix = unit_longname[suffix];
		suffixes["magnetic"] = suffix;
		values["magnetic"] = atof(numeral.c_str());
		cerr << "Need to confirm the implementation for magnetic_amplitude" << endl;
		exit(1);
	} else if (keyword == "monolayer") {
		p.monolayer = str2bool(value);
	} else if (keyword == "unscaled_contactmodel") {
		p.unscaled_contactmodel = str2bool(value);
	} else if (keyword == "repulsiveforce_length") {
		p.repulsive_length = atof(value.c_str());
	} else if (keyword == "repulsive_max_length") {
		p.repulsive_max_length = atof(value.c_str());
	} else if (keyword == "lub_reduce_parameter") {
		p.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time") {
		p.contact_relaxation_time = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time_tan"){
		p.contact_relaxation_time_tan = atof(value.c_str());
	} else if (keyword == "disp_max") {
		p.disp_max = atof(value.c_str());
	} else if (keyword == "time_end") {
		catchSuffixedValue("time", keyword, value, &p.time_end);
	} else if (keyword == "integration_method") {
		p.integration_method = atoi(value.c_str());
	} else if (keyword == "lub_max_gap") {
		p.lub_max_gap = atof(value.c_str());
	} else if (keyword == "interaction_range") {
		p.interaction_range = atof(value.c_str());
	} else if (keyword == "sd_coeff") {
		p.sd_coeff = atof(value.c_str());
	} else if (keyword == "kn") {
		catchSuffixedValue("stiffness", keyword, value, &p.kn);
	} else if (keyword == "kt") {
		catchSuffixedValue("stiffness", keyword, value, &p.kt);
	} else if (keyword == "kr") {
		catchSuffixedValue("stiffness", keyword, value, &p.kr);
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
		catchSuffixedValue("time", keyword, value, &p.time_interval_output_config);
	} else if (keyword == "time_interval_output_data") {
		catchSuffixedValue("time", keyword, value, &p.time_interval_output_data);
	} else if (keyword == "out_data_particle") {
		p.out_data_particle = str2bool(value);
	} else if (keyword == "out_data_interaction") {
		p.out_data_interaction = str2bool(value);
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
	} else if (keyword == "min_kn") {
		p.min_kn = atof(value.c_str());
	} else if (keyword == "max_kn") {
		p.max_kn = atof(value.c_str());
	} else if (keyword == "min_kt") {
		p.min_kt = atof(value.c_str());
	} else if (keyword == "max_kt") {
		p.max_kt = atof(value.c_str());
	} else if (keyword == "rest_threshold") {
		p.rest_threshold = atof(value.c_str());
	} else if (keyword == "ft_max") {
		caught_suffix = getSuffix(value, numeral, suffix);
		suffix = unit_longname[suffix];
		suffixes["ft"] = suffix;
		values["ft"] = atof(numeral.c_str());
	} else if (keyword == "fixed_dt") {
		p.fixed_dt = str2bool(value);
	} else if (keyword == "magnetic_type") {
		p.magnetic_type = atoi(value.c_str());
	} else if (keyword == "magnetic_field_type") {
		p.magnetic_field_type = atoi(value.c_str());
	} else if (keyword == "init_angle_external_magnetic_field") {
		p.init_angle_external_magnetic_field = atof(value.c_str());
	} else if (keyword == "rot_step_external_magnetic_field") {
		p.rot_step_external_magnetic_field = atof(value.c_str());
	} else if (keyword == "step_interval_external_magnetic_field") {
		p.step_interval_external_magnetic_field = atof(value.c_str());
	} else if (keyword == "magnetic_interaction_range") {
		p.magnetic_interaction_range = atof(value.c_str());
	} else if (keyword == "cross_shear") {
		sys.cross_shear = str2bool(value); // temporary, should get a parameter
	} else {
		cerr << "keyword " << keyword << " is not associated with an parameter" << endl;
		exit(1);
	}
	if (!caught_suffix) {
		errorNoSuffix(keyword);
	}
}

void Simulation::readParameterFile()
{
	ifstream fin;
	fin.open(filename_parameters.c_str());
	if (!fin) {
		cerr << " Parameter file '" << filename_parameters << "' not found." <<endl;
		exit(1);
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
			cerr << "syntax error in the parameter file." << endl;
			exit(1);
		}
		string::size_type pos_slashslash = str_parameter.find("//");
		if (pos_slashslash != string::npos) {
			cerr << " // is not the syntax to comment out. Use /* comment */" << endl;
			exit(1);
		}
		Str2KeyValue(str_parameter, keyword, value);
		autoSetParameters(keyword, value);
	}
	fin.close();
	return;
}

void Simulation::openOutputFiles(bool binary_conf)
{
	/*
	 * Set simulation name and name of output files.
	 */
	prepareSimulationName(binary_conf);
	string st_filename = "st_" +sys.simu_name + ".dat";
	fout_st.open(st_filename.c_str());
	outputDataHeader(fout_st);
	string data_filename = "data_" + sys.simu_name + ".dat";
	fout_data.open(data_filename.c_str());
	string time_filename = "t_" + sys.simu_name + ".dat";
	fout_time.open(time_filename.c_str());
	string input_filename = "input_" + sys.simu_name + ".dat";
	fout_input.open(input_filename.c_str());
	
	outputDataHeader(fout_data);
	
	if (p.out_data_particle) {
		string particle_filename = "par_" + sys.simu_name + ".dat";
		fout_particle.open(particle_filename.c_str());
		outputDataHeader(fout_particle);
		//
		string fout_par_col_def =
		"#1: number of the particle\n"
		"#2: radius\n"
		"#3: position x\n"
		"#4: position y\n"
		"#5: position z\n"
		"#6: velocity x\n"
		"#7: velocity y\n"
		"#8: velocity z\n"
		"#9: angular velocity x\n"
		"#10: angular velocity y\n"
		"#11: angular velocity z\n"
		"#12: viscosity contribution of lubrication\n"
		"#13: viscosity contributon of contact GU xz\n"
		"#14: viscosity contributon of brownian xz\n"
		"#15: angle (for 2D simulation only)\n";
		//
		fout_particle << fout_par_col_def << endl;
	}
	if (p.out_data_interaction) {
		string interaction_filename = "int_" + sys.simu_name + ".dat";
		fout_interaction.open(interaction_filename.c_str());
		outputDataHeader(fout_interaction);
		string fout_int_col_def =
		"#1: particle 1 label\n"
		"#2: particle 2 label\n"
		"#3: contact state (0 = no contact, 1 = frictionless contact, 2 = non-sliding frictional, 3 = sliding frictional)\n"
		"#4: normal vector, oriented from particle 1 to particle 2 x\n"
		"#5: normal vector, oriented from particle 1 to particle 2 y\n"
		"#6: normal vector, oriented from particle 1 to particle 2 z\n"
		"#7: dimensionless gap = s-2, s = 2r/(a1+a2)\n"
		"#8: normal part of the lubrication force\n"
		"#9: tangential part of the lubrication force x\n"
		"#10: tangential part of the lubrication force y\n"
		"#11: tangential part of the lubrication force z\n"
		"#12: norm of the normal part of the contact force\n"
		"#13: tangential part of the contact force, x\n"
		"#14: tangential part of the contact force, y\n"
		"#15: tangential part of the contact force, z\n"
		"#16: norm of the normal repulsive force\n"
		"#17: Viscosity contribution of contact xF\n";
		fout_interaction << fout_int_col_def << endl;
	}
}

void Simulation::setDefaultParameters()
{
	p.brownian_amplitude = 0;
	p.repulsion_amplitude = 0;
	p.cohesion_amplitude = 0;
	p.critical_load_amplitude = 0;
	p.Pe_switch = 5;
	p.dt = 1e-4;
	p.monolayer = false;
	p.rest_threshold = 1e-4;
	p.integration_method = 1;
	/*
	 * Stokes drag coeffient
	 */
	p.sd_coeff = 1;
	/*
	 * Lubrication model
	 * 0 no lubrication
	 * 1 1/xi lubrication (only squeeze mode)
	 * 2 log(1/xi) lubrication
	 * 3 ???
	 */
	p.lubrication_model = 2;
	/*
	 * 0 No friction
	 * 1 Linear friction law Ft < mu Fn
	 * 2 Threshold friction without repulsive force
	 * 3 Threshold friction without repulsion + mu inf
	 */
	p.friction_model = 1;
	p.time_end = 10;
	p.lub_max_gap = 0.5;
	/* This is cutoff distance (center-to-center) for interactions (repulsive force, magnetic force, etc.).
	 * If interaction_range is not indicated, this value will be set from lub_max_gap.
	 */
	p.interaction_range = -1;
	/*
	 * reduced_gap_min: gives reduced lubrication (maximum coeeffient).
	 *
	 */
	p.lub_reduce_parameter = 1e-3;
	/*
	 * contact_relaxation_factore
	 *
	 * This gives the coeffient of the resistance term for h < 0.
	 * - If the value is negative, the value of 1/lub_reduce_parameter is used.
	 *
	 */
	p.contact_relaxation_time = 1e-3;
	p.contact_relaxation_time_tan = 0;
	if (control_var == "stress") {
		p.unscaled_contactmodel = true;
		p.kn = 2000;
		p.kt = 1000;
		p.kr = 1000;
	} else {
		p.unscaled_contactmodel = false;
		p.kn = 10000;
		p.kt = 6000;
		p.kr = 6000;
	}
	p.auto_determine_knkt = false;
	p.overlap_target = 0.05;
	p.disp_tan_target = 0.05;
	p.memory_strain_avg = 0.01;
	p.memory_strain_k = 0.02;
	p.start_adjust = 0.2;
	p.min_kn = 1000;
	p.max_kn = 1000000;
	p.min_kt = 1000;
	p.max_kt = 1000000;
	p.repulsive_length = 0.05;
	p.repulsive_max_length = -1;
	p.mu_static = 1;
	p.mu_dynamic = -1;
	p.mu_rolling = 0;
	p.time_interval_output_data = 0.01;
	p.time_interval_output_config = 0.1;
	p.origin_zero_flow = true;
	p.out_data_particle = true;
	p.out_data_interaction = true;
	p.ft_max = 1;
	p.fixed_dt = false;
	/*
	 * Parameters for magnetic colloid simulation.
	 */
	p.magnetic_type = 0;
	p.magnetic_field_type = 0;
	p.magnetic_interaction_range = 20;
	p.timeinterval_update_magnetic_pair = 0.02;
}

void Simulation::importConfiguration()
{
	fstream file_import;
	file_import.open(filename_import_positions.c_str());
	if (!file_import) {
		cerr << " Position file '" << filename_import_positions << "' not found." <<endl;
		exit(1);
	}
	char buf;
	int n1, n2;
	double lx, ly, lz, vf1, vf2;
	initial_lees_edwards_disp.reset();
	getline(file_import, header_imported_configulation[0]);
	getline(file_import, header_imported_configulation[1]);
	stringstream ss(header_imported_configulation[1]);
	ss >> buf >> n1 >> n2 >> volume_or_area_fraction >> lx >> ly >> lz >> vf1 >> vf2 >> initial_lees_edwards_disp.x;
	double a;
	if(ss >> a) {
		initial_lees_edwards_disp.y = a;
	}
	else {
		initial_lees_edwards_disp.y = 0;
	}
	vector<vec3d> initial_position;
	vector <double> radius;
	if (sys.p.magnetic_type == 0) {
		double x_, y_, z_, a_;
		while (file_import >> x_ >> y_ >> z_ >> a_) {
			initial_position.push_back(vec3d(x_, y_, z_));
			radius.push_back(a_);
		}
		sys.setConfiguration(initial_position, radius, lx, ly, lz);
	} else {
		double x_, y_, z_, a_, mx_, my_, mz_, sus_;
		vector<vec3d> magnetic_moment;
		vector<double> magnetic_susceptibility;
		while (file_import >> x_ >> y_ >> z_ >> a_ >> mx_ >> my_ >> mz_ >> sus_) {
			initial_position.push_back(vec3d(x_, y_, z_));
			radius.push_back(a_);
			magnetic_moment.push_back(vec3d(mx_, my_, mz_));
			magnetic_susceptibility.push_back(sus_);
		}
		sys.setConfiguration(initial_position, radius, lx, ly, lz);
		sys.setMagneticConfiguration(magnetic_moment, magnetic_susceptibility);
	}
	file_import.close();
}


void Simulation::importConfigurationBinary()
{
	ifstream file_import;
	initial_lees_edwards_disp.reset();
	file_import.open(filename_import_positions.c_str(), ios::binary | ios::in);
	if (!file_import) {
		cerr << " Position file '" << filename_import_positions << "' not found." <<endl;
		exit(1);
	}
	int np;
	double lx, ly, lz;
	file_import.read((char*)&np, sizeof(int));
	file_import.read((char*)&volume_or_area_fraction, sizeof(double));
	file_import.read((char*)&lx, sizeof(double));
	file_import.read((char*)&ly, sizeof(double));
	file_import.read((char*)&lz, sizeof(double));
	file_import.read((char*)&initial_lees_edwards_disp.x, sizeof(double));
	file_import.read((char*)&initial_lees_edwards_disp.y, sizeof(double));
	double x_, y_, z_, r_;
	vector <vec3d> initial_position;
	vector <double> radius;
	for (int i=0; i<np; i++) {
		file_import.read((char*)&x_, sizeof(double));
		file_import.read((char*)&y_, sizeof(double));
		file_import.read((char*)&z_, sizeof(double));
		file_import.read((char*)&r_, sizeof(double));
		initial_position.push_back(vec3d(x_,y_,z_));
		radius.push_back(r_);
	}
	file_import.close();
	sys.setConfiguration(initial_position, radius, lx, ly, lz);
}

void Simulation::prepareSimulationName(bool binary_conf)
{
	ostringstream ss_simu_name;
	string::size_type pos_name_end = filename_import_positions.find_last_of(".");
	string::size_type param_name_end = filename_parameters.find_last_of(".");
	string::size_type pos_name_start;
	if (binary_conf) { // TO DO: improve name generation for binary input
		pos_name_start = filename_import_positions.find_last_of("/");
	} else {
		pos_name_start = filename_import_positions.find_last_of("/");
	}
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
	ss_simu_name << string_control_parameters.str();
	sys.simu_name = ss_simu_name.str();
	cout << "filename: " << sys.simu_name << endl;
}
