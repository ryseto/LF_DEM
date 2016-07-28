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

using namespace std;

void Simulation::contactForceParameter(string filename)
{
	/** @@@ DEPRECATED?
	 \brief Load a file containing spring constants and time steps as functions of the volume fraction.

	 Input file must be formatted as:
	 phi kn kt dt
	 */
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
		if (phi_ == volume_or_area_fraction) {
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
		error_str  << " Error: file " << filename.c_str() << " contains no data for vf = " << phi_ << endl;
		throw runtime_error(error_str.str());
	}
}

void Simulation::contactForceParameterBrownian(string filename)
{
	/** @@@ DEPRECATED?
	 \brief Load a file containing spring constants and time steps as functions of the volume fraction and Peclet number.

	 Input file must be formatted as:
	 phi peclet kn kt dt
	 */
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
	while (fin_knktdt >> phi_ >> peclet_ >> kn_ >> kt_ >> dt_) {
		if (phi_ == volume_or_area_fraction && peclet_ == dimensionless_numbers["hydro/thermal"]) {
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
		error_str  << " Error: file " << filename.c_str() << " contains no data for vf = " << volume_or_area_fraction << " and Pe = " << dimensionless_numbers["hydro/thermal"] << endl;
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
	shear_rate_expectation = shear_rate_;
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

void Simulation::resolveUnitSystem(string unit_force)
{
	/**
		\brief Check force units consistency, expresses all input forces in the unit "unit_force".

		In input, forces are given with suffixes.
		This function checks that we can make sense of these suffixes, and if so,
		it converts all the forces in the unit given as a parameter.

		It does it iteratively, ie:
		1. I know that unit_force = 1*unit_force
		2. I know the value of any force f1 that has been given in unit_force in input as
	      f1 = value*unit_force
		3. I can determine the value of any force f2 expressed as f2 = x*f1
	 4. I can then determine the value of any force f3 expressed as f3 = y*f2
		5. etc

		If there is any force undetermined at the end of this algorithm, the force unit system is inconsistent/incomplete.

		If the force unit system is consistent, this function determines the dimensionless numbers, i.e., the ratios F_A/F_B for any pair of force scales present in the system.

		\b Note: a priori, we could determine force A from force B if A is defined in B units or if B is defined in A units: for example in the step 3 of the above algorithm we could determine f2 knowing f1 if f1 is defined as f1=x^{-1}f2. However this is \b not implemented and will fail with the current implementation.
	 */

	// we keep the forces for which we can convert in unit_force units in a set.
	set <string> resolved_forces;
	resolved_forces.clear();

	// the unit_force has a value of 1*unit_force (says captain obvious)
	input_force_values[unit_force] = 1;
	input_force_units[unit_force] = unit_force;
	resolved_forces.insert(unit_force);

	// now resolve the other force units, iterativley
	auto previous_resolved = resolved_forces.size();
	do {
		previous_resolved = resolved_forces.size();
		for (const auto& f: input_force_units) {
			string force_type = f.first;
			string unit = f.second;
			if (resolved_forces.find(unit) != resolved_forces.end()) {  // the unit has a meaning
				input_force_values[force_type] *= input_force_values[unit]; // input_force_values[unit] is in unit_force units
				input_force_units[force_type] = unit_force;
				resolved_forces.insert(force_type);
			}
		}
	} while (previous_resolved < resolved_forces.size());

	// check we found everyone
	if (resolved_forces.size() < input_force_units.size()) {
		ostringstream error_str;
		for (const auto& f: input_force_units) {
			string force_type = f.first;
			string unit = f.second;
			if (resolved_forces.find(unit) == resolved_forces.end()) {
				error_str << "Error: force type \"" << force_type << "\" has an unknown unit \"" << unit << "\"" << endl;
			}
		}
		throw runtime_error(error_str.str());
	}
	// determine the dimensionless_numbers
	for (const auto& f1: input_force_units) {
		string force1_type = f1.first;
		for (const auto& f2: input_force_units) {
			string force2_type = f2.first;
			dimensionless_numbers[force2_type+'/'+force1_type] = input_force_values[force2_type]/input_force_values[force1_type];
			dimensionless_numbers[force1_type+'/'+force2_type] = 1/dimensionless_numbers[force2_type+'/'+force1_type];
		}
	}
}

void Simulation::catchForcesInStressUnits(const string &stress_unit)
{
	for (const auto& f: input_force_units) {
		if (f.second == "stress") {
			input_force_units[f.first] = stress_unit;
			input_force_values[f.first] *= abs(sys.target_stress);
		}
	}
}

void Simulation::convertInputForcesStressControlled(double dimensionlessnumber,
                                                    string stress_unit)
{
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
	string force_type = stress_unit; // our force defining the shear rate
	if (force_type == "hydro") {
		throw runtime_error(" Error: please give a stress in non-hydro units.");
		/*
		 Note:
		 Although it is in some cases possible to run under stress control with hydro units,
		 it is not always possible and as a consequence it is a bit dangerous to let the user do so.

		 With hydro units, the problem is that the target stress \f$\tilde{S}\f$ cannot take any possible value, as
		 \f$\tilde{S} = S/(\eta_0 \dot \gamma) = \eta/\eta_0\f$
		 --> It is limited to the available range of viscosity.
		 If you give a \f$\tilde{S}\f$ outside this range (for example \f$\tilde{S}=0.5\f$), you run into troubles.
		 */
	}
	if (force_type == "thermal") {
		throw runtime_error(" Error: stress controlled Brownian simulations are not yet implemented.");
	}
	sys.set_shear_rate(0);
	// we take as a unit scale the one given by the user with the stress
	// TODO: other choices may be better when several forces are used.
	internal_unit_scales = force_type;
	target_stress_input = dimensionlessnumber;
	sys.target_stress = target_stress_input/6/M_PI;

	// convert the forces expressed in "s" units, i.e. proportional to the stress
	catchForcesInStressUnits(stress_unit);

	// convert all other forces to internal_unit_scales
	resolveUnitSystem(internal_unit_scales);
}

// Command option -r indicates "rate controlled" simulation.
// -r [val]r  ---> val = F_H0/F_R0 = shear_rate/shear_rate_R0
// -r [val]b  ---> val = F_H0/F_B0 = shear_rate/shear_rate_B0
void Simulation::convertInputForcesRateControlled(double dimensionlessnumber,
												  string input_scale)
{
	/**
	 \brief Choose units for the simulation and convert the forces to this unit (rate controlled case).

	 The strategy is the following:
	 1. Convert all the forces in hydro force unit, ie the value for force f1 is f1/F_H
	 2. Decide the unit force F_W for the simulation
	 3. Convert all the forces to this unit (by multiplying every force by F_H/F_W)
	 */
	string force_type = input_scale; // the force defining the shear rate
	if (input_force_values[force_type] > 0) { // if the force defining the shear rate is redefined in the parameter file, throw an error
		ostringstream error_str;
		error_str  << "Error: redefinition of the rate (given both in the command line and in the parameter file with \"" << force_type << "\" force)" << endl;
		throw runtime_error(error_str.str());
	}
	/* Switch this force in hydro units.
	 The dimensionlessnumber in input is the shear rate in "force_type" units,
	 i.e.  is the ratio between the hydrodynamic force and the "force_type" force.
	 So in hydrodynamic force units, the value of the "force_type" force is 1/dimensionlessnumber.
	 */
	if (dimensionlessnumber == 0) {
		throw runtime_error("Vanishing rate not handled... yet! ");
	}
	input_force_values[force_type] = 1/dimensionlessnumber;
	input_force_units[force_type] = "hydro";
	// convert all other forces to hydro
	resolveUnitSystem("hydro");
	// chose simulation unit
	// the chosen unit is called internal_unit_scales
	setUnitScaleRateControlled();
	// convert from hydro scale to chosen scale
	convertForceValues(internal_unit_scales);
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
	/**
	 \brief Convert all the input forces to the unit given in parameter.

		The input forces velues can be expressed in any units,
		but these units must be known before calling this function
		(e.g. by calling Simulation::resolveUnitSystem(string unit_force) beforehand).
	 */
	for (auto& f: input_force_units) {
		string force_type = f.first;
		string old_unit = f.second;
		if (old_unit != new_unit) {
			input_force_values[force_type] *= dimensionless_numbers[old_unit+'/'+new_unit];
			f.second = new_unit;
		}
	}
}

void Simulation::setUnitScaleRateControlled()
{
	/**
	 \brief Determine the best internal force scale to run the simulation (rate controlled case).

		If the system is non-Brownian, the hydrodynamic force unit is taken (\b note: this will change in the future). If the system is Brownian, the Brownian force unit is selected at low Peclet (i.e., Peclet numbers smaller that ParameterSet::Pe_switch) and the hydrodynamic force unit is selected at high Peclet.
	 */
	bool is_brownian;
	if (dimensionless_numbers.find("hydro/thermal") != dimensionless_numbers.end()) {
		is_brownian = true;
	} else {
		is_brownian = false;
	}
	if (is_brownian) {
		if (dimensionless_numbers["hydro/thermal"] > p.Pe_switch && !sys.zero_shear) { // hydro units
			internal_unit_scales = "hydro";
		} else { // low Peclet mode
			internal_unit_scales = "thermal";
			setLowPeclet();
		}
	} else {
		internal_unit_scales = "hydro";
	}
	sys.set_shear_rate(dimensionless_numbers["hydro/"+internal_unit_scales]);
}

void Simulation::exportForceAmplitudes()
{
	/**
	 \brief Copy the input_force_values in the ForceAmplitude struct of the System class
	 */
	string indent = "  Simulation::\t";
	cout << indent+"Forces used:" << endl;
	indent += "\t";
	bool is_repulsive = input_force_values.find("repulsive") != input_force_values.end();
	if (is_repulsive) {
		sys.repulsiveforce = true;
		sys.amplitudes.repulsion = input_force_values["repulsive"];
		cout << indent+"Repulsive force (in \"" << input_force_units["repulsive"] << "\" units): " << sys.amplitudes.repulsion << endl;
	}
	bool is_critical_load = input_force_values.find("critical_load") != input_force_values.end();
	if (is_critical_load) {
		sys.critical_load = true;
		sys.amplitudes.critical_normal_force = input_force_values["critical_load"];
		cout << indent+"Critical Load (in \"" << input_force_units["critical_load"] << "\" units): " << sys.amplitudes.critical_normal_force << endl;
	}
	bool is_cohesive = input_force_values.find("cohesive") != input_force_values.end();
	if (is_cohesive) {
		sys.cohesion = true;
		sys.amplitudes.cohesion = input_force_values["cohesive"];
		cout << indent+"Cohesion (in \"" << input_force_units["cohesive"] << "\" units): " << sys.amplitudes.cohesion << endl;
	}
	bool is_ft_max = input_force_values.find("ft") != input_force_values.end();
	if (is_ft_max) {
		sys.amplitudes.ft_max = input_force_values["ft"];
		cout << indent+"Max tangential load (in \"" << input_force_units["ft"] << "\" units): " << sys.amplitudes.ft_max << endl;
	}

	bool is_brownian = input_force_values.find("thermal") != input_force_values.end();
	if (is_brownian) {
		sys.brownian = true;
		p.brownian_amplitude = input_force_values["thermal"];
		sys.amplitudes.sqrt_temperature = 1/sqrt(dimensionless_numbers[internal_unit_scales+"/thermal"]);
		cout << indent+"Brownian force (in \"" << input_force_units["thermal"] << "\" units): " << dimensionless_numbers["thermal/"+internal_unit_scales] << endl;
	}
	p.kn = input_force_values["normal_stiffness"];
	cout << indent+"Normal contact stiffness (in \"" << input_force_units["normal_stiffness"] << "\" units): " << dimensionless_numbers["normal_stiffness/"+internal_unit_scales] << endl;
	p.kt = input_force_values["tan_stiffness"];
	cout << indent+"Sliding contact stiffness (in \"" << input_force_units["tan_stiffness"] << "\" units): " << dimensionless_numbers["tan_stiffness/"+internal_unit_scales] << endl;
	p.kr = input_force_values["roll_stiffness"];
	cout << indent+"Rolling contact stiffness (in \"" << input_force_units["roll_stiffness"] << "\" units): " << dimensionless_numbers["roll_stiffness/"+internal_unit_scales] << endl;
}

void Simulation::convertInputValues(string new_unit)
{
	/**
	 \brief Convert all the non-force input parameters given with units to the unit given in parameter.

		\b Note Forces are treated with Simulation::convertForceValues(string new_unit) .
	 */
	for (auto& inval: input_values) {
		auto &inv = inval.second;
		string name = inval.first;
		string old_unit = inv.unit;
		if (old_unit == "hydro" && control_var == "stress") {
			throw runtime_error ("Simulation:: "+name+" cannot be given in hydro units for stress controlled simulations (non-constant hydro unit of time)");
		}
		if (old_unit != new_unit && old_unit!="strain") {
			if (old_unit != "hydro" && input_force_values.find(old_unit) == input_force_values.end()) {
				ostringstream error_str;
				error_str  << " Error: trying to convert " << name << " from an unknown unit \"" << inv.unit 	<< "\"" << endl;
				throw runtime_error(error_str.str());
			}
			if (inv.type == "time") {
				*(inv.value) *= dimensionless_numbers[new_unit+'/'+old_unit];
				inv.unit = new_unit;
			} else if (inv.type == "force") {
				*(inv.value) *= dimensionless_numbers[old_unit+'/'+new_unit];
				inv.unit = new_unit;
			} else {
				throw runtime_error ("Simulation:: do not know how to convert value of type "+inv.type);
			}
		}

		string indent = "  Simulation::\t";
		cout << indent << name << " (in \"" << inv.unit << "\" units): " << *(inv.value) << endl;
	}
}

void Simulation::setupNonDimensionalization(double dimensionlessnumber,
											string input_scale){
	/**
	 \brief Non-dimensionalize the simulation.

		This function determines the most appropriate unit scales to use in the System class depending on the input parameters (Brownian/non-Brownian, shear rate, stress/rate controlled), and converts all the input values in these units.
	 */
	input_scale = unit_longname[input_scale];
	if (control_var == "rate") {
		input_rate = dimensionlessnumber; // @@@ Renaming is required?
	}
	if (control_var == "rate") {
		convertInputForcesRateControlled(dimensionlessnumber, input_scale);
	} else if (control_var == "stress") {
		convertInputForcesStressControlled(dimensionlessnumber, input_scale);
	} else {
		ostringstream error_str;
		error_str  << " Error: unknown control variable \"" << control_var 	<< "\"" << endl;
		throw runtime_error(error_str.str());
	}
	exportForceAmplitudes();
	string indent = "  Simulation::\t";
	cout << indent << "internal_unit_scales = " << internal_unit_scales << endl;
	sys.ratio_unit_time = &dimensionless_numbers[input_scale+"/"+internal_unit_scales];
	convertInputValues(internal_unit_scales);
	output_unit_scales = input_scale;
}

void Simulation::assertParameterCompatibility()
{
	// test for incompatibilities
	if (sys.brownian == true) {
		if ((p.lubrication_model != "none" || p.mu_static > 0)
			&& p.integration_method != 1) {
			ostringstream error_str;
			error_str << "Brownian simulation needs to use the Predictor-Corrector method." << endl;
			error_str << "Modify the parameter file." << endl;
			throw runtime_error(error_str.str());
		}
	}
	if (control_var == "stress") {
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

void Simulation::tagStrainParameters()
{
	/**
	\brief Tag some time parameters like time_end given strain units (suffix h)

		In stress controlled simulations, hydro time units are ill-defined
		as the rate is changing in time. If times are given with these units, LF_DEM will
		throw a runtime_error if run under stress controlled conditions.
		However, there are some parameters for which giving strain units make sense,
		like time_end, etc. Here we tag those special parameters as having explicit
		"strain" units, which will later go through unaffected.

		Note that the System class must be aware that those parameters are strains
		to be compared with System::shear_strain() (as opposed to times to be compared
		with System::time()). The mechanism to declare those parameters to System is
		implemented in Simulation::resolveTimeOrStrainParameters().
	 */
	for (auto& name: {"time_interval_output_data", "time_interval_output_config", "time_end", "initial_log_time"}) {
		auto &inv =  input_values[name];
		if (inv.unit == "hydro") {
			inv.unit = "strain";
		}
	}
}

void Simulation::resolveTimeOrStrainParameters()
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
	if (input_values["time_end"].unit == "strain") {
		time_end = -1;
		strain_end = p.time_end;
	} else {
		time_end = p.time_end;
	}
	if (p.log_time_interval) {
		if (input_values["time_end"].unit != input_values["initial_log_time"].unit &&
			(input_values["time_end"].unit == "strain" || input_values["initial_log_time"].unit == "strain")) {
			throw runtime_error(" If one of time_end or initial_log_time is a strain (\"h\" unit), than both must be.\n");
		}
	}
}

void Simulation::setupSimulation(string in_args,
								 vector<string>& input_files,
								 bool binary_conf,
								 double dimensionlessnumber,
								 string input_scale,
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

	if (filename_parameters.find("init_relax", 0) != string::npos) {
		cout << "init_relax" << endl;
		sys.zero_shear = true;
	} else {
		sys.zero_shear = false;
	}
	if (sys.test_simulation > 0 && (sys.test_simulation < 20 || sys.test_simulation > 30) ) {
		sys.zero_shear = true;
		sys.mobile_fixed = true;
	}
	setDefaultParameters(input_scale);
	readParameterFile(filename_parameters);
	tagStrainParameters();
	setupNonDimensionalization(dimensionlessnumber, input_scale);

	assertParameterCompatibility();

	if (input_files[2] != "not_given") {
		if (sys.brownian && !p.auto_determine_knkt) {
			contactForceParameterBrownian(input_files[2]);
		} else {
			contactForceParameter(input_files[2]);
		}
	}
	if (input_files[3] != "not_given") {
		throw runtime_error("pre-simulation data deprecated?");
	}
	resolveTimeOrStrainParameters();

	bool is2d;
	pair<int,int> np_np_fixed;
	if (binary_conf) {
		np_np_fixed = get_np_Binary(filename_import_positions);
		is2d = isTwoDimensionBinary(filename_import_positions);
	} else {
		np_np_fixed = get_np(filename_import_positions);
		is2d = isTwoDimension(filename_import_positions);
	}
	sys.set_np(np_np_fixed.first);
	if (np_np_fixed.second > 0){
		p.np_fixed = np_np_fixed.second;
	}

	sys.setupSystemPreConfiguration(control_var, is2d);
	if (binary_conf) {
		importConfigurationBinary(filename_import_positions);
	} else {
		importConfiguration(filename_import_positions);
	}
	sys.setupSystemPostConfiguration();
	p_initial = p;
	prepareSimulationName(binary_conf, filename_import_positions, filename_parameters,
						  simu_identifier, dimensionlessnumber, input_scale);
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
	if (keyword == "lubrication_model") {
		p.lubrication_model = value;
	} else if (keyword == "friction_model") {
		if (p.friction_model == 2) {
			cerr << "!!Neglected friction_model in parameter file!!" << endl;
		} else {
			p.friction_model = atoi(value.c_str());
		}
	} else if (keyword == "repulsion_amplitude") {
		catchSuffixedForce("repulsive", value);
	} else if (keyword == "cohesion_amplitude") {
		catchSuffixedForce("cohesive", value);
	} else if (keyword == "brownian_amplitude") {
		catchSuffixedForce("thermal", value);
	} else if (keyword == "critical_load_amplitude") {
		catchSuffixedForce("critical_load", value);
	} else if (keyword == "monolayer") {
		p.monolayer = str2bool(value);
	} else if (keyword == "repulsiveforce_length") {
		p.repulsive_length = atof(value.c_str());
	} else if (keyword == "repulsive_max_length") {
		p.repulsive_max_length = atof(value.c_str());
	} else if (keyword == "lub_reduce_parameter") {
		p.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time") {
		catchSuffixedValue("time", keyword, value, &p.contact_relaxation_time);
	} else if (keyword == "contact_relaxation_time_tan"){
		catchSuffixedValue("time", keyword, value, &p.contact_relaxation_time_tan);
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
		catchSuffixedForce("normal_stiffness", value);
	} else if (keyword == "kt") {
		catchSuffixedForce("tan_stiffness", value);
	} else if (keyword == "kr") {
		catchSuffixedForce("roll_stiffness", value);
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
	} else if (keyword == "log_time_interval") {
		p.log_time_interval = str2bool(value);
	} else if (keyword == "initial_log_time") {
		catchSuffixedValue("time", keyword, value, &p.initial_log_time);
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
	} else if (keyword == "min_kn") {
		catchSuffixedValue("force", keyword, value, &p.min_kn);
	} else if (keyword == "max_kn") {
		catchSuffixedValue("force", keyword, value, &p.max_kn);
	} else if (keyword == "min_kt") {
		catchSuffixedValue("force", keyword, value, &p.min_kt);
	} else if (keyword == "max_kt") {
		catchSuffixedValue("force", keyword, value, &p.max_kt);
	} else if (keyword == "rest_threshold") {
		p.rest_threshold = atof(value.c_str());
	} else if (keyword == "ft_max") {
		catchSuffixedForce("ft", value);
	} else if (keyword == "fixed_dt") {
		p.fixed_dt = str2bool(value);
	} else if (keyword == "cross_shear") {
		p.cross_shear = str2bool(value);
	} else if (keyword == "theta_shear") {
		p.theta_shear = atof(value.c_str());
		p.theta_shear *= M_PI/180;  // convert in radians
	} else if (keyword == "strain_reversal") {
		p.strain_reversal = atof(value.c_str());
	} else if (keyword == "event_handler") {
		p.event_handler = value;
		p.event_handler.erase(remove(p.event_handler.begin(), p.event_handler.end(), '\"' ), p.event_handler.end());
	} else if (keyword == "time_init_relax") {
		catchSuffixedValue("time", keyword, value, &p.time_init_relax);
	} else if (keyword == "out_particle_stress") {
		p.out_particle_stress = value;
		p.out_particle_stress.erase(remove(p.out_particle_stress.begin(), p.out_particle_stress.end(), '\"' ), p.out_particle_stress.end());
	} else if (keyword == "out_binary_conf") {
		p.out_binary_conf = str2bool(value);
	} else if (keyword == "np_fixed") {
		p.np_fixed = atoi(value.c_str());
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

void Simulation::setDefaultParameters(string input_scale)
{
	/**
	 \brief Set default values for ParameterSet parameters.
	 */
	p.brownian_amplitude = 0;
	p.repulsion_amplitude = 0;
	p.cohesion_amplitude = 0;
	p.critical_load_amplitude = 0;
	p.Pe_switch = 5;
	p.dt = 1e-4;
	p.disp_max = 1e-3;
	p.monolayer = false;
	p.rest_threshold = 1e-4;
	p.integration_method = 1;
	p.np_fixed = 0;
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
	p.lubrication_model = "tangential";
	/*
	 * 0 No friction
	 * 1 Linear friction law Ft < mu Fn
	 * 2 Threshold friction without repulsive force
	 * 3 Threshold friction without repulsion + mu inf
	 */
	p.friction_model = 1;
	catchSuffixedValue("time", "time_end", "10h", &p.time_end);
	p.time_init_relax = 0;
	p.lub_max_gap = 0.5;
	/* This is cutoff distance (center-to-center) for interactions (repulsive force, etc.).
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
	catchSuffixedValue("time", "contact_relaxation_time", "1e-3"+input_scale, &p.contact_relaxation_time);
	catchSuffixedValue("time", "contact_relaxation_time_tan", "-1"+input_scale, &p.contact_relaxation_time_tan);
	catchSuffixedForce("kn", "2000"+input_scale);
	catchSuffixedForce("kt", "0.5kn");
	catchSuffixedForce("kr", "0.5kn");
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
	catchSuffixedValue("time", "time_interval_output_data", "1e-2h", &p.time_interval_output_data);
	catchSuffixedValue("time", "time_interval_output_config", "1e-1h", &p.time_interval_output_config);
	p.log_time_interval = false;
	catchSuffixedValue("time", "initial_log_time", "1e-4h", &p.initial_log_time);
	p.nb_output_data_log_time = 100;
	p.nb_output_config_log_time = 100;
	p.origin_zero_flow = true;
	p.out_data_particle = true;
	p.out_data_interaction = true;
	p.out_particle_stress = "";
	p.out_binary_conf = false;
	p.out_data_vel_components = false;
	p.ft_max = 1;
	p.fixed_dt = false;
	p.cross_shear = false;
	p.theta_shear = 0;
	p.event_handler = "";
}

inline string columnDefinition(int &cnb, const string &type, const string &name)
{
	stringstream defs;
	if (type == "vec3d") {
		array<string, 3> xyz = {"x", "y", "z"};
		for (auto &u : xyz) {
			stringstream col_def_complement;
			defs << "#" << cnb << ": "<< name << " " << u << "\n";
			cnb ++;
		}
	} else if (type=="scalar") {
		defs << "#" << cnb << ": "<< name << "\n";
	} else {
		throw runtime_error(" unknown type for column def\n");
	}
	return defs.str();
}

void Simulation::openOutputFiles()
{
	/**
	 \brief Set up the output files

		This function determines a simulation name from the parameters, opens the output files with the corresponding name and prints their header.
	 */

	stringstream data_header;
	createDataHeader(data_header);
	outdata.setFile("data_"+sys.simu_name+".dat", data_header.str(), force_to_run);
	outdata_st.setFile("st_"+sys.simu_name+".dat", data_header.str(), force_to_run);
	if (!p.out_particle_stress.empty()) {
		outdata_pst.setFile("pst_"+sys.simu_name+".dat", data_header.str(), force_to_run);
	}
	string time_filename = "t_"+sys.simu_name+".dat";
	fout_time.open(time_filename.c_str());
	string input_filename = "input_"+sys.simu_name+".dat";
	fout_input.open(input_filename.c_str());
	if (p.out_data_particle) {
		outdata_par.setFile("par_"+sys.simu_name+".dat", data_header.str(), force_to_run);
	}
	if (p.out_data_interaction) {
		outdata_int.setFile("int_"+sys.simu_name+".dat", data_header.str(), force_to_run);
	}
}

map<string,string> Simulation::getConfMetaData(const string &line1, const string &line2)
{
	vector<string> l1_split = splitString(line1);
	vector<string> l2_split = splitString(line2);
	if (l1_split.size() != l2_split.size()) {
		throw runtime_error("Simulation:: Ill-formed header in the configuration file.\n");
	}
	map<string,string> meta_data;
	for (unsigned int i=1; i<l1_split.size(); i++) {
		meta_data[l1_split[i]] = l2_split[i];
	}
	return meta_data;
}

string Simulation::getMetaParameter(map<string,string> &meta_data, string &key)
{
	if (meta_data.find(key) != meta_data.end()) {
		return meta_data[key];
	} else {
		ostringstream error_str;
		error_str  << " Simulation:: parameter '" << key << "' not found in the header of the configuration file." <<endl;
		throw runtime_error(error_str.str());
	}
}

string Simulation::getMetaParameter(map<string,string> &meta_data, string &key, const string &default_val)
{
	if (meta_data.find(key) != meta_data.end()) {
		return meta_data[key];
	} else {
		return default_val;
	}
}

std::pair<int,int> Simulation::get_np(const string& filename_import_positions)
{
	/**
	 \brief Read np from a text file input configuration.
	 */
	fstream file_import;
	file_import.open(filename_import_positions.c_str());
	if (!file_import) {
		ostringstream error_str;
		error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
		throw runtime_error(error_str.str());
	}
	int np, np_fixed;
	getline(file_import, header_imported_configulation[0]);
	getline(file_import, header_imported_configulation[1]);
	map<string,string> meta_data = getConfMetaData(header_imported_configulation[0], header_imported_configulation[1]);
	string key, def;
	key = "np_fixed";
	def = "0";
	np_fixed = atoi(getMetaParameter(meta_data, key, def).c_str());

	string line;
	np = 0;
	while (getline(file_import, line)) {
		np++;
	}
	file_import.close();

	return make_pair(np, np_fixed);
}

bool Simulation::isTwoDimension(const string& filename_import_positions)
{
	fstream file_import;
	file_import.open(filename_import_positions.c_str());
	if (!file_import) {
		ostringstream error_str;
		error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
		throw runtime_error(error_str.str());
	}
	double ly;
	getline(file_import, header_imported_configulation[0]);
	getline(file_import, header_imported_configulation[1]);
	map<string,string> meta_data = getConfMetaData(header_imported_configulation[0], header_imported_configulation[1]);
	string key;
	key = "ly";
	ly = atof(getMetaParameter(meta_data, key).c_str());
	file_import.close();
	return ly==0;
}

bool Simulation::isTwoDimensionBinary(const string& filename_import_positions)
{
	ifstream file_import;
	file_import.open(filename_import_positions.c_str(), ios::binary | ios::in);
	if (!file_import) {
		ostringstream error_str;
		error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
		throw runtime_error(error_str.str());
	}
	double ddumb, ly;
	int np = 0;
	file_import.read((char*)&np, sizeof(int));
	int binary_format_version;
	if (np == -1) {
		int np_fixed = 0;
		file_import.read((char*)&binary_format_version, sizeof(int));
		file_import.read((char*)&np, sizeof(int));
		if (binary_format_version == 3) {
			file_import.read((char*)&np_fixed, sizeof(int));
		}
	}
	file_import.read((char*)&ddumb, sizeof(double)); // vf
	file_import.read((char*)&ddumb, sizeof(double)); // lx
	file_import.read((char*)&ly, sizeof(double));
	return ly==0;
}

std::pair<int,int> Simulation::get_np_Binary(const string& filename_import_positions)
{
	/**
	 \brief Read np from binary file input configuration.
	 */
	ifstream file_import;
	file_import.open(filename_import_positions.c_str(), ios::binary | ios::in);
	if (!file_import) {
		ostringstream error_str;
		error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
		throw runtime_error(error_str.str());
	}
	int np = 0;
	int np_fixed = 0;
	file_import.read((char*)&np, sizeof(int));
	if (np == -1) {
		int binary_format_version;
		file_import.read((char*)&binary_format_version, sizeof(int));
		file_import.read((char*)&np, sizeof(int));
		if (binary_format_version == 3) {
			file_import.read((char*)&np_fixed, sizeof(int));
		}
	}
	file_import.close();
	return make_pair(np, np_fixed);
}

void Simulation::setMetadata(fstream &file_import)
{
	double lx, ly, lz;
	vec3d initial_lees_edwards_disp;
	initial_lees_edwards_disp.reset();
	getline(file_import, header_imported_configulation[0]);
	getline(file_import, header_imported_configulation[1]);

	map<string,string> meta_data = getConfMetaData(header_imported_configulation[0], header_imported_configulation[1]);
	string key, def;
	// key = "np_fixed";
	// def = "0";
	// sys.p.np_fixed = atoi(getMetaParameter(meta_data, key, def).c_str());
	key = "lx";
	lx = atof(getMetaParameter(meta_data, key).c_str());
	key = "ly";
	ly = atof(getMetaParameter(meta_data, key).c_str());
	key = "lz";
	lz = atof(getMetaParameter(meta_data, key).c_str());
	sys.setBoxSize(lx, ly, lz);
	key = "vf";
	volume_or_area_fraction = atof(getMetaParameter(meta_data, key).c_str());
	key = "dispx";
	def = "0";
	initial_lees_edwards_disp.x = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "dispy";
	def = "0";
	initial_lees_edwards_disp.y = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "np_wall1";
	def = "-1";
	sys.np_wall1 = atoi(getMetaParameter(meta_data, key, def).c_str());
	key = "np_wall2";
	def = "-1";
	sys.np_wall2 = atoi(getMetaParameter(meta_data, key, def).c_str());
	key = "radius_in";
	def = "0";
	sys.radius_in = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "radius_out";
	def = "0";
	sys.radius_out = atof(getMetaParameter(meta_data, key, def).c_str());
	// @@ This is temporally used for wtestA and wtestB
	key = "z_bot";
	def = "-1";
	sys.z_bot = atof(getMetaParameter(meta_data, key, def).c_str());
	key = "z_top";
	def = "-1";
	sys.z_top = atof(getMetaParameter(meta_data, key, def).c_str());

	if (sys.np_wall1 != -1) {
		sys.p.np_fixed = sys.np_wall1+sys.np_wall2;
	}
	sys.shear_disp = initial_lees_edwards_disp;
}

void Simulation::readPositions(fstream &file_import)
{
	/**
		\brief Import a regular text file configuration.

		File format:
		# header

		x y z radius
		...
	*/
	double x_, y_, z_, a_;
	vector<vec3d> initial_position;
	vector<double> radius;
	while (file_import >> x_ >> y_ >> z_ >> a_) {
		initial_position.push_back(vec3d(x_, y_, z_));
		radius.push_back(a_);
	}
	sys.setConfiguration(initial_position, radius);
}

void Simulation::readPositionsImposedVelocity(fstream &file_import)
{
	/**
		\brief Import a text file configuration with imposed velocity particles.

		File format:
		# header

		x y z radius
		...
		x y z radius vx vy vz
		...
	*/
	// http://stackoverflow.com/questions/743191/how-to-parse-lines-with-differing-number-of-fields-in-c
	double x_, y_, z_, a_, vx_, vy_, vz_;
	vector<vec3d> initial_position;
	vector<vec3d> fixed_velocities;
	vector <double> radius;
	string line;
	while(getline(file_import, line)) {
		istringstream is;
		is.str(line);
		if (!(is >> x_ >> y_ >> z_ >> a_ >> vx_ >> vy_ >> vz_)) {
			is.str(line);
			is >> x_ >> y_ >> z_ >> a_;
			initial_position.push_back(vec3d(x_, y_, z_));
			radius.push_back(a_);
		} else {
			initial_position.push_back(vec3d(x_, y_, z_));
			radius.push_back(a_);
			fixed_velocities.push_back(vec3d(vx_, vy_, vz_));
		}
	}
	if (sys.p.np_fixed != (int)fixed_velocities.size()) {
		throw runtime_error(" Simulation:: ill-formed input configuration, np_fixed != fixed_velocities.size()");
	}
	sys.setConfiguration(initial_position, radius);
	sys.setFixedVelocities(fixed_velocities);
}

void Simulation::importConfiguration(const string& filename_import_positions)
{
	/**
	 \brief Read a text file input configuration.
	 */
	fstream file_import;
	file_import.open(filename_import_positions.c_str());
	if (!file_import) {
		ostringstream error_str;
		error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
		throw runtime_error(error_str.str());
	}

	setMetadata(file_import);
	if (sys.test_simulation != 31) {
		readPositions(file_import);
	} else {
		readPositionsImposedVelocity(file_import);
	}

	file_import.close();
}

void Simulation::importConfigurationBinary(const string& filename_import_positions)
{
	/**
	 \brief Read a binary file input configuration.

	 Depending on the type of simulation, we store the data differently, defined by
	 binary format version numbers:
	 v1 : no version number field
	      metadata : np, vf, lx, ly, lz, disp_x, disp_y
	      particle data : [x, y, z, radius]*np
	      contact data : nb_interactions,
	      [p0, p1, dtx, dty, dtz, drx, dry, drz]*nb_interactions
				(with p0, p1 unsigned short)

	 v2 : no version number field
	      metadata : same as v1
	      particle data : as v1
	      contact data : as v1, except that p0 and p1 are unsigned int

	 v3 : (fixed wall particle case)
	      version nb: -1, 3  (-1 to distinguish from v1:np or v2:np)
	      metadata : np, np_fixed, vf, lx, ly, lz, disp_x, disp_y
	      particle data : [x, y, z, radius]*np, [vx, vy, vz]*np_fixed
				contact data : as v2
	 */


	// first positions
	vec3d initial_lees_edwards_disp;
	initial_lees_edwards_disp.reset();
	ifstream file_import;
	file_import.open(filename_import_positions.c_str(), ios::binary | ios::in);
	if (!file_import) {
		ostringstream error_str;
		error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
		throw runtime_error(error_str.str());
	}

	int binary_format_version;
	int np;
	int np_fixed;
	double lx, ly, lz;
	file_import.read((char*)&np, sizeof(int));
	if (np == -1) {
		file_import.read((char*)&binary_format_version, sizeof(int));
		file_import.read((char*)&np, sizeof(int));
	} else {
		binary_format_version = 2; // may also be 1, but will be determined later
	}
	if (binary_format_version == 3) {
		file_import.read((char*)&np_fixed, sizeof(int));
	}
	file_import.read((char*)&volume_or_area_fraction, sizeof(double));
	file_import.read((char*)&lx, sizeof(double));
	file_import.read((char*)&ly, sizeof(double));
	file_import.read((char*)&lz, sizeof(double));
	file_import.read((char*)&initial_lees_edwards_disp.x, sizeof(double));
	file_import.read((char*)&initial_lees_edwards_disp.y, sizeof(double));
	sys.shear_disp = initial_lees_edwards_disp;

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

	sys.setBoxSize(lx,ly,lz);
	sys.setConfiguration(initial_position, radius);

	if (binary_format_version == 3) {
		vector<vec3d> fixed_velocities;
		double vx_, vy_, vz_;
		for (int i=0; i<np_fixed; i++) {
			file_import.read((char*)&vx_, sizeof(double));
			file_import.read((char*)&vy_, sizeof(double));
			file_import.read((char*)&vz_, sizeof(double));
			fixed_velocities.push_back(vec3d(vx_, vy_, vz_));
		}
		sys.setFixedVelocities(fixed_velocities);
	}

	// now contacts
	int ncont;
	double dt_x, dt_y, dt_z, dr_x, dr_y, dr_z;
	vector <struct contact_state> cont_states;
	file_import.read((char*)&ncont, sizeof(unsigned int));
	std::iostream::pos_type file_pos = file_import.tellg();
	for (int i=0; i<ncont; i++) {
		unsigned int p0, p1;
		file_import.read((char*)&p0, sizeof(unsigned int));
		// hacky thing to guess if this is an old format with particle numbers as unsigned short
		if((int)p0>sys.get_np()){
			binary_format_version = 1;
			file_import.seekg(file_pos);
			break;
		}
		file_import.read((char*)&p1, sizeof(unsigned int));
		file_import.read((char*)&dt_x, sizeof(double));
		file_import.read((char*)&dt_y, sizeof(double));
		file_import.read((char*)&dt_z, sizeof(double));
		file_import.read((char*)&dr_x, sizeof(double));
		file_import.read((char*)&dr_y, sizeof(double));
		file_import.read((char*)&dr_z, sizeof(double));
		struct contact_state cs;
		cs.p0 = p0;
		cs.p1 = p1;
		cs.disp_tan = vec3d(dt_x, dt_y, dt_z);
		cs.disp_rolling = vec3d(dr_x, dr_y, dr_z);
		cont_states.push_back(cs);
	}
	if (binary_format_version == 1) {
		for (int i=0; i<ncont; i++) {
			unsigned short p0, p1;

			file_import.read((char*)&p0, sizeof(unsigned short));
			file_import.read((char*)&p1, sizeof(unsigned short));
			file_import.read((char*)&dt_x, sizeof(double));
			file_import.read((char*)&dt_y, sizeof(double));
			file_import.read((char*)&dt_z, sizeof(double));
			file_import.read((char*)&dr_x, sizeof(double));
			file_import.read((char*)&dr_y, sizeof(double));
			file_import.read((char*)&dr_z, sizeof(double));
			struct contact_state cs;
			cs.p0 = (int)p0;
			cs.p1 = (int)p1;
			cs.disp_tan = vec3d(dt_x, dt_y, dt_z);
			cs.disp_rolling = vec3d(dr_x, dr_y, dr_z);
			cont_states.push_back(cs);
		}
	}
	file_import.close();
	sys.setContacts(cont_states);
}

void Simulation::prepareSimulationName(bool binary_conf,
									   const string& filename_import_positions,
									   const string& filename_parameters,
									   const string& simu_identifier,
									   double dimensionlessnumber,
									   const string& input_scale)
{
	/**
	 \brief Determine simulation name.
	 */
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
	ostringstream string_control_parameters;
	if (long_file_name) {
		for (const auto& f: input_force_units) {
			if (f.first.find("stiffness") == std::string::npos) {
				string_control_parameters << "_" << f.first << input_force_values[f.first] << f.second;
			}
		}
	}
	string_control_parameters << "_" << control_var << dimensionlessnumber << input_scale;
	ss_simu_name << string_control_parameters.str();
	if (simu_identifier != "") {
		ss_simu_name << "_";
		ss_simu_name << simu_identifier;
	}
	sys.simu_name = ss_simu_name.str();
	string indent = "  Simulation::\t";
	cout << indent << "filename: " << sys.simu_name << endl;
}
