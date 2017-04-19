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

void Simulation::buildFullSetOfForceRatios(){
	// determine the complete set of force_ratios from the dimensionless forces
	std::map <std::string, DimensionalValue> input_forces;
	for (const auto& x: input_values) {
		if (x.second.type == "force") {
			input_forces[x.first] = x.second;
		}
	}
	for (const auto& f1: input_forces) {
		string force1_type = f1.first;
		double f1_val = *(f1.second.value);
		for (const auto& f2: input_forces) {
			string force2_type = f2.first;
			double f2_val = *(f2.second.value);
			force_ratios[force2_type+'/'+force1_type] = f2_val/f1_val;
			force_ratios[force1_type+'/'+force2_type] = 1/force_ratios[force2_type+'/'+force1_type];
		}
	}
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

	// the unit_force has a value of 1*unit_force (says captain obvious)
	DimensionalValue inv;
	inv.type = "force";
	inv.value = force_value_ptr[unit_force];
	inv.unit = unit_force;
	input_values[unit_force] = inv;
	*(input_values[unit_force].value) = 1;

	for (const auto& x: input_values) {
		if (x.second.type == "force") {
			string force_name = x.first;
			string unit = x.second.unit;
			force_ratios[force_name+'/'+unit] = *(x.second.value);
		}
	}

	// now resolve the other force units, iterativley

	bool unsolved_remaining;
	bool newly_solved;
	do {
		unsolved_remaining = false;
		newly_solved = false;
		for (auto& x: input_values) {
			string value_name = x.first;
			string unit = x.second.unit;
			if (unit != unit_force && unit != "strain") {
				if (force_ratios.find(unit+'/'+unit_force) != force_ratios.end()) {
					changeUnit(x.second, unit_force);
					if (x.second.type == "force") {
						force_ratios[value_name+'/'+unit_force] = *(x.second.value);
					}
					newly_solved = true;
				} else {
					unsolved_remaining = true;
				}
			}
		}
	} while (unsolved_remaining && newly_solved);

	// complain if we have not found everyone
	if (unsolved_remaining) {
		ostringstream error_str;
		for (const auto& x: input_values) {
			string value_name = x.first;
			string unit = x.second.unit;
			if (unit != unit_force && unit != "strain") {
				error_str << "Error: input value \"" << value_name << "\" has an unknown unit \"" << unit << "\"" << endl;
			}
		}
		throw runtime_error(error_str.str());
	}

	// determine the remaining force_ratios
	buildFullSetOfForceRatios();
}

void Simulation::catchForcesInStressUnits(const string &stress_unit)
{
	for (auto& x: input_values) {
		if (x.second.unit == "stress") {
			if (x.second.type != "force") {
				throw runtime_error(" Simulation:: "+x.first+" is a "+x.second.type+", which cannot be given in stress units.");
			} else {
				x.second.unit = stress_unit;
				*(x.second.value) *= abs(sys.target_stress);
			}
		}
	}
}

void Simulation::setupNonDimensionalizationStressControlled(double dimensionlessnumber,
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
	if (stress_unit == "hydro") {
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
	if (stress_unit == "thermal") {
		throw runtime_error(" Error: stress controlled Brownian simulations are not yet implemented.");
	}
	sys.set_shear_rate(0);
	// we take as a unit scale the one given by the user with the stress
	// TODO: other choices may be better when several forces are used.
	internal_unit_scales = stress_unit;
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
void Simulation::setupNonDimensionalizationRateControlled(double dimensionlessnumber,
														  string input_scale)
{
	/**
	 \brief Choose units for the simulation and convert the forces to this unit (rate controlled case).

	 The strategy is the following:
	 1. Convert all the forces in hydro force unit, ie the value for force f1 is f1/F_H
	 2. Decide the unit force F_W for the simulation
	 3. Convert all the forces to this unit (by multiplying every force by F_H/F_W)
	 */
	if (input_values.find(input_scale) != input_values.end()) { // if the force defining the shear rate is redefined in the parameter file, throw an error
		ostringstream error_str;
		error_str  << "Error: redefinition of the rate (given both in the command line and in the parameter file with \"" << input_scale << "\" force)" << endl;
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
	DimensionalValue inv;
	inv.type = "force";
	inv.value = force_value_ptr[input_scale];
	inv.unit = "hydro";
	input_values[input_scale] = inv;
	*(input_values[input_scale].value) = 1/dimensionlessnumber;
	force_ratios[input_scale+"/hydro"] = *(input_values[input_scale].value);
	// convert all other forces to hydro
	resolveUnitSystem("hydro");
	// chose simulation unit
	// the chosen unit is called internal_unit_scales
	setUnitScaleRateControlled();
	// convert from hydro scale to chosen scale
	for (auto& x: input_values) {
			changeUnit(x.second, internal_unit_scales);
	}
}

void Simulation::setLowPeclet()
{
	sys.lowPeclet = true;
	double scale_factor_SmallPe = p.Pe_switch/force_ratios["hydro/thermal"];
	p.dt *= p.Pe_switch; // to make things continuous at Pe_switch
}

void Simulation::changeUnit(DimensionalValue &x, string new_unit)
{
	/**
	 \brief Convert DimensionalValue x from unit x.unit to unit new_unit.
	 */
	if (new_unit != x.unit) {
		if (x.type == "force") {
			*(x.value) *= force_ratios[x.unit+'/'+new_unit];
			x.unit = new_unit;
		} else if (x.type == "time") {
			if (x.unit != "strain") {
				*(x.value) /= force_ratios[x.unit+'/'+new_unit];
				x.unit = new_unit;
			}
		} else {
			throw runtime_error(" Simulation:: Don't know how to change unit for DimensionalValue of type "+x.type+"\n");
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
	if (force_ratios.find("hydro/thermal") != force_ratios.end()) {
		is_brownian = true;
	} else {
		is_brownian = false;
	}
	if (is_brownian) {
		if (force_ratios["hydro/thermal"] > p.Pe_switch && !sys.zero_shear) { // hydro units
			internal_unit_scales = "hydro";
		} else { // low Peclet mode
			internal_unit_scales = "thermal";
			setLowPeclet();
		}
	} else {
		internal_unit_scales = "hydro";
	}
}

void Simulation::exportForceAmplitudes()
{
	/**
	 \brief Copy the input force alues in the ForceAmplitude struct of the System class
	 */
	string indent = "  Simulation::\t";
	cout << indent+"Forces used:" << endl;
	indent += "\t";

	sys.repulsiveforce = input_values.find("repulsion") != input_values.end();
	if (sys.repulsiveforce) {
		cout << indent+"Repulsive force (in \"" << input_values["repulsion"].unit << "\" units): " << sys.p.repulsion << endl;
	}

	sys.critical_load = input_values.find("critical_load") != input_values.end();
	if (sys.critical_load) {
		cout << indent+"Critical Load (in \"" << input_values["critical_load"].unit << "\" units): " << sys.p.critical_load << endl;
	}

	sys.cohesion = input_values.find("cohesion") != input_values.end();
	if (sys.cohesion) {
		cout << indent+"Cohesion (in \"" << input_values["cohesion"].unit << "\" units): " << sys.p.cohesion << endl;
	}

	bool is_ft_max = input_values.find("ft_max") != input_values.end();
	if (is_ft_max) {
		cout << indent+"Max tangential load (in \"" << input_values["ft_max"].unit << "\" units): " << sys.p.ft_max << endl;
	}

	sys.brownian = input_values.find("brownian") != input_values.end();
	if (sys.brownian) {
		cout << indent+"Brownian force (in \"" << input_values["brownian"].unit << "\" units): " << sys.p.brownian << endl;
	}
	cout << indent+"Normal contact stiffness (in \"" << input_values["kn"].unit << "\" units): " << p.kn << endl;
	cout << indent+"Sliding contact stiffness (in \"" << input_values["kt"].unit << "\" units): " << p.kt << endl;
	cout << indent+"Rolling contact stiffness (in \"" << input_values["kr"].unit << "\" units): " << p.kr << endl;

	if (input_values.find("hydro") != input_values.end()) { // == if rate controlled
		sys.set_shear_rate(*(input_values["hydro"].value));
	}
}

void Simulation::setupNonDimensionalization(double dimensionlessnumber,
											string input_scale){
	/**
	 \brief Non-dimensionalize the simulation.

		This function determines the most appropriate unit scales to use in the System class depending on the input parameters (Brownian/non-Brownian, shear rate, stress/rate controlled), and converts all the input values in these units.
	 */
	input_scale = unit_longname[input_scale];
	if (control_var == rate) {
		input_rate = dimensionlessnumber; // @@@ Renaming is required?
	}
	if (control_var == rate) {
		setupNonDimensionalizationRateControlled(dimensionlessnumber, input_scale);
	} else if (control_var == stress) {
		setupNonDimensionalizationStressControlled(dimensionlessnumber, input_scale);
	} else {
		ostringstream error_str;
		error_str  << " Error: unknown control variable \"" << control_var 	<< "\"" << endl;
		throw runtime_error(error_str.str());
	}
	exportForceAmplitudes();
	string indent = "  Simulation::\t";
	cout << indent << "internal_unit_scales = " << internal_unit_scales << endl;
	sys.ratio_unit_time = &force_ratios[input_scale+"/"+internal_unit_scales];
	output_unit_scales = input_scale;
}

void Simulation::assertParameterCompatibility()
{
	// test for incompatibilities
	if (sys.brownian == true) {
		if (sys.pairwise_resistance && p.integration_method != 1) {
			ostringstream error_str;
			error_str << "Brownian simulation needs to use the Predictor-Corrector method." << endl;
			error_str << "Modify the parameter file." << endl;
			throw runtime_error(error_str.str());
		}
	}
	if (control_var == stress) {
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
	setDefaultParameters(input_scale);
	readParameterFile(filename_parameters);
	setupOptionalSimulation(indent);
	tagStrainParameters();
	setupNonDimensionalization(dimensionlessnumber, input_scale);

	assertParameterCompatibility();

	if (input_files[3] != "not_given") {
		throw runtime_error("pre-simulation data deprecated?");
	}
	resolveTimeOrStrainParameters();

	if (binary_conf) {
		int format = getBinaryConfigurationFileFormat(filename_import_positions);
		ifstream file_import;
		file_import.open(filename_import_positions.c_str(), ios::binary | ios::in);
		if (!file_import) {
			ostringstream error_str;
			error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
			throw runtime_error(error_str.str());
		}

		switch(format) {
			case BIN_FORMAT_BASE_NEW:
				{
					auto conf = readBinaryBaseConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case BIN_FORMAT_FIXED_VEL:
				{
					auto conf = readBinaryFixedVeloConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
		}
	} else {
		int format = getTxtConfigurationFileFormat(filename_import_positions);
		ifstream file_import;
		file_import.open(filename_import_positions.c_str());
		if (!file_import) {
			ostringstream error_str;
			error_str  << " Position file '" << filename_import_positions << "' not found." <<endl;
			throw runtime_error(error_str.str());
		}

		switch(format) {
			case TXT_FORMAT_BASE_OLD:
				{
					auto conf = readTxtBaseConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case TXT_FORMAT_BASE_NEW:
				{
					auto conf = readTxtBaseConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case TXT_FORMAT_FIXED_VEL:
				{
					auto conf = readTxtFixedVeloConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
			case TXT_FORMAT_CIRCULAR_COUETTE:
				{
					auto conf = readTxtCircularCouetteConfiguration(file_import);
					sys.setupConfiguration(conf, control_var);
					break;
				}
		}
	}

	p_initial = p;
	simu_name = prepareSimulationName(binary_conf, filename_import_positions, filename_parameters,
									  simu_identifier, dimensionlessnumber, input_scale);
	openOutputFiles(simu_name);
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
		input_values[keyword] = str2DimensionalValue("force", keyword, value, force_value_ptr[keyword]);
	} else if (keyword == "cohesion") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, force_value_ptr[keyword]);
	} else if (keyword == "brownian") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, force_value_ptr[keyword]);
	} else if (keyword == "critical_load") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, force_value_ptr[keyword]);
	} else if (keyword == "monolayer") {
		p.monolayer = str2bool(value);
	} else if (keyword == "repulsive_length") {
		p.repulsive_length = atof(value.c_str());
	} else if (keyword == "repulsive_max_length") {
		p.repulsive_max_length = atof(value.c_str());
	} else if (keyword == "lub_reduce_parameter") {
		p.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time") {
		input_values[keyword] = str2DimensionalValue("time", keyword, value, &p.contact_relaxation_time);
	} else if (keyword == "contact_relaxation_time_tan"){
		input_values[keyword] = str2DimensionalValue("time", keyword, value, &p.contact_relaxation_time_tan);
	} else if (keyword == "disp_max") {
		p.disp_max = atof(value.c_str());
	} else if (keyword == "time_end") {
		input_values[keyword] = str2DimensionalValue("time", keyword, value, &p.time_end);
	} else if (keyword == "integration_method") {
		p.integration_method = atoi(value.c_str());
	} else if (keyword == "lub_max_gap") {
		p.lub_max_gap = atof(value.c_str());
	} else if (keyword == "interaction_range") {
		p.interaction_range = atof(value.c_str());
	} else if (keyword == "sd_coeff") {
		p.sd_coeff = atof(value.c_str());
	} else if (keyword == "kn") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, force_value_ptr[keyword]);
	} else if (keyword == "kt") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, force_value_ptr[keyword]);
	} else if (keyword == "kr") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, force_value_ptr[keyword]);
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
		input_values[keyword] = str2DimensionalValue("time", keyword, value, &p.time_interval_output_config);
	} else if (keyword == "time_interval_output_data") {
		input_values[keyword] = str2DimensionalValue("time", keyword, value, &p.time_interval_output_data);
	} else if (keyword == "log_time_interval") {
		p.log_time_interval = str2bool(value);
	} else if (keyword == "initial_log_time") {
		input_values[keyword] = str2DimensionalValue("time", keyword, value, &p.initial_log_time);
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
		input_values[keyword] = str2DimensionalValue("force", keyword, value, &p.min_kn);
	} else if (keyword == "max_kn") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, &p.max_kn);
	} else if (keyword == "min_kt") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, &p.min_kt);
	} else if (keyword == "max_kt") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, &p.max_kt);
	} else if (keyword == "min_dt") {
		p.min_dt = atof(value.c_str());
	} else if (keyword == "max_dt") {
		p.max_dt = atof(value.c_str());
	} else if (keyword == "rest_threshold") {
		p.rest_threshold = atof(value.c_str());
	} else if (keyword == "ft_max") {
		input_values[keyword] = str2DimensionalValue("force", keyword, value, force_value_ptr[keyword]);
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
	} else if (keyword == "out_binary_conf") {
		p.out_binary_conf = str2bool(value);
	} else if (keyword == "np_fixed") {
		p.np_fixed = atoi(value.c_str());
	} else if (keyword == "keep_input_strain") {
		p.keep_input_strain = str2bool(value);
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
	autoSetParameters("Pe_switch", "5");
	autoSetParameters("dt", "1e-4");
	autoSetParameters("disp_max", "1e-3");
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
	autoSetParameters("contact_relaxation_time", "1e-3"+input_scale);
	autoSetParameters("contact_relaxation_time_tan", "-1"+input_scale);
	if (input_scale != "kn") {
		autoSetParameters("kn", "2000"+input_scale);
		autoSetParameters("min_kn", "1000"+input_scale);
		autoSetParameters("max_kn", "1000000"+input_scale);
	}
	if (input_scale != "kt") {
		autoSetParameters("kt", "0.5kn");
		autoSetParameters("min_kt", "1000"+input_scale);
		autoSetParameters("max_kt", "1000000"+input_scale);
	}
	if (input_scale != "kr") {
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
	autoSetParameters("theta_shear", "0");
	autoSetParameters("event_handler", "");
	autoSetParameters("simulation_mode", "0");
	autoSetParameters("keep_input_strain", "false");
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

void Simulation::openOutputFiles(string simu_name)
{
	/**
	 \brief Set up the output files

		This function determines a simulation name from the parameters, opens the output files with the corresponding name and prints their header.
	 */

	stringstream data_header;
	createDataHeader(data_header);
	outdata.setFile("data_"+simu_name+".dat", data_header.str(), force_to_run);
	outdata_st.setFile("st_"+simu_name+".dat", data_header.str(), force_to_run);
	if (!p.out_particle_stress.empty()) {
		outdata_pst.setFile("pst_"+simu_name+".dat", data_header.str(), force_to_run);
	}
	string time_filename = "t_"+simu_name+".dat";
	fout_time.open(time_filename.c_str());
	string input_filename = "input_"+simu_name+".dat";
	fout_input.open(input_filename.c_str());
	if (p.out_data_particle) {
		outdata_par.setFile("par_"+simu_name+".dat", data_header.str(), force_to_run);
	}
	if (p.out_data_interaction) {
		outdata_int.setFile("int_"+simu_name+".dat", data_header.str(), force_to_run);
	}
}

string Simulation::prepareSimulationName(bool binary_conf,
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
		for (const auto& f: input_values) {
			if (f.first.find("stiffness") == std::string::npos) {
				string_control_parameters << "_" << f.first << *(f.second.value) << f.second.unit;
			}
		}
	}
	if (control_var==rate) {
		string_control_parameters << "_" << "rate";
	}
	if (control_var==stress) {
		string_control_parameters << "_" << "stress";
	}
	if (control_var==viscnb) {
		string_control_parameters << "_" << "viscnb";
	}
	string_control_parameters << dimensionlessnumber << input_scale;
	ss_simu_name << string_control_parameters.str();
	if (simu_identifier != "") {
		ss_simu_name << "_";
		ss_simu_name << simu_identifier;
	}
	string indent = "  Simulation::\t";
	cout << indent << "filename: " << ss_simu_name.str() << endl;

	return ss_simu_name.str();
}

TimeKeeper Simulation::initTimeKeeper() {
	TimeKeeper tk;
	if (p.log_time_interval) {
		tk.addClock("data", LogClock(p.initial_log_time,
									 p.time_end,
									 p.nb_output_data_log_time,
									 input_values["time_end"].unit == "strain"));
	} else {
		tk.addClock("data", LinearClock(p.time_interval_output_data,
										input_values["time_interval_output_data"].unit == "strain"));
	}
	if (p.log_time_interval) {
		tk.addClock("config", LogClock(p.initial_log_time,
									   p.time_end,
									   p.nb_output_config_log_time,
									   input_values["time_end"].unit == "strain"));
	} else {
		tk.addClock("config", LinearClock(p.time_interval_output_config,
										  input_values["time_interval_output_config"].unit == "strain"));
	}
	return tk;
}
