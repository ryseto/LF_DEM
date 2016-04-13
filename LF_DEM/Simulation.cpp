//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#define _USE_MATH_DEFINES
#include "Simulation.h"
#include <cmath>
#include <map>
#include <string>
#include <cstring>
#include <sstream>
#include <cctype>
#include <stdexcept>
#include "Timer.h"

using namespace std;

Simulation::Simulation():
sys(System(p, events)),
shear_rate_expectation(-1),
internal_unit_scales("hydro"),
target_stress_input(0),
diminish_output(false)
{
	unit_longname["h"] = "hydro";
	unit_longname["r"] = "repulsive";
	unit_longname["b"] = "thermal";
	unit_longname["c"] = "cohesive";
	unit_longname["cl"] = "critical_load";
	unit_longname["m"] = "magnetic";
	unit_longname["ft"] = "ft";
	unit_longname["kn"] = "normal_stiffness";
	unit_longname["kt"] = "tan_stiffness";
	unit_longname["kr"] = "roll_stiffness";
	kill = false;
};

Simulation::~Simulation()
{
	if (fout_particle.is_open()) {
		fout_particle.close();
	}
	if (fout_interaction.is_open()) {
		fout_interaction.close();
	}
};

bool Simulation::keepRunning()
{
	/** \brief Determine if we reached the end of the simulation.

		Returns true when ParameterSet::time_end is reached or if an event handler threw a kill signal.
	 */
	if (time_end == -1) {
		return (fabs(sys.get_shear_strain()) < strain_end-1e-8) && !kill;
	} else {
		return (sys.get_time() < time_end-1e-8) && !kill;
	}
}

void Simulation::setupEvents()
{
	/** \brief Set up the types of events to be watched by the System class.

		Links System::eventLookUp to a specialized function according to the value of ParameterSet::event_handler .
	 */
	if (p.event_handler == "shear_jamming") {
		sys.eventLookUp = &System::eventShearJamming;
		return;
	}
	if (p.event_handler == "fragility") {
		sys.eventLookUp = &System::eventShearJamming;
		return;
	}
	sys.eventLookUp = NULL;
}

void Simulation::handleEventsShearJamming()
{
	/** \brief Event handler to test for shear jamming

		When a negative_shear_rate event is thrown, p.disp_max is decreased.
	 If p.disp_max is below a minimal value, the shear direction is switched to y-shear.
	 */
	for (const auto& ev : events) {
		if (ev.type == "negative_shear_rate") {
			cout << " negative rate " << endl;
			p.disp_max /= 1.1;
		}
	}
	if (p.disp_max < 1e-6) {
		cout << "jammed" << endl;
		evaluateData();
		outputData(); // new
		outputConfigurationBinary();
		outputConfigurationData();
		kill = true;
	}
}

void Simulation::handleEventsFragility()
{
	/** \brief Event handler to test for shear jamming

		When a negative_shear_rate event is thrown, p.disp_max is decreased.
	 If p.disp_max is below a minimal value, the shear direction is switched to y-shear.
	 */
	for (const auto& ev : events) {
		if (ev.type == "negative_shear_rate") {
			cout << " negative rate " << endl;
			p.disp_max /= 1.1;
		}
	}
	if (p.disp_max < 1e-6 || sys.get_shear_strain() > 3.) {
		p.cross_shear = true; //!p.cross_shear;
		p.disp_max = p_initial.disp_max;
		cout << "Event Fragility : starting cross shear" << endl;
	}
}

void Simulation::handleEvents()
{
	/** \brief Handle the list of events that appears in the previous time step

		This function dispatches to specialized handlers according to the value of ParameterSet::event_handler .
	 */
	if (p.event_handler == "shear_jamming") {
		handleEventsShearJamming();
	}
	if (p.event_handler == "fragility") {
		handleEventsFragility();
	}
	events.clear();
}

void Simulation::generateOutput(const set<string> &output_events, int& binconf_counter)
{
	outputConfigurationBinary(); // generic, for recovery if crash
	if(output_events.find("data") != output_events.end()) {
		evaluateData();
		outputData();
	}

	if(output_events.find("config") != output_events.end()) {
		if(p.out_binary_conf){
			string binconf_filename = "conf_" + sys.simu_name + "_" + to_string(++binconf_counter) + ".bin";
			outputConfigurationBinary(binconf_filename);
		} else {
			outputConfigurationData();
		}
	}
}

/*
 * Main simulation
 */
void Simulation::simulationSteadyShear(string in_args,
									   vector<string>& input_files,
									   bool binary_conf,
									   double dimensionless_number,
									   string input_scale,
									   string control_variable,
									   string simu_identifier)
{
	string indent = "  Simulation::\t";
	control_var = control_variable;
	/***************  This part is temporal ********************/
	// TODO
	if (simu_identifier == "mtest1") {
		cout << indent << "Test simulation for reversibility in mixed problem" << endl;
		sys.test_simulation = 1;//mtest1
	} else if (simu_identifier == "mtest2") {
		cout << indent << "Test simulation for a mixed problem" << endl;
		sys.test_simulation = 2;//mtest2
	} else if (simu_identifier == "mtest3") {
		cout << indent << "Test simulation for a mixed problem" << endl;
		sys.test_simulation = 3;//mtest2
	} else if (simu_identifier == "ctest1") {
		cout << indent << "Test simulation with co-axial cylinders (rotate outer clynder)" << endl;
		sys.wall_rheology = true;
		sys.test_simulation = 11;//ctest1
	} else if (simu_identifier == "ctest2") {
		cout << indent << "Test simulation with co-axial cylinders (rotate inner clynder)" << endl;
		sys.wall_rheology = true;
		sys.test_simulation = 12;//ctest2
	} else if (simu_identifier == "ctest3") {
		cout << indent << "Test simulation with co-axial cylinders (rotate both inner and outer clynder)" << endl;
		sys.wall_rheology = true;
		sys.test_simulation = 13;//ctest3
	} else if (simu_identifier == "ctest0") {
		cout << indent << "Test simulation with co-axial cylinders (rotate both inner and outer clynder)" << endl;
		sys.wall_rheology = true;
		sys.test_simulation = 10;//ctest3
	} else if (simu_identifier == "rtest1") {
		cout << indent << "Test simulation for shear reversibility" << endl;
		sys.test_simulation = 21;//rtest1
	} else if (simu_identifier == "wtest1") {
		cout << indent << "Test simulation (wtest1), simple shear with walls" << endl;
		sys.test_simulation = 31;//wtest1
	} else if (simu_identifier == "wtestA") {
		cout << indent << "Test simulation (wtestA), simple shear with walls" << endl;
		sys.wall_rheology = true;
		sys.test_simulation = 41;//wtestA
	} else if (simu_identifier == "wtestB") {
		cout << indent << "Test simulation (wtestB), simple shear with walls" << endl;
		sys.wall_rheology = true;
		sys.test_simulation = 42;//wtestB
	} else if (simu_identifier == "stest1") {
		cout << indent << "Test simulation (wtestB), simple shear with walls" << endl;
		sys.wall_rheology = true;
		sys.test_simulation = 51;//stest1
	}

	/*************************************************************/
	setupSimulation(in_args, input_files, binary_conf, dimensionless_number, input_scale, simu_identifier);
	time_t now;
	time_strain_1 = 0;
	now = time(NULL);
	time_strain_0 = now;

	setupEvents();
	cout << indent << "Time evolution started" << endl << endl;
	TimeKeeper tk;
	if (p.log_time_interval) {
		tk.addClock("data", LogClock(p.initial_log_time,
									 p.time_end,
									 p.nb_output_data_log_time,
									 input_values["time_end"].unit == "strain"));
	} else {
		tk.addClock("data", LinearClock(time_end,
										p.time_interval_output_data,
										input_values["time_end"].unit == "strain"));
	}
	if (p.log_time_interval) {
		tk.addClock("config", LogClock(p.initial_log_time,
									   p.time_end,
									   p.nb_output_config_log_time,
									   input_values["time_end"].unit == "strain"));
	} else {
		tk.addClock("config", LinearClock(time_end,
										  p.time_interval_output_config,
										  input_values["time_end"].unit == "strain"));
	}
	int binconf_counter = 0;
	while (keepRunning()) {
		pair<double, string> t = tk.nextTime();
		pair<double, string> s = tk.nextStrain();
		if (t.second.empty()) { // no next time
			sys.timeEvolution(-1, s.first);
		} else if (s.second.empty()) { // no next strain
			sys.timeEvolution(t.first, -1);
		} else { // either next time or next strain
			sys.timeEvolution(t.first, s.first);
		}
		handleEvents();

		set<string> output_events = tk.getElapsedClocks(sys.get_time(), fabs(sys.get_shear_strain()));
		generateOutput(output_events, binconf_counter);

		if (time_end != -1) {
			cout << "time: " << sys.get_time_in_simulation_units() << " , " << sys.get_time() << " / " << time_end << " , strain: " << sys.get_shear_strain() << endl;
		} else {
			cout << "time: " << sys.get_time_in_simulation_units() << " , strain: " << sys.get_shear_strain() << " / " << strain_end << endl;
		}
		if (time_strain_1 == 0 && fabs(sys.get_shear_strain()) > 1) {
			now = time(NULL);
			time_strain_1 = now;
			timestep_1 = sys.get_total_num_timesteps();
		}
	}
	now = time(NULL);
	time_strain_end = now;
	timestep_end = sys.get_total_num_timesteps();
	outputComputationTime();
	string filename_parameters = input_files[1];
	if (filename_parameters.find("init_relax", 0)) {
		/* To prepare relaxed initial configuration,
		 * we can use Brownian simulation for a short interval.
		 * Here is just to export the position data.
		 */
		string filename_configuration = input_files[0];
		outputFinalConfiguration(filename_configuration);
	}
	cout << indent << "Time evolution done" << endl << endl;
}

void Simulation::simulationInverseYield(string in_args,
										vector<string>& input_files,
										bool binary_conf,
										double dimensionless_number,
										string input_scale,
										string control_variable,
										string simu_identifier)
{
	control_var = control_variable;
	setupSimulation(in_args, input_files, binary_conf, dimensionless_number, input_scale, simu_identifier);

	int jammed = 0;
	time_t now;
	time_strain_1 = 0;
	now = time(NULL);
	time_strain_0 = now;
	/******************** OUTPUT INITIAL DATA ********************/
	evaluateData();
	outputData(); // new
	outputConfigurationBinary();
	outputConfigurationData();
	/*************************************************************/

	TimeKeeper tk;
	tk.addClock("data", LinearClock(time_end,
									p.time_interval_output_data,
									input_values["time_end"].unit == "strain"));
	tk.addClock("config", LinearClock(time_end,
										p.time_interval_output_config,
										input_values["time_end"].unit == "strain"));
	int binconf_counter = 0;
	while (keepRunning()) {
		pair<double, string> t = tk.nextTime();
		pair<double, string> s = tk.nextStrain();

		if (t.second.empty()) { // no next time
			sys.timeEvolution(-1, s.first);
		} else if (s.second.empty()) { // no next strain
			sys.timeEvolution(t.first, -1);
		} else { // either next time or next strain
			sys.timeEvolution(t.first, s.first);
		}
		set<string> output_events = tk.getElapsedClocks(sys.get_time(), fabs(sys.get_shear_strain()));
		generateOutput(output_events, binconf_counter);


		cout << "time: " << sys.get_time() << " / " << p.time_end << endl;
		if (!sys.zero_shear
			&& abs(sys.get_shear_rate()) < p.rest_threshold) {
			cout << "shear jamming " << jammed << endl;
			jammed ++;
			if (jammed > 20) {
				sys.set_shear_rate(1);
				cout << "target_stress = " << target_stress_input << endl;
				target_stress_input *= 0.95;
				sys.target_stress = target_stress_input/6/M_PI;
				sys.updateUnscaledContactmodel();
				cout << "new target_stress = " << target_stress_input << endl;
				jammed = 0;
			}
		} else {
			jammed = 0;
		}
		if (time_strain_1 == 0 && sys.get_shear_strain() > 1) {
			now = time(NULL);
			time_strain_1 = now;
			timestep_1 = sys.get_total_num_timesteps();
		}
	}

	now = time(NULL);
	time_strain_end = now;
	timestep_end = sys.get_total_num_timesteps();
	outputComputationTime();

	string	filename_parameters = input_files[1];
	if (filename_parameters.find("init_relax", 0)) {
		/* To prepare relaxed initial configuration,
		 * we can use Brownian simulation for a short interval.
		 * Here is just to export the position data.
		 */
		string filename_configuration = input_files[0];
		outputFinalConfiguration(filename_configuration);
	}
}

void Simulation::simulationMagnetic(string in_args,
									vector<string>& input_files,
									bool binary_conf,
									double dimensionless_number,
									string input_scale,
									string control_variable,
									string simu_identifier)
{
	/* Monolayer: Particles are confined in y = 0 plane.
	 *
	 *
	 */
	control_var = control_variable;
	setupSimulation(in_args, input_files, binary_conf,
					dimensionless_number, input_scale, simu_identifier);
	int cnt_simu_loop = 1;
	int cnt_config_out = 1;
	double time_output_data = 0;
	double time_output_config = 0;
	/******************** OUTPUT INITIAL DATA ********************/
	evaluateData();
	outputDataMagnetic();
	outputConfigurationBinary();
	outputConfigurationData();
	/*************************************************************/
	if (sys.p.magnetic_field_type == 0) {
		// Field direction is fixed
		sys.external_magnetic_field.set(0, 1, 0);
	} else if (sys.p.magnetic_field_type == 1) {
		sys.external_magnetic_field.set(sin(sys.angle_external_magnetic_field),
										cos(sys.angle_external_magnetic_field),
										0);
		throw runtime_error("magnetic_field_type == 1 not yet implemented");// @not yet
	} else if (sys.p.magnetic_field_type == 2) {
		sys.external_magnetic_field.set(cos(sys.angle_external_magnetic_field),
										0,
										sin(sys.angle_external_magnetic_field));
		throw runtime_error("magnetic_field_type == 2 not yet implemented");// @not yet
	}
	sys.setMagneticMomentZero();
	bool initial_relax = true;
	// Main simulation loop
	double initial_time = sys.get_time();
	string indent = "  Simulation::\t";
	cout << indent << "Time evolution started" << endl << endl;
	while (keepRunning()) {
		if (initial_relax && sys.get_time() >= 0) {
			sys.setInducedMagneticMoment();
			initial_relax = false;
		}
		time_output_data = initial_time+cnt_simu_loop*p.time_interval_output_data;
		time_output_config = initial_time+cnt_config_out*p.time_interval_output_config;
		sys.timeEvolution(time_output_data, -1); // @@@ I changed to new timeEvolution method, is that ok? The old one is not as flexible so I would like to deprecate it
		cnt_simu_loop ++;
		/******************** OUTPUT DATA ********************/
		evaluateData();
		outputDataMagnetic();
		outputConfigurationBinary();
		if (sys.get_time() >= time_output_config-1e-8) {
			outputConfigurationData();
			cnt_config_out ++;
		}
		cout << "time: " << sys.get_time() << " / " << time_end << endl;
	}
	outputComputationTime();
	string	filename_parameters = input_files[1];
	if (filename_parameters.find("init_relax", 0) < filename_parameters.size()) {
		/* To prepare relaxed initial configuration,
		 * we can use Brownian simulation for a short interval.
		 * Here is just to export the position data.
		 */
		string filename_configuration = input_files[0];
		outputFinalConfiguration(filename_configuration);
	}
}

void Simulation::outputComputationTime()
{
	int time_from_1 = time_strain_end-time_strain_1;
	int time_from_0 = time_strain_end-time_strain_0;
	int timestep_from_1 = timestep_end-timestep_1;
	fout_time << "# np time_from_0 time_from_1 timestep_end timestep_from_1" << endl;
	fout_time << sys.get_np() << ' ';
	fout_time << time_from_0 << ' ';
	fout_time << time_from_1 << ' ';
	fout_time << timestep_end << ' ';
	fout_time << timestep_from_1 << endl;
}

void Simulation::catchSuffixedForce(const string& keyword,
									const string& value)
{
	string numeral, suffix;
	bool caught_suffix = getSuffix(value, numeral, suffix);
	suffix = unit_longname[suffix];
	input_force_units[keyword] = suffix;
	input_force_values[keyword] = atof(numeral.c_str());

	if (!caught_suffix) {
		errorNoSuffix(keyword);
	}
}

void Simulation::catchSuffixedValue(string type, string keyword,
									string value_str, double *value_ptr)
{
	InputValue inv;
	inv.type = type;
	inv.value = value_ptr;

	string numeral, suffix;
	bool caught_suffix = true;
	caught_suffix = getSuffix(value_str, numeral, suffix);
	if (!caught_suffix) {
		errorNoSuffix(keyword);
	}
	suffix = unit_longname[suffix];
	*(inv.value) = atof(numeral.c_str());
	inv.unit = suffix;
	input_values[keyword] = inv;
}

void Simulation::outputConfigurationBinary()
{
	string conf_filename;
	//	conf_filename =  "conf_" + sys.simu_name + "_strain" + to_string(sys.get_shear_strain()) + ".dat";
	conf_filename = "conf_" + sys.simu_name + ".dat";
	outputConfigurationBinary(conf_filename);
}

void Simulation::outputConfigurationBinary(string conf_filename)
{
	/**
		\brief Saves the current configuration of the system in a binary file.

		In the current implementation, it stores particle positions, x and y strain, and contact states.
	 */
	int np = sys.get_np();
	vector< vector<double> > pos(np);
	int dims = 4;
	for (int i=0; i<np; i++) {
		pos[i].resize(dims);
		pos[i][0] = sys.position[i].x;
		pos[i][1] = sys.position[i].y;
		pos[i][2] = sys.position[i].z;
		pos[i][3] = sys.radius[i];
	}
	ofstream conf_export;
	double lx = sys.get_lx();
	double ly = sys.get_ly();
	double lz = sys.get_lz();
	conf_export.open(conf_filename.c_str(), ios::binary | ios::out);
	conf_export.write((char*)&np, sizeof(int));
	conf_export.write((char*)&volume_or_area_fraction, sizeof(double));
	conf_export.write((char*)&lx, sizeof(double));
	conf_export.write((char*)&ly, sizeof(double));
	conf_export.write((char*)&lz, sizeof(double));
	conf_export.write((char*)&(sys.shear_disp.x), sizeof(double));
	conf_export.write((char*)&(sys.shear_disp.y), sizeof(double));
	for (int i=0; i<np; i++) {
		conf_export.write((char*)&pos[i][0], dims*sizeof(double));
	}
	vector <struct contact_state> cs;
	sys.getContacts(cs);
	int ncont = cs.size();
	conf_export.write((char*)&ncont, sizeof(unsigned int));
	for (int i=0; i<ncont; i++) {
		conf_export.write((char*)&(cs[i].p0), sizeof(unsigned int));
		conf_export.write((char*)&(cs[i].p1), sizeof(unsigned int));
		conf_export.write((char*)&(cs[i].disp_tan.x), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_tan.y), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_tan.z), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_rolling.x), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_rolling.y), sizeof(double));
		conf_export.write((char*)&(cs[i].disp_rolling.z), sizeof(double));
	}
	//	conf_export.write((char*)&(sys.dt), sizeof(double));
	conf_export.close();
}

double Simulation::getRate()
{
	/**
	 \brief The shear rate in the input units
	 */
	if (control_var == "rate") {
		return input_rate;
	} else if (control_var == "stress") {
		return sys.get_shear_rate();
	} else {
		return 1;
	}
}

void Simulation::evaluateData()
{
	/**
	 \brief Get rheological data from the System class.

	 In this method we keep the internal units. There is no conversion to output units at this stage

	 */
	sys.analyzeState();
	sys.calcStress();
	//	if (sys.p.lubrication_model > 0) {
	//		sys.calcLubricationForce();
	//	}
}

void Simulation::outputData()
{
	/**
	 \brief Output data.

	 Stress data are converted in units of
	 \f$\eta_0\dot\gamma\f$. Other data are output in the units used
	 in the System class (these can be hydrodynamic, Brownian or
	 repulsive force units).

	 \b NOTE: this behavior should be changed
	 and made more consistent in the future.
	 */

	string dimless_nb_label = internal_unit_scales+"/"+output_unit_scales;
	if (dimensionless_numbers.find(dimless_nb_label) == dimensionless_numbers.end()) {
		ostringstream error_str;
		error_str << " Error : don't manage to convert from \"" << internal_unit_scales << "\" units to \"" << output_unit_scales << "\" units to output data." << endl;
		throw runtime_error(error_str.str());
	}
	outdata.setDimensionlessNumber(dimensionless_numbers[dimless_nb_label]);

	int number_of_data = 37;
	if (sys.wall_rheology) {
		number_of_data = 40;
	}
	if (sys.test_simulation==31) {
		number_of_data = 39;
	}
	outdata.init(number_of_data, output_unit_scales);
	double sr = sys.get_shear_rate();
	double shear_stress = shearStressComponent(sys.total_stress, p.theta_shear);
	outdata.entryData(1, "time", "time", sys.get_time());
	outdata.entryData(2, "shear strain", "none", sys.get_shear_strain());
	outdata.entryData(3, "shear rate", "rate", sys.get_shear_rate());

	outdata.entryData(5, "viscosity", "viscosity", shear_stress/sr);
	outdata.entryData(6, "Viscosity(lub)", "viscosity", shearStressComponent(sys.total_hydro_stress, p.theta_shear)/sr);
	outdata.entryData(7, "Viscosity(xF_contact part)", "viscosity", shearStressComponent(sys.total_contact_stressXF, p.theta_shear)/sr);
	outdata.entryData(8, "Viscosity(GU_contact part)", "viscosity", shearStressComponent(sys.total_contact_stressGU, p.theta_shear)/sr);
	if (sys.repulsiveforce) {
		outdata.entryData(9, "Viscosity(repulsive force XF)", "viscosity", shearStressComponent(sys.total_repulsive_stressXF, p.theta_shear)/sr);
		outdata.entryData(10, "Viscosity(repulsive force GU)", "viscosity", shearStressComponent(sys.total_repulsive_stressGU, p.theta_shear)/sr);
	}
	if (sys.brownian) {
		outdata.entryData(11, "Viscosity(brownian)", "viscosity", shearStressComponent(sys.total_brownian_stressGU, p.theta_shear)/sr);
	}
	/*
	 * Stress
	 */
	outdata.entryData(14, "shear stress", "stress", shear_stress);
	outdata.entryData(15, "N1 viscosity", "viscosity", sys.total_stress.getNormalStress1()/sr);
	outdata.entryData(16, "N2 viscosity", "viscosity", sys.total_stress.getNormalStress2()/sr);
	outdata.entryData(17, "particle pressure", "stress", sys.total_stress.getParticlePressure());
	outdata.entryData(18, "particle pressure contact", "stress", sys.total_contact_stressXF.getParticlePressure());
	/* energy
	 */
	outdata.entryData(21, "energy", "none", sys.get_total_energy());
	/* maximum deformation of contact bond
	 */
	outdata.entryData(22, "min gap", "none", sys.min_reduced_gap);
	outdata.entryData(23, "max gap(cohesion)", "none", sys.max_contact_gap);
	outdata.entryData(24, "max tangential displacement", "none", sys.max_disp_tan);
	outdata.entryData(25, "max rolling displacement", "none", sys.max_disp_rolling);
	/* contact number
	 */
	outdata.entryData(26, "contact number", "none", sys.getContactNumber());
	outdata.entryData(27, "frictional contact number", "none", sys.getFrictionalContactNumber());
	outdata.entryData(28, "number of interaction", "none", sys.get_nb_of_active_interactions());
	/* maximum velocity
	 */
	outdata.entryData(29, "max velocity", "velocity", sys.max_velocity);
	outdata.entryData(30, "max angular velocity", "velocity", sys.max_ang_velocity);
	/* simulation parameter
	 */
	outdata.entryData(31, "dt", "time", sys.avg_dt);
	outdata.entryData(32, "kn", "none", p.kn);
	outdata.entryData(33, "kt", "none", p.kt);
	outdata.entryData(34, "kr", "none", p.kr);
	outdata.entryData(35, "shear displacement x", "none", sys.shear_disp.x);
	outdata.entryData(36, "shear displacement y", "none", sys.shear_disp.y);
	if (sys.wall_rheology) {
		outdata.entryData(37, "shear viscosity wall 1", "viscosity", sys.shearstress_wall1/sr);
		outdata.entryData(38, "shear viscosity wall 2", "viscosity", sys.shearstress_wall2/sr);
		outdata.entryData(39, "normal stress/rate wall 1", "viscosity", sys.normalstress_wall1/sr);
		outdata.entryData(40, "normal stress/rate wall 2", "viscosity", sys.normalstress_wall2/sr);
	}
	if (sys.test_simulation == 31) {
		pair<double, double> fw = sys.checkForceOnWalls();
		outdata.entryData(37, "force top wall", "force", fw.first);
		outdata.entryData(38, "force bottom wall", "force", fw.second);
	}
	outdata.writeToFile();
	/****************************   Stress Tensor Output *****************/
	outdata_st.setDimensionlessNumber(dimensionless_numbers[dimless_nb_label]);
	outdata_st.init(8, output_unit_scales);
	outdata_st.entryData(1, "time", "time", sys.get_time());
	outdata_st.entryData(2, "shear strain", "none", sys.get_shear_strain());
	outdata_st.entryData(3, "shear rate", "rate", sys.get_shear_rate());
	outdata_st.entryData(4, "total stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_stress);
	outdata_st.entryData(5, "hydro stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_hydro_stress);
	outdata_st.entryData(6, "contact stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_contact_stressXF+sys.total_contact_stressGU);
	outdata_st.entryData(7, "repulsive stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_repulsive_stressXF+sys.total_repulsive_stressGU);
	outdata_st.entryData(8, "brownian stress tensor (xx, xy, xz, yz, yy, zz)", "stress", sys.total_brownian_stressGU);
	outdata_st.writeToFile();

	if (!p.out_particle_stress.empty()) {
		outdata_pst.setDimensionlessNumber(dimensionless_numbers[dimless_nb_label]);
		int nb_of_fields = strlen(p.out_particle_stress.c_str());
		outdata_pst.init(nb_of_fields, output_unit_scales);
		for (int i=0; i<sys.get_np(); i++) {
			int field_index = 0;
			if (p.out_particle_stress.find('t') != string::npos) {
				StressTensor s = sys.lubstress[i]+sys.contactstressXF[i]+sys.contactstressGU[i];
				if (sys.brownian) {
					s += sys.brownianstressGU[i];
				}
				if (sys.repulsiveforce) {
					s += sys.repulsivestressGU[i]+sys.repulsivestressXF[i];
				}
				outdata_pst.entryData(++field_index, "total stress (xx, xy, xz, yz, yy, zz), excluding magnetic stress", "stress", s);
			}
			if (p.out_particle_stress.find('l') != string::npos) {
				outdata_pst.entryData(++field_index, "lubrication stress (xx, xy, xz, yz, yy, zz)", "stress", sys.lubstress[i]);
			}
			if (p.out_particle_stress.find('c') != string::npos) {
				outdata_pst.entryData(++field_index, "contact stress (xx, xy, xz, yz, yy, zz)", "stress", sys.contactstressXF[i] + sys.contactstressGU[i]);
			}
			if (sys.brownian && p.out_particle_stress.find('b') != string::npos ) {
				outdata_pst.entryData(++field_index, "Brownian stress (xx, xy, xz, yz, yy, zz)", "stress", sys.brownianstressGU[i]);
			}
			if (sys.repulsiveforce && p.out_particle_stress.find('r') != string::npos ) {
				outdata_pst.entryData(++field_index, "repulsive stress (xx, xy, xz, yz, yy, zz)", "stress", sys.repulsivestressGU[i] + sys.repulsivestressXF[i]);
			}
		}
		stringstream snapshot_header;
		getSnapshotHeader(snapshot_header);
		outdata_pst.writeToFile(snapshot_header.str());
	}
}

void Simulation::getSnapshotHeader(stringstream& snapshot_header)
{
	snapshot_header << "# " << sys.get_shear_strain() << ' ';
	snapshot_header << sys.shear_disp.x << ' ';
	snapshot_header << getRate() << ' ';
	snapshot_header << target_stress_input << ' ';
	snapshot_header << sys.get_time() << ' ';
	if (p.magnetic_type != 0) {
		snapshot_header << sys.angle_external_magnetic_field;
	}
	snapshot_header << endl;
}

void Simulation::outputDataMagnetic()
{
	/**
	 \brief Output data for magnetic colloid work
	 */
	string dimless_nb_label = internal_unit_scales+"/"+output_unit_scales;
	if (dimensionless_numbers.find(dimless_nb_label) == dimensionless_numbers.end()) {
		ostringstream error_str;
		error_str << " Error : don't manage to convert from \"" << internal_unit_scales << "\" units to \"" << output_unit_scales << "\" units to output data." << endl;
		throw runtime_error(error_str.str());
	}
	outdata.setDimensionlessNumber(dimensionless_numbers[dimless_nb_label]);
	outdata.init(27, output_unit_scales);
	outdata.entryData(1, "time", "time", sys.get_time());
	/* energy */
	outdata.entryData(3, "energy", "none", sys.get_total_energy());
	outdata.entryData(4, "magnetic energy", "none", sys.magnetic_dd_energy);
	/* contact number */
	outdata.entryData(5, "contact number", "none", sys.getContactNumber());
	outdata.entryData(6, "frictional contact number", "none", sys.getFrictionalContactNumber());
	outdata.entryData(7, "number of interaction", "none", sys.get_nb_of_active_interactions());
	/* maximum deformation of contact bond */
	outdata.entryData(8, "min gap", "none", sys.min_reduced_gap);
	outdata.entryData(9, "max gap(cohesion)", "none", sys.max_contact_gap);
	outdata.entryData(10, "magnetic field strength", "none", sys.p.external_magnetic_field_norm);
	outdata.entryData(11, "magnetic field angle", "none", sys.p.external_magnetic_field_ang_theta);
	outdata.entryData(12, "magnetic field angle", "none", sys.p.external_magnetic_field_ang_phi);
	/* pressure */
	outdata.entryData(13, "particle pressure", "stress", sys.total_stress.getParticlePressure());
	outdata.entryData(14, "particle pressure contact XF", "stress", sys.total_contact_stressXF.getParticlePressure());
	outdata.entryData(15, "particle pressure contact GU", "stress", sys.total_contact_stressGU.getParticlePressure());
	outdata.entryData(16, "particle pressure brownian", "stress", sys.total_brownian_stressGU.getParticlePressure());
	outdata.entryData(17, "particle pressure magnetic XF", "stress", sys.total_magnetic_stressXF.getParticlePressure());
	outdata.entryData(18, "particle pressure magnetic GU", "stress", sys.total_magnetic_stressGU.getParticlePressure());
	/* maximum velocity */
	outdata.entryData(20, "max velocity", "velocity", sys.max_velocity);
	outdata.entryData(21, "max angular velocity", "velocity", sys.max_ang_velocity);
	/* simulation parameter */
	outdata.entryData(22, "dt", "none", sys.dt);
	outdata.entryData(23, "kn", "none", sys.p.kn);
	outdata.entryData(24, "kt", "none", sys.p.kt);
	outdata.entryData(25, "kr", "none", sys.p.kr);
	/* misc */
	outdata.writeToFile();
}

vec3d Simulation::shiftUpCoordinate(double x, double y, double z)
{
	if (p.origin_zero_flow) {
		z += sys.Lz_half();
		if (z > sys.Lz_half()) {
			x -= sys.shear_disp.x;
			y -= sys.shear_disp.y;
			if (x < -sys.Lx_half()) {
				x += sys.get_lx();
			}
			if (y < -sys.Ly_half()) {
				y += sys.get_ly();
			}
			z -= sys.get_lz();
		}
	}
	return vec3d(x,y,z);
}

void Simulation::createDataHeader(stringstream& data_header)
{
	data_header << "# LF_DEM version " << GIT_VERSION << endl;
	data_header << "# np " << sys.get_np() << endl;
	data_header << "# VF " << volume_or_area_fraction << endl;
	data_header << "# Lx " << sys.get_lx() << endl;
	data_header << "# Ly " << sys.get_ly() << endl;
	data_header << "# Lz " << sys.get_lz() << endl;
}

void Simulation::outputDataHeader(ofstream& fout)
{
	stringstream data_header;
	createDataHeader(data_header);
	fout << data_header.str();
}

void Simulation::outputConfigurationData()
{
	int np = sys.get_np();
	int output_precision = 6;
	if (diminish_output) {
		output_precision = 4;
	}

	vector<vec3d> pos(np);
	vector<vec3d> vel(np);
	for (int i=0; i<np; i++) {
		pos[i] = shiftUpCoordinate(sys.position[i].x-sys.Lx_half(),
								   sys.position[i].y-sys.Ly_half(),
								   sys.position[i].z-sys.Lz_half());
	}
	/* If the origin is shifted,
	 * we need to change the velocities of particles as well.
	 */
	for (int i=0; i<np; i++) {
		vel[i] = sys.velocity[i];
		if (p.origin_zero_flow) {
			if (pos[i].z < 0) {
				vel[i].x -= sys.get_shear_rate()*sys.get_lz();
			}
		}
	}
	/*
	 * shear_disp = sys.strain() - (int)(sys.strain()/Lx)*Lx
	 */
	if (p.out_data_particle) {
		cout << "   out config: " << sys.get_shear_strain() << endl;
		fout_particle << "# " << sys.get_shear_strain() << ' ';
		fout_particle << sys.shear_disp.x << ' ';
		fout_particle << getRate() << ' ';
		fout_particle << target_stress_input << ' ';
		fout_particle << sys.get_time() << ' ';
		if (p.magnetic_type != 0) {
			fout_particle << sys.angle_external_magnetic_field;
		}
		fout_particle << endl;

		for (int i=0; i<np; i++) {
			const vec3d& r = pos[i];
			const vec3d& v = vel[i];
			const vec3d& o = sys.ang_velocity[i];
			fout_particle << i; //1: number
			fout_particle << ' ' << sys.radius[i]; //2: radius
			fout_particle << setprecision(6);
			fout_particle << ' ' << r.x << ' ' << r.y << ' ' << r.z; //3, 4, 5: position
			fout_particle << setprecision(output_precision) << ' ' << v.x << ' ' << v.y << ' ' << v.z; //6, 7, 8: velocity
			if (control_var != "magnetic") {
				fout_particle << ' ' << o.x << ' ' << o.y << ' ' << o.z; //9, 10, 11: angular velocity
				if (sys.couette_stress) {
					double stress_rr, stress_thetatheta, stress_rtheta;
					sys.getStressCouette(i, stress_rr, stress_thetatheta, stress_rtheta);
					fout_particle << ' ' << stress_rr << ' ' << stress_thetatheta << ' ' << stress_rtheta;
				} else {
					if (diminish_output == false) {
						double lub_xzstress = shearStressComponent(sys.lubstress[i], p.theta_shear);
						double contact_xzstressGU = shearStressComponent(sys.contactstressGU[i], p.theta_shear);
						double brownian_xzstressGU = 0;
						if (sys.brownian) {
							brownian_xzstressGU = shearStressComponent(sys.brownianstressGU[i], p.theta_shear);
						}
						fout_particle << ' ' << 6*M_PI*lub_xzstress; //12: xz stress contributions //@@@ remove?
						fout_particle << ' ' << 6*M_PI*contact_xzstressGU; //13: xz stress contributions //@@@ remove?
						fout_particle << ' ' << 6*M_PI*brownian_xzstressGU; //14: xz stress contributions //@@@ remove?
					} else {
						fout_particle << " d d d";
					}
				}
				if (sys.twodimension) {
					fout_particle << ' ' << sys.angle[i]; // 15
				} else {
					fout_particle << ' ' << 0; // 15
				}
			} else {
				if (sys.p.magnetic_type == 1) {
					fout_particle << ' ' << o.x << ' ' << o.y << ' ' << o.z; //9, 10, 11: angular velocity
					fout_particle << ' ' << sys.magnetic_moment[i].x;
					fout_particle << ' ' << sys.magnetic_moment[i].y;
					fout_particle << ' ' << sys.magnetic_moment[i].z;
				} else {
					fout_particle << ' ' << sys.magnetic_susceptibility[i]; // 1: magnetic 0: non-magnetic
				}
			}
			if (p.out_data_vel_components) {
				fout_particle << setprecision(output_precision);
				fout_particle << ' ' << sys.vel_hydro[i];
				fout_particle << ' ' << sys.ang_vel_hydro[i];
				fout_particle << ' ' << sys.vel_contact[i];
				fout_particle << ' ' << sys.ang_vel_contact[i];
				if (sys.repulsiveforce) {
					fout_particle << ' ' << sys.vel_repulsive[i];
					fout_particle << ' ' << sys.ang_vel_repulsive[i];
				}
				if (sys.brownian) {
					fout_particle << ' ' << sys.vel_brownian[i];
					fout_particle << ' ' << sys.ang_vel_brownian[i];
				}
				if (sys.magnetic) {
					fout_particle << ' ' << sys.vel_magnetic[i];
					fout_particle << ' ' << sys.ang_vel_magnetic[i];
				}
				if (sys.mobile_fixed) {
					fout_particle << ' ' << sys.vel_hydro_from_fixed[i];
					fout_particle << ' ' << sys.ang_vel_hydro_from_fixed[i];
				}
			}
			fout_particle << endl;
		}
	}
	int cnt_interaction = 0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active()) {
			cnt_interaction ++;
		}
	}
	if (p.out_data_interaction) {
		fout_interaction << "# " << sys.get_shear_strain();
		fout_interaction << ' ' << cnt_interaction;
		fout_interaction << ' ' << sys.get_time();
		fout_interaction << endl;
		for (int k=0; k<sys.nb_interaction; k++) {
			if (sys.interaction[k].is_active()) {
				unsigned int i, j;
				sys.interaction[k].get_par_num(i, j);
				vec3d& nr_vec = sys.interaction[k].nvec;
				StressTensor stress_contact = sys.interaction[k].contact.getContactStressXF();
				fout_interaction << setprecision(6);
				fout_interaction << i << ' ' << j << ' '; // 1, 2
				/* contact.state:
				 * 0 no contact
				 * 1 Friction is not activated (critical load model)
				 * 2 Static friction
				 * 3 Sliding
				 */
				fout_interaction << sys.interaction[k].contact.state << ' '; //3
				if (diminish_output == false) {
					fout_interaction << nr_vec.x << ' '; // 4
					fout_interaction << nr_vec.y << ' '; // 5
					fout_interaction << nr_vec.z << ' '; // 6
					fout_interaction << sys.interaction[k].get_reduced_gap() << ' '; // 7
				} else {
					fout_interaction << "d d d d ";
				}
				/* [NOTE]
				 * Lubrication forces are reference values
				 * in the Brownian case. The force balancing
				 * velocities are recalculated without
				 * including the Brownian forces.
				 * It seems there is no better way to visualize
				 * the lubrication forces.
				 */
				fout_interaction << setprecision(output_precision) << sys.interaction[k].lubrication.get_lubforce_normal() << ' '; // 8
				fout_interaction << setprecision(output_precision) << sys.interaction[k].lubrication.get_lubforce_tan() << ' '; // 9, 10, 11
				/*
				 * Contact forces include only spring forces.
				 */
				fout_interaction << setprecision(output_precision) << sys.interaction[k].contact.get_f_contact_normal_norm() << ' '; // 12
				fout_interaction << setprecision(output_precision) << sys.interaction[k].contact.get_f_contact_tan() << ' '; // 13, 14, 15
				fout_interaction << setprecision(output_precision) << sys.interaction[k].repulsion.getForceNorm() << ' '; // 16

				if (diminish_output == false) {
					fout_interaction << 6*M_PI*shearStressComponent(stress_contact, p.theta_shear) << ' '; // 17
				} else {
					fout_interaction << "d";
				}
				fout_interaction << endl;
			}
		}
	}
}

void Simulation::outputFinalConfiguration(const string& filename_import_positions)
{
	cout << "Output final configuration" << endl;
	ofstream fout_finalconfig;
	string filename_final_configuration = "relaxed_"+filename_import_positions;
	fout_finalconfig.open(filename_final_configuration.c_str());
	fout_finalconfig << header_imported_configulation[0] << endl;
	fout_finalconfig << header_imported_configulation[1] << endl;
	int np = sys.get_np();
	if (control_var != "magnetic") {
		for (int i=0; i<np; i++) {
			fout_finalconfig << sys.position[i].x << ' ';
			fout_finalconfig << sys.position[i].y << ' ';
			fout_finalconfig << sys.position[i].z << ' ';
			fout_finalconfig << sys.radius[i] << endl;
		}
	} else {
		for (int i=0; i<np; i++) {
			fout_finalconfig << sys.position[i].x << ' ';
			fout_finalconfig << sys.position[i].y << ' ';
			fout_finalconfig << sys.position[i].z << ' ';
			fout_finalconfig << sys.radius[i] << ' ';
			if (sys.p.magnetic_type == 1) {
				fout_finalconfig << sys.magnetic_moment[i].x << ' ';
				fout_finalconfig << sys.magnetic_moment[i].y << ' ';
				fout_finalconfig << sys.magnetic_moment[i].z << ' ';
				fout_finalconfig << sys.magnetic_susceptibility[i] << endl;
			} else if (sys.p.magnetic_type == 2) {
				fout_finalconfig << "0 0 0 ";
				fout_finalconfig << sys.magnetic_susceptibility[i] << endl;
			}
		}
	}
	string filename_bin = filename_final_configuration;
	string ext = ".dat";
	size_t start_pos = filename_bin.find(ext);
	if (start_pos == string::npos) {
		cerr << " WARNING, no binary output generated " << endl;
		return;
	}
	filename_bin.replace(start_pos, ext.length(), ".bin");
	outputConfigurationBinary(filename_bin);
}
