//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012-2015 Ryohei Seto and Romain Mari. All rights reserved.
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <map>
#include <string>
#include <cstring>
#include <sstream>
#include <cctype>
#include <stdexcept>
#include "Simulation.h"
#include "SystemHelperFunctions.h"

using namespace std;

Simulation::Simulation():
sys(System(p, events)),
shear_rate_expectation(-1),
internal_unit_scales("hydro"),
target_stress_input(0),
diminish_output(false)
{
	unit_longname["h"] = "hydro";
	unit_longname["r"] = "repulsion";
	unit_longname["b"] = "brownian";
	unit_longname["c"] = "cohesion";
	unit_longname["cl"] = "critical_load";
	unit_longname["ft"] = "ft_max";
	unit_longname["kn"] = "kn";
	unit_longname["kt"] = "kt";
	unit_longname["kr"] = "kr";
	unit_longname["s"] = "stress";

	force_value_ptr["hydro"] = &dimensionless_rate; // the dimensionless hydrodynamic force is also the dimensionless shear rate
	force_value_ptr["repulsion"] = &sys.amplitudes.repulsion;
	force_value_ptr["critical_load"] = &sys.amplitudes.critical_normal_force;
	force_value_ptr["cohesion"] = &sys.amplitudes.cohesion;
	force_value_ptr["ft_max"] = &sys.amplitudes.ft_max;
	force_value_ptr["brownian"] = &sys.amplitudes.temperature;
	force_value_ptr["kn"] = &sys.p.kn;
	force_value_ptr["kt"] = &sys.p.kt;
	force_value_ptr["kr"] = &sys.p.kr;
	kill = false;
};

Simulation::~Simulation(){};

string Simulation::gitVersion(){
	return GIT_VERSION;
}

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
		sys.calcStress();
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
	if (output_events.find("data") != output_events.end()) {
		sys.calcStress();
		outputData();
	}

	if (output_events.find("config") != output_events.end()) {
		if (p.out_binary_conf) {
			string binconf_filename = "conf_" + simu_name + "_" + to_string(++binconf_counter) + ".bin";
			outputConfigurationBinary(binconf_filename);
		} else {
			outputConfigurationData();
		}
	}
}


void Simulation::setupOptionalSimulation(string indent)
{
	cerr << indent << "simulation_mode = " << sys.p.simulation_mode << endl;
	switch (sys.p.simulation_mode) {
		case 0:
			cout << indent << "basic simulation" << endl;
			break;
		case 1:
			cout << indent << "Test simulation for reversibility in mixed problem" << endl;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			break;
		case 2:
			cout << indent << "Test simulation for a mixed problem" << endl;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			break;
		case 3:
			cout << indent << "Test simulation for a mixed problem" << endl;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			break;
		case 4:
			cout << indent << "Test simulation for relax" << endl;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			sys.p.cross_shear = false;
			break;
		case 11: //ctest1
			cout << indent << "Test simulation with co-axial cylinders (rotate outer clynder)" << endl;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			sys.wall_rheology = true;
			break;
		case 12: // ctest2
			cout << indent << "Test simulation with co-axial cylinders (rotate inner clynder)" << endl;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			sys.wall_rheology = true;
			break;
		case 13: // ctest3
			cout << indent << "Test simulation with co-axial cylinders (rotate both inner and outer clynder)" << endl;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			sys.wall_rheology = true;
			break;
		case 10:
			cout << indent << "Test simulation with co-axial cylinders (rotate both inner and outer clynder)" << endl;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			sys.wall_rheology = true;
			break;
		case 21:
			cout << indent << "Test simulation for shear reversibility" << endl;
			sys.p.cross_shear = true;
			break;
		case 31:
			cout << indent << "Test simulation (wtest1), simple shear with walls" << endl;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			break;
		case 41:
			cout << indent << "Test simulation (wtestA), simple shear with walls" << endl;
			sys.wall_rheology = true;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			break;
		case 42:
			cout << indent << "Test simulation (wtestB), simple shear with walls" << endl;
			sys.wall_rheology = true;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			break;
		case 51:
			cout << indent << "Test simulation (wtestB), simple shear with walls" << endl;
			sys.wall_rheology = true;
			sys.zero_shear = true;
			sys.mobile_fixed = true;
			break;
		default:
			break;
	}
}

void Simulation::timeEvolutionUntilNextOutput(const TimeKeeper &tk)
{
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
}

void Simulation::printProgress()
{
	if (time_end != -1) {
				cout << "time: " << sys.get_time_in_simulation_units() << " , " << sys.get_time() << " / " << time_end << " , strain: " << sys.get_shear_strain() << endl;
	} else {
				cout << "time: " << sys.get_time_in_simulation_units() << " , strain: " << sys.get_shear_strain() << " / " << strain_end << endl;
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
	setupSimulation(in_args, input_files, binary_conf, dimensionless_number, input_scale, simu_identifier);
	time_t now;
	time_strain_1 = 0;
	now = time(NULL);
	time_strain_0 = now;

	setupEvents();
	cout << indent << "Time evolution started" << endl << endl;
	TimeKeeper tk = initTimeKeeper();

	int binconf_counter = 0;
	while (keepRunning()) {
		timeEvolutionUntilNextOutput(tk);

		set<string> output_events = tk.getElapsedClocks(sys.get_time(), fabs(sys.get_shear_strain()));
		generateOutput(output_events, binconf_counter);

		printProgress();

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
	if (filename_parameters.find("init_relax", 0) != string::npos) {
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
	sys.calcStress();
	outputData(); // new
	outputConfigurationBinary();
	outputConfigurationData();
	/*************************************************************/

	TimeKeeper tk;
	tk.addClock("data", LinearClock(p.time_interval_output_data,
	                                input_values["time_interval_output_data"].unit == "strain"));
	tk.addClock("config", LinearClock(p.time_interval_output_config,
	                                  input_values["time_interval_output_config"].unit == "strain"));
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
				throw runtime_error("use of updateUnscaledContactmodel() invalid");
				// sys.updateUnscaledContactmodel();
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
	if (filename_parameters.find("init_relax", 0) != string::npos) {
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
	time_t time_from_1 = time_strain_end-time_strain_1;
	time_t time_from_0 = time_strain_end-time_strain_0;
	int timestep_from_1 = timestep_end-timestep_1;
	fout_time << "# np time_from_0 time_from_1 timestep_end timestep_from_1" << endl;
	fout_time << sys.get_np() << ' ';
	fout_time << time_from_0 << ' ';
	fout_time << time_from_1 << ' ';
	fout_time << timestep_end << ' ';
	fout_time << timestep_from_1 << endl;
}

DimensionalValue Simulation::str2DimensionalValue(string type,
                                                  string name,
                                                  string value_str,
                                                  double *value_ptr)
{
	DimensionalValue inv;
	inv.type = type;
	inv.value = value_ptr;

	string numeral, suffix;
	bool caught_suffix = true;
	caught_suffix = getSuffix(value_str, numeral, suffix);
	if (!caught_suffix) {
		errorNoSuffix(name);
	}
	suffix = unit_longname[suffix];
	*(inv.value) = atof(numeral.c_str());
	inv.unit = suffix;
	return inv;
}


void Simulation::outputConfigurationBinary()
{
	string conf_filename;
	conf_filename = "conf_" + simu_name + ".dat";
	outputConfigurationBinary(conf_filename);
}

void Simulation::outputConfigurationBinary(string conf_filename)
{
	/**
		\brief Saves the current configuration of the system in a binary file.

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

	int conf_switch = -1; // older formats did not have labels, -1 signs for a labeled binary
	int binary_conf_format = 2; // v2 as default. v1 deprecated.
	if (sys.p.simulation_mode == 31) {
		binary_conf_format = 3;
	}
	conf_export.write((char*)&conf_switch, sizeof(int));
	conf_export.write((char*)&binary_conf_format, sizeof(int));

	conf_export.write((char*)&np, sizeof(int));
	if (binary_conf_format == 3) {
		int np_fixed = sys.get_np() - sys.np_mobile;
		conf_export.write((char*)&np_fixed, sizeof(int));
	}
	conf_export.write((char*)&volume_or_area_fraction, sizeof(double));
	conf_export.write((char*)&lx, sizeof(double));
	conf_export.write((char*)&ly, sizeof(double));
	conf_export.write((char*)&lz, sizeof(double));
	conf_export.write((char*)&(sys.shear_disp.x), sizeof(double));
	conf_export.write((char*)&(sys.shear_disp.y), sizeof(double));
	for (int i=0; i<np; i++) {
		conf_export.write((char*)&pos[i][0], dims*sizeof(double));
	}
	if (binary_conf_format == 3) {
		int np_fixed = sys.get_np() - sys.np_mobile;
		vector< vector<double> > vel(np_fixed);
		for (int i=0; i<np_fixed; i++) {
			vel[i].resize(3);
			vel[i][0] = sys.fixed_velocities[i].x;
			vel[i][1] = sys.fixed_velocities[i].y;
			vel[i][2] = sys.fixed_velocities[i].z;
		}
		for (int i=0; i<np_fixed; i++) {
			conf_export.write((char*)&vel[i][0], 3*sizeof(double));
		}
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
	if (force_ratios.find(dimless_nb_label) == force_ratios.end()) {
		ostringstream error_str;
		error_str << " Error : don't manage to convert from \"" << internal_unit_scales << "\" units to \"" << output_unit_scales << "\" units to output data." << endl;
		throw runtime_error(error_str.str());
	}
	outdata.setDimensionlessNumber(force_ratios[dimless_nb_label]);

	outdata.setUnit(output_unit_scales);
	double sr = sys.get_shear_rate();
	double shear_stress = shearStressComponent(sys.total_stress, p.theta_shear);
	outdata.entryData("time", "time", 1, sys.get_time());
	if (sys.get_omega_wheel() == 0 || sys.wall_rheology == false) {
		// Simple shear geometry
		outdata.entryData("shear strain", "none", 1, sys.get_shear_strain());
		outdata.entryData("shear rate", "rate", 1, sys.get_shear_rate());
	} else {
		// Rotary Couette geometry
		outdata.entryData("rotation angle", "none", 1, sys.get_angle_wheel());
		outdata.entryData("omega wheel", "rate", 1, sys.get_omega_wheel());
	}
	outdata.entryData("viscosity", "viscosity", 1, shear_stress/sr);
	for (const auto &stress_comp: sys.total_stress_groups) {
		string entry_name = "Viscosity("+stress_comp.first+")";
		outdata.entryData(entry_name, "viscosity", 1, shearStressComponent(stress_comp.second, p.theta_shear)/sr);
	}
	/*
	 * Stress
	 */
	outdata.entryData("shear stress", "stress", 1, shear_stress);
	outdata.entryData("N1 viscosity", "viscosity", 1, sys.total_stress.getNormalStress1()/sr);
	outdata.entryData("N2 viscosity", "viscosity", 1, sys.total_stress.getNormalStress2()/sr);
	outdata.entryData("particle pressure", "stress", 1, sys.total_stress.getParticlePressure());
	outdata.entryData("particle pressure contact", "stress", 1, sys.total_stress_groups["contact"].getParticlePressure());
	/* energy
	 */
	outdata.entryData("energy", "none", 1, getPotentialEnergy(sys));
	/* maximum deformation of contact bond
	 */
	outdata.entryData("min gap", "none", 1, evaluateMinGap(sys));
	if (sys.cohesion) {
		outdata.entryData("max gap(cohesion)", "none", 1, evaluateMaxContactGap(sys));
	}
	outdata.entryData("max tangential displacement", "none", 1, evaluateMaxDispTan(sys));
	outdata.entryData("max rolling displacement", "none", 1, evaluateMaxDispRolling(sys));
	/* contact number
	 */
	unsigned int contact_nb, frictional_contact_nb;
	std::tie(contact_nb, frictional_contact_nb) = countNumberOfContact(sys);
	double contact_nb_per_particle = (double)2*contact_nb/sys.get_np();
	double frictional_contact_nb_per_particle = (double)2*frictional_contact_nb/sys.get_np();
	outdata.entryData("contact number", "none", 1, contact_nb_per_particle);
	outdata.entryData("frictional contact number", "none", 1, frictional_contact_nb_per_particle);
	outdata.entryData("number of interaction", "none", 1, sys.get_nb_of_active_interactions());
	/* maximum velocity
	 */
	outdata.entryData("max velocity", "velocity", 1, sys.max_velocity);
	outdata.entryData("max angular velocity", "velocity", 1, evaluateMaxAngVelocity(sys));
	/* simulation parameter
	 */
	outdata.entryData("dt", "time", 1, sys.avg_dt);
	outdata.entryData("kn", "none", 1, p.kn);
	outdata.entryData("kt", "none", 1, p.kt);
	outdata.entryData("kr", "none", 1, p.kr);
	outdata.entryData("shear displacement x", "none", 1, sys.shear_disp.x);
	outdata.entryData("shear displacement y", "none", 1, sys.shear_disp.y);
	if (sys.wall_rheology) {
		outdata.entryData("shear viscosity wall 1", "viscosity", 1, sys.shearstress_wall1/sr);
		outdata.entryData("shear viscosity wall 2", "viscosity", 1, sys.shearstress_wall2/sr);
		outdata.entryData("normal stress/rate wall 1", "viscosity", 1, sys.normalstress_wall1/sr);
		outdata.entryData("normal stress/rate wall 2", "viscosity", 1, sys.normalstress_wall2/sr);
	}
	if (sys.p.simulation_mode == 31) {
		outdata.entryData("force top wall", "force", 3, sys.force_upwall);
		outdata.entryData("force bottom wall", "force", 3, sys.force_downwall);
	}
	outdata.writeToFile();
	/****************************   Stress Tensor Output *****************/
	outdata_st.setDimensionlessNumber(force_ratios[dimless_nb_label]);
	outdata_st.setUnit(output_unit_scales);
	outdata_st.entryData("time", "time", 1, sys.get_time());
	outdata_st.entryData("shear strain", "none", 1, sys.get_shear_strain());
	outdata_st.entryData("shear rate", "rate", 1, sys.get_shear_rate());
	outdata_st.entryData("total stress tensor (xx, xy, xz, yz, yy, zz)", "stress", 6, sys.total_stress);
	for (const auto &stress_comp: sys.total_stress_groups) {
		string entry_name = stress_comp.first+" stress tensor (xx, xy, xz, yz, yy, zz)";
		outdata_st.entryData(entry_name, "stress", 6, stress_comp.second);
	}
	outdata_st.writeToFile();

	if (!p.out_particle_stress.empty()) {
		outdata_pst.setDimensionlessNumber(force_ratios[dimless_nb_label]);
		outdata_pst.setUnit(output_unit_scales);

		cerr << "Warning, particle stress data output temporarily disabled " << endl;
		// for (int i=0; i<sys.get_np(); i++) {
		// 	if (p.out_particle_stress.find('t') != string::npos) {
		// 		StressTensor s = sys.lubstress[i]+sys.contactstressXF[i]+sys.contactstressGU[i];
		// 		if (sys.brownian) {
		// 			s += sys.brownianstressGU[i];
		// 		}
		// 		if (sys.repulsiveforce) {
		// 			s += sys.repulsivestressGU[i]+sys.repulsivestressXF[i];
		// 		}
		// 		outdata_pst.entryData("total stress (xx, xy, xz, yz, yy, zz), excluding magnetic stress", "stress", 6, s);
		// 	}
		// 	if (p.out_particle_stress.find('l') != string::npos) {
		// 		outdata_pst.entryData("lubrication stress (xx, xy, xz, yz, yy, zz)", "stress", 6, sys.lubstress[i]);
		// 	}
		// 	if (p.out_particle_stress.find('c') != string::npos) {
		// 		outdata_pst.entryData("contact stress (xx, xy, xz, yz, yy, zz)", "stress", 6, sys.contactstressXF[i] + sys.contactstressGU[i]);
		// 	}
		// 	if (sys.brownian && p.out_particle_stress.find('b') != string::npos ) {
		// 		outdata_pst.entryData("Brownian stress (xx, xy, xz, yz, yy, zz)", "stress", 6, sys.brownianstressGU[i]);
		// 	}
		// 	if (sys.repulsiveforce && p.out_particle_stress.find('r') != string::npos ) {
		// 		outdata_pst.entryData("repulsive stress (xx, xy, xz, yz, yy, zz)", "stress", 6, sys.repulsivestressGU[i] + sys.repulsivestressXF[i]);
		// 	}
		// }
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
	snapshot_header << endl;
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

void Simulation::outputParFileTxt()
{
	int np = sys.get_np();

	int output_precision = 6;
	if (diminish_output) {
		output_precision = 4;
	}
	outdata_par.setDefaultPrecision(output_precision);
	outdata_int.setDefaultPrecision(output_precision);

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
				vel[i] -= sys.vel_difference;
			}
		}
	}

	string dimless_nb_label = internal_unit_scales+"/"+output_unit_scales;
	if (force_ratios.find(dimless_nb_label) == force_ratios.end()) {
		ostringstream error_str;
		error_str << " Error : don't manage to convert from \"" << internal_unit_scales << "\" units to \"" << output_unit_scales << "\" units to output data." << endl;
		throw runtime_error(error_str.str());
	}
	cout << "   out config: " << sys.get_shear_strain() << endl;
	outdata_par.setDimensionlessNumber(force_ratios[dimless_nb_label]);
	outdata_par.setUnit(output_unit_scales);
	for (int i=0; i<sys.get_np(); i++) {
		outdata_par.entryData("particle index", "none", 1, i);
		outdata_par.entryData("radius", "none", 1, sys.radius[i]);
		outdata_par.entryData("position (x, y, z)", "none", 3, pos[i], 6);
		outdata_par.entryData("velocity (x, y, z)", "velocity", 3, vel[i]);


		outdata_par.entryData("angular velocity (x, y, z)", "velocity", 3, sys.ang_velocity[i]);
		if (sys.couette_stress) {
			double stress_rr, stress_thetatheta, stress_rtheta;
			sys.getStressCouette(i, stress_rr, stress_thetatheta, stress_rtheta);
			double sr = sys.get_shear_rate();
			outdata_par.entryData("stress_rr", "viscosity", 1, stress_rr/sr);
			outdata_par.entryData("stress_thetatheta", "viscosity", 1, stress_thetatheta/sr);
			outdata_par.entryData("stress_rtheta", "viscosity", 1, stress_rtheta/sr);
		}
		if (sys.twodimension) {
			outdata_par.entryData("angle", "none", 1, sys.angle[i]);
		}

		if (p.out_data_vel_components) {
			for (const auto &vc: sys.velocity_components) {
				string entry_name_vel = "non-affine "+vc.first+" velocity (x, y, z)";
				string entry_name_ang_vel = "non-affine angular "+vc.first+" velocity (x, y, z)";
				outdata_par.entryData(entry_name_vel, "velocity", 3, vc.second.vel[i]);
				outdata_par.entryData(entry_name_ang_vel, "velocity", 3, vc.second.ang_vel[i]);
			}
		}
	}
	stringstream snapshot_header;
	getSnapshotHeader(snapshot_header);
	outdata_par.writeToFile(snapshot_header.str());
}

void Simulation::outputIntFileTxt()
{
	string dimless_nb_label = internal_unit_scales+"/"+output_unit_scales;

	outdata_int.setDimensionlessNumber(force_ratios[dimless_nb_label]);
	outdata_int.setUnit(output_unit_scales);
	stringstream snapshot_header;
	getSnapshotHeader(snapshot_header);
	for (const auto &inter: sys.interaction) {
		unsigned int i, j;
		std::tie(i, j) = inter.get_par_num();
		StressTensor stress_contact = inter.contact.getContactStressXF();
		outdata_int.entryData("particle 1 label", "none", 1, i);
		outdata_int.entryData("particle 2 label", "none", 1, j);
		outdata_int.entryData("contact state "
		                      "(0 = no contact, "
		                      "1 = frictionless contact, "
		                      "2 = non-sliding frictional, "
		                      "3 = sliding frictional)",
		                      "none", 1, inter.contact.getFrictionState());
		if (diminish_output == false) {
			outdata_int.entryData("normal vector, oriented from particle 1 to particle 2", \
			                      "none", 3, inter.nvec);
			outdata_int.entryData("dimensionless gap = s-2, s = 2r/(a1+a2)", \
			                      "none", 1,  inter.get_reduced_gap());
		}
		/* [NOTE]
		 * Lubrication forces are reference values
		 * in the Brownian case. The force balancing
		 * velocities are recalculated without
		 * including the Brownian forces.
		 * It seems there is no better way to visualize
		 * the lubrication forces.
		 */
		if (sys.lubrication) {
			double normal_part = -dot(inter.lubrication.getTotalForce(), inter.nvec);
			outdata_int.entryData("normal part of the lubrication force (positive for compression)", "force", 1, \
			                      normal_part);
			outdata_int.entryData("tangential part of the lubrication force", "force", 3, \
			                      inter.lubrication.getTangentialForce());
		}
		/*
		 * Contact forces include only spring forces.
		 */
		outdata_int.entryData("norm of the normal part of the contact force", "force", 1, \
		                      inter.contact.getNormalForce().norm());
		outdata_int.entryData("tangential part of the contact force", "force", 3, \
		                      inter.contact.getTangentialForce());
		outdata_int.entryData("norm of the normal repulsive force", "force", 1, \
		                      inter.repulsion.getForceNorm());
		if (diminish_output == false) {
			outdata_int.entryData("Viscosity contribution of contact xF", "stress", 1, \
			                      shearStressComponent(stress_contact, p.theta_shear));
		}
	}
	outdata_int.writeToFile(snapshot_header.str());

}
void Simulation::outputConfigurationData()
{
	if (p.out_data_particle) {
		outputParFileTxt();
	}
	if (p.out_data_interaction) {
		outputIntFileTxt();
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
	for (int i=0; i<np; i++) {
		fout_finalconfig << sys.position[i].x << ' ';
		fout_finalconfig << sys.position[i].y << ' ';
		fout_finalconfig << sys.position[i].z << ' ';
		fout_finalconfig << sys.radius[i] << endl;
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
