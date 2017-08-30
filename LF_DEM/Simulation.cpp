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
#include <complex>
#include "Simulation.h"
#include "SystemHelperFunctions.h"

using namespace std;

Simulation::Simulation():
sys(System(p, events)),
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
	force_value_ptr["repulsion"] = &sys.p.repulsion;
	force_value_ptr["critical_load"] = &sys.p.critical_load;
	force_value_ptr["cohesion"] = &sys.p.cohesion;
	force_value_ptr["ft_max"] = &sys.p.ft_max;
	force_value_ptr["brownian"] = &sys.p.brownian;
	force_value_ptr["kn"] = &sys.p.kn;
	force_value_ptr["kt"] = &sys.p.kt;
	force_value_ptr["kr"] = &sys.p.kr;
	kill = false;
};

Simulation::~Simulation(){};

string Simulation::gitVersion()
{
	return GIT_VERSION;
}

bool Simulation::keepRunning()
{
	/** \brief Determine if we reached the end of the simulation.

		Returns true when ParameterSet::time_end is reached or if an event handler threw a kill signal.
	 */
	if (time_end == -1) {
		return (sys.get_cumulated_strain() < strain_end-1e-8) && !kill;
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
	if (p.disp_max < 1e-6 || sys.get_cumulated_strain() > 3.) {
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
	if (sys.p.out_bond_order_parameter6) {
		sys.calcOrderParameter();
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
	cout << indent << "simulation_mode = " << sys.p.simulation_mode << endl;
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
		cout << "time: " << sys.get_time_in_simulation_units() << " , "\
		     << sys.get_time() << " / " << time_end\
		     << " , strain: " << sys.get_cumulated_strain() << endl;
	} else {
		cout << "time: " << sys.get_time_in_simulation_units()\
		     << " , strain: " << sys.get_cumulated_strain() << " / " << strain_end << endl;
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
									   ControlVariable control_variable,
									   string flow_type,
									   string simu_identifier)
{
	string indent = "  Simulation::\t";
	if (flow_type == "extension") {
		sys.ext_flow = true;
	} else {
		sys.ext_flow = false;
	}
	control_var = control_variable;
	setupSimulation(in_args, input_files, binary_conf, dimensionless_number, input_scale,
					flow_type, simu_identifier);
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

		set<string> output_events = tk.getElapsedClocks(sys.get_time(), sys.get_cumulated_strain());
		if (sys.retrim_ext_flow) {
			output_events.insert("data");
			output_events.insert("config");
		}
		generateOutput(output_events, binconf_counter);

		printProgress();

		if (time_strain_1 == 0 && sys.get_cumulated_strain() > 1) {
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
										ControlVariable control_variable,
										string flow_type,
										string simu_identifier)
{
	control_var = control_variable;
	setupSimulation(in_args, input_files, binary_conf, dimensionless_number, input_scale, flow_type, simu_identifier);

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
		set<string> output_events = tk.getElapsedClocks(sys.get_time(), sys.get_cumulated_strain());
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
		if (time_strain_1 == 0 && sys.get_cumulated_strain() > 1) {
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
	auto conf = sys.getConfiguration();
	int np = conf.position.size();
	vector< vector<double> > pos(np);
	int dims = 4;
	for (int i=0; i<np; i++) {
		pos[i].resize(dims);
		pos[i][0] = conf.position[i].x;
		pos[i][1] = conf.position[i].y;
		pos[i][2] = conf.position[i].z;
		pos[i][3] = conf.radius[i];
	}
	ofstream conf_export;
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
		int np_fixed = sys.fixed_velocities.size();
		conf_export.write((char*)&np_fixed, sizeof(int));
	}
	conf_export.write((char*)&conf.volume_or_area_fraction, sizeof(double));
	conf_export.write((char*)&conf.lx, sizeof(double));
	conf_export.write((char*)&conf.ly, sizeof(double));
	conf_export.write((char*)&conf.lz, sizeof(double));
	conf_export.write((char*)&(conf.lees_edwards_disp.x), sizeof(double));
	conf_export.write((char*)&(conf.lees_edwards_disp.y), sizeof(double));
	for (int i=0; i<np; i++) {
		conf_export.write((char*)&pos[i][0], dims*sizeof(double));
	}
	if (binary_conf_format == 3) {
		int np_fixed = sys.fixed_velocities.size();
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
	const auto &cs = conf.contact_states;
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
	if (control_var == rate) {
		return input_rate;
	} else if (control_var == stress) {
		return sys.get_shear_rate();
	} else {
		return 1;
	}
}

vector<Sym2Tensor> Simulation::getParticleStressGroup(string group)
{
	vector<Sym2Tensor> S (sys.get_np(), 0);
	for (const auto &sc: sys.stress_components) {
		if (sc.second.group == group || group=="total") {
			for (size_t i=0; i<S.size(); i++) {
				S[i] += sc.second.particle_stress[i];
			}
		}
	}
	return S;
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
	double viscous_material_function, inviscid_material_function0, inviscid_material_function3;
	if (sr != 0) {
		// generalized viscosity kappa (= 2*eta)
		viscous_material_function   = doubledot(sys.total_stress, sys.getEinfty())/ sys.getEinfty().selfdoubledot();
		inviscid_material_function0 = doubledot(sys.total_stress, stress_basis_0) / stress_basis_0.selfdoubledot();
		inviscid_material_function3 = doubledot(sys.total_stress, stress_basis_3) / stress_basis_3.selfdoubledot();
	} else {
		// @@@ tentative ouptut for Pe = 0 simulation
		// output xz component of stress tensor
		viscous_material_function = sys.total_stress.elm[2];
	}
	outdata.entryData("time", "time", 1, sys.get_time());
	if (sys.get_omega_wheel() == 0 || sys.wall_rheology == false) {
		// Simple shear geometry
		outdata.entryData("cumulated shear strain", "none", 1, sys.get_cumulated_strain());
		outdata.entryData("shear rate", "rate", 1, sys.get_shear_rate());
	} else {
		// Rotary Couette geometry
		outdata.entryData("rotation angle", "none", 1, sys.get_angle_wheel());
		outdata.entryData("omega wheel", "rate", 1, sys.get_omega_wheel());
	}
	outdata.entryData("viscosity", "viscosity", 1, 0.5*viscous_material_function);
	for (const auto &stress_comp: sys.total_stress_groups) {
		string entry_name = "Viscosity("+stress_comp.first+")";
		double viscous_mf_component = doubledot(stress_comp.second, sys.getEinfty())/sys.getEinfty().selfdoubledot();
		outdata.entryData(entry_name, "viscosity", 1, 0.5*viscous_mf_component);
	}
	/*
	 * Stress
	 */
	//outdata.entryData("shear stress", "stress", 1, shear_stress);
	auto stress_diag = sys.total_stress.diag();
	outdata.entryData("inviscid function 0th", "viscosity", 1, inviscid_material_function0);
	outdata.entryData("inviscid function 3rd", "viscosity", 1, inviscid_material_function3);
	outdata.entryData("N1 viscosity", "viscosity", 1, (stress_diag.x-stress_diag.z)/sr);
	outdata.entryData("N2 viscosity", "viscosity", 1, (stress_diag.z-stress_diag.y)/sr);
	outdata.entryData("particle pressure", "stress", 1, -sys.total_stress.trace()/3);
	outdata.entryData("particle pressure contact", "stress", 1, -sys.total_stress_groups["contact"].trace()/3);
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
	outdata.entryData("number of interaction", "none", 1, sys.get_nb_interactions());
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
	vec3d shear_strain = sys.get_shear_strain();
	outdata.entryData("shear strain", "none", 3, shear_strain);
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
	if (sys.brownian) {
		outdata.entryData("max_velocity_brownian", "velocity", 1, sys.max_velocity_brownian);
		outdata.entryData("max_velocity_contact", "velocity", 1, sys.max_velocity_contact);
	}
	outdata.writeToFile();
	/****************************   Stress Tensor Output *****************/
	outdata_st.setDimensionlessNumber(force_ratios[dimless_nb_label]);
	outdata_st.setUnit(output_unit_scales);
	outdata_st.entryData("time", "time", 1, sys.get_time());
	outdata_st.entryData("cumulated strain", "none", 1, sys.get_cumulated_strain());
	outdata_st.entryData("shear rate", "rate", 1, sys.get_shear_rate());
	outdata_st.entryData("total stress tensor (xx, xy, xz, yz, yy, zz)", "stress", 6, sys.total_stress);
	for (const auto &stress_comp: sys.total_stress_groups) {
		string entry_name = stress_comp.first+" stress tensor (xx, xy, xz, yz, yy, zz)";
		outdata_st.entryData(entry_name, "stress", 6, stress_comp.second);
	}
	outdata_st.writeToFile();
}

void Simulation::getSnapshotHeader(stringstream& snapshot_header)
{
	snapshot_header << "# "; //1
	snapshot_header << sys.get_cumulated_strain() << ' ';//2
	snapshot_header << sys.shear_disp.x << ' ';//3
	snapshot_header << getRate() << ' ';//4
	snapshot_header << target_stress_input << ' ';//5
	snapshot_header << sys.get_cumulated_strain() << ' ';//6
	snapshot_header << sys.get_cumulated_strain()-sys.strain_retrim+sys.strain_retrim_interval << ' ';//7
	snapshot_header << sys.get_shear_rate() << ' '; //8
	if (sys.ext_flow) {
		if (sys.retrim_ext_flow) {
			snapshot_header << "retrim" << ' '; //9
		}
	}
	snapshot_header << endl;
}

vec3d Simulation::shiftUpCoordinate(double x, double y, double z)
{
	if (p.origin_zero_flow) {
		z += 0.5*sys.get_lz();
		if (z > 0.5*sys.get_lz()) {
			x -= sys.shear_disp.x;
			y -= sys.shear_disp.y;
			if (x < -0.5*sys.get_lx()) {
				x += sys.get_lx();
			}
			if (y < -0.5*sys.get_ly()) {
				y += sys.get_ly();
			}
			z -= sys.get_lz();
		}
	}
	return vec3d(x,y,z);
}

void Simulation::createDataHeader(stringstream& data_header)
{
	auto conf = sys.getConfiguration();
	data_header << "# LF_DEM version " << GIT_VERSION << endl;
	data_header << "# np " << conf.position.size() << endl;
	data_header << "# VF " << conf.volume_or_area_fraction << endl;
	data_header << "# Lx " << conf.lx << endl;
	data_header << "# Ly " << conf.ly << endl;
	data_header << "# Lz " << conf.lz << endl;
	data_header << "# flow_type " << sys.p.flow_type << endl;
}

void Simulation::outputDataHeader(ofstream& fout)
{
	stringstream data_header;
	createDataHeader(data_header);
	fout << data_header.str();
}

void Simulation::outputPstFileTxt()
{
	string dimless_nb_label = internal_unit_scales+"/"+output_unit_scales;
	outdata_pst.setDimensionlessNumber(force_ratios[dimless_nb_label]);
	outdata_pst.setUnit(output_unit_scales);

	map<string, string> group_shorts;
	group_shorts["l"] = "hydro";
	group_shorts["c"] = "contact";
	group_shorts["r"] = "repulsion";
	group_shorts["b"] = "brownian";
	group_shorts["d"] = "dashpot";
	group_shorts["t"] = "total";
	map<string, vector<Sym2Tensor>> particle_stress;
	for (auto &type: p.out_particle_stress) {
		auto group_name = group_shorts[string(1, type)];
		particle_stress[group_name] = getParticleStressGroup(group_name);
	}
	for (int i=0; i<sys.get_np(); i++) {
		for (const auto &pst: particle_stress) {
			outdata_pst.entryData(pst.first + " stress (xx, xy, xz, yz, yy, zz)", "stress", 6, pst.second[i]);
		}
	}
	stringstream snapshot_header;
	getSnapshotHeader(snapshot_header);
	outdata_pst.writeToFile(snapshot_header.str());
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
	if (!sys.ext_flow) {
		for (int i=0; i<np; i++) {
			pos[i] = shiftUpCoordinate(sys.position[i].x-0.5*sys.get_lx(),
									   sys.position[i].y-0.5*sys.get_ly(),
									   sys.position[i].z-0.5*sys.get_lz());
		}
	} else {
		for (int i=0; i<np; i++) {
			pos[i] = shiftUpCoordinate(sys.position[i].x,
									   sys.position[i].y,
									   sys.position[i].z);
		}
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
	cout << "   out config: " << sys.get_cumulated_strain() << endl;
	outdata_par.setDimensionlessNumber(force_ratios[dimless_nb_label]);
	outdata_par.setUnit(output_unit_scales);
	auto na_disp = sys.getNonAffineDisp();
	for (int i=0; i<sys.get_np(); i++) {
		outdata_par.entryData("particle index", "none", 1, i);
		outdata_par.entryData("radius", "none", 1, sys.radius[i]);
		outdata_par.entryData("position (x, y, z)", "none", 3, pos[i], 6);
		outdata_par.entryData("velocity (x, y, z)", "velocity", 3, vel[i]);
		outdata_par.entryData("angular velocity (x, y, z)", "none", 3, sys.ang_velocity[i]);
		if (!sys.ext_flow) {
			outdata_par.entryData("non affine displacement (x, y, z)", "none", 3, na_disp[i]);
		}
		if (sys.twodimension) {
			outdata_par.entryData("angle", "none", 1, sys.angle[i]);
		}
		//		if (sys.couette_stress) {
		//			double stress_rr, stress_thetatheta, stress_rtheta;
		//			sys.getStressCouette(i, stress_rr, stress_thetatheta, stress_rtheta);
		//			double sr = sys.get_shear_rate();
		//			outdata_par.entryData("stress_rr", "viscosity", 1, stress_rr/sr);
		//			outdata_par.entryData("stress_thetatheta", "viscosity", 1, stress_thetatheta/sr);
		//			outdata_par.entryData("stress_rtheta", "viscosity", 1, stress_rtheta/sr);
		//		}
		if (p.out_data_vel_components) {
			for (const auto &vc: sys.na_velo_components) {
				string entry_name_vel = "non-affine "+vc.first+" velocity (x, y, z)";
				string entry_name_ang_vel = "non-affine angular "+vc.first+" velocity (x, y, z)";
				outdata_par.entryData(entry_name_vel, "velocity", 3, vc.second.vel[i]);
				outdata_par.entryData(entry_name_ang_vel, "velocity", 3, vc.second.ang_vel[i]);
			}
		}
		if (p.out_bond_order_parameter6) {
			outdata_par.entryData("abs_phi6", "none", 1, abs(sys.phi6[i]));
			outdata_par.entryData("arg_phi6", "none", 1, arg(sys.phi6[i]));
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
	double sr = sys.get_shear_rate();
	for (const auto &inter: sys.interaction) {
		unsigned int i, j;
		std::tie(i, j) = inter.get_par_num();
		Sym2Tensor stress_contact = inter.contact.getContactStressXF();
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
			if (inter.get_reduced_gap() > 0) {
				double normal_part = -dot(inter.lubrication.getTotalForce(), inter.nvec);
				outdata_int.entryData("normal part of the lubrication force (positive for compression)", "force", 1, \
									  normal_part);
				outdata_int.entryData("tangential part of the lubrication force", "force", 3, \
									  inter.lubrication.getTangentialForce());
			} else {
				outdata_int.entryData("normal part of the lubrication force (positive for compression)", "force", 1, 0);
				outdata_int.entryData("tangential part of the lubrication force", "force", 3, vec3d(0,0,0));
			}
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
								  doubledot(stress_contact, sys.getEinfty()/sr)/sr);
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
	if (!p.out_particle_stress.empty()) {
		outputPstFileTxt();
	}
	//if (sys.ext_flow) {
	//	sys.yaplotBoxing(fout_boxing); // for debugging.
	//}
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
