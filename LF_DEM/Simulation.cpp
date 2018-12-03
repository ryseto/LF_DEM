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

Simulation::Simulation(State::BasicCheckpoint chkp):
sys(System(p, events, chkp)),
target_stress_input(0),
restart_from_chkp(false),
stress_reversal(false),
timestep_1(0),
diminish_output(false)
{
	kill = false;
	restart_from_chkp = !isZeroTimeChkp(chkp);
};

string Simulation::gitVersion()
{
	return GIT_VERSION;
}

bool Simulation::keepRunning()
{
	/** \brief Determine if we reached the end of the simulation.

		Returns true when ParameterSet::time_end is reached or if an event handler threw a kill signal.
	 */
	if (sys.p.time_end.dimension == Dimensional::Dimension::Strain) {
		return (sys.get_cumulated_strain() < sys.p.time_end.value-1e-8) && !kill;
	} else {
		return (sys.get_time() < sys.p.time_end.value-1e-8) && !kill;
	}
}

void Simulation::setupEvents()
{
	/** \brief Set up the types of events to be watched by the System class.

		Links System::eventLookUp to a specialized function according to the value of ParameterSet::event_handler .
	 */
	if (sys.p.event_handler == "shear_jamming") {
		sys.eventLookUp = &System::eventShearJamming;
	} else if (sys.p.event_handler == "fragility") {
		sys.eventLookUp = &System::eventShearJamming;
	} else if (sys.p.event_handler == "jamming_stress_reversal") {
		sys.eventLookUp = &System::eventShearJamming;
		ifstream sj_program_file;
		sj_program_file.open(sys.p.sj_program_file.c_str());
		int sj_stress;
		double sj_duration;
		while(sj_program_file >> sj_stress >> sj_duration) {
			cerr << "stress program, sj_duration : ";
			cerr << sj_stress << " , " << sj_duration << endl;
			sj_program_stress.push_back(sj_stress);
			sj_program_duration.push_back(sj_duration);
		}
		sj_program_file.close();
	} else {
		sys.eventLookUp = NULL;
	}
}

void Simulation::handleEventsShearJamming()
{
	/** \brief Event handler to test for shear jamming
	 
		When a negative_shear_rate event is thrown, p.disp_max is decreased.
	 If p.disp_max is below a minimal value, the shear direction is switched to y-shear.
	 */
	bool ending_simulation = false;
	if (sys.p.fixed_dt == false) {
		for (const auto& ev : events) {
			if (ev.type == "jammed_shear_rate") {
				cout << " jammed rate " << endl;
				sys.p.disp_max /= sys.p.sj_disp_max_shrink_factor;
			}
		}
		if (sys.p.disp_max < sys.p.sj_disp_max_goal) {
			ending_simulation = true;
		}
	} else {
		exit(1);
//		static int shear_jam_counter = 0;
//		double sr = sqrt(2*sys.getEinfty().selfdoubledot()); // shear rate for simple shear.
//		if (abs(sr) < sys.p.shear_jamming_rate) {
//			shear_jam_counter ++;
//			cerr << "shear_jam_counter = " << shear_jam_counter << endl;
//		} else {
//			shear_jam_counter = 0;
//		}
//		if (shear_jam_counter == sys.p.shear_jamming_max_count) {
//			ending_simulation = true;
//		}
	}
	if (ending_simulation == true) {
		cout << "jammed" << endl;
		sys.calcStress();
		outputData();
		outputConfigurationData();
		checkpoint();
		kill = true;
	}
}

void Simulation::handleEventsFragility()
{
	/** \brief Event handler to test for shear jamming

		When a jammed_shear_rate event is thrown, p.disp_max is decreased.
	 If p.disp_max is below a minimal value, the shear direction is switched to y-shear.
	 */
	for (const auto& ev : events) {
		if (ev.type == "jammed_shear_rate") {
			cout << " jammed rate " << endl;
			sys.p.disp_max /= sys.p.sj_disp_max_shrink_factor;
		}
	}
	if (sys.p.disp_max < sys.p.sj_disp_max_goal || sys.get_cumulated_strain() > 3.) {
		sys.p.disp_max = p_initial.disp_max;
		cout << "Event Fragility : starting cross shear" << endl;
	}
}

void Simulation::handleEventsJammingStressReversal()
{
	//	double sr = sqrt(2*sys.getEinfty().selfdoubledot()); // shear rate for simple shear.
	stress_reversal = false;
	for (const auto& ev : events) {
		if (ev.type == "jammed_shear_rate") {
			if (sys.p.fixed_dt) {
				stress_reversal = true;
			} else {
				/* @@ I think this is not effective way to handle jamming.
				 */
				sys.p.disp_max /= sys.p.sj_disp_max_shrink_factor;
				if (sys.p.disp_max < sys.p.sj_disp_max_goal) {
					stress_reversal = true;
				}
			}
		}
	}
}

void Simulation::handleEvents()
{
	/** \brief Handle the list of events that appears in the previous time step

		This function dispatches to specialized handlers according to the value of ParameterSet::event_handler .
	 */
	if (sys.p.event_handler == "shear_jamming") {
		handleEventsShearJamming();
	}
	if (sys.p.event_handler == "fragility") {
		handleEventsFragility();
	}
	if (sys.p.event_handler == "jamming_stress_reversal") {
		handleEventsJammingStressReversal();
	}
	events.clear();
}

void Simulation::generateOutput(const set<string> &output_events, int& binconf_counter)
{
	checkpoint(); // generic, for recovery if crash
	if (output_events.find("data") != output_events.end()) {
		sys.calcStress();
		outputData();
	}
	if (output_events.find("config") != output_events.end()) {
		if (sys.p.output.out_binary_conf) {
			string binconf_filename = "conf_" + simu_name + "_" + to_string(++binconf_counter) + ".bin";
			outputConfigurationBinary(binconf_filename);
			outputConfigurationData();
		} else {
			outputConfigurationData();
		}
		if (sys.p.output.out_gsd) {
			outputGSD();
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
			cout << indent << "Stress reversal test (fragility of shear jamming)" << endl;
			cout << indent << "stress is reversed once jammed" << endl;
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
		case 22:
			cout << indent << "Test simulation for shear stop" << endl;
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
#ifdef SIGINT_CATCH
	if (sig_caught == SIGINT) {
		checkpoint();
		std::cerr << "exiting after SIGINT " << std::endl;
		exit(99);
	}
#endif
	handleEvents();
}

void Simulation::printProgress()
{
	Dimensional::DimensionalQty<double> current_time = {Dimensional::Dimension::Time, sys.get_time(), system_of_units.getInternalUnit()};
	system_of_units.convertFromInternalUnit(current_time, output_unit);
	if (sys.p.time_end.dimension == Dimensional::Dimension::Time) {
		cout << "time: " << current_time.value << " / " << sys.p.time_end.value << " , strain: " << sys.get_cumulated_strain() << endl;
	} else {
		cout << "time: " << current_time.value << " , strain: " << sys.get_cumulated_strain() << " / " << sys.p.time_end.value << ' ';
		if (sys.p.sj_program_file != "") {
			cout << "sjp " << sj_program_stress.front();
		}
		cout << endl;
	}
}

/*
 * Main simulation
 * @@@ Discrepancy between the name and functions.
 * @@@ (Flow rates in stress control simulations are not steady)
 */
void Simulation::simulationSteadyShear(string in_args,
									   vector<string>& input_files,
									   bool binary_conf,
									   Parameters::ControlVariable control_variable,
									   Dimensional::DimensionalQty<double> control_value,
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
	setupSimulation(in_args, input_files, binary_conf, control_value,
	                flow_type, simu_identifier);
	time_t now;
	time_strain_1 = 0;
	now = time(NULL);
	time_strain_0 = now;
	setupEvents();
	cout << indent << "Time evolution started" << endl << endl;
	TimeKeeper tk = initTimeKeeper();
	if (restart_from_chkp) {
		std::set<std::string> elapsed;
		do {
			elapsed = tk.getElapsedClocks(sys.get_time(), sys.get_cumulated_strain());
		}
		while (!elapsed.empty()); // flush tk to not output on first time step
	}
	int binconf_counter = 0;
	if (sys.p.sj_program_file != "") {
		stressProgram();
	}
	while (keepRunning()) {
		if (sys.p.simulation_mode == 22) {
			stopShearing(tk);
			if (sys.get_time() > 20) {
				break;
			}
		}
		timeEvolutionUntilNextOutput(tk);
		set<string> output_events = tk.getElapsedClocks(sys.get_time(), sys.get_cumulated_strain());
		if (sys.retrim_ext_flow) {
			output_events.insert("data");
			output_events.insert("config");
		}
		if (sys.get_time()-time_last_sj_program > sj_duration_min) {
			if (sys.p.event_handler == "jamming_stress_reversal") {
				operateJammingStressReversal(output_events);
			}
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

void Simulation::operateJammingStressReversal(std::set<std::string> &output_events)
{
	if (stress_reversal) {
		if (sys.p.sj_program_file == "") {
			jamming_strain = sys.get_cumulated_strain();
			static int cnt_shear_jamming_repetation = 0;
			cnt_shear_jamming_repetation ++;
			if (cnt_shear_jamming_repetation > sys.p.sj_reversal_repetition) {
				kill = true;
			}
			stressReversal();
			sys.p.disp_max = p_initial.disp_max;
			sys.dt = sys.p.dt;
		} else {
			stressProgram();
			jamming_strain = sys.get_cumulated_strain();
		}
		output_events.insert("data");
		output_events.insert("config");
	} else {
		jamming_strain = 0;
		if (sys.get_shear_rate() < 100*sys.p.sj_shear_rate) {
			if (sys.dt > sys.p.dt_jamming) {
				sys.dt = sys.p.dt_jamming;
			}
		} else {
			if (sys.p.sj_program_file == "") {
				sys.dt = sys.p.dt;
			}
		}
	}
}

void Simulation::stopShearing(TimeKeeper &tk)
{
	static bool initial_shearing = true;
	double strain_to_stop;
	if (!sys.ext_flow) {
		strain_to_stop = 2;
	} else {
		strain_to_stop = sys.strain_retrim_interval;
	}
	if (sys.get_cumulated_strain() >= strain_to_stop-1e-8) {
		sys.zero_shear = true;
		sys.set_shear_rate(0);
		if (!sys.ext_flow) {
			// simple shear
			Sym2Tensor Einf_common = {0, 0, 0, 0, 0, 0};
			vec3d Omegainf(0, 0, 0);
			sys.setImposedFlow(Einf_common, Omegainf);
			sys.setVelocityDifference();
		} else {
			// extensional flow
			sys.vel_difference.reset();
			sys.grad_u.set_zero();
		}
		if (initial_shearing) {
			cerr << "Stop shear at " << sys.get_cumulated_strain() << endl;
			tk.removeClock();
			tk.addClock("data", LogClock(sys.get_time()+sys.dt, sys.get_time()+1, 100, false));
			tk.addClock("config", LogClock(sys.get_time()+sys.dt, sys.get_time()+1, 100, false));
			initial_shearing = false;
		}
	}
}

void Simulation::stressReversal()
{
	jamming_strain = sys.get_cumulated_strain();
	static int shear_direction = 0;
	double theta_shear = (shear_direction % 2) ? 0 : M_PI;
	sys.setShearDirection(theta_shear);
	shear_direction ++;
	sys.reset_cumulated_strain();
}

void Simulation::stressProgram()
{
	/* Return true when the program is accepted.
	 */
	static bool first_time = true;
	static double stress_original;
	static double sj_rate_original;
	static double kn_original;
	static double kt_original;
	if (first_time) {
		first_time = false;
		stress_original = sys.target_stress;
		sj_rate_original = sys.p.sj_shear_rate;
		kn_original = sys.p.kn;
		kt_original = sys.p.kt;
	}
	cerr << " shear jamming stress program " << sj_program_stress.front() << endl;
	if (sj_program_stress.front() == 0) {
		double infinitesimal_stress = 1e-4;
		sys.target_stress = infinitesimal_stress*stress_original;
		sys.p.kn = kn_original*infinitesimal_stress;
		sys.p.kt = kt_original*infinitesimal_stress;
	} else if (sj_program_stress.front() == 1) {
		sys.setShearDirection(0);
		sys.target_stress = stress_original;
		sys.p.kn = kn_original;
		sys.p.kt = kt_original;
	} else if (sj_program_stress.front() == -1) {
		sys.setShearDirection(M_PI);
		sys.target_stress = stress_original;
		sys.p.kn = kn_original;
		sys.p.kt = kt_original;
	} else if (sj_program_stress.front() == 999) {
		kill = true;
	} else {
		cerr << "Only 1, 0, -1, 999 in the program file\n";
		exit(1);
	}
	sj_duration_min = sj_program_duration.front();
	sys.resetContactModelParameer();
	sys.calculateForces();
	time_last_sj_program = sys.get_time();
	sj_program_stress.pop_front();
	sj_program_duration.pop_front();
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

void Simulation::outputConfigurationBinary(string conf_filename)
{
	/**
		\brief Saves the current configuration of the system in a binary file.

	 */

	ConfFileFormat binary_conf_format = ConfFileFormat::bin_format_base_shear; // v2 as default. v1 deprecated.
	if (sys.p.simulation_mode == 31) {
		binary_conf_format = ConfFileFormat::bin_format_fixed_vel_shear;
	}
	if (sys.delayed_adhesion) {
		binary_conf_format = ConfFileFormat::bin_delayed_adhesion;
	}
	outputBinaryConfiguration(sys, conf_filename, binary_conf_format);
}

void Simulation::checkpoint()
{
	string conf_filename;
	conf_filename = "conf_" + simu_name + ".dat";
	outputConfigurationBinary(conf_filename);
	string state_filename;
	state_filename = "chk_" + simu_name + ".dat";
	State::outputStateBinary(state_filename, sys);
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
	outdata.setUnits(system_of_units, output_unit);
	double sr = sqrt(2*sys.getEinfty().selfdoubledot()); // shear rate for simple shear.
	if (sys.p.output.effective_coordination_number) {
		sys.countContactNumber();
	}
	outdata.entryData("time", Dimensional::Dimension::Time, 1, sys.get_time());
	if (sys.wall_rheology == false || sys.get_omega_wheel() == 0) {
		// Simple shear geometry
		outdata.entryData("cumulated shear strain", Dimensional::Dimension::none, 1, sys.get_cumulated_strain());
		/* Note: shear rate
		 * shear rate = 2*sqrt(sys.getEinfty().selfdoubledot()/2);
		 * - shear rate (\dot{\gamma}) in simple shear flow
		 * - twice of extensional rate (\dot{\varepsilon})in extensional flow
		 */
		outdata.entryData("shear rate", Dimensional::Dimension::Rate, 1, sys.get_shear_rate());
	} else {
		// Rotary Couette geometry
		outdata.entryData("rotation angle", Dimensional::Dimension::none, 1, sys.get_angle_wheel());
		outdata.entryData("omega wheel", Dimensional::Dimension::Rate, 1, sys.get_omega_wheel());
	}
	/************** viscosity **********************************************************************/
	if (sr != 0) {
		viscosity = 0.5*doubledot(sys.total_stress, sys.getEinfty())/sys.getEinfty().selfdoubledot();
	} else {
		// @@@ tentative ouptut for Pe = 0 simulation
		// output xz component of stress tensor
		//viscous_material_function = sys.total_stress.elm[2];
		viscosity = 0.5*doubledot(sys.total_stress, Einf_base)/ Einf_base.selfdoubledot();
		/* D = ((0, 0, 1/2), (0, 0, 0), (1/2, 0, 0))
		 */
	}
	outdata.entryData("viscosity", Dimensional::Dimension::Viscosity, 1, viscosity);
	for (const auto &stress_comp: sys.total_stress_groups) {
		string entry_name = "Viscosity("+stress_comp.first+")";
		double viscosity_component = 0.5*doubledot(stress_comp.second, sys.getEinfty())/sys.getEinfty().selfdoubledot();
		outdata.entryData(entry_name, Dimensional::Dimension::Viscosity, 1, viscosity_component);
	}
	//outdata.entryData("shear stress", Dimensional::Dimension::Stress, 1, shear_stress);
	auto stress_diag = sys.total_stress.diag();
	/************** isotropic stress (particle pressure) **************************************/
	outdata.entryData("particle pressure", Dimensional::Dimension::Stress, 1, -sys.total_stress.trace()/3);
	outdata.entryData("particle pressure contact", Dimensional::Dimension::Stress, 1, -sys.total_stress_groups["contact"].trace()/3);
	/************** normal stress anisotropy  *************************************************/
	if (sys.p.output.new_material_functions) {
		/************** material function lambda0 *********************************************
		 * Anisotropy between isotropic stress in the flow plane and the out-of-plane normal stress.
		 * lambda0 has a better physical meaning than N2.
		 * lambda0 = - (2/3)*(N2 + 0.5*N1)/shear_rate = (2/3)*N0/shear_rate N0 = 1.5 lambda0
		 * N2 = rate*(-1.5*lambda0 + lambda3)
		 * lambda0*ep_dot / p = lambda0*gamma_dot / 2p = - (N2 + 0.5*N1)* / 3 p
		 **************************************************************************************/
		double mf_inplane_pressure; // lambda_0
		if (sr != 0) {
			/* E = ((-1/4, 0, 0), (0, 1/2, 0), (0, 0, -1/4))
			 */
			mf_inplane_pressure = 0.5*doubledot(sys.total_stress, stress_basis_0)/stress_basis_0.selfdoubledot();
		} else {
			mf_inplane_pressure = 0;
		}
		outdata.entryData("inviscid function 0th", Dimensional::Dimension::Viscosity, 1, mf_inplane_pressure);
		for (const auto &stress_comp: sys.total_stress_groups) {
			string entry_name = "inviscid function 0th("+stress_comp.first+")";
			double mf_inplane_pressure_component = 0.5*doubledot(stress_comp.second, stress_basis_0)/stress_basis_0.selfdoubledot();
			outdata.entryData(entry_name, Dimensional::Dimension::Viscosity, 1, mf_inplane_pressure_component);
		}
		/************** material function lambda3 ****************************************
		 * lambda3 induces a reoreientation of the stress eigenvectors.
		 * lambda3 is equivalent to N1.
		 * lambda3 = - (1/2)*N1/shear_rate
		 * N1 = -2*shear_rate*lambda_3
		 *********************************************************************************/
		double mf_reorientation; // lambda_3
		if (sr != 0) {
			mf_reorientation = 0.5*doubledot(sys.total_stress, stress_basis_3)/stress_basis_3.selfdoubledot();
			/* G = ((-1/2, 0, 0), (0, 0, 0), (0, 0, 1/2)
			 */
		} else {
			mf_reorientation = 0;
		}
		outdata.entryData("inviscid function 3rd", Dimensional::Dimension::Viscosity, 1, mf_reorientation);
		for (const auto &stress_comp: sys.total_stress_groups) {
			string entry_name = "inviscid function 3rd("+stress_comp.first+")";
			double mf_reorientation_component = 0.5*doubledot(stress_comp.second, stress_basis_3)/stress_basis_3.selfdoubledot();
			outdata.entryData(entry_name, Dimensional::Dimension::Viscosity, 1, mf_reorientation_component);
		}
	} else {
		/************** Normal stress differences **************************************
		 * N1 = sigma11 - sigma22
		 * N2 = sigma22 - sigma33
		 *******************************************************************************/
		normal_stress_diff1 = (stress_diag.x-stress_diag.z)/sr;
		normal_stress_diff2 = (stress_diag.z-stress_diag.y)/sr;
		outdata.entryData("N1 viscosity", Dimensional::Dimension::Viscosity, 1, normal_stress_diff1);
		outdata.entryData("N2 viscosity", Dimensional::Dimension::Viscosity, 1, normal_stress_diff2);
	}
	/***************************************************************************************************************/
	/* energy
	 */
	outdata.entryData("energy", Dimensional::Dimension::none, 1, getPotentialEnergy(sys));
	/* maximum deformation of contact bond
	 */
	outdata.entryData("min gap", Dimensional::Dimension::none, 1, evaluateMinGap(sys));
	if (sys.cohesion) {
		outdata.entryData("max gap(cohesion)", Dimensional::Dimension::none, 1, evaluateMaxContactGap(sys));
	}
	if (sys.friction) {
		outdata.entryData("max tangential displacement", Dimensional::Dimension::none, 1, evaluateMaxDispTan(sys));
	}
	if (sys.rolling_friction) {
		outdata.entryData("max rolling displacement", Dimensional::Dimension::none, 1, evaluateMaxDispRolling(sys));
	}
	/* contact number
	 */
	unsigned int contact_nb, frictional_contact_nb;
	std::tie(contact_nb, frictional_contact_nb) = countNumberOfContact(sys);
	double contact_nb_per_particle = (double)2*contact_nb/sys.get_np();
	double frictional_contact_nb_per_particle = (double)2*frictional_contact_nb/sys.get_np();
	outdata.entryData("contact number", Dimensional::Dimension::none, 1, contact_nb_per_particle);
	outdata.entryData("frictional contact number", Dimensional::Dimension::none, 1, frictional_contact_nb_per_particle);
	outdata.entryData("number of interaction", Dimensional::Dimension::none, 1, sys.get_nb_interactions());
	if (sys.delayed_adhesion) {
		unsigned active_nb;
		double active_ratio;
		std::tie(active_nb, active_ratio) = getTAAdhesionActivityStatistics(sys);
		outdata.entryData("delayed adhesion active nb", Dimensional::Dimension::none, 1, active_nb);
		outdata.entryData("delayed adhesion active ratio", Dimensional::Dimension::none, 1, active_ratio);
	}
	/* maximum velocity
	 */
	outdata.entryData("max velocity", Dimensional::Dimension::Velocity, 1, sys.max_velocity);
	outdata.entryData("max angular velocity", Dimensional::Dimension::Velocity, 1, evaluateMaxAngVelocity(sys));
	/* simulation parameter
	 */
	outdata.entryData("dt", Dimensional::Dimension::Time, 1, sys.avg_dt);
	outdata.entryData("kn", Dimensional::Dimension::none, 1, sys.p.kn);
	outdata.entryData("kt", Dimensional::Dimension::none, 1, sys.p.kt);
	outdata.entryData("kr", Dimensional::Dimension::none, 1, sys.p.kr);
	vec3d shear_strain = sys.get_shear_strain();
	outdata.entryData("shear strain", Dimensional::Dimension::none, 3, shear_strain);
	if (sys.wall_rheology) {
		outdata.entryData("shear viscosity wall 1", Dimensional::Dimension::Viscosity, 1, sys.shearstress_wall1/sr);
		outdata.entryData("shear viscosity wall 2", Dimensional::Dimension::Viscosity, 1, sys.shearstress_wall2/sr);
		outdata.entryData("normal stress/rate wall 1", Dimensional::Dimension::Viscosity, 1, sys.normalstress_wall1/sr);
		outdata.entryData("normal stress/rate wall 2", Dimensional::Dimension::Viscosity, 1, sys.normalstress_wall2/sr);
	}
	if (sys.p.simulation_mode == 31) {
		outdata.entryData("force top wall", Dimensional::Dimension::Force, 3, sys.force_upwall);
		outdata.entryData("force bottom wall", Dimensional::Dimension::Force, 3, sys.force_downwall);
	}
	if (sys.brownian) {
		outdata.entryData("max_velocity_brownian", Dimensional::Dimension::Velocity, 1, sys.max_velocity_brownian);
		outdata.entryData("max_velocity_contact", Dimensional::Dimension::Velocity, 1, sys.max_velocity_contact);
	}
	if (sys.p.output.effective_coordination_number) {
		outdata.entryData("eff_coordination_number", Dimensional::Dimension::none, 1, sys.effective_coordination_number);
	}
	outdata.entryData("shear stress", Dimensional::Dimension::Stress, 1, sys.target_stress);
	outdata.entryData("theta shear", Dimensional::Dimension::none, 1, sys.p.theta_shear);
	if (sys.p.event_handler == "jamming_stress_reversal") {
		outdata.entryData("jamming strain", Dimensional::Dimension::none, 1, jamming_strain);
	}
	outdata.writeToFile();
	/****************************   Stress Tensor Output *****************/
	outdata_st.setUnits(system_of_units, output_unit);
	outdata_st.entryData("time", Dimensional::Dimension::Time, 1, sys.get_time());
	outdata_st.entryData("cumulated strain", Dimensional::Dimension::none, 1, sys.get_cumulated_strain());
	outdata_st.entryData("shear rate", Dimensional::Dimension::Rate, 1, sys.get_shear_rate());
	outdata_st.entryData("total stress tensor (xx, xy, xz, yz, yy, zz)", Dimensional::Dimension::Stress, 6, sys.total_stress);
	for (const auto &stress_comp: sys.total_stress_groups) {
		string entry_name = stress_comp.first+" stress tensor (xx, xy, xz, yz, yy, zz)";
		outdata_st.entryData(entry_name, Dimensional::Dimension::Stress, 6, stress_comp.second);
	}
	outdata_st.writeToFile();
}

void Simulation::getSnapshotHeader(stringstream& snapshot_header)
{
	string sep = " : ";
	snapshot_header << "# cumulated strain" << sep << sys.get_cumulated_strain() << endl;
	snapshot_header << "# shear disp" << sep << sys.shear_disp.x << endl;
	Dimensional::DimensionalQty<double> rate = {Dimensional::Dimension::Rate, sys.get_shear_rate(), system_of_units.getInternalUnit()};
	system_of_units.convertFromInternalUnit(rate, output_unit);
	Dimensional::DimensionalQty<double> stress = {Dimensional::Dimension::Stress, sys.target_stress, system_of_units.getInternalUnit()};
	system_of_units.convertFromInternalUnit(stress, output_unit);
	snapshot_header << "# shear rate" << sep << rate.value << endl;
	if (control_var == Parameters::ControlVariable::stress) {
		Dimensional::DimensionalQty<double> time = {Dimensional::Dimension::Time, sys.get_time(), system_of_units.getInternalUnit()};
		snapshot_header << "# target stress" << sep << stress.value << endl;
		snapshot_header << "# time" << sep << time.value << endl;
		snapshot_header << "# theta" << sep << sys.p.theta_shear << endl;
	}
	if (sys.ext_flow) {
		/* The following snapshot data is required to
		 * construct visualization file for extensional flow simulation in the script
		 * generateYaplotFile_extflow.pl
		 */
		double strain_retrimed = sys.strain_retrim-sys.strain_retrim_interval;
		snapshot_header << "# cumulated strain - strain_retrim" << sep << sys.get_cumulated_strain()-strain_retrimed << endl;
		if (sys.retrim_ext_flow) {
			snapshot_header << "# retrim ext flow " << sep << 1 << endl;
		} else {
			snapshot_header << "# retrim ext flow " << sep << 0 << endl;
		}
	}
	//snapshot_header << "# viscosity" << sep << viscosity << endl;
	//snapshot_header << "# n1" << sep << normal_stress_diff1 << endl;
}

vec3d Simulation::shiftUpCoordinate(double x, double y, double z)
{
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
	return vec3d(x,y,z);
}

void Simulation::createDataHeader(stringstream& data_header)
{
	if (!restart_from_chkp) {
		auto conf = sys.getBaseConfiguration();
		data_header << "# LF_DEM version " << GIT_VERSION << endl;
		data_header << "# np " << conf.position.size() << endl;
		data_header << "# VF " << conf.volume_or_area_fraction << endl;
		data_header << "# Lx " << conf.lx << endl;
		data_header << "# Ly " << conf.ly << endl;
		data_header << "# Lz " << conf.lz << endl;
		data_header << "# flow_type " << sys.p.flow_type << endl;
	}
}

void Simulation::outputPstFileTxt()
{
	outdata_pst.setUnits(system_of_units, output_unit);

	map<string, string> group_shorts;
	group_shorts["l"] = "hydro";
	group_shorts["c"] = "contact";
	group_shorts["r"] = "repulsion";
	group_shorts["b"] = "brownian";
	group_shorts["d"] = "dashpot";
	group_shorts["t"] = "total";
	map<string, vector<Sym2Tensor>> particle_stress;
	for (auto &type: sys.p.output.out_particle_stress) {
		auto group_name = group_shorts[string(1, type)];
		particle_stress[group_name] = getParticleStressGroup(group_name);
	}
	for (int i=0; i<sys.get_np(); i++) {
		for (const auto &pst: particle_stress) {
			outdata_pst.entryData(pst.first + " stress (xx, xy, xz, yz, yy, zz)", Dimensional::Dimension::Stress, 6, pst.second[i]);
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
	auto pos = sys.position;
	auto vel = sys.velocity;
	if (sys.p.output.origin_zero_flow) {
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
			if (pos[i].z < 0) {
				vel[i] -= sys.vel_difference;
			}
		}
	} else if (sys.p.output.relative_position_view) {
		relativePositionView(pos, vel);
	} else {
		for (int i=0; i<np; i++) {
			pos[i].x = sys.position[i].x-0.5*sys.get_lx();
			pos[i].y = sys.position[i].y-0.5*sys.get_ly();
			pos[i].z = sys.position[i].z-0.5*sys.get_lz();
		}
	}
	cout << "   out config: " << sys.get_cumulated_strain() << endl;
	outdata_par.setUnits(system_of_units, output_unit);

	for (int i=0; i<sys.get_np(); i++) {
		outdata_par.entryData("particle index", Dimensional::Dimension::none, 1, i);
		outdata_par.entryData("radius", Dimensional::Dimension::none, 1, sys.radius[i]);
		if (!sys.twodimension) {
			outdata_par.entryData("position (x, y, z)", Dimensional::Dimension::none, 3, pos[i], 6);
		} else {
			outdata_par.entryData("position x", Dimensional::Dimension::none, 1, pos[i].x, 6);
			outdata_par.entryData("position z", Dimensional::Dimension::none, 1, pos[i].z, 6);
		}
		if (diminish_output == false) {
			outdata_par.entryData("velocity (x, y, z)", Dimensional::Dimension::Velocity, 3, vel[i]);
			outdata_par.entryData("angular velocity (x, y, z)", Dimensional::Dimension::none, 3, sys.ang_velocity[i]);
			if (sys.twodimension) {
				outdata_par.entryData("angle", Dimensional::Dimension::none, 1, sys.angle[i]);
			}
		}
		//		if (sys.couette_stress) {
		//			double stress_rr, stress_thetatheta, stress_rtheta;
		//			sys.getStressCouette(i, stress_rr, stress_thetatheta, stress_rtheta);
		//			double sr = sys.get_shear_rate();
		//			outdata_par.entryData("stress_rr", Dimensional::Dimension::Viscosity, 1, stress_rr/sr);
		//			outdata_par.entryData("stress_thetatheta", Dimensional::Dimension::Viscosity, 1, stress_thetatheta/sr);
		//			outdata_par.entryData("stress_rtheta", Dimensional::Dimension::Viscosity, 1, stress_rtheta/sr);
		//		}
		if (sys.p.output.out_na_vel) {
			if (sys.twodimension) {
				outdata_par.entryData("non-affine velocity x", Dimensional::Dimension::Velocity, 1, sys.na_velocity[i].x);
				outdata_par.entryData("non-affine velocity z", Dimensional::Dimension::Velocity, 1, sys.na_velocity[i].z);
			} else {
				outdata_par.entryData("non-affine velocity (x, y, z)", Dimensional::Dimension::Velocity, 3, sys.na_velocity[i]);
			}
		}
		if (sys.p.output.out_na_disp) {
			outdata_par.entryData("non affine displacement (x, y, z)", Dimensional::Dimension::none, 3, sys.getNonAffineDisp()[i]);
		}
		if (sys.p.output.out_data_vel_components) {
			for (const auto &vc: sys.na_velo_components) {
				string entry_name_vel = "non-affine "+vc.first+" velocity (x, y, z)";
				string entry_name_ang_vel = "non-affine angular "+vc.first+" velocity (x, y, z)";
				outdata_par.entryData(entry_name_vel, Dimensional::Dimension::Velocity, 3, vc.second.vel[i]);
				outdata_par.entryData(entry_name_ang_vel, Dimensional::Dimension::Velocity, 3, vc.second.ang_vel[i]);
			}
		}
		if (sys.p.output.effective_coordination_number) {
			outdata_par.entryData("contact_number", Dimensional::Dimension::none, 1, sys.n_contact[i]);
		}
		
	}
	sys.resetNonAffineDispData();
	stringstream snapshot_header;
	getSnapshotHeader(snapshot_header);
	outdata_par.writeToFile(snapshot_header.str());
}

void Simulation::relativePositionView(std::vector<vec3d> &pos, std::vector<vec3d> &vel)
{
	int np = sys.get_np();
	for (int i=0; i<np; i++) {
		pos[i].set(sys.position[i].x-sys.position[0].x,
				   sys.position[i].y-sys.position[0].y,
				   sys.position[i].z-sys.position[0].z);
		if (pos[i].z > 0.5*sys.get_lz()) {
			pos[i].x -= sys.shear_disp.x;
			pos[i].y -= sys.shear_disp.y;
			if (pos[i].x < -0.5*sys.get_lx()) {
				pos[i].x += sys.get_lx();
			}
			if (pos[i].y < -0.5*sys.get_ly()) {
				pos[i].y += sys.get_ly();
			}
			pos[i].z -= sys.get_lz();
			for (int ii=0; i<np; i++) {
				if (pos[ii].z < 0) {
					vel[ii] -= sys.vel_difference;
				}
			}
		} else if (pos[i].z < -0.5*sys.get_lz()) {
			pos[i].x += sys.shear_disp.x;
			pos[i].y += sys.shear_disp.y;
			if (pos[i].x > 0.5*sys.get_lx()) {
				pos[i].x -= sys.get_lx();
			}
			if (pos[i].y > 0.5*sys.get_ly()) {
				pos[i].y -= sys.get_ly();
			}
			pos[i].z += sys.get_lz();
			for (int ii=0; i<np; i++) {
				if (pos[ii].z < 0) {
					vel[ii] += sys.vel_difference;
				}
			}
		}
		while (pos[i].x < -0.5*sys.get_lx()) {
			pos[i].x += sys.get_lx();
		}
		while (pos[i].x > 0.5*sys.get_lx()) {
			pos[i].x -= sys.get_lx();
		}
		while (pos[i].y < -0.5*sys.get_ly()) {
			pos[i].y += sys.get_ly();
		}
		while (pos[i].y > 0.5*sys.get_ly()) {
			pos[i].y -= sys.get_ly();
		}
	}
}

void Simulation::outputIntFileTxt()
{
	outdata_int.setUnits(system_of_units, output_unit);
	stringstream snapshot_header;
	getSnapshotHeader(snapshot_header);
	double sr = sys.get_shear_rate();
	for (const auto &inter: sys.interaction) {
		unsigned int i, j;
		std::tie(i, j) = inter.get_par_num();
		Sym2Tensor stress_contact = inter.contact.getContactStressXF();
		outdata_int.entryData("particle 1 label", Dimensional::Dimension::none, 1, i);
		outdata_int.entryData("particle 2 label", Dimensional::Dimension::none, 1, j);
		if (diminish_output == false) {
			outdata_int.entryData("contact state "
								  "(0 = no contact, "
								  "1 = frictionless contact, "
								  "2 = non-sliding frictional, "
								  "3 = sliding frictional)",
								  Dimensional::Dimension::none, 1, inter.contact.getFrictionState());
			outdata_int.entryData("normal vector, oriented from particle 1 to particle 2", \
								  Dimensional::Dimension::none, 3, inter.nvec);
			outdata_int.entryData("dimensionless gap = s-2, s = 2r/(a1+a2)", \
								  Dimensional::Dimension::none, 1,  inter.get_reduced_gap());
		}
		/* [NOTE]
		 * Lubrication forces are reference values
		 * in the Brownian case. The force balancing
		 * velocities are recalculated without
		 * including the Brownian forces.
		 * It seems there is no better way to visualize
		 * the lubrication forces.
		 */
		if (diminish_output == false) {
			if (sys.lubrication) {
				if (inter.get_reduced_gap() > 0) {
					double normal_part = -dot(inter.lubrication.getTotalForce(), inter.nvec);
					outdata_int.entryData("normal part of the lubrication force (positive for compression)", Dimensional::Dimension::Force, 1, \
										  normal_part);
					outdata_int.entryData("tangential part of the lubrication force", Dimensional::Dimension::Force, 3, \
										  inter.lubrication.getTangentialForce());
				} else {
					outdata_int.entryData("normal part of the lubrication force (positive for compression)", Dimensional::Dimension::Force, 1, 0);
					outdata_int.entryData("tangential part of the lubrication force", Dimensional::Dimension::Force, 3, vec3d(0,0,0));
				}
			}
		}
		/*
		 * Contact forces are the sums of spring forces and dashpot forces.
		 * (It can be negative even repulsive contact force).
		 */
		outdata_int.entryData("norm of the normal part of the contact force", Dimensional::Dimension::Force, 1, \
							  -inter.contact.getNormalForceValue());
		
		if (diminish_output == false) {
			outdata_int.entryData("tangential part of the contact force", Dimensional::Dimension::Force, 3, \
								  inter.contact.getTangentialForce());
		}
		if (sys.repulsiveforce) {
			outdata_int.entryData("norm of the normal repulsive force", Dimensional::Dimension::Force, 1, \
								  inter.repulsion.getForceNorm());
		}
		if (diminish_output == false) {
			outdata_int.entryData("Viscosity contribution of contact xF", Dimensional::Dimension::Stress, 1, \
								  doubledot(stress_contact, sys.getEinfty()/sr)/sr);
		}
		if (sys.delayed_adhesion) {
			outdata_int.entryData("norm of the normal adhesion force", Dimensional::Dimension::Force, 1, \
								  inter.delayed_adhesion->getForceNorm());
			if (diminish_output == false) {
				outdata_int.entryData("adhesion ratio uptime to activation time", Dimensional::Dimension::none, 1, \
									  inter.delayed_adhesion->ratioUptimeToActivation());
			}
		}
	}
	if (sys.interaction.size() > 0) {
		outdata_int.writeToFile(snapshot_header.str());
	}
}

void Simulation::outputConfigurationData()
{
	if (sys.p.output.out_data_particle) {
		outputParFileTxt();
	}
	if (sys.p.output.out_data_interaction) {
		outputIntFileTxt();
	}
	if (!sys.p.output.out_particle_stress.empty()) {
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

void Simulation::outputGSD()
{
	static std::vector<float> vectorBuffer;    // DIM * bufferSize
	static std::vector<float> scalarBuffer;    // bufferSize
	static std::vector<float> quaternionBuffer;    // bufferSize
	static bool first_time = true;
	static int ts = 0;
	int np = sys.get_np();
	if (first_time) {
		first_time = false;
		vectorBuffer.resize(3*np, 0);
		quaternionBuffer.resize(4*np, 0);
		scalarBuffer.resize(np, 0);
	}
	auto pos = sys.position;
	auto vel = sys.velocity;
	double lx = sys.get_lx();
	double ly = sys.get_ly();
	double lz = sys.get_lz();
	vec3d shear_strain = sys.get_shear_strain();
	if (sys.eventLookUp == NULL) {
		/* modulate for 0 < strain < 1
		 */
		shear_strain = sys.shear_disp/lz;
	} else {
		/* no modulation for deformed simulation cell.
		 * This is useful to visualize shear jamming.
		 */
		shear_strain = sys.get_shear_strain();
	}
	if (sys.twodimension) {
		ly = 2*sys.radius[np-1];
	}
	for (int i=0; i<np; i++) {
		if (-pos[i].x+shear_strain.x*pos[i].z > 0) {
			pos[i].x += lx;
		}
		if (-pos[i].x+shear_strain.x*pos[i].z < -lx) {
			pos[i].x -= lx;
		}
		pos[i].x -= (lx+shear_strain.x*lz)/2;
		pos[i].z -= lz/2;
		if (!sys.twodimension) {
			pos[i].y -= ly/2;
		}
	}
	uint64_t _ts = ts;
	uint8_t dim = 3;
	float box[6] = {static_cast<float>(lx), static_cast<float>(lz), static_cast<float>(ly),
		static_cast<float>(shear_strain.x),
		static_cast<float>(shear_strain.y),
		static_cast<float>(shear_strain.z)};
	
	gsd_write_chunk(&gsdOut, "confix1guration/step", GSD_TYPE_UINT64, 1, 1, 0, &_ts);
	gsd_write_chunk(&gsdOut, "configuration/dimensions", GSD_TYPE_UINT8, 1, 1, 0, &dim);
	gsd_write_chunk(&gsdOut, "configuration/box", GSD_TYPE_FLOAT, 6, 1, 0, &box);
	
	// Total number of elements / particles
	uint32_t n = np;
	gsd_write_chunk(&gsdOut, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n);

	if (ts == 0) { // Write type/particle names
		const int max_size = 63;
		int  n_types;
		char *types;
		if (dispersion_type == DispersionType::mono) {
			n_types = 1;
			types = new char [n_types*max_size];
			snprintf(types, max_size, "colloid");
		} else if (dispersion_type == DispersionType::bi) {
			n_types = 2;
			types = new char [n_types*max_size];
			snprintf(types, max_size, "colloid1");
			snprintf(types+max_size, max_size, "colloid2");
		} else {
			exit(1);
		}
		gsd_write_chunk(&gsdOut, "particles/types", GSD_TYPE_INT8, n_types, max_size, 0, types);
		delete[] types;
	}
	
	// particle diameters
	// particle IDs
	{
		float* fptr = scalarBuffer.data();
		unsigned int*  uptr = new unsigned int [np];
		// particle radius
		for (int i=0; i<np; i++) {
			fptr[i] = 2*sys.radius[i];
			if (i < np1) {
				uptr[i] = 0;
			} else {
				uptr[i] = 1;
			}
		}
		gsd_write_chunk(&gsdOut, "particles/diameter", GSD_TYPE_FLOAT, np, 1, 0, fptr);
		gsd_write_chunk(&gsdOut, "particles/typeid", GSD_TYPE_UINT32, np, 1, 0, uptr);
		delete[] uptr;
	}
	
	// particle positions
	{
		float* fptr = vectorBuffer.data();
		int j = 0;
		for (int i=0; i<np; i++) {
			fptr[j++] = pos[i].x;
			fptr[j++] = pos[i].z;
			fptr[j++] = pos[i].y;
		}
		gsd_write_chunk(&gsdOut, "particles/position", GSD_TYPE_FLOAT, np, 3, 0, fptr);
	}

	// particle velocities
	{
		float* fptr = vectorBuffer.data();
		int j = 0;
		for (int i=0; i<np; i++) {
			fptr[j++] = vel[i].x;
			fptr[j++] = vel[i].z;
			fptr[j++] = vel[i].y;
		}
		gsd_write_chunk(&gsdOut, "particles/velocity", GSD_TYPE_FLOAT, np, 3, 0, fptr);
	}

	{
		//particles/charge
		float* fptr = scalarBuffer.data();
		sys.calcTotalStressPerParticle();
		//float* fptr = scalarBuffer.data();
		for (int i=0; i<np; i++) {
			fptr[i] = -sys.total_stress_pp[i].trace()/3;
		}
		gsd_write_chunk(&gsdOut, "particles/charge", GSD_TYPE_FLOAT, np, 1, 0, fptr);
	}
	
	//	particles/orientation
	{
		if (sys.twodimension) {
			float* fptr = quaternionBuffer.data();
			int j = 0;
			// outdata_par.entryData("angle", Dimensional::Dimension::none, 1, sys.angle[i]);
			
			/* q = cos(theta/2) + rot_vec * sin(theta/2)
			 *   = (cos(theta/2), 0, 0 , sin(theta/2)
			 */
			for (int i=0; i<np; i++) {
				fptr[j++] = cos(sys.angle[i]/2);
				fptr[j++] = 0;
				fptr[j++] = 0;
				fptr[j++] = sin(sys.angle[i]/2);
			}
			gsd_write_chunk(&gsdOut, "particles/orientation", GSD_TYPE_FLOAT, np, 4, 0, fptr);
		}
	}

	gsd_end_frame(&gsdOut);
	ts++;
}
