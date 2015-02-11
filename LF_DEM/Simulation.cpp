//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#define _USE_MATH_DEFINES
#include "Simulation.h"
#ifndef GIT_VERSION
#include "VersionInfo.h"
#endif
/*
 * VersionInfo.h is automatically generated
 * before compiling source codes.
 * In Xcode, the following script is run in the Pre-Action of Build.
 * -----------------------------------
 * git=/usr/bin/git
 * cd ${PROJECT_DIR}/LF_DEM
 * version=`$git describe --dirty`
 * echo "#define GIT_VERSION \"$version\"" > VersionInfo.h
 * -----------------------------------
 */
#include <cmath>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <cctype>

Simulation::Simulation(): shear_rate_expectation(-1) {};
Simulation::~Simulation(){
	if (fout_rheo.is_open()) {
		fout_rheo.close();
	}
	if (fout_particle.is_open()) {
		fout_particle.close();
	}
	if (fout_interaction.is_open()) {
		fout_interaction.close();
	}
};

void
Simulation::contactForceParameter(string filename){
	ifstream fin_knktdt;
	fin_knktdt.open(filename.c_str());
	if (!fin_knktdt) {
		cerr << " Contact parameter file '" << filename << "' not found." <<endl;
		exit(1);
	}

	double phi_;
	double kn_;
	double kt_;
	double dt_;
	while (fin_knktdt >> phi_ >> kn_ >> kt_ >> dt_) {
		if (phi_ == volume_or_area_fraction) {
			break;
		}
	}
	fin_knktdt.close();
	p.kn = kn_;
	p.kt = kt_;
    p.dt = dt_;
	cerr << phi_ << ' ' << kn_ << ' ' << kt_ << ' ' << dt_ << endl;
}

void
Simulation::contactForceParameterBrownian(string filename){
	ifstream fin_knktdt;
	fin_knktdt.open(filename.c_str());
	if (!fin_knktdt) {
		cerr << " Contact parameter file '" << filename << "' not found." <<endl;
		exit(1);
	}

	double phi_;
	double peclet_;
	double kn_;
	double kt_;
	double dt_;
	bool found=false;
	while (fin_knktdt >> phi_ >> peclet_ >> kn_ >> kt_ >> dt_) {
		if (phi_ == volume_or_area_fraction && peclet_ == sys.dimensionless_shear_rate) {
			found = true;
			break;
		}
	}
	fin_knktdt.close();
	if(found){
		p.kn = kn_;
		p.kt = kt_;
		p.dt = dt_;
		cout << "Input for vf = " << phi_ << " and Pe = " << peclet_ << " : kn = " << kn_ << ", kt = " << kt_ << " and dt = " << dt_ << endl;
	}
	else{
		cerr << " Error: file " << filename.c_str() << " contains no data for vf = " << volume_or_area_fraction << " and Pe = " << sys.dimensionless_shear_rate << endl;
		exit(1);
	}
}

void
Simulation::importPreSimulationData(string filename){
	ifstream fin_PreSimulationData;
	fin_PreSimulationData.open(filename.c_str());
	if (!fin_PreSimulationData) {
		cerr << " Pre-simulation data file '" << filename << "' not found." <<endl;
		exit(1);
	}

	double stress_;
	double shear_rate_;
	while (fin_PreSimulationData >> stress_ >> shear_rate_) {
		if (stress_ == sys.target_stress_input) {
			break;
		}
	}
	shear_rate_expectation = shear_rate_;
}

void
Simulation::setupSimulationSteadyShear(vector<string> &input_files,
									   bool binary_conf,
									   double peclet_num,
									   double ratio_repulsion,
									   double ratio_cohesion,
									   double ratio_critical_load,
									   string control_variable){
	control_var = control_variable;
	filename_import_positions = input_files[0];
	filename_parameters = input_files[1];
	if (control_var == "rate") {
		if (peclet_num > 0) {
			cerr << "Brownian" << endl;
			sys.brownian = true;
			sys.dimensionless_shear_rate = peclet_num;
			if (ratio_repulsion > 0) {
				cerr << "Repulsive force" << endl;
				/* When both Brownian and repulsive forces exist
				 * `ratio_repulsion' = F_rep(0)/(kT/a)
				 * Filename includes "RepXXXX_PeXXXX".
				 */
				sys.set_ratio_repulsion(ratio_repulsion);
				sys.repulsiveforce_amplitude = ratio_repulsion/peclet_num;
				sys.repulsiveforce = true;
				string_control_parameters << "_r" << ratio_repulsion << "_p" << peclet_num;
			} else if (ratio_critical_load > 0) {
				cerr << "Critical load" << endl;
				sys.critical_normal_force = ratio_critical_load/peclet_num;
				p.friction_model = 2;
				string_control_parameters << "_c" << ratio_critical_load << "_p" << peclet_num;
			} else {
				cerr << "Only Brownian" << endl;
				string_control_parameters << "_p" << peclet_num;
			}
		} else {
			cerr << "non-Brownian" << endl;
			if (ratio_repulsion > 0 && ratio_cohesion == 0) {
				cerr << "Repulsive force" << endl;
				sys.dimensionless_shear_rate = ratio_repulsion;
				sys.repulsiveforce_amplitude = 1/ratio_repulsion;
				sys.repulsiveforce = true;
				string_control_parameters << "_r" << ratio_repulsion;
			} else if (ratio_repulsion == -1 && ratio_cohesion == 0) {
				cerr << "Infinite shear rate" << endl;
				sys.dimensionless_shear_rate = -1;
				sys.repulsiveforce_amplitude = 0;
				sys.repulsiveforce = true;
				string_control_parameters << "_rINF";
			} else if (ratio_critical_load > 0 && ratio_cohesion == 0) {
				cerr << "Critical load" << endl;
				sys.dimensionless_shear_rate = ratio_critical_load;
				p.friction_model = 2;
				sys.critical_normal_force = 1/ratio_critical_load;
				string_control_parameters << "_c" << ratio_critical_load;
			} else if (ratio_repulsion == 0 && ratio_cohesion > 0) {
				cerr << "Cohesive force" << endl;
				sys.dimensionless_shear_rate = ratio_cohesion;
				sys.cohesive_force = 1/ratio_cohesion;
				string_control_parameters << "_a" << ratio_cohesion;
			} else if (ratio_repulsion > 0 && ratio_cohesion > 0) {
				cerr << "Repulsive force + Cohesive force" << endl;
				sys.dimensionless_shear_rate = ratio_repulsion;
				sys.repulsiveforce_amplitude = 1/ratio_repulsion;
				sys.repulsiveforce = true;
				sys.cohesive_force = ratio_cohesion/ratio_repulsion;
				string_control_parameters << "_a" << ratio_cohesion << "_r" << ratio_repulsion;
			} else {
				cerr << "strain -> non-Brownian -> ???" << endl;
			}
		}
	} else if (control_var == "stress") {
		p.unscaled_contactmodel = true;
		sys.brownian = false;
		if (ratio_critical_load > 0) {
			cerr << " Stress controlled simulations for CLM not implemented ! " << endl;
			exit(1);
		}
		if (ratio_repulsion == 0 && ratio_cohesion == 0) {
			cerr << " Stress controlled simulations need a repulsive force ! " << endl;
			exit(1);
		} else {
			if (ratio_repulsion != 0) {
				cerr << "Repulsive force" << endl;
				sys.repulsiveforce = true;
				sys.repulsiveforce_amplitude = 1;
				sys.target_stress_input = ratio_repulsion;
				sys.target_stress = ratio_repulsion/6/M_PI;
				string_control_parameters << "_s" << ratio_repulsion;
			} else if (ratio_cohesion != 0) {
				cerr << "Cohesive force" << endl;
				p.unscaled_contactmodel = false;
				sys.cohesion = true;
				//sys.dimensionless_shear_rate = ratio_cohesion;
				sys.cohesive_force = 1;
				sys.target_stress_input = ratio_cohesion;
				sys.target_stress = ratio_cohesion/6/M_PI;
				/* Initial relaxation for stress control simulation.
				 * (To avoid breaking bonds due to startup flows.)
				 */
				sys.init_strain_shear_rate_limit = -9999;
				sys.init_shear_rate_limit = 9999;
				string_control_parameters << "_b" << ratio_cohesion;
			}
			sys.dimensionless_shear_rate = 1; // needed for 1st time step
			/* The target stress (``ratio_repulsion'') is given trough the command argument
			 * with an unit stres: eta_0*gammmadot_0.
			 * However, in the code, sys.target_stress is computed as an unit F_rep/a^2.
			 */
		}
	}
	setDefaultParameters();
	readParameterFile();

 	if(binary_conf){
		importConfigurationBinary();
	}
	else{
		importInitialPositionFile();
	}
	if (initial_lees_edwards_disp > 0){
		sys.shear_disp = initial_lees_edwards_disp;
	} else {
		sys.shear_disp = 0;
	}


	if (input_files[2] != "not_given") {
		if(sys.brownian&&!p.auto_determine_knkt){
			contactForceParameterBrownian(input_files[2]);
		}
		else{
			contactForceParameter(input_files[2]);
		}
	}
	if (input_files[3] != "not_given") {
		importPreSimulationData(input_files[3]);
		time_interval_output_data = p.strain_interval_output_data/shear_rate_expectation;
		time_interval_output_config = p.strain_interval_output_config/shear_rate_expectation;
	} else {
		time_interval_output_data = -1;
		time_interval_output_config = -1;
	}
	exportParameterSet();
	if (sys.brownian) {
		sys.setupBrownian();
	}
	sys.setupSystem(control_var);
	openOutputFiles(binary_conf);
	if (filename_parameters == "init_relax.txt") {
		sys.zero_shear = true;
	}
	sys.setupShearFlow(true);
	if (control_var == "stress") {
		p.integration_method = 0;
	}
}

/*
 * Main simulation
 */
void
Simulation::simulationSteadyShear(vector<string> &input_files,
								  bool binary_conf,
								  double peclet_num, double ratio_repulsion, double ratio_cohesion,
								  double ratio_critical_load, string control_variable){
	user_sequence = false;
	control_var = control_variable;
	setupSimulationSteadyShear(input_files, binary_conf, peclet_num,
							   ratio_repulsion, ratio_cohesion, ratio_critical_load, control_var);
	int cnt_simu_loop = 1;
	int cnt_config_out = 1;
	double strain_output_data = 0;
	double strain_output_config = 0;
	double time_output_data = 0;
	double time_output_config = 0;
	sys.new_contact_gap = 0.02;
	int jammed = 0;
	while (sys.get_shear_strain() < p.shear_strain_end-1e-8) {
		if (time_interval_output_data == -1) {
			strain_output_data = cnt_simu_loop*p.strain_interval_output_data;
			strain_output_config = cnt_config_out*p.strain_interval_output_config;
		} else {
			time_output_data = cnt_simu_loop*time_interval_output_data;
			time_output_config = cnt_config_out*time_interval_output_config;
		}

		sys.timeEvolution(strain_output_data, time_output_data);
		cnt_simu_loop ++;
		evaluateData();
		outputRheologyData();
		outputStressTensorData();
		outputConfigurationBinary();
		if (time_interval_output_data == -1) {
			if (sys.get_shear_strain() >= strain_output_config-1e-8) {
				cerr << "   out config: " << sys.get_shear_strain() << endl;
				outputConfigurationData();
				cnt_config_out ++;
				
			}
		} else {
			if (sys.get_time() >= time_output_config-1e-8) {
				cerr << "   out config: " << sys.get_shear_strain() << endl;
				outputConfigurationData();
				cnt_config_out ++;
			}
		}
		cerr << "strain: " << sys.get_shear_strain() << " / " << p.shear_strain_end << endl;
		if (abs(sys.dimensionless_shear_rate) < 1e-4 ){
			cerr << "shear jamming " << jammed << endl;
			jammed ++;
			if (jammed > 10) {
				jammed = 0;
				cerr << "shear jamming";
				break;
			}
		} else {
			jammed = 0;
		}
		sys.new_contact_gap = 0;
	}
	if (filename_parameters == "init_relax.txt") {
		/* To prepare relaxed initial configuration,
		 * we can use Brownian simulation for a short interval.
		 * Here is just to export the position data.
		 */
		outputFinalConfiguration();
	}
}

/*
 * Main simulation
 */
void
Simulation::simulationUserDefinedSequence(string seq_type, vector<string> &input_files, bool binary_conf, string control_variable){
	user_sequence = true;
	control_var = control_variable;
	filename_import_positions = input_files[0];
	filename_parameters = input_files[1];
	filename_sequence = input_files[4];
	string::size_type pos_ext_sequence = filename_sequence.find(".dat");
	sys.brownian = false;
	sys.target_stress_input = 0;
	sys.target_stress = 0;
	cerr << seq_type << endl;
	if (seq_type == "S") {
		p.unscaled_contactmodel = true;
		cerr << "Repulsive force" << endl;
		sys.repulsiveforce = true;
		sys.repulsiveforce_amplitude = 1;
		string_control_parameters << "_S" << filename_sequence.substr(0, pos_ext_sequence);
	} else if (seq_type == "R") {
		//p.unscaled_contactmodel
		sys.repulsiveforce = true;
		cerr << "Repulsive force" << endl;
		cerr << " User Defined Sequence only implemented for ....\n";
		exit(1);
	} else if (seq_type == "B") {
		cerr << "Cohesive force" << endl;
		sys.repulsiveforce = false;
		sys.cohesion = true;
		sys.cohesive_force = 1;
		string_control_parameters << "_B" << filename_sequence.substr(0, pos_ext_sequence);
	} else {
		cerr << " User Defined Sequence only implemented for ....\n";
		exit(1);
	}
	setDefaultParameters();
	readParameterFile();
 	if(binary_conf){
		importConfigurationBinary();
	}
	else{
		importInitialPositionFile();
	}
	if (input_files[3] != "not_given") {
		importPreSimulationData(input_files[3]);
		// strain_interval_out
		time_interval_output_data = p.strain_interval_output_data/shear_rate_expectation;
		time_interval_output_config = p.strain_interval_output_config/shear_rate_expectation;
	} else {
		time_interval_output_data = -1;
		time_interval_output_config = -1;
	}
	
	p.integration_method = 0;
	sys.importParameterSet(p);
	sys.setupSystem(control_var);
	openOutputFiles(binary_conf);
	outputConfigurationData();
	sys.setupShearFlow(true);
	vector <double> strain_sequence;
	vector <double> rsequence;
	ifstream fin_seq;
	fin_seq.open(filename_sequence.c_str());
	if (!fin_seq) {
		cerr << " Sequence file '" << filename_sequence << "' not found." <<endl;
		exit(1);
	}
		
	double strain;
	double targ_st;
	while (fin_seq >> strain >> targ_st){
		strain_sequence.push_back(strain);
		rsequence.push_back(targ_st);
	}
	int cnt_simu_loop = 1;
	int cnt_config_out = 1;
	double next_strain = 0;
	double strain_output_data = 0;
	double strain_output_config = 0;
	double time_output_data = 0;
	double time_output_config = 0;
	int jammed = 0;
	for (unsigned int step = 0; step<strain_sequence.size(); step++){
		/* The target stress (``rsequence'') is given trough the command argument
		 * with an unit stres: eta_0*gammmadot_0.
		 * However, in the code, sys.target_stress is computed as an unit F_rep/a^2.
		 */
		sys.target_stress_input = rsequence[step];
		sys.target_stress = rsequence[step]/6/M_PI;
		cerr << "Target stress " << sys.target_stress_input << endl;
		sys.updateUnscaledContactmodel();
		sys.dimensionless_shear_rate = 1; // needed for 1st time step
		next_strain = sys.get_shear_strain()+strain_sequence[step];
		while (sys.get_shear_strain() < next_strain-1e-8) {
			if (time_interval_output_data == -1) {
				strain_output_data = cnt_simu_loop*p.strain_interval_output_data;
				strain_output_config = cnt_config_out*p.strain_interval_output_config;
			} else {
				time_output_data = cnt_simu_loop*time_interval_output_data;
				time_output_config = cnt_config_out*time_interval_output_config;
			}
			sys.timeEvolution(strain_output_data, time_output_data);
			cnt_simu_loop ++;
			evaluateData();
			outputRheologyData();
			outputStressTensorData();
			outputConfigurationBinary();
			if (time_interval_output_data == -1) {
				if (sys.get_shear_strain() >= strain_output_config-1e-8) {
					cerr << "   out config: " << sys.get_shear_strain() << endl;
					outputConfigurationData();
					cnt_config_out ++;
				}
			} else {
				if (sys.get_time() >= time_output_config-1e-8) {
					cerr << "   out config: " << sys.get_shear_strain() << endl;
					outputConfigurationData();
					cnt_config_out ++;
				}
			}
 			if (abs(sys.dimensionless_shear_rate) < 1e-4 ){
				cerr << "shear jamming " << jammed << endl;
				jammed ++;
				if (jammed > 10) {
					jammed = 0;
					cerr << "shear jamming";
					break;
				}
			} else {
				jammed = 0;
			}
			cerr << "strain: " << sys.get_shear_strain() << " / " << p.shear_strain_end;
			cerr << "      stress = " << sys.target_stress_input << endl;
		}
	}
}

bool
str2bool(string value){
	if (value == "true") {
		return true;
	} else if (value == "false") {
		return false;
	} else {
		cerr << "The value should be true or false" << endl;
		exit(1);
	}
}

void
Str2KeyValue(string &str_parameter,
			 string &keyword,
			 string &value){
	string::size_type pos_equal = str_parameter.find("=");

	keyword = str_parameter.substr(0, pos_equal);
	value = str_parameter.substr(pos_equal+1);
	return;
}

void
removeBlank(string &str){
	str.erase(std::remove_if(str.begin(), str.end(), (int(*)(int))isspace), str.end());
}

void
Simulation::autoSetParameters(const string &keyword, const string &value){
	
	if (keyword == "lubrication_model") {
		p.lubrication_model = atoi(value.c_str());
	} else if (keyword == "friction_model") {
		if (p.friction_model == 2) {
			cerr << "!!Neglected friction_model in parameter file!!" << endl;
	} else {
			p.friction_model = atoi(value.c_str());
		}
	} else if (keyword == "rolling_friction") {
		p.rolling_friction = str2bool(value);
	} else if (keyword == "unscaled_contactmodel") {
		p.unscaled_contactmodel = str2bool(value);
	} else if (keyword == "repulsiveforce_length") {
		if (sys.repulsiveforce) {
			p.repulsive_length = atof(value.c_str());
		}
	} else if (keyword == "lub_reduce_parameter") {
		p.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time") {
		p.contact_relaxation_time = atof(value.c_str());
	} else if (keyword == "contact_relaxation_time_tan"){
		p.contact_relaxation_time_tan =  atof(value.c_str());
	} else if (keyword == "disp_max") {
		p.disp_max = atof(value.c_str());
	} else if (keyword == "shear_strain_end") {
		p.shear_strain_end = atof(value.c_str());
	} else if (keyword == "integration_method") {
		p.integration_method = atoi(value.c_str());
	} else if (keyword == "lub_max") {
		p.lub_max = atof(value.c_str());
	} else if (keyword == "sd_coeff") {
		p.sd_coeff = atof(value.c_str());
	} else if (keyword == "kn") {
		p.kn = atof(value.c_str());
	} else if (keyword == "kt") {
		p.kt = atof(value.c_str());
	} else if (keyword == "kr") {
		p.kr = atof(value.c_str());
	} else if (keyword == "dt") {
		p.dt = atof(value.c_str());
		//	} else if (keyword == "kn_lowPeclet") {
		//		p.kn_lowPeclet = atof(value.c_str());
		//	} else if (keyword == "kt_lowPeclet") {
		//		p.kt_lowPeclet = atof(value.c_str());
		//	} else if (keyword == "dt_lowPeclet") {
		//		p.dt_lowPeclet = atof(value.c_str());
	} else if (keyword == "Pe_switch") {
		p.Pe_switch = atof(value.c_str());
	} else if (keyword == "mu_static") {
		p.mu_static = atof(value.c_str());
	} else if (keyword == "strain_interval_output_config") {
		p.strain_interval_output_config = atof(value.c_str());
	} else if (keyword == "strain_interval_output_data") {
		p.strain_interval_output_data = atof(value.c_str());
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
	} else {
		cerr << "keyword " << keyword << " is not associated with an parameter" << endl;
		exit(1);
	}
}

void
Simulation::readParameterFile(){
	ifstream fin;
	fin.open(filename_parameters.c_str());
	if (!fin) {
		cerr << " Parameter file '" << filename_parameters << "' not found." <<endl;
		exit(1);
	}

	string keyword;
	string value;
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
		if (begin_comment > end_comment ) {
			cerr << str_parameter.find("/*") << endl;
			cerr << str_parameter.find("*/") << endl;
			cerr << "syntax error in the parameter file." << endl;
			exit(1);
		}
		string::size_type pos_slashslash = str_parameter.find("//");
		if( pos_slashslash != string::npos) {
			cerr << " // is not the syntax to comment out. Use /* comment */" << endl;
			exit(1);
		}
		Str2KeyValue(str_parameter, keyword, value);
		autoSetParameters(keyword, value);
	}
	fin.close();
	return;
}

void
Simulation::openOutputFiles(bool binary_conf){
	/*
	 * Set simulation name and name of output files.
	 */
	prepareSimulationName(binary_conf);


	string st_filename = "st_" +sys.simu_name + ".dat";
	fout_st.open(st_filename.c_str());
	outputDataHeader(fout_st);

	string rheo_filename = "rheo_" + sys.simu_name + ".dat";
	fout_rheo.open(rheo_filename.c_str());
	outputDataHeader(fout_rheo);
	//
	string fout_rheo_col_def =
	"#1: shear strain\n"
	"#2: Viscosity\n"
	"#3: N1\n"
	"#4: N2\n"
	"#5: Viscosity(lub)\n"
	"#6: N1(lub)\n"
	"#7: N2(lub)\n"
	"#8: Viscosity(xF_contact part)\n"
	"#9: N1(xF_contact part)\n"
	"#10: N2(xF_contact part)\n"
	"#11: Viscosity(GU_contact part)\n"
	"#12: N1(GU_contact part)\n"
	"#13: N2(GU_contact part)\n"
	"#14: Viscosity(friction)\n"
	"#15: N1(friction)\n"
	"#16: N2(friction)\n"
	"#17: Viscosity(repulsive force XF)\n"
	"#18: N1(repulsive force XF)\n"
	"#19: N2(repulsive force XF)\n"
	"#20: Viscosity(repulsive force GU)\n"
	"#21: N1(repulsive force GU)\n"
	"#22: N2(repulsive force GU)\n"
	"#23: Viscosity(brownian)\n"
	"#24: N1(brownian)\n"
	"#25: N2(brownian)\n"
	"#26: particle pressure\n"
	"#27: particle pressure contact\n"
	"#28: min gap (non-dim)\n"
	"#29: max tangential displacement\n"
	"#30: max Fc_normal\n"
	"#31: max Fc_tan\n"
	"#32: max velocity\n"
	"#33: max angular velocity\n"
	"#34: ave contact normal velocity\n"
	"#35: max contact normal velocity\n"
	"#36: ave contact tangential velocity\n"
	"#37: max contact tangential velocity\n"
	"#38: ave sliding velocity\n"
	"#39: max sliding velocity\n"
	"#40: ave contact number per particle\n"
	"#41: num of interaction\n"
	"#42: num of contacts\n"
	"#43: num of frictional contacts\n"
	"#44: kn\n"
	"#45: kt\n"
	"#46: dt\n"
	"#47: time\n"
	"#48: dimensionless_shear_rate\n"
	"#49: stress\n"
	"#50: shear_disp\n";
	//
	fout_rheo << fout_rheo_col_def << endl;

	
	if(p.out_data_particle){
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

	if(p.out_data_interaction){
		string interaction_filename = "int_" + sys.simu_name + ".dat";
		fout_interaction.open(interaction_filename.c_str());
		outputDataHeader(fout_interaction);
		string fout_int_col_def =
			"#1: particle 1 label\n"
			"#2: particle 2 label\n"
			"#3: contact state (0 = no contact, 1 = frictionless contact, 1 = non-sliding frictional, 2 = sliding frictional)\n"
			"#4: normal vector, oriented from particle 1 to particle 2 x\n"
			"#5: normal vector, oriented from particle 1 to particle 2 y\n"
			"#6: normal vector, oriented from particle 1 to particle 2 z\n"
			"#7: dimensionless gap = s-2, s = 2r/(a1+a2)\n"
			"#8: norm of the normal part of the lubrication force\n"
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


void
Simulation::exportParameterSet(){
	sys.importParameterSet(p);
}

void
Simulation::setDefaultParameters(){
	p.Pe_switch = 5;
	p.dt = 1e-4;
	//	p.dt_lowPeclet = 1e-4;
	p.disp_max = 2e-3;
	
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

	p.rolling_friction = false;
	p.shear_strain_end = 10;
	p.lub_max = 2.5;
	/*
	 * reduced_gap_min: gives reduced lubrication (maximum coeeffient).
	 *
	 */
	p.lub_reduce_parameter = 1e-3;
	/*
	 * contact_relaxation_factor:
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

	p.mu_static = 1;

	p.strain_interval_output_data = 0.01;
	p.strain_interval_output_config = 0.1;
	p.origin_zero_flow = true;

	p.out_data_particle = true;
	p.out_data_interaction = true;

}

void
Simulation::importInitialPositionFile(){
	fstream file_import;
	file_import.open(filename_import_positions.c_str());
	if (!file_import) {
		cerr << " Position file '" << filename_import_positions << "' not found." <<endl;
		exit(1);
	}
	int n1, n2;
	double lx, ly, lz;
	double vf1, vf2;
	
	char buf;
	getline(file_import, import_line[0]);
	getline(file_import, import_line[1]);
	stringstream ss(import_line[1]);
	ss >> buf >> n1 >> n2 >> volume_or_area_fraction >> lx >> ly >> lz >> vf1 >> vf2 >> initial_lees_edwards_disp;
	double x_, y_, z_, a_;
	vector<vec3d> initial_position;
	vector <double> radius;
	while (file_import >> x_ >> y_ >> z_ >> a_) {
		initial_position.push_back(vec3d(x_, y_, z_));
		radius.push_back(a_);
	}
	file_import.close();
	sys.setConfiguration(initial_position, radius, lx, ly, lz);
}

void
Simulation::outputConfigurationBinary(){
	string conf_filename;
	//	conf_filename =  "conf_" + sys.simu_name + "_strain" + to_string(sys.get_shear_strain()) + ".dat";
	conf_filename =  "conf_" + sys.simu_name + ".dat";
	outputConfigurationBinary(conf_filename);
}
void
Simulation::outputConfigurationBinary(string conf_filename){

	vector < vector <double> > pos;
	int np = sys.get_np();
	int dims = 4;
	pos.resize(np);
	
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
	double shear_disp = sys.shear_disp;

	conf_export.open(conf_filename.c_str(), ios::binary | ios::out);
	conf_export.write((char*)&np, sizeof(int));
	conf_export.write((char*)&volume_or_area_fraction, sizeof(double));
	conf_export.write((char*)&lx, sizeof(double));
	conf_export.write((char*)&ly, sizeof(double));
	conf_export.write((char*)&lz, sizeof(double));
	conf_export.write((char*)&shear_disp, sizeof(double));
	for (int i=0; i<np; i++) {
		conf_export.write((char*)&pos[i][0], dims*sizeof(double));
	}
	
	conf_export.close();
}

void
Simulation::importConfigurationBinary(){
	ifstream file_import;
	file_import.open(filename_import_positions.c_str(), ios::binary | ios::in);
	if (!file_import) {
		cerr << " Position file '" << filename_import_positions << "' not found." <<endl;
		exit(1);
	}

	int np;
	double lx;
	double ly;
	double lz;
	file_import.read((char*)&np, sizeof(int));
	file_import.read((char*)&volume_or_area_fraction, sizeof(double));
	file_import.read((char*)&lx, sizeof(double));
	file_import.read((char*)&ly, sizeof(double));
	file_import.read((char*)&lz, sizeof(double));
	file_import.read((char*)&initial_lees_edwards_disp, sizeof(double));

	double x_;
	double y_;
	double z_;
	double r_;
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

void
Simulation::prepareSimulationName(bool binary_conf){
	ostringstream ss_simu_name;
	string::size_type pos_name_end = filename_import_positions.find_last_of(".");
	string::size_type param_name_end = filename_parameters.find_last_of(".");
	string::size_type pos_name_start;
	if(binary_conf){ // TO DO: improve name generation for binary input
		pos_name_start = filename_import_positions.find_last_of("/");
	}else{
		pos_name_start = filename_import_positions.find_last_of("/");
	}
	string::size_type param_name_start = filename_parameters.find_last_of("/");
	if(pos_name_start == std::string::npos){
		pos_name_start = -1;
	}
	if(param_name_start == std::string::npos){
		param_name_start = -1;
	}
	pos_name_start += 1;
	param_name_start += 1;
	ss_simu_name << filename_import_positions.substr(pos_name_start, pos_name_end-pos_name_start);
	ss_simu_name << "_";
	ss_simu_name << filename_parameters.substr(param_name_start, param_name_end-param_name_start);
	ss_simu_name << string_control_parameters.str();
	sys.simu_name = ss_simu_name.str();
	cerr << "filename: " << sys.simu_name << endl;	
}

void
Simulation::evaluateData(){
	sys.analyzeState();
	sys.calcStress();
	sys.calcLubricationForce();
		
	viscosity = sys.einstein_viscosity+sys.total_stress.getStressXZ();
	normalstress_diff_1 = sys.total_stress.getNormalStress1();
	normalstress_diff_2 = sys.total_stress.getNormalStress2();
	particle_pressure = sys.total_stress.getParticlePressure();
	viscosity_hydro = sys.total_hydro_stress.getStressXZ();
	normalstress_diff_1_hydro = sys.total_hydro_stress.getNormalStress1();
	normalstress_diff_2_hydro = sys.total_hydro_stress.getNormalStress2();
	viscosity_cont_XF = sys.total_contact_stressXF.getStressXZ();
	normalstress_diff_1_cont_XF = sys.total_contact_stressXF.getNormalStress1();
	normalstress_diff_2_cont_XF = sys.total_contact_stressXF.getNormalStress2();
	particle_pressure_cont = sys.total_contact_stressXF.getParticlePressure();
	viscosity_friction = sys.total_contact_stressXF_tan.getStressXZ();
	normalstress_diff_1_friction = sys.total_contact_stressXF_tan.getNormalStress1();
	normalstress_diff_2_friction = sys.total_contact_stressXF_tan.getNormalStress2();
	viscosity_cont_GU = sys.total_contact_stressGU.getStressXZ();
	normalstress_diff_1_cont_GU = sys.total_contact_stressGU.getNormalStress1();
	normalstress_diff_2_cont_GU = sys.total_contact_stressGU.getNormalStress2();
	if (sys.repulsiveforce) {
		viscosity_repulsive_XF = sys.total_repulsive_stressXF.getStressXZ();
		normalstress_diff_1_repulsive_XF = sys.total_repulsive_stressXF.getNormalStress1();
		normalstress_diff_2_repulsive_XF = sys.total_repulsive_stressXF.getNormalStress2();
		particle_pressure_repulsive = sys.total_repulsive_stressXF.getParticlePressure();
		viscosity_repulsive_GU = sys.total_repulsive_stressGU.getStressXZ();
		normalstress_diff_1_repulsive_GU = sys.total_repulsive_stressGU.getNormalStress1();
		normalstress_diff_2_repulsive_GU = sys.total_repulsive_stressGU.getNormalStress2();
	} else {
		viscosity_repulsive_XF = 0;
		normalstress_diff_1_repulsive_XF = 0;
		normalstress_diff_2_repulsive_XF = 0;
		particle_pressure_repulsive = 0;
		viscosity_repulsive_GU = 0;
		normalstress_diff_1_repulsive_GU = 0;
		normalstress_diff_2_repulsive_GU = 0;
	}
	if (sys.brownian) {
		viscosity_brownian = sys.total_brownian_stressGU.getStressXZ();
		normalstress_diff_1_brownian = sys.total_brownian_stressGU.getNormalStress1();
		normalstress_diff_2_brownian = sys.total_brownian_stressGU.getNormalStress2();
	} else {
		viscosity_brownian = 0;
		normalstress_diff_1_brownian = 0;
		normalstress_diff_2_brownian = 0;
	}
}

void
Simulation::outputStressTensorData(){
	fout_st << sys.get_shear_strain() << ' ';
	fout_st << 6*M_PI*viscosity << ' ';
	/* total_stress = sys.total_hydro_stress;
	 * + total_contact_stressXF + total_repulsive_stress;
	 */
	// As it is, the output stress lacks a 6pi factor (as the viscosity)
	sys.total_stress.outputStressTensor(fout_st); // (3,4,5,6,7,8)
	sys.total_hydro_stress.outputStressTensor(fout_st); // (9,10,11,12,13,14)
	sys.total_contact_stressXF.outputStressTensor(fout_st); // (15,16,17,18,19,20)
	sys.total_contact_stressGU.outputStressTensor(fout_st); // (21,22,23,24,25,26)
	sys.total_repulsive_stress.outputStressTensor(fout_st); // (27,28,29,30,31,32)
	sys.total_brownian_stressGU.outputStressTensor(fout_st); // (33,34,35,36,37,38)
	fout_st << sys.dimensionless_shear_rate << ' '; // 39
	fout_st << endl;
}

void
Simulation::outputRheologyData(){
	/*
	 * Output the sum of the normal forces.
	 *
	 *  Viscosity = S_{xz} / shear_rate
	 *  N1 = S_{xx}-S_{zz}
	 *  N2 = S_{zz}-S_{yy} = S_zz-(-S_xx-S_zz) = S_xx+2*S_zz
	 *
	 * Relative viscosity = Viscosity / viscosity_solvent
	 */
	
	/*
	 * hat(...) indicates dimensionless quantities.
	 * (1) relative viscosity = Sxz/(eta0*shear_rate) = 6*pi*hat(Sxz)
	 * (2) N1/(eta0*shear_rate) = 6*pi*hat(N1)
	 * (3) N2/(eta0*shear_rate) = 6*pi*hat(N2)
	 *
	 * In simulation, we use the force unit where Stokes drag is F = -(U-U^inf)
	 *
	 * [note] In stress controlled simulation,
	 * Averaged viscosity need to be calculated with dimensionless_shear_rate_averaged,
	 * i.e. <viscosity> = taget_stress / dimensionless_shear_rate_averaged.
	 */
	fout_rheo << sys.get_shear_strain() << ' '; //1
	fout_rheo << 6*M_PI*viscosity << ' '; //2
	fout_rheo << 6*M_PI*normalstress_diff_1 << ' '; //3
	fout_rheo << 6*M_PI*normalstress_diff_2 << ' '; //4
	/*
	 * Hydrodynamic contribution means
	 * stresslet_hydro_GU_i+stresslet_ME_i from vel_hydro
	 * vel_hydro is obtained with GE for the rhs.
	 *
	 * "_hydro" might be bit confusing.
	 * Something indicating "E_inf" would be better.
	 */
	fout_rheo << 6*M_PI*viscosity_hydro << ' '; //5
	fout_rheo << 6*M_PI*normalstress_diff_1_hydro << ' '; //6
	fout_rheo << 6*M_PI*normalstress_diff_2_hydro << ' '; //7
	/*
	 * Contact force contribution seems to be
	 * the sum of viscosity_cont_XF and viscosity_cont_GU.
	 */
	fout_rheo << 6*M_PI*viscosity_cont_XF << ' '; //8
	fout_rheo << 6*M_PI*normalstress_diff_1_cont_XF << ' '; //9
	fout_rheo << 6*M_PI*normalstress_diff_2_cont_XF << ' '; //10
	fout_rheo << 6*M_PI*viscosity_cont_GU << ' ' ; //11
	fout_rheo << 6*M_PI*normalstress_diff_1_cont_GU << ' ' ; //12
	fout_rheo << 6*M_PI*normalstress_diff_2_cont_GU << ' ' ; //13
	/*
	 *
	 */
	fout_rheo << 6*M_PI*viscosity_friction << ' '; //14
	fout_rheo << 6*M_PI*normalstress_diff_1_friction << ' '; //15
	fout_rheo << 6*M_PI*normalstress_diff_2_friction  << ' '; //16
	fout_rheo << 6*M_PI*viscosity_repulsive_XF << ' '; //17
	fout_rheo << 6*M_PI*normalstress_diff_1_repulsive_XF << ' '; //18
	fout_rheo << 6*M_PI*normalstress_diff_2_repulsive_XF << ' '; //19
	fout_rheo << 6*M_PI*viscosity_repulsive_GU << ' '; //20
	fout_rheo << 6*M_PI*normalstress_diff_1_repulsive_GU << ' '; //21
	fout_rheo << 6*M_PI*normalstress_diff_2_repulsive_GU << ' '; //22
	fout_rheo << 6*M_PI*viscosity_brownian << ' ' ; //23
	fout_rheo << 6*M_PI*normalstress_diff_1_brownian << ' ' ; //24
	fout_rheo << 6*M_PI*normalstress_diff_2_brownian << ' ' ; //25
	fout_rheo << 6*M_PI*particle_pressure << ' ';//26
	fout_rheo << 6*M_PI*particle_pressure_cont << ' ';//27
	fout_rheo << sys.min_reduced_gap << ' '; //28
	fout_rheo << sys.max_disp_tan << ' '; //29
	fout_rheo << sys.max_fc_normal << ' '; //30
	fout_rheo << sys.max_fc_tan << ' ';//31
	fout_rheo << sys.max_velocity << ' '; //32
	fout_rheo << sys.max_ang_velocity << ' '; //33
	fout_rheo << sys.ave_contact_velo_normal << ' '; //34
	fout_rheo << sys.max_contact_velo_normal << ' '; //35
	fout_rheo << sys.ave_contact_velo_tan << ' '; //36
	fout_rheo << sys.max_contact_velo_tan << ' '; //37
	fout_rheo << sys.ave_sliding_velocity << ' ' ; //38
	fout_rheo << sys.max_sliding_velocity << ' ' ; //39
	fout_rheo << sys.getParticleContactNumber() << ' ';//40
	fout_rheo << sys.get_nb_of_active_interactions() << ' ';//41
	fout_rheo << sys.contact_nb << ' '; //42
	fout_rheo << sys.fric_contact_nb << ' '; //43
	fout_rheo << sys.kn << ' '; //44
	fout_rheo << sys.kt << ' '; //45
	fout_rheo << sys.dt << ' '; //46
	fout_rheo << sys.get_time() << ' ' ; //47
	/* In stress control simulation,
	 * shear jammed state may cause oscilation of dimensionless_shear_rate around 0.
	 * Then, time step also oscilate.
	 * This is why we need to take time average to have correct value of dimensionless_shear_rate.
	 */
	fout_rheo << sys.dimensionless_shear_rate << ' '; // 48
	if(control_var == "stress"){
		fout_rheo << sys.target_stress_input << ' '; // 49
	}
	else{
		fout_rheo << 6*M_PI*viscosity*sys.dimensionless_shear_rate << ' '; // 49
	}
	fout_rheo << sys.shear_disp << ' '; // 50
	fout_rheo << endl;
}

vec3d
Simulation::shiftUpCoordinate(double x, double y, double z){
	if (p.origin_zero_flow) {
		z += sys.Lz_half();
		if (z > sys.Lz_half()) {
			x -= sys.shear_disp;
			if (x < -sys.Lx_half()) {
				x += sys.get_lx();
			}
			z -= sys.get_lz();
		}
	}
	return vec3d(x,y,z);
}

void
Simulation::outputDataHeader(ofstream &fout){
	fout << "# LF_DEM version " << GIT_VERSION << endl;
	fout << "# np " << sys.get_np() << endl;
	fout << "# VF " << sys.volume_fraction << endl;
	fout << "# Lx " << sys.get_lx() << endl;
	fout << "# Ly " << sys.get_ly() << endl;
	fout << "# Lz " << sys.get_lz() << endl;
}

void
Simulation::outputConfigurationData(){
	vector<vec3d> pos;
	vector<vec3d> vel;
	int np = sys.get_np();
	pos.resize(np);
	vel.resize(np);
	for (int i=0; i<np; i++) {
		pos[i] = shiftUpCoordinate(sys.position[i].x-sys.Lx_half(),
								   sys.position[i].y-sys.Ly_half(),
								   sys.position[i].z-sys.Lz_half());
	}
	/* If the origin is shifted,
	 * we need to change the velocities of particles as well.
	 */
	if (p.origin_zero_flow) {
		for (int i=0; i<np; i++) {
			vel[i] = sys.velocity[i];
			if (pos[i].z < 0) {
				vel[i].x -= sys.get_lz();
			}
		}
	}
	/*
	 * shear_disp = sys.strain() - (int)(sys.strain()/Lx)*Lx
	 */
	if (p.out_data_particle) {
		fout_particle << "# " << sys.get_shear_strain() << ' ';
		fout_particle << sys.shear_disp << ' ';
		fout_particle << sys.dimensionless_shear_rate << ' ';
		fout_particle << sys.target_stress_input << endl;
		for (int i=0; i<np; i++) {
			vec3d &r = pos[i];
			vec3d &v = vel[i];
			vec3d &o = sys.ang_velocity[i];
			double lub_xzstress = sys.lubstress[i].getStressXZ();
			double contact_xzstressGU = sys.contactstressGU[i].getStressXZ();
			double brownian_xzstressGU = 0;
			if (sys.brownian) {
				brownian_xzstressGU = sys.brownianstressGU[i].getStressXZ();
			}
			fout_particle << i; //1: number
			fout_particle << ' ' << sys.radius[i]; //2: radius
			fout_particle << ' ' << r.x << ' ' << r.y << ' ' << r.z; //3, 4, 5: position
			fout_particle << ' ' << v.x << ' ' << v.y << ' ' << v.z; //6, 7, 8: velocity
			fout_particle << ' ' << o.x << ' ' << o.y << ' ' << o.z; //9, 10, 11: angular velocity
			fout_particle << ' ' << 6*M_PI*lub_xzstress; //12: xz stress contributions
			fout_particle << ' ' << 6*M_PI*contact_xzstressGU; //13: xz stress contributions
			fout_particle << ' ' << 6*M_PI*brownian_xzstressGU; //14: xz stress contributions
			if (sys.twodimension) {
				fout_particle << ' ' << sys.angle[i]; // 15
			}
			fout_particle << endl;
		}
	}
	int cnt_interaction = 0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_active()) {
			cnt_interaction++;
		}
	}
	if (p.out_data_interaction) {
		fout_interaction << "# " << sys.get_shear_strain();
		fout_interaction << ' ' << cnt_interaction << endl;
		for (int k=0; k<sys.nb_interaction; k++) {
			if (sys.interaction[k].is_active()) {
				unsigned short i, j;
				sys.interaction[k].get_par_num(i, j);
				vec3d nr_vec = sys.interaction[k].get_nvec();
				StressTensor stress_contact = sys.interaction[k].contact.getContactStressXF();
				fout_interaction << i << ' ' << j << ' '; // 1, 2
				/* contact.state:
				 * 0 no contact
				 * 1 Friction is not activated (critical load model)
				 * 2 Static friction
				 * 3 Sliding
				 */
				fout_interaction << sys.interaction[k].contact.state << ' '; //3
				fout_interaction << nr_vec.x << ' '; // 4
				fout_interaction << nr_vec.y << ' '; // 5
				fout_interaction << nr_vec.z << ' '; // 6
				fout_interaction << sys.interaction[k].get_reduced_gap() << ' '; // 7
				/* [NOTE]
				 * Lubrication forces are reference values
				 * in the Brownian case. The force balancing
				 * velocities are recalculated without
				 * including the Brownian forces.
				 * It seems there is no better way to visualize
				 * the lubrication forces.
				 */
				fout_interaction << sys.interaction[k].lubrication.get_lubforce_normal() << ' '; // 8
				fout_interaction << sys.interaction[k].lubrication.get_lubforce_tan() << ' '; // 9, 10, 11
				/*
				 * Contact forces include only spring forces.
				 */
				fout_interaction << sys.interaction[k].contact.get_f_contact_normal_norm() << ' '; // 12
				fout_interaction << sys.interaction[k].contact.get_f_contact_tan() << ' '; // 13, 14, 15
				fout_interaction << sys.interaction[k].repulsion.getForceNorm() << ' '; // 16
				fout_interaction << 6*M_PI*stress_contact.getStressXZ() << ' '; // 17
				sys.interaction[k].contact.addUpContactForceTorque();
				//fout_interaction << 6*M_PI*stress_contact.getNormalStress1() << ' ';
				//fout_interaction << 6*M_PI*stress_contact.getNormalStress2() << ' ';
				fout_interaction << endl;
			}
		}
	}
}

void
Simulation::outputFinalConfiguration(){
	ofstream fout_finalconfig;
	string filename_final_configuration = "./after_relax/"+filename_import_positions;
	fout_finalconfig.open(filename_final_configuration.c_str());
	fout_finalconfig << import_line[0] << endl;
	fout_finalconfig << import_line[1] << endl;
	int np = sys.get_np();
	for (int i=0; i<np; i++) {
		fout_finalconfig << sys.position[i].x << ' ';
		fout_finalconfig << sys.position[i].y << ' ';
		fout_finalconfig << sys.position[i].z << ' ';
		fout_finalconfig << sys.radius[i] << endl;
	}

	string filename_bin = filename_final_configuration;
	string ext=".dat";
	size_t start_pos = filename_bin.find(ext);
    if(start_pos == string::npos){
		cerr << " WARNING, no binary output generated " << endl;
        return;
	}
    filename_bin.replace(start_pos, ext.length(), ".bin");
	outputConfigurationBinary(filename_bin);
}


