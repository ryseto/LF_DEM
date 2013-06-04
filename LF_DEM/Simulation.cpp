//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//
#define _USE_MATH_DEFINES
#define VERSION "1.1"
#include "Simulation.h"
#include <cmath>
#include <map>
#include <string>
#include <algorithm>
#include <cctype>
Simulation::Simulation(){};
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
	double phi_;
	double kn_;
	double kt_;
	double dt_max_;
	while (fin_knktdt >> phi_ >> kn_ >> kt_ >> dt_max_) {
		if (phi_ == volume_fraction){
			break;
		}
	}
	fin_knktdt.clear();
	sys.Kn(kn_);
	sys.Kt(kt_);
	sys.set_dt_max(dt_max_);
	cerr << phi_ << ' ' << kn_ << ' ' << kt_ << ' ' << dt_max_ << endl;
}

/*
 * Main simulation
 */
void
Simulation::simulationMain(int argc, const char * argv[]){
	filename_import_positions = argv[1];
	filename_parameters = argv[2];
	sys.dimensionless_shear_rate = atof(argv[3]);
	setDefaultParameters();
	readParameterFile();
	importInitialPositionFile();
	if (argc == 5) {
		contactForceParameter(argv[4]);
	}
	openOutputFiles();
	outputDataHeader(fout_particle);
	outputDataHeader(fout_interaction);
	outputDataHeader(fout_rheo);
	outputDataHeader(fout_st);
	sys.setupSystem();
	outputConfigurationData();
	sys.setupShearFlow(true);
	int cnt_simu_loop = 1;
	int cnt_knkt_adjustment = 1;
	int cnt_config_out = 1;
	double strain_next_config_out = 0;
	while (sys.Shear_strain() < shear_strain_end-1e-8) {
		double strain_knkt_adjustment = cnt_knkt_adjustment*strain_interval_knkt_adjustment;
		strain_next_config_out = cnt_config_out*strain_interval_output;
		double strain_next = (cnt_simu_loop)*strain_interval_output_data;
		sys.timeEvolution(strain_next);
		evaluateData();
		outputRheologyData();
		outputStressTensorData();
		if (sys.Shear_strain() >= strain_next_config_out-1e-8) {
			outputConfigurationData();
			cnt_config_out ++;
		}
		if (kn_kt_adjustment) {
			if (sys.Shear_strain() >= strain_knkt_adjustment-1e-8) {
				sys.adjustContactModelParameters();
				cnt_knkt_adjustment ++;
			}
		}
		cnt_simu_loop ++;
		cerr << "strain: " << sys.Shear_strain() << " / " << shear_strain_end << endl;
	}
}

void
Simulation::relaxationZeroShear(vector<vec3d> &position_,
								vector<double> &radius_,
								double lx_, double ly_, double lz_){
	sys.setConfiguration(position_, radius_, lx_, ly_, lz_);
	setDefaultParameters();
	sys.Integration_method(0);
	sys.Dt(1e-4);
	sys.Kn(2000);
	sys.kb_T = 0;
	sys.Mu_static(0);
	sys.dimensionless_shear_rate = 1;
	sys.Colloidalforce_length(0); // dimensionless
	sys.setupSystem();
	sys.Colloidalforce_amplitude(0);
	sys.setupShearFlow(false);
	double energy_previous = 0;
	while (true) {
		int i_time_interval = 1000;
		sys.timeEvolutionRelax(i_time_interval);
		evaluateData();
		//sys.analyzeState();
		sys.calcTotalPotentialEnergy();
		cout << sys.min_gap_nondim << ' ' << sys.total_energy << endl;
		cerr << energy_previous-sys.total_energy << endl;
		if (sys.min_gap_nondim > 0 &&
			energy_previous-sys.total_energy < 0.01) {
			cerr << "finish" << endl;
			break;
		}
		energy_previous = sys.total_energy;
	}
	for (int i=0; i<sys.Np(); i++) {
		position_[i] = sys.position[i];
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
Simulation::autoSetParameters(const string &keyword,
							  const string &value){
	if (keyword == "bgf_factor") {
		sys.Bgf_factor(atof(value.c_str()));
	} else if (keyword == "kn_kt_adjustment") {
		kn_kt_adjustment = str2bool(value);
	} else if (keyword == "strain_interval_knkt_adjustment") {
		strain_interval_knkt_adjustment = atof(value.c_str());
	} else if (keyword == "colloidalforce_length") {
		sys.Colloidalforce_length(atof(value.c_str()));
	} else if (keyword == "lub_reduce_parameter") {
		sys.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxzation_time") {
		sys.contact_relaxzation_time = atof(value.c_str());
	} else if (keyword == "kb_T") {
		sys.kb_T = atof(value.c_str());
	} else if (keyword == "dt_max") {
		sys.set_dt_max(atof(value.c_str()));
	} else if (keyword == "disp_max") {
		sys.set_disp_max(atof(value.c_str()));
	} else if (keyword == "shear_strain_end") {
		shear_strain_end = atof(value.c_str());
	} else if (keyword == "integration_method") {
		sys.Integration_method(atoi(value.c_str()));
	} else if (keyword == "lub_max") {
		sys.Lub_max(atof(value.c_str()));
	} else if (keyword == "kn") {
		sys.Kn(atof(value.c_str()));
	} else if (keyword == "kt") {
		sys.Kt(atof(value.c_str()));
	} else if (keyword == "mu_static") {
		sys.Mu_static(atof(value.c_str()));
	} else if (keyword == "strain_interval_out") {
		strain_interval_output = atof(value.c_str());
	} else if (keyword == "strain_interval_out_data") {
		strain_interval_output_data = atof(value.c_str());
	} else if (keyword == "out_data_particle") {
		out_data_particle = str2bool(value);
	} else if (keyword == "out_data_interaction") {
		out_data_interaction = str2bool(value);
	} else if (keyword == "origin_zero_flow") {
		origin_zero_flow = str2bool(value);
	} else if (keyword == "overlap_target") {
		sys.overlap_target = atof(value.c_str());
	} else if (keyword == "disp_tan_target") {
		sys.disp_tan_target = atof(value.c_str());
	} else {
		cerr << "keyword " << keyword << " is not associated with an parameter" << endl;
		exit(1);
	}
}

void
Simulation::readParameterFile(){
	ifstream fin;
	fin.open(filename_parameters.c_str());
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
		} while(true);
		if (begin_comment > end_comment ) {
			cerr << str_parameter.find("/*") << endl;
			cerr << str_parameter.find("*/") << endl;
			cerr << "syntax error in the parameter file." << endl;
			exit(1);
		}
		string::size_type pos_slashslash = str_parameter.find("//");
		if( pos_slashslash != string::npos) {
			cerr << " // is not syntax for comment out." << endl;
			exit(1);
		}
		Str2KeyValue(str_parameter, keyword, value);
		autoSetParameters(keyword, value);
	}
	fin.close();
	return;
}

void
Simulation::openOutputFiles(){
	/*
	 * Set simulation name and name of output files.
	 */
	prepareSimulationName();
	string particle_filename = "par_" + sys.simu_name + ".dat";
	string interaction_filename = "int_" + sys.simu_name + ".dat";
	string vel_filename = "rheo_" + sys.simu_name + ".dat";
	string st_filename = "st_" +sys.simu_name + ".dat";
	fout_particle.open(particle_filename.c_str());
	fout_interaction.open(interaction_filename.c_str());
	fout_rheo.open(vel_filename.c_str());
	fout_st.open(st_filename.c_str());
	sys.openFileInteractionData();
}

void
Simulation::setDefaultParameters(){
	/*
	 * Simulation
	 *
	 * dt: the time step to integrate the equation of motion.
	 *     We need to give a good criterion to give.
	 * dt_mid: the intermediate time step for the mid-point
	 *     algortithm. dt/dt_mid = dt_ratio
	 *     Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 *    ASD code from Brady has dt_ratio=150
	 *
	 */
	double _dt = 1e-4;
	/*
	 * integration_method:
	 * 0 Euler's Method,
	 * 1 predictor-corrector,
	 * 2 Brownian (if kT > 0).
	 */
	int _integration_method = 1;
	/*
	 * Shear flow
	 *  shear_rate: shear rate
	 *  strain(): total strain (length of simulation)
	 *
	 */
	shear_strain_end = 10;
	/*
	 * Lubrication force
	 * lub_max: reduced large cutoff distance for lubrication
	 * I think lub_max = 2.5 and 3 generate different results.
	 * We should give suffiently larger value.
	 * The value 3 or 3.5 should be better (To be checked.)
	 */
	double _lub_max = 2.5;
	/*
	 * gap_nondim_min: gives reduced lubrication (maximum coeeffient).
	 *
	 */
	sys.lub_reduce_parameter = 1e-3;
	/*
	 * contact_relaxzation_factor:
	 *
	 * This gives the coeffient of the resistance term for h < 0.
	 * - If the value is negative, the value of 1/lub_reduce_parameter is used.
	 *
	 */
	sys.contact_relaxzation_time = 1e-3;
	/*
	 *  bgf_factor: background flow factor gives the weight between the one-body force and two-body force.
	 *   bgf_factor = 1.0 means full drag forces from undisturbed shear flow, that should be overestimate.
	 *   The optimal value of bgf_factor (< 1.0) may exist.
	 *
	 */
	double _bgf_factor = 1;
	/*
	 * Brownian force
	 * kb_T: Thermal energy kb*T
	 * kb_T = 0 ---> non-brownian
	 * kb_T > 0 ---> brownian
	 */
	sys.kb_T = 0;
	/*
	 * Contact force parameters
	 * kn: normal spring constant
	 * kt: tangential spring constant
	 */
	double _kn = 5000;
	double _kt = 1000;
	kn_kt_adjustment = false;
	strain_interval_knkt_adjustment = 5;
	sys.overlap_target = 0.03;
	sys.disp_tan_target = 0.03;
	
	/*
	 * Colloidal force parameter
	 * Short range repulsion is assumed.
	 * cf_amp_dl0: cf_amp_dl at shearrate = 1
	 */
	double _colloidalforce_length = 0.1;
	/*
	 * mu_static: static friction coeffient
	 * mu_dynamic: dynamic friction coeffient
	 */
	double _mu_static = 1;
	/*
	 * Output interval:
	 * strain_interval_output_data is for outputing rheo_...
	 * strain_interval_output is for outputing int_... and par_...
	 */
	strain_interval_output_data = 0.01;
	strain_interval_output = 0.05;
	/*
	 *  Data output
	 */
	/*
	 * The middle height of the simulation box is set to the flow zero level.
	 */
	origin_zero_flow = true;
	/*
	 * position and interaction data
	 */
	out_data_particle = true;
	out_data_interaction = true;
	
	sys.Integration_method(_integration_method);
	sys.Bgf_factor(_bgf_factor);
	sys.Lub_max(_lub_max);
	sys.set_dt_max(_dt);
	sys.Kn(_kn);
	sys.Kt(_kt);
	sys.Mu_static(_mu_static);
	sys.Colloidalforce_length(_colloidalforce_length);
}

void
Simulation::importInitialPositionFile(){
	fstream file_import;
	file_import.open(filename_import_positions.c_str());
	if(!file_import) {
		cerr << " Position file '" << filename_import_positions << "' not found." <<endl;
		exit(1);
	}
	string line;
	getline(file_import, line);
	int n1, n2;
	double volume_fraction_;
	double lx, ly, lz;
	double vf1, vf2;
	char buf;
	file_import >> buf >> n1 >> n2 >> volume_fraction_ >> lx >> ly >> lz >> vf1 >> vf2;
	volume_fraction = volume_fraction_;
	double x_, y_, z_, a_;
	while (file_import >> x_ >> y_ >> z_ >> a_) {
		//initial_position.push_back(vec3d (x_,y_,z_));
		vec3d tmp(x_, y_, z_);
		initial_position.push_back(tmp);
		radius.push_back(a_);
	}
	file_import.close();
	sys.setConfiguration(initial_position, radius, lx, ly, lz);
}

void
Simulation::prepareSimulationName(){
	ostringstream ss_simu_name;
	string::size_type pos_ext_position = filename_import_positions.find(".dat");
	string::size_type pos_ext_parameter = filename_parameters.find(".txt");
	ss_simu_name << filename_import_positions.substr(0, pos_ext_position);
	ss_simu_name << "_";
	ss_simu_name << filename_parameters.substr(0, pos_ext_parameter);
	ss_simu_name << "_sr" << sys.dimensionless_shear_rate;
	sys.simu_name = ss_simu_name.str();
	cerr << sys.simu_name << endl;
}

void
Simulation::evaluateData(){
	sys.calcStress();
	sys.analyzeState();
	/* NOTE:
	 * The total stress does not include contact GU term
	 * due to the asumption of hard-sphere model.
	 */
	total_contact_stressXF = sys.total_contact_stressXF_normal+sys.total_contact_stressXF_tan;
	total_colloidal_stress = sys.total_colloidal_stressXF+sys.total_colloidal_stressGU;
	
	total_stress = sys.total_hydro_stress;
	total_stress += total_contact_stressXF;
	total_stress += total_colloidal_stress;
	
	if (sys.brownian) {
		total_stress += sys.total_brownian_stress;
	}
	
	viscosity = total_stress.getStressXZ()+5*volume_fraction/(12*M_PI);
	normalstress_diff_1 = total_stress.getNormalStress1();
	normalstress_diff_2 = total_stress.getNormalStress2();
	particle_pressure = total_stress.getParticlePressure();

	viscosity_hydro = sys.total_hydro_stress.getStressXZ();
	normalstress_diff_1_hydro = sys.total_hydro_stress.getNormalStress1();
	normalstress_diff_2_hydro = sys.total_hydro_stress.getNormalStress2();
	viscosity_cont_XF = total_contact_stressXF.getStressXZ();
	normalstress_diff_1_cont_XF = total_contact_stressXF.getNormalStress1();
	normalstress_diff_2_cont_XF = total_contact_stressXF.getNormalStress2();
	particle_pressure_cont = total_contact_stressXF.getParticlePressure();
	
	viscosity_friction = sys.total_contact_stressXF_tan.getStressXZ();
	normalstress_diff_1_friction = sys.total_contact_stressXF_tan.getNormalStress1();
	normalstress_diff_2_friction = sys.total_contact_stressXF_tan.getNormalStress2();
	viscosity_cont_GU = sys.total_contact_stressGU.getStressXZ();
	normalstress_diff_1_cont_GU = sys.total_contact_stressGU.getNormalStress1();
	normalstress_diff_2_cont_GU = sys.total_contact_stressGU.getNormalStress2();
	viscosity_col_XF = sys.total_colloidal_stressXF.getStressXZ();
	normalstress_diff_1_col_XF = sys.total_colloidal_stressXF.getNormalStress1();
	normalstress_diff_2_col_XF = sys.total_colloidal_stressXF.getNormalStress2();
	particle_pressure_col = sys.total_colloidal_stressXF.getParticlePressure();
	viscosity_col_GU = sys.total_colloidal_stressGU.getStressXZ();
	normalstress_diff_1_col_GU = sys.total_colloidal_stressGU.getNormalStress1();
	normalstress_diff_2_col_GU = sys.total_colloidal_stressGU.getNormalStress2();
	viscosity_brownian = sys.total_brownian_stress.getStressXZ();
	normalstress_diff_1_brownian = sys.total_brownian_stress.getNormalStress1();
	normalstress_diff_2_brownian = sys.total_brownian_stress.getNormalStress2();
}

void
Simulation::outputStressTensorData(){
	fout_st << sys.Shear_strain() << ' ';
	fout_st << 6*M_PI*viscosity << ' ';
	/* total_stress = sys.total_hydro_stress;
	 * + total_contact_stressXF + total_colloidal_stress;
	 */
	total_stress.outputStressTensor(fout_st);
	sys.total_hydro_stress.outputStressTensor(fout_st);
	total_contact_stressXF.outputStressTensor(fout_st);
	sys.total_contact_stressGU.outputStressTensor(fout_st);
	total_colloidal_stress.outputStressTensor(fout_st);
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
	static bool firsttime = true;
	if (firsttime) {
		firsttime = false;
		fout_rheo << "#1: shear strain" << endl;
		fout_rheo << "#2: Viscosity" << endl;
		fout_rheo << "#3: N1" << endl;
		fout_rheo << "#4: N2" << endl;
		fout_rheo << "#5: Viscosity(lub)" << endl;
		fout_rheo << "#6: N1(lub)" << endl;
		fout_rheo << "#7: N2(lub)" << endl;
		fout_rheo << "#8: Viscosity(xF_contact part)" << endl;
		fout_rheo << "#9: N1(xF_contact part)" << endl;
		fout_rheo << "#10: N2(xF_contact part)" << endl;
		fout_rheo << "#11: Viscosity(GU_contact part)" << endl;
		fout_rheo << "#12: N1(GU_contact part)" << endl;
		fout_rheo << "#13: N2(GU_contact part)" << endl;
		fout_rheo << "#14: Viscosity(friction)" << endl;
		fout_rheo << "#15: N1(friction)" << endl;
		fout_rheo << "#16: N2(friction)" << endl;
		fout_rheo << "#17: Viscosity(Colloidal force XF)" << endl;
		fout_rheo << "#18: N1(Colloidal force XF)" << endl;
		fout_rheo << "#19: N2(Colloidal force XF)" << endl;
		fout_rheo << "#20: Viscosity(Colloidal force GU)" << endl;
		fout_rheo << "#21: N1(Colloidal force GU)" << endl;
		fout_rheo << "#22: N2(Colloidal force GU)" << endl;
		fout_rheo << "#23: Viscosity(brownian)" << endl;
		fout_rheo << "#24: N1(brownian)" << endl;
		fout_rheo << "#25: N2(brownian)" << endl;
		fout_rheo << "#26: min gap (non-dim)" << endl;
		fout_rheo << "#27: max tangential displacement" << endl;
		fout_rheo << "#28: Average normal contact force" << endl;
		fout_rheo << "#29: max Fc_normal" << endl;
		fout_rheo << "#30: max velocity" << endl;
		fout_rheo << "#31: max angular velocity" << endl;
		fout_rheo << "#32: max contact normal velocity" << endl;
		fout_rheo << "#33: max contact tangential velocity" << endl;
		fout_rheo << "#34: average contact number per particle" << endl;
		fout_rheo << "#35: kn" << endl;
		fout_rheo << "#36: kt" << endl;
		fout_rheo << "#37: dt" << endl;
		fout_rheo << "#38: particle pressure" << endl;
		fout_rheo << "#39: particle pressure contact" << endl;
		fout_rheo << "#40: particle pressure colloidal force" << endl;
		
	}
	/*
	 * hat(...) indicates dimensionless quantities.
	 * (1) relative viscosity = Sxz/(eta0*shear_rate) = 6*pi*hat(Sxz)
	 * (2) N1/(eta0*shear_rate) = 6*pi*hat(N1)
	 * (3) N2/(eta0*shear_rate) = 6*pi*hat(N2)
	 *
	 * In simulation, we use the force unit where Stokes drag is F = -(U-U^inf)
	 *
	 */
	fout_rheo << sys.Shear_strain() << ' '; //1
	fout_rheo << 6*M_PI*viscosity << ' '; //2
	fout_rheo << 6*M_PI*normalstress_diff_1 << ' '; //3
	fout_rheo << 6*M_PI*normalstress_diff_2 << ' '; //4
	fout_rheo << 6*M_PI*viscosity_hydro << ' '; //5
	fout_rheo << 6*M_PI*normalstress_diff_1_hydro << ' '; //6
	fout_rheo << 6*M_PI*normalstress_diff_2_hydro << ' '; //7
	fout_rheo << 6*M_PI*viscosity_cont_XF << ' '; //8
	fout_rheo << 6*M_PI*normalstress_diff_1_cont_XF << ' '; //9
	fout_rheo << 6*M_PI*normalstress_diff_2_cont_XF << ' '; //10
	fout_rheo << 6*M_PI*viscosity_cont_GU << ' ' ; //11
	fout_rheo << 6*M_PI*normalstress_diff_1_cont_GU << ' ' ; //12
	fout_rheo << 6*M_PI*normalstress_diff_2_cont_GU << ' ' ; //13
	fout_rheo << 6*M_PI*viscosity_friction << ' '; //14
	fout_rheo << 6*M_PI*normalstress_diff_1_friction << ' '; //15
	fout_rheo << 6*M_PI*normalstress_diff_2_friction  << ' '; //16
	fout_rheo << 6*M_PI*viscosity_col_XF << ' '; //17
	fout_rheo << 6*M_PI*normalstress_diff_1_col_XF << ' '; //18
	fout_rheo << 6*M_PI*normalstress_diff_2_col_XF << ' '; //19
	fout_rheo << 6*M_PI*viscosity_col_GU << ' '; //20
	fout_rheo << 6*M_PI*normalstress_diff_1_col_GU << ' '; //21
	fout_rheo << 6*M_PI*normalstress_diff_2_col_GU << ' '; //22
	fout_rheo << 6*M_PI*viscosity_brownian << ' ' ; //23
	fout_rheo << 6*M_PI*normalstress_diff_1_brownian << ' ' ; //24
	fout_rheo << 6*M_PI*normalstress_diff_2_brownian << ' ' ; //25
	fout_rheo << sys.min_gap_nondim << ' '; //26
	fout_rheo << sys.max_disp_tan << ' '; //27
	fout_rheo << sys.average_fc_normal << ' '; //28
	fout_rheo << sys.max_fc_normal << ' '; //29
	fout_rheo << sys.max_velocity << ' '; //30
	fout_rheo << sys.max_ang_velocity << ' '; //31
	fout_rheo << sys.max_contact_velo_normal << ' '; //32
	fout_rheo << sys.max_contact_velo_tan << ' '; //33
	fout_rheo << sys.getParticleContactNumber() << ' '; //34
	fout_rheo << sys.Kn() << ' '; //35
	fout_rheo << sys.Kt() << ' '; //36
	fout_rheo << sys.Dt() << ' '; //37
	fout_rheo << 6*M_PI*particle_pressure << ' ';//38
	fout_rheo << 6*M_PI*particle_pressure_cont << ' ';//39
	fout_rheo << 6*M_PI*particle_pressure_col << ' ';//40
	fout_rheo << endl;
}

vec3d
Simulation::shiftUpCoordinate(double x, double y, double z){
	if (origin_zero_flow) {
		z += sys.Lz_half();
		if (z > sys.Lz_half()) {
			x -= sys.shear_disp;
			if (x < -sys.Lx_half()) {
				x += sys.Lx();
			}
			z -= sys.Lz();
		}
	}
	return vec3d(x,y,z);
}

void
Simulation::outputDataHeader(ofstream &fout){
	fout << "# LF_DEM version " << VERSION << endl;
	fout << "# np " << sys.Np() << endl;
	fout << "# VF " << volume_fraction << endl;
	fout << "# Lx " << sys.Lx() << endl;
	fout << "# Ly " << sys.Ly() << endl;
	fout << "# Lz " << sys.Lz() << endl;
}

void
Simulation::outputConfigurationData(){
	vector<vec3d> pos;
	vector<vec3d> vel;
	int np = sys.Np();
	pos.resize(np);
	vel.resize(np);
	for (int i=0; i < np; i++) {
		pos[i] = shiftUpCoordinate(sys.position[i].x-sys.Lx_half(),
								   sys.position[i].y-sys.Ly_half(),
								   sys.position[i].z-sys.Lz_half());
	}
	/* If the origin is shifted,
	 * we need to change the velocities of particles as well.
	 */
	if (origin_zero_flow) {
		for (int i=0; i<np; i++) {
			vel[i] = sys.velocity[i];
			if (pos[i].z < 0) {
				vel[i].x -= sys.Lz();
			}
		}
	}
	/*
	 * shear_disp = sys.strain() - (int)(sys.strain()/Lx)*Lx
	 */
	if (out_data_particle) {
		fout_particle << "# " << sys.Shear_strain() << ' ';
		fout_particle << sys.shear_disp << ' ';
		fout_particle << sys.dimensionless_shear_rate << endl;
		for (int i=0; i < np; i++) {
			vec3d &p = pos[i];
			vec3d &v = vel[i];
			vec3d &o = sys.ang_velocity[i];
			double lub_xzstress = sys.lubstress[i].getStressXZ();
			double contact_xzstressGU = sys.contactstressGU[i].getStressXZ();
			double brownian_xzstress = sys.brownianstress[i].getStressXZ();
			/* 1: number of the particle
			 * 2: radius
			 * 3, 4, 5: position
			 * 6, 7, 8: velocity
			 * 9, 10, 11: angular velocity
			 * 12: viscosity contribution of lubrication
			 * 13: viscosity contributon of contact GU xz
			 * 14: viscosity contributon of brownian xz
			 * (15: angle for 2D simulation)
			 */
			fout_particle << i; //1: number
			fout_particle << ' ' << sys.radius[i]; //2: radius
			fout_particle << ' ' << p.x << ' ' << p.y << ' ' << p.z; //3, 4, 5: position
			fout_particle << ' ' << v.x << ' ' << v.y << ' ' << v.z; //6, 7, 8: velocity
			fout_particle << ' ' << o.x << ' ' << o.y << ' ' << o.z; //9, 10, 11: angular velocity
			fout_particle << ' ' << 6*M_PI*lub_xzstress; //12: xz stress contributions
			fout_particle << ' ' << 6*M_PI*contact_xzstressGU; //13: xz stress contributions
			fout_particle << ' ' << 6*M_PI*brownian_xzstress; //14: xz stress contributions
			if (sys.dimension == 2) {
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
	if (out_data_interaction) {
		fout_interaction << "# " << sys.Shear_strain();
		fout_interaction << ' ' << cnt_interaction << endl;
		for (int k=0; k<sys.nb_interaction; k++) {
			if (sys.interaction[k].is_active()) {
				/* 1, 2: numbers of the interacting particles
				 * 3: 1=contact, 0=apart
				 * 4, 5, 6: normal vector
				 * 7: dimensionless gap = s - 2, s = 2r/(a1+a2)
				 * 8: lubrication force
				 * 9: Normal part of contact force
				 * 10: Tangential part of contact force
				 * 11: Colloidal force
				 * 12: Viscosity contribution of contact xF
				 * 13: N1 contribution of contact xF
				 * 14: N2 contribution of contact xF
				 */
				unsigned int i, j;
				sys.interaction[k].get_par_num(i, j);
				vec3d nr_vec = sys.interaction[k].Nr_vec();
				StressTensor stress_contact = sys.interaction[k].getContactStressXF();
				fout_interaction << i << ' ' << j << ' '; // 1, 2
				fout_interaction << sys.interaction[k].is_contact() << ' '; // 3
				fout_interaction << nr_vec.x << ' '; // 4
				fout_interaction << nr_vec.y << ' '; // 5
				fout_interaction << nr_vec.z << ' '; // 6
				fout_interaction << sys.interaction[k].Gap_nondim() << ' '; // 7
				fout_interaction << sys.interaction[k].getLubForce() << ' '; // 8
				fout_interaction << sys.interaction[k].getFcNormal() << ' '; // 9
				fout_interaction << sys.interaction[k].getFcTan() << ' '; // 10
				fout_interaction << sys.interaction[k].getColloidalForce() << ' '; // 11
				fout_interaction << 6*M_PI*stress_contact.getStressXZ() << ' '; // 12
				fout_interaction << 6*M_PI*stress_contact.getNormalStress1() << ' '; // 13
				fout_interaction << 6*M_PI*stress_contact.getNormalStress2() << endl; // 14
			}
		}
	}
}
