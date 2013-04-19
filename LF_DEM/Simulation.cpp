//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

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

/*
 * Main simulation
 */
void
Simulation::SimulationMain(int argc, const char * argv[]){
	filename_import_positions = argv[1];
	importInitialPositionFile();
	setDefaultParameters();
	sys.np(num_of_particle);
	filename_parameters = argv[2];
	readParameterFile();
	sys.dimensionless_shear_rate = atof(argv[3]);
	openOutputFiles();
	sys.setupSystem(initial_position, radius);
	outputDataHeader(fout_particle);
	outputConfigurationData();
	sys.setupShearFlow(true);
	double strain_next_config_out = strain_interval_output;
	do {
		sys.timeEvolution(strain_interval_output_data);
		evaluateData();
		outputRheologyData();
		if (sys.strain() >= strain_next_config_out) {
			outputConfigurationData();
			strain_next_config_out = sys.strain()+strain_interval_output-1e-6;
		}
		if (kn_kt_adjustment) {
			sys.adjustContactModelParameters(10);
		}
		cerr << "strain: " << sys.strain() << endl;
	} while (strain_next_config_out < shear_strain_end);
}

void
Simulation::RelaxationZeroShear(vector<vec3d> &positions,
								 vector<double> &radii,
								 double lx,
								 double ly,
								 double lz){
	num_of_particle = positions.size();
	sys.np(num_of_particle);
	sys.lx(lx);
	sys.ly(ly);
	sys.lz(lz);
	if (ly == 0) {
		sys.dimension = 2;
	} else {
		sys.dimension = 3;
	}
	double max_radius = 0;
	for (int i=0; i<num_of_particle; i++) {
		if (max_radius < radii[i]) {
			max_radius = radii[i];
		}
	}
	sys.setRadiusMax(max_radius);
	setDefaultParameters();
	sys.integration_method = 0;
	sys.dt = 1e-4;
	sys.kn = 2000;
	sys.mu_static = 0;
	sys.dimensionless_shear_rate = 1;
	sys.colloidalforce_length = 0.05; // dimensionless
	sys.colloidalforce_amplitude = 10;
	sys.setupSystem(positions, radii);
	sys.setupShearFlow(false);
	double energy_previous = 0;
	while (true) {
		int i_time_interval = 1000;
		sys.timeEvolutionRelax(i_time_interval);
		evaluateData();
		sys.calcTotalPotentialEnergy();
		cout << sys.strain() << ' ' << sys.min_gap_nondim << ' ' << sys.total_energy << endl;
		cerr << energy_previous-sys.total_energy << endl;
		if (sys.min_gap_nondim > 0 &&
			energy_previous-sys.total_energy < 0.01) {
			cerr << "finish" << endl;
			break;
		}
		energy_previous = sys.total_energy;
		//		outputRheologyData();
		//		outputConfigurationData();
	}
	for (int i=0; i<num_of_particle; i++) {
		positions[i] = sys.position[i];
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
		sys.bgf_factor = atof(value.c_str());
	} else if (keyword == "kn_kt_adjustment") {
		kn_kt_adjustment = str2bool(value);
	} else if (keyword == "colloidalforce_length") {
		sys.colloidalforce_length = atof(value.c_str());
	} else if (keyword == "lub_reduce_parameter") {
		sys.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxzation_time") {
		sys.contact_relaxzation_time = atof(value.c_str());
	} else if (keyword == "kb_T") {
		sys.kb_T = atof(value.c_str());
	} else if (keyword == "dt") {
		sys.dt = atof(value.c_str());
	} else if (keyword == "shear_strain_end") {
		shear_strain_end = atof(value.c_str());
	} else if (keyword == "integration_method") {
		sys.integration_method = atof(value.c_str());
	} else if (keyword == "lub_max") {
		sys.lub_max = atof(value.c_str());
	} else if (keyword == "kn") {
		sys.kn = atof(value.c_str());
	} else if (keyword == "kt") {
		sys.kt = atof(value.c_str());
	} else if (keyword == "mu_static") {
		sys.mu_static = atof(value.c_str()) ;
	} else if (keyword == "strain_interval_out") {
		strain_interval_output = atof(value.c_str());
	} else if (keyword == "strain_interval_out_data") {
		strain_interval_output_data = atof(value.c_str());
	} else if (keyword == "draw_rotation_2d") {
		sys.draw_rotation_2d = str2bool(value);
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
	fout_particle.open(particle_filename.c_str());
	fout_interaction.open(interaction_filename.c_str());
	fout_rheo.open(vel_filename.c_str());
	sys.openFileInteractionData();
}

void
Simulation::setDefaultParameters(){
	kn_kt_adjustment = false;
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
	sys.dt = 1e-4;
	/*
	 * integration_method:
	 * 0 Euler's Method,
	 * 1 predictor-corrector,
	 * 2 Brownian (if kT > 0).
	 */
	sys.integration_method = 1;
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
	sys.lub_max = 2.5;
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
	sys.bgf_factor = 1;
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
	sys.kn = 5000;
	sys.kt = 1000;
	sys.overlap_target = 0.03;
	sys.disp_tan_target = 0.03;
	/*
	 * Colloidal force parameter
	 * Short range repulsion is assumed.
	 * cf_amp_dl0: cf_amp_dl at shearrate = 1
	 */
	sys.colloidalforce_length = 0.1;
	/*
	 * mu_static: static friction coeffient
	 * mu_dynamic: dynamic friction coeffient
	 */
	sys.mu_static = 1;	
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
	double lx_, ly_, lz_;
	double vf1, vf2;
	char buf;
	file_import >> buf >> n1 >> n2 >> volume_fraction_ >> lx_ >> ly_ >> lz_ >> vf1 >> vf2;
	num_of_particle = n1+n2;
	if (ly_ == 0) {
		sys.dimension = 2;
	} else {
		sys.dimension = 3;
	}
	sys.lx(lx_);
	sys.ly(ly_);
	sys.lz(lz_);
	cerr << "box: " << lx_ << ' ' <<  ly_ << ' ' << lz_ << endl;
	sys.volume_fraction = volume_fraction_;
	initial_position.resize(num_of_particle);
	radius.resize(num_of_particle);
	double _radius;
	vec3d _pos;
	for (int i=0; i<num_of_particle ; i++) {
		file_import >> _pos.x >> _pos.y >> _pos.z >> _radius;
		initial_position[i] = _pos;
		radius[i] = _radius;
	}
	file_import.close();
	double max_radius = 0;
	for (int i=0; i<num_of_particle; i++) {
		if (max_radius < radius[i]) {
			max_radius = radius[i];
		}
	}
	sys.setRadiusMax(max_radius);
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
	StressTensor total_stress;
	total_stress = sys.total_hydro_stress;
	total_stress += sys.total_contact_stressXF_normal;
	total_stress += sys.total_contact_stressXF_tan;
	total_stress += sys.total_colloidal_stressXF;
	total_stress += sys.total_colloidal_stressGU;
	if (sys.brownian) {
		total_stress += sys.total_brownian_stress;
	}
	StressTensor total_contact_stressXF = sys.total_contact_stressXF_normal+sys.total_contact_stressXF_tan;
	Viscosity = total_stress.getStressXZ();
	N1 = total_stress.getNormalStress1();
	N2 = total_stress.getNormalStress2();
	Viscosity_h = sys.total_hydro_stress.getStressXZ();
	N1_h = sys.total_hydro_stress.getNormalStress1();
	N2_h = sys.total_hydro_stress.getNormalStress2();
	Viscosity_cont_XF = total_contact_stressXF.getStressXZ();
	N1_cont_XF = total_contact_stressXF.getNormalStress1();
	N2_cont_XF = total_contact_stressXF.getNormalStress2();
	Viscosity_friction = sys.total_contact_stressXF_tan.getStressXZ();
	N1_friction = sys.total_contact_stressXF_tan.getNormalStress1();
	N2_friction = sys.total_contact_stressXF_tan.getNormalStress2();
	Viscosity_cont_GU = sys.total_contact_stressGU.getStressXZ();
	N1_cont_GU = sys.total_contact_stressGU.getNormalStress1();
	N2_cont_GU = sys.total_contact_stressGU.getNormalStress2();
	Viscosity_col_XF = sys.total_colloidal_stressXF.getStressXZ();
	N1_col_XF = sys.total_colloidal_stressXF.getNormalStress1();
	N2_col_XF = sys.total_colloidal_stressXF.getNormalStress2();
	Viscosity_col_GU = sys.total_colloidal_stressGU.getStressXZ();
	N1_col_GU = sys.total_colloidal_stressGU.getNormalStress1();
	N2_col_GU = sys.total_colloidal_stressGU.getNormalStress2();
	Viscosity_brownian = sys.total_brownian_stress.getStressXZ();
	N1_brownian = sys.total_brownian_stress.getNormalStress1();
	N2_brownian = sys.total_brownian_stress.getNormalStress2();
}

void
Simulation::outputRheologyData(){
	static bool firsttime = true;
	/*
	 * Output the sum of the normal forces.
	 *
	 *  Viscosity = S_{xz} / shear_rate
	 *  N1 = S_{xx}-S_{zz}
	 *  N2 = S_{zz}-S_{yy} = S_zz-(-S_xx-S_zz) = S_xx+2*S_zz
	 * 
	 * Relative viscosity = Viscosity / viscosity_solvent
	 */
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
		fout_rheo << "#29: max Fc_normal_norm" << endl;
		fout_rheo << "#30: max velocity" << endl;
		fout_rheo << "#31: max angular velocity" << endl;
		fout_rheo << "#32: max contact normal velocity" << endl;
		fout_rheo << "#33: max contact tangential velocity" << endl;
		fout_rheo << "#34: average contact number per particle" << endl;
		fout_rheo << "#35: kn" << endl;
		fout_rheo << "#36: kt" << endl;
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
	fout_rheo << sys.strain() << ' '; //1
	fout_rheo << 6*M_PI*Viscosity << ' '; //2
	fout_rheo << 6*M_PI*N1 << ' '; //3
	fout_rheo << 6*M_PI*N2 << ' '; //4
	fout_rheo << 6*M_PI*Viscosity_h << ' '; //5
	fout_rheo << 6*M_PI*N1_h << ' '; //6
	fout_rheo << 6*M_PI*N2_h << ' '; //7
	fout_rheo << 6*M_PI*Viscosity_cont_XF << ' '; //8
	fout_rheo << 6*M_PI*N1_cont_XF << ' '; //9
	fout_rheo << 6*M_PI*N2_cont_XF << ' '; //10
	fout_rheo << 6*M_PI*Viscosity_cont_GU << ' ' ; //11
	fout_rheo << 6*M_PI*N1_cont_GU << ' ' ; //12
	fout_rheo << 6*M_PI*N2_cont_GU << ' ' ; //13
	fout_rheo << 6*M_PI*Viscosity_friction << ' '; //14
	fout_rheo << 6*M_PI*N1_friction << ' '; //15
	fout_rheo << 6*M_PI*N2_friction  << ' '; //16
	fout_rheo << 6*M_PI*Viscosity_col_XF << ' '; //17
	fout_rheo << 6*M_PI*N1_col_XF << ' '; //18
	fout_rheo << 6*M_PI*N2_col_XF << ' '; //19
	fout_rheo << 6*M_PI*Viscosity_col_GU << ' '; //20
	fout_rheo << 6*M_PI*N1_col_GU << ' '; //21
	fout_rheo << 6*M_PI*N2_col_GU << ' '; //22
	fout_rheo << 6*M_PI*Viscosity_brownian << ' ' ; //23
	fout_rheo << 6*M_PI*N1_brownian << ' ' ; //24
	fout_rheo << 6*M_PI*N2_brownian << ' ' ; //25
	fout_rheo << sys.min_gap_nondim << ' '; //26
	fout_rheo << sys.max_disp_tan << ' '; //27
	fout_rheo << sys.average_Fc_normal_norm << ' '; //28
	fout_rheo << sys.max_Fc_normal_norm << ' '; //29
	fout_rheo << sys.max_velocity << ' '; //30
	fout_rheo << sys.max_ang_velocity << ' '; //31
	fout_rheo << sys.max_contact_velo_normal << ' '; //32
	fout_rheo << sys.max_contact_velo_tan << ' '; //33
	fout_rheo << sys.getParticleContactNumber() << ' '; //34
	fout_rheo << sys.kn << ' '; //35
	fout_rheo << sys.kt << ' '; //36
	fout_rheo << endl;
}

vec3d
Simulation::shiftUpCoordinate(double x, double y, double z){
	if (origin_zero_flow) {
		z += sys.lz_half();
		if (z > sys.lz_half()) {
			x -= sys.shear_disp;
			if (x < -sys.lx_half()) {
				x += sys.lx();
			}
			z -= sys.lz();
		}
	}
	return vec3d(x,y,z);
}

void
Simulation::outputDataHeader(ofstream &fout){
	char sp = ' ';
	fout << "np" << sp << sys.np() << endl;
	fout << "VF" << sp << sys.volume_fraction << endl;
	fout << "Lx" << sp << sys.lx() << endl;
	fout << "Ly" << sp << sys.ly() << endl;
	fout << "Lz" << sp << sys.lz() << endl;
}

void
Simulation::outputConfigurationData(){
	vector<vec3d> pos;
	vector<vec3d> vel;
	char sp = ' ';
	int np = sys.np();
	pos.resize(np);
	vel.resize(np);
	for (int i=0; i < np; i++) {
		pos[i] = shiftUpCoordinate(sys.position[i].x-sys.lx_half(),
								   sys.position[i].y-sys.ly_half(),
								   sys.position[i].z-sys.lz_half());
	}
	/* If the origin is shifted,
	 * we need to change the velocities of particles as well.
	 */
	if (origin_zero_flow) {
		for (int i=0; i < np; i++) {
			vel[i] = sys.velocity[i];
			if (pos[i].z < 0) {
				vel[i].x -= sys.lz();
			}
		}
	}
	/*
	 * shear_disp = sys.strain() - (int)(sys.strain()/Lx)*Lx
	 */
	fout_particle << "#" << sp << sys.strain() << sp;
	fout_particle << sys.shear_disp << sp << sys.dimensionless_shear_rate << endl;
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
		fout_particle << i << sp; //1: number
		fout_particle << sys.radius[i] << sp; //2: radius
		fout_particle << p.x << sp << p.y << sp << p.z << sp; //3, 4, 5: position
		fout_particle << v.x << sp << v.y << sp << v.z << sp; //6, 7, 8: velocity
		fout_particle << o.x << sp << o.y << sp << o.z << sp; //9, 10, 11: angular velocity
		fout_particle << 6*M_PI*lub_xzstress << sp; //12: xz stress contributions
		fout_particle << 6*M_PI*contact_xzstressGU << sp; //13: xz stress contributions
		fout_particle << 6*M_PI*brownian_xzstress << sp; //14: xz stress contributions
		if (sys.dimension == 2) {
			fout_particle << sys.angle[i] << sp; // 15
		}
		fout_particle << endl;
	}
	int cnt_interaction = 0;
	for (int k=0; k<sys.num_interaction; k++) {
		if (sys.interaction[k].active) {
			cnt_interaction++;
		}
	}
	fout_interaction << "#" << sp << sys.strain() << sp ;
	fout_interaction << cnt_interaction << endl;
	for (int k=0; k<sys.num_interaction; k++) {
		if (sys.interaction[k].active) {
			vec3d fc_tan = sys.interaction[k].getFcTan();
			StressTensor stress_contact = sys.interaction[k].getContactStressXF();
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
			fout_interaction << sys.interaction[k].par_num[0] << sp; // 1
			fout_interaction << sys.interaction[k].par_num[1] << sp; // 2
			fout_interaction << sys.interaction[k].contact << sp; // 3
			fout_interaction << sys.interaction[k].nr_vec.x << sp; // 4
			fout_interaction << sys.interaction[k].nr_vec.y << sp; // 5
			fout_interaction << sys.interaction[k].nr_vec.z << sp; // 6
			fout_interaction << sys.interaction[k].gap_nondim() << sp; // 7
			fout_interaction << sys.interaction[k].getLubForce() << sp; // 8
			fout_interaction << sys.interaction[k].getFcNormal() << sp; // 9
			fout_interaction << sys.interaction[k].getFcTan_norm() << sp; // 10
			fout_interaction << sys.interaction[k].getColloidalForce() << sp; // 11
			fout_interaction << 6*M_PI*stress_contact.getStressXZ() << sp; // 12
			fout_interaction << 6*M_PI*stress_contact.getNormalStress1() << sp; // 13
			fout_interaction << 6*M_PI*stress_contact.getNormalStress2() << sp; // 14
			fout_interaction << endl;
		}
	}
}
