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
	if (fout_rheo.is_open()){
		fout_rheo.close();
	}
	if (fout_particle.is_open()){
		fout_particle.close();
	}
	if (fout_interaction.is_open()){
		fout_interaction.close();
	}
	if (fout_vpy.is_open()){
		fout_vpy.close();
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
	filename_parameters = argv[2];
	readParameterFile();
	if ( argc == 3){
		sys.shear_rate = 1.0;
		filename_addition = "sr1";
	} else {
		sys.shear_rate = atof(argv[3]);
		filename_addition = "sr";
		filename_addition += argv[3];
	}
	setUnits();
	
	openOutputFiles();
	sys.setupSystem(initial_positions, radii);
	outputDataHeader(fout_particle);
	int i_time_interval = strain_interval_out/sys.dt;

	outputConfigurationData();
	while(sys.strain() <= shear_strain_end){
		cerr << "strain: " << sys.strain() << endl;
		sys.timeEvolution(i_time_interval);
		evaluateData();
		outputRheologyData();
		outputConfigurationData();
	}
}

void
Simulation::setUnits(){
	unit_of_length = radius_of_particle; // = radius of smaller particle (a0)
	unit_of_velocity = sys.shear_rate*unit_of_length;
	unit_of_force = 6*M_PI*viscosity_solvent*radius_of_particle*unit_of_velocity;
	sys.cf_amp_dl = sys.cf_amp/unit_of_force;
	cerr << "cf_amp_dl = " << sys.cf_amp_dl << endl;
}

bool
str2bool(string value){
	if (value == "true"){
		return true;
	} else if (value == "false"){
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
	if (keyword == "lubrication"){
		sys.lubrication = str2bool(value);
	} else if (keyword == "friction"){
		sys.friction = str2bool(value);
	} else if (keyword == "bgf_factor"){
		sys.bgf_factor = atof(value.c_str());
	} else if (keyword == "cf_amp"){
		sys.cf_amp = atof(value.c_str());
	} else if (keyword == "cf_range"){
		sys.cf_range = atof(value.c_str());
	} else if (keyword == "lub_reduce_parameter"){
		sys.lub_reduce_parameter = atof(value.c_str());
	} else if (keyword == "contact_relaxzation_time"){
		sys.contact_relaxzation_time = atof(value.c_str());
//	} else if (keyword == "shear_rate"){
//		sys.shear_rate = atof(value.c_str());
	} else if (keyword == "kb_T"){
		sys.kb_T = atof(value.c_str());
	} else if (keyword == "dt"){
		sys.dt = atof(value.c_str());
	} else if (keyword == "shear_strain_end"){
		shear_strain_end = atof(value.c_str());
	} else if (keyword == "dt_ratio"){
		sys.dt_ratio = atof(value.c_str());
	} else if (keyword == "integration_method"){
		sys.integration_method = atof(value.c_str());
	} else if (keyword == "lub_max"){
		sys.lub_max = atof(value.c_str());
	} else if (keyword == "shearrate_scale_Fc_normal"){
		sys.shearrate_scale_Fc_normal = str2bool(value);
	} else if (keyword == "kn"){
		sys.kn = atof(value.c_str());
	} else if (keyword == "kt"){
		sys.kt = atof(value.c_str());
	} else if (keyword == "mu_static"){
		sys.mu_static = atof(value.c_str());
	} else if (keyword == "strain_interval_out") {
		strain_interval_out = atof(value.c_str());
	} else if (keyword == "draw_rotation_2d"){
		sys.draw_rotation_2d = str2bool(value);
	} else if (keyword == "out_data_particle"){
		out_data_particle = str2bool(value);
	} else if (keyword == "out_data_interaction"){
		out_data_interaction = str2bool(value);
	} else if (keyword == "origin_zero_flow"){
		origin_zero_flow = str2bool(value);
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
	while (!fin.eof()){
		string line;
		if (!getline(fin, line, ';'))
			break;
		if (fin.eof())
			break;
		string str_parameter;
		removeBlank(line);
		str_parameter = line;
		string::size_type begin_comment;
		string::size_type end_comment;
		do {
			begin_comment = str_parameter.find("/*");
			end_comment = str_parameter.find("*/");
			if (begin_comment > 10000 )
				break;
			str_parameter = str_parameter.substr(end_comment+2);
		}while (true);
		if (begin_comment > end_comment ){
			cerr << str_parameter.find("/*") << endl;
			cerr << str_parameter.find("*/") << endl;
			cerr << "syntax error in the parameter file." << endl;
			exit(1);
		}
		string::size_type pos_slashslash = str_parameter.find("//");
		if( pos_slashslash != string::npos){
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
	string vpy_filename = "vpy_" + sys.simu_name + ".dat";
	fout_particle.open(particle_filename.c_str());
	fout_interaction.open(interaction_filename.c_str());
	fout_rheo.open(vel_filename.c_str());
}

void
Simulation::setDefaultParameters(){
	/*
	 * If lubrication is false, it should be free-draining approximation.
	 * So far, we do not test this case well.
	 */
	sys.lubrication = true;
	/*
	 * If friction is false, the contact forces are only normal repulsion.
	 * It is important to find the difference between with and without friction.
	 */
	sys.friction = true;
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
	sys.dt_ratio = 2;
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
//	sys.shear_rate = 1.0;
	shear_strain_end = 10.;

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
	 * - The value
	 *
	 *
	 */
	sys.contact_relaxzation_time = 0.001;
	/*
	 *  bgf_factor: background flow factor gives the weight between the one-body force and two-body force.
	 *   bgf_factor = 1.0 means full drag forces from undisturbed shear flow, that should be overestimate.
	 *   The optimal value of bgf_factor (< 1.0) may exist.
	 *
	 */
	sys.bgf_factor = 1.0;
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
	sys.shearrate_scale_Fc_normal = true;
	sys.kn = 5000;
	sys.kt = 1000;
	/*
	 * mu_static: static friction coeffient
	 * mu_dynamic: dynamic friction coeffient
	 */
	sys.mu_static = 10;	
	/* Colloidal force
	 * Short range repulsion is assumed.
	 */
	sys.cf_amp = 1.0;
	sys.cf_range = 0.1;
	/*
	 * Output interval
	 */
	strain_interval_out = 0.01;
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

	/*
	 *
	 *
	 */
	viscosity_solvent = 1.0;
	radius_of_particle = 1.0;
}

void
Simulation::importInitialPositionFile(){
	fstream file_import;
	file_import.open( filename_import_positions.c_str());
	if(!file_import){
		cerr << " Position file '" << filename_import_positions << "' not found." <<endl;
		exit(1);
	}
	string line;
	getline(file_import, line);
	
	vec3d pos;
	double radius;
	int np_a_, np_b_;
	double volume_fraction_;
	double lx_, ly_, lz_;
	char buf;
	file_import >> buf >> np_a_ >> np_b_ >> volume_fraction_ >> lx_ >> ly_ >> lz_ ;
	np_a = np_a_;
	np_b = np_b_;

	int num_of_particle = np_a_ + np_b_;
	sys.np(num_of_particle);
	if (np_b_ > 0){
		sys.poly = true;
	}else{
		sys.poly = false;
	}
	cerr << "np = " << num_of_particle << endl;
	if (ly_ == 0){
		sys.dimension = 2;
	} else {
		sys.dimension = 3;
	}
	cerr << "dimension = " << sys.dimension << endl;
	sys.lx(lx_);
	sys.ly(ly_);
	sys.lz(lz_);
	cerr << "box: " << lx_ << ' ' <<  ly_ << ' ' << lz_ << endl;
	sys.volume_fraction = volume_fraction_;
	initial_positions.resize(num_of_particle);
	radii.resize(num_of_particle);
	double radius_max=1.0;
	for (int i = 0; i < num_of_particle ; i++){
		file_import >> pos.x >> pos.y >> pos.z >> radius;
		initial_positions[i] = pos;
		radii[i] = radius;
		if (radius > radius_max){
			radius_max = radius;
		}
	}
	sys.setRadiusMax(radius_max);
	sys.setSystemVolume();
	radius_a = radii[0];
	if (sys.poly ){
		radius_b = radii[np_a];
	} else {
		radius_b = 0;
	}
	file_import.close();
}

void
Simulation::prepareSimulationName(){
	ostringstream ss_simu_name;
	string::size_type pos_ext_position = filename_import_positions.find(".dat");
	string::size_type pos_ext_parameter = filename_parameters.find(".txt");
	ss_simu_name << filename_import_positions.substr(0, pos_ext_position);
	ss_simu_name << "_";
	ss_simu_name << filename_parameters.substr(0, pos_ext_parameter);
	ss_simu_name << filename_addition;
	sys.simu_name = ss_simu_name.str();
	cerr << sys.simu_name << endl;
}

void
Simulation::evaluateData(){
	sys.calcStress();
	sys.analyzeState();
	double total_stress[5];
	for (int u=0; u<5; u++){
		total_stress[u] = sys.total_hydro_stress[u];
		total_stress[u] += sys.total_contact_stressXF[u];
		total_stress[u] += sys.total_colloidal_stressXF[u];
		total_stress[u] += sys.total_colloidal_stressGU[u];
		// + sys.total_contact_stressGU[u];
	}
	if (sys.brownian){
		for (int u=0; u<5; u++){
			total_stress[u] += sys.total_brownian_stress[u];
		}
	}
	Viscosity = total_stress[2];
	Viscosity_h = sys.total_hydro_stress[2];
	Viscosity_c_XF = sys.total_contact_stressXF[2];
	Viscosity_c_GU = sys.total_contact_stressGU[2];
	Viscosity_col_XF = sys.total_colloidal_stressXF[2];
	Viscosity_col_GU = sys.total_colloidal_stressGU[2];
	
	/* N1 = tau_xx-tau_zz = tau_xx-(-tau_xx-tau_yy) = 2tau_xx+tau_yy
	 * N2 = tau_zz-tau_yy = (-tau_xx-tau_yy)-tau_yy = -tau_xx-2tau_yy
	 */
	N1 = 2*total_stress[0]+total_stress[4];
	N2 = -total_stress[0]-2*total_stress[4];
	N1_h = 2*sys.total_hydro_stress[0]+sys.total_hydro_stress[4];
	N2_h = -sys.total_hydro_stress[0]-2*sys.total_hydro_stress[4];
	N1_c_XF = 2*sys.total_contact_stressXF[0]+sys.total_contact_stressXF[4];
	N2_c_XF = -sys.total_contact_stressXF[0]-2*sys.total_contact_stressXF[4];
	N1_c_GU = 2*sys.total_contact_stressGU[0]+sys.total_contact_stressGU[4];
	N2_c_GU = -sys.total_contact_stressGU[0]-2*sys.total_contact_stressGU[4];
	if (sys.brownian){
		Viscosity_b = sys.total_brownian_stress[2];
		N1_b = 2*sys.total_brownian_stress[0]+sys.total_brownian_stress[4];
		N2_b = -sys.total_brownian_stress[0]-2*sys.total_brownian_stress[4];
	} else {
		Viscosity_b = 0;
		N1_b = 0;
		N2_b = 0;
	}
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
	 */
	if ( firsttime ){
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
		fout_rheo << "#14: Viscosity(brownian)" << endl;
		fout_rheo << "#15: N1(brownian)" << endl;
		fout_rheo << "#16: N2(brownian)" << endl;
		fout_rheo << "#17: min gap (non-dim)" << endl;
		fout_rheo << "#18: Average normal contact force" << endl;
		fout_rheo << "#19: max Fc_normal_norm" << endl;
	}
	double unit_of_viscosity = unit_of_force/(unit_of_velocity*unit_of_length);
	double unit_of_stress = unit_of_force/(unit_of_length*unit_of_length);
	fout_rheo << sys.strain() << ' '; //1
	fout_rheo << Viscosity*unit_of_viscosity << ' ' ; //2
	fout_rheo << N1*unit_of_stress << ' ' ; //3
	fout_rheo << N2*unit_of_stress << ' ' ; //4
	fout_rheo << Viscosity_h*unit_of_viscosity << ' ' ; //5
	fout_rheo << N1_h*unit_of_stress << ' ' ; //6
	fout_rheo << N2_h*unit_of_stress << ' ' ; //7
	fout_rheo << Viscosity_c_XF*unit_of_viscosity << ' ' ; //8
	fout_rheo << N1_c_XF*unit_of_stress << ' ' ; //9
	fout_rheo << N2_c_XF*unit_of_stress << ' ' ; //10
	fout_rheo << Viscosity_c_GU*unit_of_viscosity << ' ' ; //11
	fout_rheo << N1_c_GU*unit_of_stress << ' ' ; //12
	fout_rheo << N2_c_GU*unit_of_stress << ' ' ; //13
	fout_rheo << Viscosity_b*unit_of_viscosity << ' ' ; //14
	fout_rheo << N1_b*unit_of_stress << ' ' ; //15
	fout_rheo << N2_b*unit_of_stress << ' ' ; //16
	fout_rheo << Viscosity_col_GU*unit_of_viscosity << ' '; //17
	fout_rheo << Viscosity_col_XF*unit_of_viscosity << ' '; //18
	fout_rheo << sys.minvalue_gap_nondim << ' '; // 19
	fout_rheo << sys.average_Fc_normal_norm << ' '; // 20
	fout_rheo << sys.max_Fc_normal_norm << ' '; // 21
	fout_rheo << endl;
}

vec3d
Simulation::shiftUpCoordinate(double x, double y, double z){
	if (origin_zero_flow){
		z += sys.lz2();
		if (z > sys.lz2()){
			x -= sys.shear_disp;
			if (x < - sys.lx2()){
				x += sys.lx();
			}
			z -=  sys.lz();
		}
	}
	return vec3d(x,y,z);
}

void
Simulation::drawLine(char type , const vec3d &pos, const vec3d &vec, ofstream &fout){
	fout << type << ' ';
	fout << pos.x << ' '<< pos.y << ' '<< pos.z << ' ';
	fout << pos.x + vec.x << ' '<< pos.y + vec.y << ' '<< pos.z + vec.z << endl;
}

void
Simulation::drawLine2(char type , const vec3d &pos1, const vec3d &pos2, ofstream &fout){
	vec3d seg = pos2 - pos1;
	vec3d pos2_ = pos2;
	fout << type << ' ';
	if (seg.z > sys.lz2()){
		pos2_.z -= sys.lz();
		pos2_.x -= sys.shear_disp;
		seg = pos2 - pos1;
	} else if (seg.z < -sys.lz2()){
		pos2_.z += sys.lz();
		pos2_.x += sys.shear_disp;
		seg = pos2 - pos1;
	}
	while (seg.x > sys.lx2()){
		pos2_.x -= sys.lx();
		seg = pos2 - pos1;
	}
	while (seg.x < -sys.lx2()){
		pos2_.x += sys.lx();
		seg = pos2 - pos1;
	}
	if (seg.y > sys.ly2()){
		pos2_.y -= sys.ly();
	} else if (seg.y < -sys.ly2()){
		pos2_.y += sys.ly();
	}
	fout << pos1.x << ' '<< pos1.y << ' '<< pos1.z << ' ';
	fout << pos2_.x << ' '<< pos2_.y << ' '<< pos2_.z << endl;
}

void
Simulation::drawLine(double x0, double y0, double z0,
			  double x1, double y1, double z1,
			  ofstream &fout){
	fout << 'l' << ' ';
	fout << x0 << ' ' << y0 << ' ' << z0 << ' ';
	fout << x1 << ' ' << y1 << ' ' << z1 << endl;
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
	int np = np_a+np_b;
	pos.resize(np);
	vel.resize(np);
	for (int i=0; i < np; i++){
		pos[i] = shiftUpCoordinate(sys.position[i].x-sys.lx2(),
								   sys.position[i].y-sys.ly2(),
								   sys.position[i].z-sys.lz2());
	}
	/* If the origin is shifted,
	 * we need to change the velocities of particles as well.
	 */
	if (origin_zero_flow){
		for (int i=0; i < np; i++){
			vel[i] = sys.velocity[i];
			if (pos[i].z < 0){
				vel[i].x -= sys.lz();
			}
		}
	}
	/*
	 * shear_disp = sys.strain() - (int)(sys.strain()/Lx)*Lx
	 */
	fout_particle << "#" << sp << sys.strain() << sp;
	fout_particle << sys.shear_disp << sp << sys.shear_rate << endl;
	for (int i=0; i < np; i++){
		vec3d &p = pos[i];
		vec3d &v = vel[i];
		vec3d &o = sys.ang_velocity[i];
		double h_xzstress = sys.lubstress[i].elm[2] + sys.bgfstress[i].elm[2];
		double c_xzstressXF = sys.contactstressXF[i].elm[2];
		double c_xzstressGU = sys.contactstressGU[i].elm[2];
		double b_xzstress = sys.brownianstress[i].elm[2];
		fout_particle << i << sp; //1: number
		fout_particle << sys.radius[i] << sp; //2: radius
		fout_particle << p.x << sp << p.y << sp << p.z << sp; //3,4,5: position
		fout_particle << v.x << sp << v.y << sp << v.z << sp; //6,7,8: velocity
		fout_particle << o.x << sp << o.y << sp << o.z << sp; //9,10,11: angular velocity
		fout_particle << h_xzstress << sp; //12: xz stress contributions
		fout_particle << c_xzstressXF << sp; //13: xz stress contributions
		fout_particle << c_xzstressGU << sp; //14: xz stress contributions
		fout_particle << b_xzstress << sp; //15: xz stress contributions
		if (sys.dimension == 2){
			fout_particle << sys.angle[i] << sp;
		}
		fout_particle << endl;
	}
	int cnt_interaction = 0;
	for (int k=0; k<sys.num_interaction; k++){
		if (sys.interaction[k].active){
			cnt_interaction++;
		}
	}
	fout_interaction << "#" << sp << sys.strain() << sp ;
	fout_interaction << cnt_interaction << " avg_contact_time: " << sys.average_contact_time;
	fout_interaction << " avg_nearing_time: " << sys.average_nearing_time << endl;
	for (int k=0; k<sys.num_interaction; k++){
		if (sys.interaction[k].active){
			fout_interaction << sys.interaction[k].par_num[0] << sp; // 1
			fout_interaction << sys.interaction[k].par_num[1] << sp; // 2
			fout_interaction << sys.interaction[k].valLubForce() << sp; // 3
			fout_interaction << sys.interaction[k].normal_force() << sp; // 5
			fout_interaction << sys.interaction[k].tangential_force().x << sp; // 6
			fout_interaction << sys.interaction[k].tangential_force().y << sp; // 7
			fout_interaction << sys.interaction[k].tangential_force().z << sp; // 8
			fout_interaction << sys.interaction[k].colloidal_force() << sp; //
			fout_interaction << sys.interaction[k].nr_vec.x << sp;
			fout_interaction << sys.interaction[k].nr_vec.y << sp;
			fout_interaction << sys.interaction[k].nr_vec.z << sp;
			fout_interaction << sys.interaction[k].gap_nondim() << sp; //7
			fout_interaction << sys.interaction[k].lubStresslet(2) << sp; //
			fout_interaction << sys.interaction[k].contact << sp; //8
			fout_interaction << endl;
		}
	}
}
