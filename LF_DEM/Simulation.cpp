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
	if (fout_yap.is_open()){
		fout_yap.close();
	}
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
Simulation::AutoSetParameters(const string &keyword,
							  const string &value){
	if (keyword == "lubrication"){
		sys.lubrication = str2bool(value);
	} else if (keyword == "friction"){
		sys.friction = str2bool(value);
	} else if (keyword == "bgf_factor"){
		sys.bgf_factor = atof(value.c_str());
	} else if (keyword == "poly"){
		sys.poly = str2bool(value);
	} else if (keyword == "gap_cutoff"){
		sys.gap_cutoff = atof(value.c_str());
	} else if (keyword == "shear_rate"){
		sys.shear_rate = atof(value.c_str());
	} else if (keyword == "kb_T"){
		sys.kb_T = atof(value.c_str());
	} else if (keyword == "dt"){
		sys.dt = atof(value.c_str());
	} else if (keyword == "shear_strain_end"){
		shear_strain_end = atof(value.c_str());
	} else if (keyword == "dt_ratio"){
		sys.dt_ratio = atof(value.c_str());
	} else if (keyword == "lub_max"){
		sys.lub_max = atof(value.c_str());
	} else if (keyword == "kn"){
		sys.kn = atof(value.c_str());
	} else if (keyword == "kt"){
		sys.kt = atof(value.c_str());
	} else if (keyword == "mu_static"){
		sys.mu_static = atof(value.c_str());
	} else if (keyword == "mu_dynamic"){
		sys.mu_dynamic = atof(value.c_str());
	} else if (keyword == "dynamic_friction_critical_velocity"){
		sys.dynamic_friction_critical_velocity = atof(value.c_str());
	} else if (keyword == "strain_interval_out") {
		strain_interval_out = atof(value.c_str());
	} else if (keyword == "yap_force_factor"){
		yap_force_factor = atof(value.c_str());
	} else if (keyword == "draw_rotation_2d"){
		sys.draw_rotation_2d = str2bool(value);
	} else if (keyword == "dist_near"){
		sys.dist_near = atof(value.c_str());
	} else if (keyword == "out_data_particle"){
		out_data_particle = str2bool(value);
	} else if (keyword == "out_data_interaction"){
		out_data_interaction = str2bool(value);
	} else if (keyword == "out_vpython"){
		out_vpython = str2bool(value);
	} else if (keyword == "out_yaplot"){
		out_yaplot = str2bool(value);
	} else if (keyword == "out_pairtrajectory"){
		sys.out_pairtrajectory = str2bool(value);
	} else if (keyword == "origin_zero_flow"){
		origin_zero_flow = str2bool(value);
	} else {
		cerr << "keyword " << keyword << " is'nt associated with an parameter" << endl;
		exit(1);
	}
}

void
Simulation::ReadParameterFile(){
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
		AutoSetParameters(keyword, value);
	}
	fin.close();
	return;
}

void
Simulation::SetParametersPostProcess(){
	if (sys.kb_T == 0){
		sys.brownian = false;
	} else {
		sys.brownian = true;
	}
	
	/*
	 * Set simulation name and name of output files.
	 */
	prepareSimulationName();
	

	string particle_filename = "par_" + sys.simu_name + ".dat";
	string interaction_filename = "int_" + sys.simu_name + ".dat";
	string yap_filename = "yap_" + sys.simu_name + ".yap";
	string vel_filename = "rheo_" + sys.simu_name + ".dat";
	string vpy_filename = "vpy_" + sys.simu_name + ".dat";
	fout_particle.open(particle_filename.c_str());
	fout_interaction.open(interaction_filename.c_str());
	
	if (out_yaplot)
	fout_yap.open(yap_filename.c_str());
	fout_rheo.open(vel_filename.c_str());
	fout_vpy.open(vpy_filename.c_str());
	if (sys.out_pairtrajectory){
		string trj_filename = "trj_" + sys.simu_name + ".dat";
		sys.fout_trajectory.open(trj_filename.c_str());
	}

	sys.sq_critical_velocity = sys.dynamic_friction_critical_velocity * sys.dynamic_friction_critical_velocity;	
	sys.sq_lub_max = sys.lub_max*sys.lub_max; // square of lubrication cutoff length.
	sys.ts = 0;
	sys.shear_disp = 0;
	sys.vel_difference = sys.shear_rate*sys.lz();
	/*
	 * dt_mid: the intermediate time step for the mid-point
	 * algortithm. dt/dt_mid = dt_ratio
	 * Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 * ASD code from Brady has dt_ratio=150
	 */
	sys.dt_mid = sys.dt/sys.dt_ratio;
	/*
	 * The time steps finishing simulation.
	 */
	//ts_max = (int)(shear_strain / sys.dt);

}

void
Simulation::SetDefaultParameters(){
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
	sys.dt_ratio = 100;
	/*
	 * Shear flow
	 *  shear_rate: shear rate
	 *  shear_strain: total strain (length of simulation)
	 *
	 */
	sys.shear_rate = 1.0;
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
	 *gap_cutoff: reduced small cutoff distance for gap diverging coeffient:
	 *
	 * 1/ksi when ksi > gap_cutoff*(a0+a1)/2. = ksi_cutoff
	 * 1/ksi_cutoff when h <= ksi_cutoff
	 */
	sys.gap_cutoff = 0.001;
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
	 *
	 *
	 */
	sys.kn = 300;
	sys.kt = 300;
	/*
	 * mu_static: static friction coeffient
	 * mu_dynamic: dynamic friction coeffient
	 */
	sys.mu_static = 10;
	sys.mu_dynamic = 8;
	/*
	 * dynamic_friction_critical_velocity:
	 * This is a threshold velocity to swich from dynamic friction to
	 * static friction. But this is a temporal provísional.
	 * There is no reference to give this value.
	 */
	sys.dynamic_friction_critical_velocity = 0.01;
	/*
	 * Output interval
	 */
	strain_interval_out = 0.01;
	/*
	 * Consider particles are `near' if the gap is less than the following value
	 */
	sys.dist_near = 1e-2;
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
	 * output data for vpython
	 */
	out_vpython = true;
	/*
	 * output data for python
	 */
	out_yaplot = false;
	/*
	 * The bond width indicates the force strength.
	 */
	yap_force_factor = 0.02;
	/*
	 * Visualize rotations by crosses.
	 * (for yaplot data)
	 */
	if (sys.dimension == 2)
		sys.draw_rotation_2d = true;
	else
		sys.draw_rotation_2d = false;
	/*
	 * output pair trajeectory
	 */
	sys.out_pairtrajectory = false;
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
	int _np1_, _np2_;
	double _volume_fraction_;
	double _lx_, _ly_, _lz_;
	char buf;
	file_import >> buf >> _np1_ >> _np2_ >> _volume_fraction_ >> _lx_ >> _ly_ >> _lz_ ;
	np1 = _np1_;
	np2 = _np2_;

	int num_of_particle = _np1_ + _np2_;
	sys.set_np(num_of_particle);
	
	if (np2 == 0)
		sys.poly = false;
	else
		sys.poly = true;
	cerr << "np = " << num_of_particle << endl;
	if (_ly_ == 0){
		sys.dimension = 2;
	} else {
		sys.dimension = 3;
	}
	cerr << "dimension = " << sys.dimension << endl;
	sys.lx( _lx_);
	sys.ly( _ly_);
	sys.lz( _lz_);
	cerr << "box: " << _lx_ << ' ' <<  _ly_ << ' ' << _lz_ << endl;
	sys.volume_fraction = _volume_fraction_;
	initial_positions.resize(num_of_particle);
	radii.resize(num_of_particle);
	for (int i = 0; i < num_of_particle ; i++){
		file_import >> pos.x >> pos.y >> pos.z >> radius;
		initial_positions[i] = pos;
		radii[i] = radius;
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
	sys.simu_name = ss_simu_name.str();
	cerr << sys.simu_name << endl;
}


/*
 * Main simulation
 */
void
Simulation::SimulationMain(int argc, const char * argv[]){
	filename_import_positions = argv[1];
	importInitialPositionFile();
	SetDefaultParameters();
	if (argc == 3){
		cerr << "Read Parameter File" << endl;
		filename_parameters = argv[2];
		ReadParameterFile();
	}
	
	SetParametersPostProcess();
	sys.allocateRessources();
	for (int i=0; i < sys.np; i++){
		sys.position[i] = initial_positions[i];
		sys.radius[i] = radii[i];
		sys.angle[i] = 0;
	}
	sys.initializeBoxing();
	sys.checkNewInteraction();
	sys.shear_strain = 0;
	double time_interval = strain_interval_out/sys.shear_rate;
	int i_time_interval = time_interval / sys.dt;
	while(sys.shear_strain <= shear_strain_end){
		sys.timeEvolution(i_time_interval);
		outputRheologyData();
		outputData();
		if(out_yaplot)
			output_yap();
		if(out_vpython)
			output_vpython(sys.ts);

	}
}

/*
 *
 */
void
Simulation::outputRheologyData(){
	/*
	 * Output the sum of the normal forces.
	 */
	fout_rheo << sys.dt * sys.ts << ' ';// 1
	fout_rheo << sys.mean_hydro_stress[2] + sys.mean_contact_stress[2]  << ' ' ; //2
	fout_rheo << sys.mean_hydro_stress[0] << ' ' ; //3
	fout_rheo << sys.mean_hydro_stress[1] << ' ' ; //4
	fout_rheo << sys.mean_hydro_stress[2] << ' ' ; //5
	fout_rheo << sys.mean_hydro_stress[3] << ' ' ; //6
	fout_rheo << sys.mean_hydro_stress[4] << ' ' ; //7
	fout_rheo << sys.mean_contact_stress[0] << ' ' ; //8
	fout_rheo << sys.mean_contact_stress[1] << ' ' ; //9
	fout_rheo << sys.mean_contact_stress[2] << ' ' ; //10
	fout_rheo << sys.mean_contact_stress[3] << ' ' ; //11
	fout_rheo << sys.mean_contact_stress[4] << ' ' ; //12
	fout_rheo << sys.gap_min << ' '; // 13
	fout_rheo << sys.ave_overlap << ' '; //14
	fout_rheo << sys.shear_strain << ' '; //15
	fout_rheo << sys.max_age << ' '; //16
	fout_rheo << sys.ave_age << ' '; //17
	if (sys.n_vec_longcontact.z > 0){
		fout_rheo << sys.n_vec_longcontact.x << ' '; //18
		fout_rheo << sys.n_vec_longcontact.y << ' '; //19
		fout_rheo << sys.n_vec_longcontact.z << ' '; //20
	} else {
		fout_rheo << -sys.n_vec_longcontact.x << ' ';
		fout_rheo << -sys.n_vec_longcontact.y << ' ';
		fout_rheo << sys.n_vec_longcontact.z << ' ';
	}
	fout_rheo << sys.total_contact << ' '; //21
	fout_rheo << sys.cnt_contact_number[0] << ' ';
	fout_rheo << sys.cnt_contact_number[1] << ' ';
	fout_rheo << sys.cnt_contact_number[2] << ' ';
	fout_rheo << sys.cnt_contact_number[3] << ' ';
	fout_rheo << sys.cnt_contact_number[4] << ' ';
	fout_rheo << sys.cnt_contact_number[5] << ' ';
	fout_rheo << sys.cnt_contact_number[6] << ' ';
	fout_rheo << sys.cnt_contact_number[7] << ' ';
	fout_rheo << sys.cnt_contact_number[8] << ' ';
	fout_rheo << sys.cnt_contact_number[9] << ' ';

	fout_rheo << endl;
}

vec3d
Simulation::shiftUpCoordinate(double x, double y, double z){
	if (origin_zero_flow){
		z += sys.lz2();
		if (z > sys.lz2()){
			x += - sys.shear_disp;
			if ( x < - sys.lx2())
				x += sys.lx();
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
Simulation::outputData(){
	
	vector<vec3d> pos;
	char sp = ' ';
	int np = np1 + np2;
	pos.resize(np);
	for (int i=0; i < np; i++){
		pos[i] = shiftUpCoordinate(sys.position[i].x - sys.lx2(),
								   sys.position[i].y - sys.ly2(),
								   sys.position[i].z - sys.lz2());
	}
	
	fout_particle << "#" << sp << sys.shear_strain  << sp << np1 << sp << np2 << sp;
	fout_particle << sys.radius[0] << sp << sys.radius[np1] << endl;
	for (int i=0; i < np; i++){
		vec3d &p = pos[i];
		vec3d &v = sys.velocity[i];
		vec3d &o = sys.ang_velocity[i];
		fout_particle << i << sp; //1: number
		fout_particle << p.x << sp << p.y << sp << p.z << sp; //2,3,4: position
		fout_particle << v.x << sp << v.y << sp << v.z << sp; //5,6,7: velocity
		fout_particle << o.x << sp << o.y << sp << o.z << sp; //5,6,7: velocity
		fout_particle << endl;
	}
	fout_interaction << "#" << sp << sys.shear_strain  << sp << np1 << sp << np2 << sp;
	fout_interaction << sys.radius[0] << sp << sys.radius[np1] << endl;

	for (int k=0; k < sys.num_interaction; k++){
		if (sys.interaction[k].active){
			fout_interaction << sys.interaction[k].particle_num[0] << sp; //1
			fout_interaction << sys.interaction[k].particle_num[1] << sp; //2
			fout_interaction << sys.interaction[k].r() << sp; //3
			fout_interaction << sys.interaction[k].valLubForce() << sp; //4
			fout_interaction << sys.interaction[k].Fc_normal << sp; // 5
			fout_interaction << sys.interaction[k].Fc_tangent.norm() << sp; // 6
			fout_interaction << sys.interaction[k].static_friction << sp; //7
			/// fout_interaction << ??? 
			fout_interaction << endl;
		}
	}

	
}


/* Output data for yaplot visualization.
 *
 */
void
Simulation::output_yap(){
	static bool fasttime = true;
	if (fasttime){
		fasttime = false;
	}else{
		fout_yap << endl;
	}
	
	double y_trimming = 2;
	/*
	 * yaplot color
	 *
	 * int color_black = 0;
	 * int color_gray = 1;
	 */
	int color_white = 2;
	int color_green = 3;
	int color_yellow = 4;
	int color_orange = 5;
	int color_blue = 6;
	
	vector<vec3d> pos;
	int np = np1 + np2 ;
	pos.resize(np);
	for (int i=0; i < np; i++){
		pos[i] = shiftUpCoordinate(sys.position[i].x - sys.lx2(),
								   sys.position[i].y - sys.ly2(),
								   sys.position[i].z - sys.lz2());
	}

	/* Layer 1: Circles for particles
	 */
	fout_yap << "y 1\n";
	fout_yap << "@ " << color_white << endl;
	//vec3d pos;
	fout_yap << "r " << sys.radius[0] << endl;
	
	for (int i=0; i < np1; i++){
		if (abs(pos[i].y) < y_trimming ){
			fout_yap << "c " << pos[i].x << ' ' << pos[i].y << ' ' << pos[i].z << endl;
		}
	}
	fout_yap << "r " << sys.radius[np1+1] << endl;
	for (int i = np1; i < sys.np ; i++){
		if (abs(pos[i].y) < y_trimming ){
			fout_yap << "c " << pos[i].x << ' ' << pos[i].y << ' ' << pos[i].z << endl;
		}
	}
	
	/* Layer 4: Orientation of particle (2D simulation)
	 * Only for small system.
	 */
	if (sys.np <= 1000 && sys.draw_rotation_2d){
		fout_yap << "y 5\n";
		fout_yap << "@ " << color_white << endl;
		for (int i=0; i < sys.np; i++){

				vec3d u(cos(-sys.angle[i]),0,sin(-sys.angle[i]));
				u = sys.radius[i]*u;
				drawLine('l', pos[i]-u, 2*u, fout_yap);
				u.set(-sin(-sys.angle[i]), 0, cos(-sys.angle[i]));
				u = sys.radius[i]*u;
				drawLine('l', pos[i]-u, 2*u, fout_yap);
			
		}
	}
	/* Layer 2: Friction
	 */
	if (sys.friction){
		fout_yap << "y 2\n";
		for (int k=0; k < sys.num_interaction; k++){
			if ( sys.interaction[k].contact){
				if (sys.interaction[k].static_friction)
					fout_yap << "@ " << color_green << endl;
				else
					fout_yap << "@ " << color_orange << endl;
				int i = sys.interaction[k].particle_num[0];
				fout_yap << "r " << yap_force_factor*sys.interaction[k].Fc_tangent.norm()  << endl;
				drawLine('s', pos[i], sys.radius[i]*sys.interaction[k].nr_vec, fout_yap);
				int j = sys.interaction[k].particle_num[1];
				drawLine('s', pos[j], -sys.radius[j]*sys.interaction[k].nr_vec, fout_yap);
			}
		}
	}
	/* Layer 3: lubrication (Compression)
	 */
	// nr_vec i ---> j
	fout_yap << "y 3\n";
	fout_yap << "@ " << color_yellow << endl;
	for (int k=0; k < sys.num_interaction; k++){
		if ( sys.interaction[k].active ){
			double lub_force = sys.interaction[k].valLubForce();
			if (lub_force >= 0){
				vec3d nvec = sys.interaction[k].nr_vec;
				int i = sys.interaction[k].particle_num[0];
				int j = sys.interaction[k].particle_num[1];
				if (abs(pos[i].y) < y_trimming || abs(pos[j].y) < y_trimming ){
					fout_yap << "r " << yap_force_factor*lub_force << endl;
					drawLine('s', pos[i],  sys.radius[i]*nvec, fout_yap);
					drawLine('s', pos[j], -sys.radius[j]*nvec, fout_yap);
				}
			}
		}
	}
	/* Layer 4: lubrication (Extension)
	 */
	fout_yap << "y 4\n";
	fout_yap << "@ " << color_green << endl;
	for (int k=0; k < sys.num_interaction; k++){
		if ( sys.interaction[k].active ){
			double lub_force = sys.interaction[k].valLubForce();
			if (lub_force < 0){
				vec3d nvec = sys.interaction[k].nr_vec;
				int i = sys.interaction[k].particle_num[0];
				int j = sys.interaction[k].particle_num[1];
				if (abs(pos[i].y) < y_trimming || abs(pos[j].y) < y_trimming ){
					fout_yap << "r " << -yap_force_factor*lub_force << endl;
					drawLine('s', pos[i],  sys.radius[i]*nvec, fout_yap);
					drawLine('s', pos[j], -sys.radius[j]*nvec, fout_yap);
				}
			}
		}
	}

	/* Layer 6: Box and guide lines
	 */
	fout_yap << "y 6\n";
	fout_yap << "@ " << color_blue << endl;
	if (sys.dimension == 2){
		drawLine(-sys.lx2(), 0, -sys.lz2(),  sys.lx2(), 0, -sys.lz2(), fout_yap);
		drawLine(-sys.lx2(), 0,  sys.lz2(), -sys.lx2(), 0, -sys.lz2(), fout_yap);
		drawLine(-sys.lx2(), 0,  sys.lz2(),  sys.lx2(), 0,  sys.lz2(), fout_yap);
		drawLine( sys.lx2(), 0, -sys.lz2(),  sys.lx2(), 0,  sys.lz2(), fout_yap);
		/*
		for (int i = 0 ; i < 10; i ++){
			drawLine(-sys.lx2()                , 0, -sys.lz2() + i*0.2*sys.lz2(),
					 -sys.lx2() + i*0.2*sys.lz2(), 0,                 -sys.lz2(), fout_yap);
			drawLine( sys.lx2()                   ,        0, sys.lz2() - i*0.2*sys.lz2(),
					  sys.lx2() - i*0.2*sys.lz2(), 0,  sys.lz2(), fout_yap);
		}
		drawLine(-sys.lx2(), 0, sys.lz2(), sys.lx2() , 0, -sys.lz2(), fout_yap);
		 */
	} else {
		drawLine(-sys.lx2(), -sys.ly2(), -sys.lz2(),  sys.lx2(), -sys.ly2(), -sys.lz2(), fout_yap);
		drawLine(-sys.lx2(),  sys.ly2(), -sys.lz2(),  sys.lx2(),  sys.ly2(), -sys.lz2(), fout_yap);

		drawLine(-sys.lx2(),  sys.ly2(),  sys.lz2(),  sys.lx2(),  sys.ly2(),  sys.lz2(), fout_yap);
		drawLine(-sys.lx2(), -sys.ly2(),  sys.lz2(),  sys.lx2(), -sys.ly2(),  sys.lz2(), fout_yap);
		
		drawLine(-sys.lx2(), -sys.ly2(), -sys.lz2(), -sys.lx2(), sys.ly2(), -sys.lz2(), fout_yap);
		drawLine(-sys.lx2(), -sys.ly2(),  sys.lz2(), -sys.lx2(), sys.ly2(),  sys.lz2(), fout_yap);
		drawLine( sys.lx2(), -sys.ly2(),  sys.lz2(),  sys.lx2(), sys.ly2(),  sys.lz2(), fout_yap);
		drawLine( sys.lx2(), -sys.ly2(), -sys.lz2(),  sys.lx2(), sys.ly2(), -sys.lz2(), fout_yap);
		
		drawLine( sys.lx2(),  sys.ly2(),  sys.lz2(),  sys.lx2(), sys.ly2(), -sys.lz2(),  fout_yap);
		drawLine(-sys.lx2(),  sys.ly2(),  sys.lz2(), -sys.lx2(), sys.ly2(), -sys.lz2(),  fout_yap);
		drawLine(-sys.lx2(), -sys.ly2(),  sys.lz2(), -sys.lx2(), -sys.ly2(), -sys.lz2(),  fout_yap);
		drawLine( sys.lx2(), -sys.ly2(),  sys.lz2(),  sys.lx2(), -sys.ly2(), -sys.lz2(),  fout_yap);
	}
}

/* Output data for vpython visualization.
 *
 */
void
Simulation::output_vpython(double time){
	vec3d pos;
	fout_vpy << "time: " << time << endl;
	for (int i=0; i < sys.np; i++){
		pos = shiftUpCoordinate(sys.position[i].x - sys.lx2(),
								sys.position[i].y - sys.ly2(),
								sys.position[i].z - sys.lz2());
		fout_vpy << i << ' ' << pos.x << ' ' << pos.y << ' ' << pos.z << ' ' << sys.radius[i] << endl;
	}
	fout_vpy << endl;
}
