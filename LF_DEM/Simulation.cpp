//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "Simulation.h"
#include <cmath>

Simulation::Simulation(){};

Simulation::~Simulation(){
	fout_yap.close();
	fout_force.close();
	delete [] fc;

	for (int i=0; i < num_particle; i++){
		delete [] contact_pair[i];
	}
	delete [] contact_pair;
};

void Simulation::SetParameters(int argc, const char * argv[]){
	sys.friction = true;
	sys.lub = true;
	filename_import_positions = argv[1];
	sys.lubcore = atof(argv[2]);
	/* take parameters from import file name.
	 *
	 */ 
	int i_D = (int)filename_import_positions.find( "D") + 1;
	sys.dimension = atoi( filename_import_positions.substr(i_D, 1).c_str() );
	if (sys.dimension == 2 ){
		// example: D2L10_10vf0.8.dat
		int i_lx = (int)filename_import_positions.find( "L") + 1;
		int j_lx = (int)filename_import_positions.find( "_" );
		int j_lz = (int)filename_import_positions.find( "vf", j_lx);
		int j_vf = (int)filename_import_positions.find( ".dat", j_lz);
		cerr << i_lx << ' ' << j_lx << ' ' << j_lz << endl;
		sys.lx = atoi( filename_import_positions.substr(i_lx, j_lx - i_lx).c_str() );
		sys.ly = 0;
		sys.lz = atoi( filename_import_positions.substr(j_lx+1, j_lz - j_lx-1).c_str() );
		sys.volume_fraction = atof( filename_import_positions.substr(j_lz + 2, j_vf - j_lz-2).c_str() );
	} else {
		// example: D3L10_10_10vf0.5.dat
		int i_lx = (int)filename_import_positions.find( "L") + 1;
		int j_lx = (int)filename_import_positions.find( "_", i_lx);
		int j_ly = (int)filename_import_positions.find( "_", j_lx+1);
		int j_lz = (int)filename_import_positions.find( "vf", j_ly+1);
		int j_vf = (int)filename_import_positions.find( ".dat", j_lz);
		sys.lx = atoi( filename_import_positions.substr(i_lx  , j_lx - i_lx).c_str() );
		sys.ly = atoi( filename_import_positions.substr(j_lx+1, j_ly - j_lx-1).c_str() );
		sys.lz = atoi( filename_import_positions.substr(j_ly+1, j_lz - j_ly-1).c_str() );
		sys.volume_fraction = atof( filename_import_positions.substr(j_lz + 2, j_vf-j_lz-2).c_str() );
	}
	cerr << sys.lx << ' ' << sys.ly << ' ' << sys.lz << ' ' << sys.volume_fraction << endl;
	/*
	 * Simulation parameters
	 */
	sys.eta = 1.0; // viscosity
	sys.shear_rate = 1; // shear rate
	shear_strain = 0.02; // shear strain
	sys.kb_T=1.;
	cutoff_distance = 2.5; // to delete possible neighbor
	sys.sq_lub_max = 2.5*2.5; // square of lubrication cutoff length.
	sys.dt = 1e-4 / sys.shear_rate; //time step.
	ts_max = (int)(shear_strain / sys.dt); // time step max

	/*
	 * Contact force parameters
	 */
	sys.kn = 200; // normal spring constant
	sys.kt = 200; // tangential spring constant
//	sys.mu_static = 0.3; // static friction coeffient
//	sys.mu_dynamic = 0.2; // dynamic friction coeffient
	sys.mu_static = 0.6; // static friction coeffient
	sys.mu_dynamic = 0.4; // dynamic friction coeffient
	sys.dynamic_friction_critical_velocity = 0.01;
	/*
	 * Visualization
	 */
	if (sys.dimension ==2)
		draw_rotation_2d = true;
	else
		draw_rotation_2d = false;
	interval_snapshot = 100;
	yap_force_factor = 0.01;
	origin_zero_flow = false;
	/********************************************************************************************/
	if (sys.dimension == 2){
		num_particle = (int)(sys.lx*sys.lz*sys.volume_fraction/M_PI);
	} else {
		num_particle = (int)(sys.lx*sys.ly*sys.lz*sys.volume_fraction/(4.0*M_PI/3.0));
	}
	cerr << "N = " << num_particle << endl;
	max_num_interaction = 20 * num_particle;
	sys.init();
	string yap_filename = "yap_" + sys.simu_name + ".yap";
	string vel_filename = "force_" + sys.simu_name + ".dat";
	fout_yap.open(yap_filename.c_str());
	fout_force.open(vel_filename.c_str());
}

void Simulation::importInitialPositionFile(){
	fstream file_import;
	file_import.open( filename_import_positions.c_str());
	vec3d pos;
	while ( !file_import.eof() ){
		file_import >> pos.x >> pos.y >> pos.z;
		cerr << pos.x << ' '<< pos.y << ' '<< pos.z << endl;
		initial_positions.push_back(pos);
	}
	file_import.close();
}



/*
 * Main simulation
 */
void Simulation::SimulationMain(int argc, const char * argv[]){
	SetParameters(argc, argv);
	importInitialPositionFile();
	unsigned long num_of_particles = initial_positions.size();
	sys.prepareSimulation(num_of_particles);
	for (int i=0; i < initial_positions.size(); i++){
		sys.position[i] = initial_positions[i];
	}
	initContactPair();
	fc = new ContactForce [max_num_interaction];

	
	for (int i = 0; i < max_num_interaction; i++){
		fc[i].init( &sys );
	}
	cerr << "set initial positions" << endl;
	//sys.setRandomPosition();
	importInitialPositionFile();
	cerr << "start simulation" << endl;
	timeEvolution();
	cerr << "finished" << endl;
}

/*
 * only elements of j > i will be used.
 */
void Simulation::initContactPair(){
	contact_pair = new int * [num_particle];
	for (int i=0; i < num_particle; i++){
		contact_pair[i] = new int [num_particle];
	}
	for (int i=0; i < num_particle-1; i++){
		for (int j=i+1; j < num_particle; j++){
			contact_pair[i][j] = -1;
		}
	}
}

/* Interaction is gone when two particles are separated enough.
 *
 */
void Simulation::checkBreak(){
	for (int k = 0; k < num_interaction; k++){
		if ( fc[k].active
			&& fc[k].r > cutoff_distance){
			fc[k].active = false;
			contact_pair[fc[k].particle_num[0]][fc[k].particle_num[1]] = -1;
			deactivated_interaction.push(k);
		}
	}
}

/* Check the distance between separating particles.
 * i < j 
 *
 * A patch-up prescription to aboid 
 * contact_pair[i][j] < 0 indicates separating particles to be checked.
 * contact_pair[i][j] = -1, the particles are near contact. So every time step, distance should be checked.a
 * contact_pair[i][j] < -1, the particles have some distance.
 */
void Simulation::checkContact(){
	for (int i=0; i < num_particle-1; i++){
		for (int j=i+1; j < num_particle; j++){
			if ( contact_pair[i][j] == -1){
				//double sq_distance = sys.sq_distance(i, j);
				double sq_distance = sys.checkContact(i, j);
				if ( sq_distance < 4){
					int new_num_interaction;
					if (deactivated_interaction.empty()){
						// add an interaction object.
						new_num_interaction = num_interaction;
						num_interaction ++;
					} else {
						// fill a deactivated interaction object.
						new_num_interaction = deactivated_interaction.front();
						deactivated_interaction.pop();
					}
					fc[new_num_interaction].create(i,j);
					contact_pair[i][j] = new_num_interaction;
				}
			}
		}
	}
}

/* Simulation for the time evolution.
 *
 */
void Simulation::timeEvolution(){
	sys.x_shift = 0;
	sys.ts = 0;
	while (sys.ts < ts_max){
		checkContact();
		sys.forceReset();
		sys.torqueReset();
		if (sys.friction){
			for (int k=0; k < num_interaction; k++){
				fc[k].calcInteraction();
			}
		} else {
			for (int k=0; k < num_interaction; k++){
				fc[k].calcInteractionNoFriction();
			}
		}
		if (sys.lub){
			// Lubrication dynamics
			sys.updateVelocityLubrication();
		} else {
			// Free-draining approximation
			sys.updateVelocity();
		}
		if (sys.ts % interval_snapshot == 0){
			output_yap();
		}
		sys.deltaTimeEvolution();
		if (sys.friction){
			for (int k = 0; k < num_interaction; k++){
				fc[k].incrementTangentialDisplacement();
			}
		}
		checkBreak();

		sys.ts ++;
	}
}

vec3d Simulation::shiftUpCoordinate(double x, double y, double z){
	if (origin_zero_flow){
		z += sys.lz2;
		if (z > sys.lz2){
			x += - sys.x_shift;
		if ( x < - sys.lx2)
			x += sys.lx;
			z -=  sys.lz;
		}
	}
	return vec3d(x,y,z);
}


void Simulation::drawLine(char type , vec3d pos, vec3d vec, ofstream &fout){
	fout << type << ' ';
	fout << pos.x << ' '<< pos.y << ' '<< pos.z << ' ';
	fout << pos.x + vec.x << ' '<< pos.y + vec.y << ' '<< pos.z + vec.z << endl;
}

void Simulation::drawLine2(char type , vec3d pos1, vec3d pos2, ofstream &fout){
	vec3d seg = pos2 - pos1;
	fout << type << ' ';
	if (seg.z > sys.lz2){
		pos2.z -= sys.lz;
		pos2.x -= sys.x_shift;
		seg = pos2 - pos1;
	} else if (seg.z < -sys.lz2){
		pos2.z += sys.lz;
		pos2.x += sys.x_shift;
		seg = pos2 - pos1;
	}
		
	while (seg.x > sys.lx2){
		pos2.x -= sys.lx;
		seg = pos2 - pos1;
	}
	while (seg.x < -sys.lx2){
		pos2.x += sys.lx;
		seg = pos2 - pos1;
	}
	
	if (seg.y > sys.ly2){
		pos2.y -= sys.ly;
	} else if (seg.y < -sys.ly2){
		pos2.y += sys.ly;
	}
	
	fout << pos1.x << ' '<< pos1.y << ' '<< pos1.z << ' ';
	fout << pos2.x << ' '<< pos2.y << ' '<< pos2.z << endl;
}

void Simulation::drawLine(double x0, double y0, double z0,
			  double x1, double y1, double z1,
			  ofstream &fout){
	fout << 'l' << ' ';
	fout << x0 << ' ' << y0 << ' ' << z0 << ' ';
	fout << x1 << ' ' << y1 << ' ' << z1 << endl;
}

/* Output data for yaplot visualization.
 *
 */
void Simulation::output_yap(){
	static bool fasttime = true;
	if (fasttime){
		fasttime = false;
	}else{
		fout_yap << endl;
	}
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
	/* Layer 1: Circles for particles
	 */
	fout_yap << "y 1\n";
	fout_yap << "@ " << color_white << endl;
	vec3d pos;
	for (int i=0; i < num_particle; i++){
		pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
								sys.position[i].y - sys.ly2,
								sys.position[i].z - sys.lz2);
		fout_yap << "c " << pos.x << ' ' << pos.y << ' ' << pos.z << endl;
	}

	/* Layer 4: Orientation of particle (2D simulation)
	 */
	if (draw_rotation_2d){
		fout_yap << "y 4\n";
		fout_yap << "@ " << color_white << endl;
		for (int i=0; i < num_particle; i++){
			vec3d u(cos(-sys.angle[i]),0,sin(-sys.angle[i]));
			pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
									sys.position[i].y - sys.ly2,
									sys.position[i].z - sys.lz2);
			drawLine('l', pos-u, 2*u, fout_yap);
			u.set(-sin(-sys.angle[i]), 0, cos(-sys.angle[i]));
			drawLine('l', pos-u, 2*u, fout_yap);
		}
	}
	/* Layer 2: Friction
	 */
	if (sys.friction){
		fout_yap << "y 2\n";
		for (int k=0; k < num_interaction; k++){
			if ( fc[k].active && fc[k].r < 2 ){
				if (fc[k].static_friction)
					fout_yap << "@ " << color_green << endl;
				else
					fout_yap << "@ " << color_orange << endl;
				
				int i = fc[k].particle_num[0];
				pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
										sys.position[i].y - sys.ly2,
										sys.position[i].z - sys.lz2);
				fout_yap << "r " << yap_force_factor*fc[k].f_tangent.norm()  << endl;
				drawLine('s', pos, fc[k].nr_vec, fout_yap);
				int j = fc[k].particle_num[1];
				pos = shiftUpCoordinate(sys.position[j].x - sys.lx2,
										sys.position[j].y - sys.ly2,
										sys.position[j].z - sys.lz2);
				drawLine('s', pos, -fc[k].nr_vec, fout_yap);
			}
		}
	}
	/* Layer 3: Normal
	 * Lubrication + contact force
	 */
	fout_yap << "y 3\n";
	fout_yap << "@ " << color_yellow << endl;
	double total_normal_force = 0;
	for (int i=0; i < num_particle; i++){
		for (int j=i+1; j < num_particle; j++){
			double f_ij = sys.lubricationForceFactor(i, j);
			if (f_ij != 0){
				if (contact_pair[i][j] != 0){
					f_ij += -fc[contact_pair[i][j]].f_normal;
				}
				fout_yap << "r " << yap_force_factor*f_ij << endl;
				total_normal_force += f_ij;
				vec3d pos1 = shiftUpCoordinate(sys.position[i].x - sys.lx2,
											   sys.position[i].y - sys.ly2,
											   sys.position[i].z - sys.lz2);
				vec3d pos2 = shiftUpCoordinate(sys.position[j].x - sys.lx2,
											   sys.position[j].y - sys.ly2,
											   sys.position[j].z - sys.lz2);
				
				drawLine2('s', pos1, pos2, fout_yap);
			}
		}
	}
	/*
	 * Output the sum of the normal forces.
	 */
	fout_force << sys.dt * sys.ts << ' ' << total_normal_force << endl;
	
	/* Layer 6: Box and guide lines
	 */
	
	fout_yap << "y 6\n";
	fout_yap << "@ " << color_blue << endl;
	if (sys.dimension == 2){
		drawLine(-sys.lx2, 0, -sys.lz2,  sys.lx2, 0, -sys.lz2, fout_yap);
		drawLine(-sys.lx2, 0,  sys.lz2, -sys.lx2, 0, -sys.lz2, fout_yap);
		drawLine(-sys.lx2, 0,  sys.lz2,  sys.lx2, 0,  sys.lz2, fout_yap);
		drawLine( sys.lx2, 0, -sys.lz2,  sys.lx2, 0,  sys.lz2, fout_yap);
		/*
		for (int i = 0 ; i < 10; i ++){
			drawLine(-sys.lx2                , 0, -sys.lz2 + i*0.2*sys.lz2,
					 -sys.lx2 + i*0.2*sys.lz2, 0,                 -sys.lz2, fout_yap);
			drawLine( sys.lx2                   ,        0, sys.lz2 - i*0.2*sys.lz2,
					  sys.lx2 - i*0.2*sys.lz2, 0,  sys.lz2, fout_yap);
		}
		drawLine(-sys.lx2, 0, sys.lz2, sys.lx2 , 0, -sys.lz2, fout_yap);
		 */
	} else {
		drawLine(-sys.lx2, -sys.ly2, -sys.lz2,  sys.lx2, -sys.ly2, -sys.lz2, fout_yap);
		drawLine(-sys.lx2,  sys.ly2, -sys.lz2,  sys.lx2,  sys.ly2, -sys.lz2, fout_yap);

		drawLine(-sys.lx2,  sys.ly2,  sys.lz2,  sys.lx2,  sys.ly2,  sys.lz2, fout_yap);
		drawLine(-sys.lx2, -sys.ly2,  sys.lz2,  sys.lx2, -sys.ly2,  sys.lz2, fout_yap);
		
		drawLine(-sys.lx2, -sys.ly2, -sys.lz2, -sys.lx2, sys.ly2, -sys.lz2, fout_yap);
		drawLine(-sys.lx2, -sys.ly2,  sys.lz2, -sys.lx2, sys.ly2,  sys.lz2, fout_yap);
		drawLine( sys.lx2, -sys.ly2,  sys.lz2,  sys.lx2, sys.ly2,  sys.lz2, fout_yap);
		drawLine( sys.lx2, -sys.ly2, -sys.lz2,  sys.lx2, sys.ly2, -sys.lz2, fout_yap);
		
		drawLine( sys.lx2,  sys.ly2,  sys.lz2,  sys.lx2, sys.ly2, -sys.lz2,  fout_yap);
		drawLine(-sys.lx2,  sys.ly2,  sys.lz2, -sys.lx2, sys.ly2, -sys.lz2,  fout_yap);
		drawLine(-sys.lx2, -sys.ly2,  sys.lz2, -sys.lx2, -sys.ly2, -sys.lz2,  fout_yap);
		drawLine( sys.lx2, -sys.ly2,  sys.lz2,  sys.lx2, -sys.ly2, -sys.lz2,  fout_yap);
	}
}

