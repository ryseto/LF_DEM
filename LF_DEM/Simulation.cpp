//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "Simulation.h"
#include <cmath>

Simulation::Simulation(){};
Simulation::~Simulation(){
	fout_yap.close();
	fout_vel.close();
	delete [] interaction;
	for (int i=0; i < num_particle; i++){
		delete [] interacting_pair[i];
	}
	delete [] interacting_pair;
};


void Simulation::SetParameters(int argc, const char * argv[]){
	sys.friction = true;
	sys.lub = true;
	
	sys.dimension = 2;
	sys.lx = 20;
	sys.ly = 20;
	sys.lz = 20;
	sys.volume_fraction = 0.70;
	sys.lubcore = 2.0;
	if ( argc == 2){
		sys.volume_fraction = atof(argv[1]);
	} else if ( argc == 3){
		sys.volume_fraction = atof(argv[1]);
		sys.lubcore = atof(argv[2]);
	} else if ( argc == 4){
		sys.volume_fraction = atof(argv[1]);
		sys.lubcore = atof(argv[2]);
		double l = atof(argv[3]);
		sys.lx = l;
		sys.ly = l;
		sys.lz = l;
	}

	/*
	 * Simulation parameters
	 */
	sys.eta = 1.0; // viscosity
	sys.shear_rate = 1; // shear rate
	shear_strain = 100; // shear strain
	cutoff_distance = 2.5; // to delete possible neighbor
	sys.sq_lub_max = 2.5*2.5; // square of lubrication cutoff length.
	sys.dt = 1e-4 / sys.shear_rate; //time step.
	ts_max = (int)(shear_strain / sys.dt); // time step max

	/*
	 * Contact force parameters
	 */
	sys.kn = 200; // normal spring constant
	sys.kt = 200; // tangential spring constant
	sys.mu_static = 0.3; // static friction coeffient
	sys.mu_dynamic = 0.2; // dynamic friction coeffient
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
	max_num_interaction = 6 * num_particle;
	sys.init();
	string yap_filename = "yap_" + sys.simu_name + ".yap";
	string vel_filename = "vel_" + sys.simu_name + ".dat";
	fout_yap.open(yap_filename.c_str());
	fout_vel.open(vel_filename.c_str());
}

/*
 * Main simulation
 */
void Simulation::SimulationMain(int argc, const char * argv[]){
	SetParameters(argc, argv);
	sys.setNumberParticle(num_particle);
	initInteractingPair();
	interaction = new Interaction [max_num_interaction];
	for (int i = 0; i < max_num_interaction; i++){
		interaction[i].init( &sys );
	}
	cerr << "set initial positions" << endl;
	sys.setRandomPosition();
	cerr << "start simulation" << endl;
	timeEvolution();
	cerr << "finished" << endl;
}

/*
 * Initialize interacting_pair object.
 * only elements of j > i will be used.
 */
void Simulation::initInteractingPair(){
	interacting_pair = new int * [num_particle];
	for (int i=0; i < num_particle; i++){
		interacting_pair[i] = new int [num_particle];
	}
	for (int i=0; i < num_particle-1; i++){
		for (int j=i+1; j < num_particle; j++){
			interacting_pair[i][j] = -1;
		}
	}
}

/* Interaction is gone when two particles are separated enough.
 *
 */
void Simulation::checkBreak(){
	for (int k = 0; k < num_interaction; k++){
		if ( interaction[k].active
			&& interaction[k].r > cutoff_distance){
			interaction[k].active = false;
			interacting_pair[interaction[k].particle_num[0]][interaction[k].particle_num[1]] = -1;
			deactivated_interaction.push(k);
		}
	}
}

/* Check the distance between separating particles.
 * i < j 
 *
 * A patch-up prescription to aboid 
 * interacting_pair[i][j] < 0 indicates separating particles to be checked.
 * interacting_pair[i][j] = -1, the particles are near contact. So every time step, distance should be checked.a
 * interacting_pair[i][j] < -1, the particles have some distance.
 */
void Simulation::checkContact(){
	for (int i=0; i < num_particle-1; i++){
		for (int j=i+1; j < num_particle; j++){
			if ( interacting_pair[i][j] == -1){
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
					interaction[new_num_interaction].create(i,j);
					interacting_pair[i][j] = new_num_interaction;
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
	for (int ts = 0 ; ts < ts_max ; ts ++){
		if (ts % interval_snapshot == 0){
			output_yap();
		}
		checkContact();
		sys.forceReset();
		sys.torqueReset();
		if (sys.friction){
			for (int k=0; k < num_interaction; k++){
				interaction[k].calcInteraction();
			}
		} else {
			for (int k=0; k < num_interaction; k++){
				interaction[k].calcInteractionNoFriction();
			}
		}
		
		// Free-draining approximation
		//sys.updateVelocity();
		sys.updateVelocityLubrication();
		sys.deltaTimeEvolution();
		if (sys.friction){
			for (int k = 0; k < num_interaction; k++){
				interaction[k].incrementTangentialDisplacement();
			}
		}
		checkBreak();
		
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
			if ( interaction[k].active && interaction[k].r < 2 ){
				if (interaction[k].static_friction)
					fout_yap << "@ " << color_green << endl;
				else
					fout_yap << "@ " << color_orange << endl;
				
				int i = interaction[k].particle_num[0];
				pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
										sys.position[i].y - sys.ly2,
										sys.position[i].z - sys.lz2);
				fout_yap << "r " << yap_force_factor*interaction[k].f_tangent.norm()  << endl;
				drawLine('s', pos, interaction[k].nr_vec, fout_yap);
				int j = interaction[k].particle_num[1];
				pos = shiftUpCoordinate(sys.position[j].x - sys.lx2,
										sys.position[j].y - sys.ly2,
										sys.position[j].z - sys.lz2);
				drawLine('s', pos, -interaction[k].nr_vec, fout_yap);
			}
		}
	}
	/* Layer 3: Normal
	 */
	fout_yap << "y 3\n";
	fout_yap << "@ " << color_yellow << endl;
	for (int k=0; k < num_interaction; k++){
		if ( interaction[k].active && interaction[k].r < 2 ){
			int i = interaction[k].particle_num[0];
			pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
									sys.position[i].y - sys.ly2,
									sys.position[i].z - sys.lz2);
			fout_yap << "r " << yap_force_factor*abs(interaction[k].f_normal) << endl;
			drawLine('s', pos, interaction[k].nr_vec, fout_yap);
			int j = interaction[k].particle_num[1];
			pos = shiftUpCoordinate(sys.position[j].x - sys.lx2,
									sys.position[j].y - sys.ly2,
									sys.position[j].z - sys.lz2);
			drawLine('s', pos, -interaction[k].nr_vec, fout_yap);
		}
	}

	fout_yap << "y 8\n";
	fout_yap << "@ " << color_orange << endl;
	for (int i=0; i < num_particle; i++){
		for (int j=i+1; j < num_particle; j++){
			double f_ij = sys.lubricationForceFactor(i, j);
			if (f_ij != 0){
				fout_yap << "r " << yap_force_factor*f_ij << endl;
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

void Simulation::output_vel(){
	for (int k=0; k < num_interaction; k++){
		if ( interaction[k].active && interaction[k].static_friction == false ){
			fout_vel << sqrt(interaction[k].sqnorm_contact_velocity) << ' ';
		}
	}
	fout_vel << endl;
}
