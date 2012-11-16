//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "Simulation.h"

Simulation::Simulation(){
	dimension = 3;
	num_particle = 80;
	sys.lx = 10;
	sys.ly = 10;
	sys.lz = 10;
	eta = 1;
	dt = 0.002;
	ts_max = 100000;
	max_num_interaction = 6* num_particle;
	cutoff_distance = 2.5;
	sys.shear_rate = 0.1;
	sys.kn = 100;
	sys.kt = 100;
	sys.mu_static = 10;
	sys.mu_dynamic = 9;

	sys.init();
};

Simulation::~Simulation(){
	delete [] interaction;
};

/*
 * Main simulation
 */
void Simulation::run(){
	sys.setNumberParticle(num_particle);
	initInteractingPair();
	interaction = new Interaction [max_num_interaction];
	for (int i = 0; i < max_num_interaction; i++){
		interaction[i].init( &sys );
	}
	cerr << "set initial positions" << endl;
	sys.setRandomPosition(dimension);
	cerr << "done" << endl;
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
		if ( interaction[k].active ){
			if (interaction[k].r > cutoff_distance ){
				interaction[k].active = false;
				interacting_pair[ interaction[k].particle_num[0]][ interaction[k].particle_num[1]] = -1;
				deactivated_interaction.push(k);
			}
		}
	}
}


/* Check the distance between separating particles.
 * i < j 
 */
void Simulation::checkContact(){
	for (int i=0; i < num_particle-1; i++){
		for (int j=i+1; j < num_particle; j++){
			if ( interacting_pair[i][j] == -1){
				double sq_distance = sys.sq_distance(i, j);
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


/* Periodic boundary condition for simple shear.
 * x : x <--> x+lx 
 * y : y <--> y+ly
 * z : (x,z) <--> (x + gammadot*t, z+lz)
 */
void Simulation::checkPeriodicBoundary(){
	if (sys.x_shift > sys.lx){
		sys.x_shift -= sys.lx;
	}
	
	for (int i=0; i < num_particle; i++){
		if ( sys.position[i].z > sys.lz ){
			sys.position[i].z -= sys.lz;
			sys.position[i].x -= sys.x_shift;
			
		} else if ( sys.position[i].z < 0 ){
			sys.position[i].z += sys.lz;
			sys.position[i].x += sys.x_shift;
		}
		
		if ( sys.position[i].x > sys.lx ){
			sys.position[i].x -= sys.lx;
		} else if (sys.position[i].x < 0 ){
			sys.position[i].x += sys.lx;
		}
		
		if ( sys.position[i].y > sys.ly ){
			sys.position[i].y -= sys.ly;
		} else if (sys.position[i].y < 0 ){
			sys.position[i].y += sys.ly;
		}
	}
}



/* Simulation for the time evolution.
 *
 */
void Simulation::timeEvolution(){
	fout_tmp.open("tmp.yap");
	sys.x_shift = 0;
	for (int ts = 0 ; ts < ts_max ; ts ++){
		if (ts % 30 == 0)
			output_tmp();
		
		checkContact();
		for (int i=0; i < num_particle; i++){
			sys.force[i].reset();
			sys.torque[i].reset();
		}
		for (int k=0; k < num_interaction; k++){
			interaction[k].calcInteraction();
		}
		
		for (int i=0; i < num_particle; i++){
			vec3d U_inf(sys.shear_rate*sys.position[i].z, 0, 0);
			vec3d O_inf(0, 0.5*sys.shear_rate, 0);
			sys.velocity[i] = (1.0/eta)*sys.force[i] + U_inf;
			sys.ang_velocity[i] = (1.33333/eta)*sys.torque[i] + O_inf;
		}
		
		sys.x_shift += sys.shear_rate*sys.lz*dt;
		for (int i=0; i < num_particle; i++){
			sys.position[i] += sys.velocity[i]*dt;
		}
		for (int k=0; k < num_interaction; k++){
			interaction[k].incrementTangentialDisplacement(dt);
		}
		checkBreak();
		checkPeriodicBoundary();
	}
	fout_tmp.close();
}

/* Output data for yaplot visualization.
 *
 */
void Simulation::output_tmp(){
	fout_tmp << "y 1\n";
	fout_tmp << "@ 2\n";
	for (int i=0; i < num_particle; i++){
		fout_tmp << "c " << sys.position[i].x - sys.lx2<< ' ';
		fout_tmp << sys.position[i].y - sys.ly2<< ' ';
		fout_tmp << sys.position[i].z - sys.lz2<< endl;
	}
	fout_tmp << "y 2\n";

	for (int k=0; k < num_interaction; k++){
		
		if ( interaction[k].active && interaction[k].r < 2 ){
			if (interaction[k].static_friction)
				fout_tmp << "@ 3\n";
			else
				fout_tmp << "@ 4\n";
			fout_tmp << "r " << interaction[k].f_tangent.norm() << endl;
			int i = interaction[k].particle_num[0];
			fout_tmp << "s " << sys.position[i].x - sys.lx2<< ' ';
			fout_tmp << sys.position[i].y - sys.ly2 -0.1<< ' ';
			fout_tmp << sys.position[i].z - sys.lz2 << ' ';
			fout_tmp << sys.position[i].x - sys.lx2 + interaction[k].nr_vec.x << ' ';
			fout_tmp << sys.position[i].y - sys.ly2 + interaction[k].nr_vec.y -0.1<< ' ';
			fout_tmp << sys.position[i].z - sys.lz2 + interaction[k].nr_vec.z << endl;
			int j = interaction[k].particle_num[1];
			fout_tmp << "s " << sys.position[j].x - sys.lx2<< ' ';
			fout_tmp << sys.position[j].y - sys.ly2 -0.1<< ' ';
			fout_tmp << sys.position[j].z - sys.lz2<< ' ';
			fout_tmp << sys.position[j].x - sys.lx2 - interaction[k].nr_vec.x << ' ';
			fout_tmp << sys.position[j].y - sys.ly2 - interaction[k].nr_vec.y -0.1<< ' ';
			fout_tmp << sys.position[j].z - sys.lz2 - interaction[k].nr_vec.z << endl;
		}
	}
	
	fout_tmp << "y 3\n";
	fout_tmp << "@ 6\n";
	for (int k=0; k < num_interaction; k++){
		if ( interaction[k].active && interaction[k].r < 2 ){
			fout_tmp << "r " << abs(interaction[k].f_normal) << endl;
			int i = interaction[k].particle_num[0];
			fout_tmp << "s " << sys.position[i].x - sys.lx2<< ' ';
			fout_tmp << sys.position[i].y - sys.ly2 -0.1<< ' ';
			fout_tmp << sys.position[i].z - sys.lz2 << ' ';
			fout_tmp << sys.position[i].x - sys.lx2 + interaction[k].nr_vec.x << ' ';
			fout_tmp << sys.position[i].y - sys.ly2 + interaction[k].nr_vec.y -0.1<< ' ';
			fout_tmp << sys.position[i].z - sys.lz2 + interaction[k].nr_vec.z << endl;
			int j = interaction[k].particle_num[1];
			fout_tmp << "s " << sys.position[j].x - sys.lx2<< ' ';
			fout_tmp << sys.position[j].y - sys.ly2 -0.1<< ' ';
			fout_tmp << sys.position[j].z - sys.lz2<< ' ';
			fout_tmp << sys.position[j].x - sys.lx2 - interaction[k].nr_vec.x << ' ';
			fout_tmp << sys.position[j].y - sys.ly2 - interaction[k].nr_vec.y -0.1<< ' ';
			fout_tmp << sys.position[j].z - sys.lz2 - interaction[k].nr_vec.z << endl;
		}
	}

//	fout_tmp << "@ 5" << endl;
//	for (int k=0; k < num_interaction; k++){
//		if ( interaction[k].active && interaction[k].r < 2){
//			int i = interaction[k].particle_num[0];
//			int j = interaction[k].particle_num[1];
//			vec3d p_mid = 0.5*(sys.position[i] +sys.position[j]);
//			
//			p_mid.x += interaction[k].pd_x*sys.lx2 + 0.5*interaction[k].pd_z*sys.x_shift;
//			p_mid.y += interaction[k].pd_y*sys.ly2;
//			
//			p_mid.z += interaction[k].pd_z*sys.lz2;
//			
//			
//			fout_tmp << "r " << interaction[k].xi.norm() << endl;
//			fout_tmp << "c " << p_mid.x - sys.lx2 << ' ' << p_mid.y - sys.ly2 << ' ' << p_mid.z - sys.lz2 <<endl;
//		}
//	}
	fout_tmp << "y 3\n";
	fout_tmp << "@ 5\n";
	fout_tmp << "l " << sys.x_shift - sys.lx2 << ' ' << 0 << ' ' << sys.lz2 ;
	fout_tmp << ' '  << sys.x_shift - sys.lx2 << ' ' << 0 << ' ' << sys.lz2 +2 << endl;
	fout_tmp << "l " << - sys.lx2 << ' ' << 0 << ' ' << - sys.lz2 ;
	fout_tmp << ' '  << sys.lx2 << ' ' << 0 << ' ' << - sys.lz2 << endl;
	fout_tmp << "l " << - sys.lx2 << ' ' << 0 << ' ' << sys.lz2 ;
	fout_tmp << ' '  << sys.lx2 << ' ' << 0 << ' ' << sys.lz2 << endl;
	fout_tmp << "l " << - sys.lx2 << ' ' << 0 << ' ' << - sys.lz2 ;
	fout_tmp << ' '  << - sys.lx2 << ' ' << 0 << ' ' << + sys.lz2 << endl;
	fout_tmp << "l " << sys.lx2 << ' ' << 0 << ' ' << - sys.lz2 ;
	fout_tmp << ' '  << sys.lx2 << ' ' << 0 << ' ' << + sys.lz2 << endl;
	fout_tmp << endl;
}
