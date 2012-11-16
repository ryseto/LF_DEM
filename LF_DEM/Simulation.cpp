//
//  Simulation.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "Simulation.h"

/*
 *
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

void Simulation::checkBreak(){
	for (int i=0; i < num_particle; i++){
		if ( interaction[i].active ){
			if (interaction[i].r > cutoff_distance ){
				interaction[i].active = false;
				interacting_pair[ interaction[i].particle_num[0]][ interaction[i].particle_num[1]] = -1;
				deactivated_interaction.push(i);
			}
		}
	}
}



void Simulation::output_tmp(){
	fout_tmp << "@ 2\n";
	for (int i=0; i < num_particle; i++){
		fout_tmp << "c " << sys.position[i].x - sys.lx2<< ' ';
		fout_tmp << sys.position[i].y - sys.ly2<< ' ';
		fout_tmp << sys.position[i].z - sys.lz2<< endl;
	}
	fout_tmp << "@ 3\n";
	for (int k=0; k < num_interaction; k++){
		if ( interaction[k].active){
			int i = interaction[k].particle_num[0];
			fout_tmp << "l " << sys.position[i].x - sys.lx2<< ' ';
			fout_tmp << sys.position[i].y - sys.ly2<< ' ';
			fout_tmp << sys.position[i].z - sys.lz2<< ' ';
			fout_tmp << sys.position[i].x - sys.lx2 + interaction[k].nr_vec.x << ' ';
			fout_tmp << sys.position[i].y - sys.ly2 + interaction[k].nr_vec.y << ' ';
			fout_tmp << sys.position[i].z - sys.lz2 + interaction[k].nr_vec.z << endl;
			i = interaction[k].particle_num[1];
			fout_tmp << "l " << sys.position[i].x - sys.lx2<< ' ';
			fout_tmp << sys.position[i].y - sys.ly2<< ' ';
			fout_tmp << sys.position[i].z - sys.lz2<< ' ';
			fout_tmp << sys.position[i].x - sys.lx2 - interaction[k].nr_vec.x << ' ';
			fout_tmp << sys.position[i].y - sys.ly2 - interaction[k].nr_vec.y << ' ';
			fout_tmp << sys.position[i].z - sys.lz2 - interaction[k].nr_vec.z << endl;
			
		}
	}
	fout_tmp << "@ 4\n";
	
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