//
//  System.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__System__
#define __LF_DEM__System__

#include <iostream>
#include <fstream>
#include <queue>
#include "Interaction.h"
#include "State.h"

class System{
private:
	State state;
	int num_particle;
	double eta;
	double dt;
	
	/*
	 * Interparticle interactions
	 */
	Interaction *interaction;
	int **init_interacting_pair;
	int max_num_interaction; // Allowed length of interaction array.
	int num_interaction; // Length of used interaction array.
	queue<int> deactivated_interaction;

	/*
	 *  Simulation parameters
	 */
	double cutoff_distance;
	
	/*
	 * For output data.
	 */
	
	
	ofstream fout_tmp;

protected:
	
	void output_tmp(){
		for (int i=0; i < num_particle; i++){
			fout_tmp << "c " << state.position[i].x - state.lx2<< ' ';
			fout_tmp << state.position[i].y - state.ly2<< ' ';
			fout_tmp << state.position[i].z - state.lz2<< endl;
		}
		for (int k=0; k < num_interaction; k++){
			if ( interaction[k].active){
				int i = interaction[k].particle_num[0];
				fout_tmp << "l " << state.position[i].x - state.lx2<< ' ';
				fout_tmp << state.position[i].y - state.ly2<< ' ';
				fout_tmp << state.position[i].z - state.lz2<< ' ';
				fout_tmp << state.position[i].x - state.lx2 + interaction[k].nr_vec.x << ' ';
				fout_tmp << state.position[i].y - state.ly2 + interaction[k].nr_vec.y << ' ';
				fout_tmp << state.position[i].z - state.lz2 + interaction[k].nr_vec.z << endl;
				i = interaction[k].particle_num[1];
				fout_tmp << "l " << state.position[i].x - state.lx2<< ' ';
				fout_tmp << state.position[i].y - state.ly2<< ' ';
				fout_tmp << state.position[i].z - state.lz2<< ' ';
				fout_tmp << state.position[i].x - state.lx2 - interaction[k].nr_vec.x << ' ';
				fout_tmp << state.position[i].y - state.ly2 - interaction[k].nr_vec.y << ' ';
				fout_tmp << state.position[i].z - state.lz2 - interaction[k].nr_vec.z << endl;

			}
		}
		
			 
		fout_tmp << "l " << state.x_shift - state.lx2 << ' ' << 0 << ' ' << state.lz2 ;
		fout_tmp << ' '  << state.x_shift - state.lx2 << ' ' << 0 << ' ' << state.lz2 +2 << endl;
		fout_tmp << "l " << - state.lx2 << ' ' << 0 << ' ' << - state.lz2 ;
		fout_tmp << ' '  << state.lx2 << ' ' << 0 << ' ' << - state.lz2 << endl;
		fout_tmp << "l " << - state.lx2 << ' ' << 0 << ' ' << state.lz2 ;
		fout_tmp << ' '  << state.lx2 << ' ' << 0 << ' ' << state.lz2 << endl;
		fout_tmp << "l " << - state.lx2 << ' ' << 0 << ' ' << - state.lz2 ;
		fout_tmp << ' '  << - state.lx2 << ' ' << 0 << ' ' << + state.lz2 << endl;
		fout_tmp << "l " << state.lx2 << ' ' << 0 << ' ' << - state.lz2 ;
		fout_tmp << ' '  << state.lx2 << ' ' << 0 << ' ' << + state.lz2 << endl;
		fout_tmp << endl;
	}
	
	void initInteractingPair(){
		init_interacting_pair = new int * [num_particle];
		for (int i=0; i < num_particle; i++){
			init_interacting_pair[i] = new int [num_particle];
		}
		for (int i=0; i < num_particle-1; i++){
			for (int j=i+1; j < num_particle; j++){
				init_interacting_pair[i][j] = -1;
			}
		}
	}
	
	void checkBreak(){
		for (int i=0; i < num_particle; i++){
			if ( interaction[i].active ){
				if (interaction[i].r > cutoff_distance ){
					interaction[i].active = false;
					init_interacting_pair[ interaction[i].particle_num[0]][ interaction[i].particle_num[1]] = -1;
					deactivated_interaction.push(i);
				}
			}
		}
	}
	
	void timeEvolution(){
	
		fout_tmp.open("tmp.yap");
		state.x_shift = 0;
		for (int tt = 0 ; tt < 10000 ; tt ++){
			checkContact();
			if (true)
				output_tmp();
			for (int i=0; i < num_particle; i++){
				state.force[i].reset();
			}
						
			for (int k=0; k < num_interaction; k++){
				interaction[k].calcInteraction();
			}
			
			for (int i=0; i < num_particle; i++){
				vec3d U_inf(0.1*state.position[i].z, 0, 0);
				state.velocity[i] = (1.0/eta)*state.force[i] + U_inf;
			}

			state.x_shift += 0.1*state.lz*dt;
			for (int i=0; i < num_particle; i++){
				state.position[i] += state.velocity[i]*dt;
			}
			
			checkBreak();
//			double ave_velocity= sum_velocity / num_particle;
//			for (int i=0; i < num_particle; i++){
//				state.velocity[i] =  state.velocity[i]/ ave_velocity;
//			}
			checkPeriodicBoundary();
		}
		fout_tmp.close();
	}
	
	void checkContact(){
		for (int i=0; i < num_particle-1; i++){
			for (int j=i+1; j < num_particle; j++){
				if ( init_interacting_pair[i][j] == -1){
					double sq_distance = state.sq_distance(i, j);
					if ( sq_distance < 4){
						int new_num_interaction;
						if (deactivated_interaction.empty()){
							new_num_interaction = num_interaction;
							num_interaction ++;
						} else {
							cerr << deactivated_interaction.size();
							new_num_interaction = deactivated_interaction.front();
							deactivated_interaction.pop();
							cerr << "==> " << deactivated_interaction.size() << endl;
						}
						init_interacting_pair[i][j] = new_num_interaction; // (i > j)
						interaction[new_num_interaction].create(i,j);
					}
				}
			}
		}
	}
	
	void checkPeriodicBoundary(){
		if (state.x_shift > state.lx){
			state.x_shift -= state.lx;
		}

		for (int i=0; i < num_particle; i++){
			if ( state.position[i].z > state.lz ){
				state.position[i].z -= state.lz;
				state.position[i].x -= state.x_shift;
				
			} else if ( state.position[i].z < 0 ){
				state.position[i].z += state.lz;
				state.position[i].x += state.x_shift;
			}
			
			if ( state.position[i].x > state.lx ){
				state.position[i].x -= state.lx;
			} else if (state.position[i].x < 0 ){
				state.position[i].x += state.lx;
			}
			
			if ( state.position[i].y > state.ly ){
				state.position[i].y -= state.ly;
			} else if (state.position[i].y < 0 ){
				state.position[i].y += state.ly;
			}
			
		
		}
	}

	
public:
    /* For DEMsystem
     */
	System(){
		eta = 1;
		dt = 0.01;
		
		
		num_particle = 50;
		max_num_interaction = 10* num_particle;
		cutoff_distance = 2.5;
	};
	~System(){
		delete [] interaction;
	};
	
	void simulation(){
		state.setNumberParticle(num_particle);
		initInteractingPair();
		interaction = new Interaction [max_num_interaction];
		for (int i = 0; i < max_num_interaction; i++){
			interaction[i].init( &state );
		}
		state.setRandomPosition();
		timeEvolution();
	}
};


#endif /* defined(__LF_DEM__System__) */











