//
//  Simulation.h
//  LF_DEM
//
//  Created by Ryohei Seto on 11/15/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef __LF_DEM__Simulation__
#define __LF_DEM__Simulation__

#include <iostream>
#include <fstream>
#include <queue>
#include "Interaction.h"
#include "System.h"

class Simulation{
private:
	System sys;
	int num_particle;
	double eta;
	double dt;
	
	/*
	 * Interparticle interactions
	 */
	Interaction *interaction;
	int **interacting_pair; // Table
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

	void output_tmp();
	void initInteractingPair();
	void checkBreak();
	


	
	void timeEvolution(){
		fout_tmp.open("tmp.yap");
		sys.x_shift = 0;
		for (int tt = 0 ; tt < 10000 ; tt ++){
			if (tt % 10 == 0 )
				output_tmp();

			checkContact();
			for (int i=0; i < num_particle; i++){
				sys.force[i].reset();
			}
			for (int k=0; k < num_interaction; k++){
				interaction[k].calcInteraction();
			}
			
			for (int i=0; i < num_particle; i++){
				vec3d U_inf(0.1*sys.position[i].z, 0, 0);
				sys.velocity[i] = (1.0/eta)*sys.force[i] + U_inf;
			}

			sys.x_shift += 0.1*sys.lz*dt;
			for (int i=0; i < num_particle; i++){
				sys.position[i] += sys.velocity[i]*dt;
			}
			checkBreak();
			checkPeriodicBoundary();
		}
		fout_tmp.close();
	}
	
	void checkContact(){
		for (int i=0; i < num_particle-1; i++){
			for (int j=i+1; j < num_particle; j++){
				if ( interacting_pair[i][j] == -1){
					double sq_distance = sys.sq_distance(i, j);
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
						interacting_pair[i][j] = new_num_interaction; // (i > j)
						interaction[new_num_interaction].create(i,j);
					}
				}
			}
		}
	}
	
	void checkPeriodicBoundary(){
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

	
public:
    /* For DEMsystem
     */
	Simulation(){
		eta = 1;
		dt = 0.01;
		num_particle = 50;
		max_num_interaction = 10* num_particle;
		cutoff_distance = 2.5;
	};
	~Simulation(){
		delete [] interaction;
	};
	
	void run(){
		sys.setNumberParticle(num_particle);
		initInteractingPair();
		interaction = new Interaction [max_num_interaction];
		for (int i = 0; i < max_num_interaction; i++){
			interaction[i].init( &sys );
		}
		sys.setRandomPosition();
		timeEvolution();
	}
};


#endif /* defined(__LF_DEM__System__) */











