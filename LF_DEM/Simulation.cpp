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
	fout_rheo.close();
};

void
Simulation::SetParameters(int argc, const char * argv[]){
	sys.lub = true;
	sys.brownian = false;
	filename_import_positions = argv[1];
	sys.lubcore = atof(argv[2]);
	string arg3 =argv[3];
	if (arg3 == "f" ){
		sys.friction = true;
		cerr <<  sys.friction << endl;
	} else if (arg3 == "nf"){
		sys.friction = false;
	} else {
		cerr << "arg3 'f' for friction / 'nf' for no friction" << endl;
		exit(1);
	}

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
	cerr << "L = " << sys.lx << ' ' << sys.ly << ' ' << sys.lz << endl;
	cerr << "VF = " << sys.volume_fraction << endl;
	/*
	 * Simulation parameters
	 */

	/*
	 * Viscosity is not implemented in the simulation.
	 * eta = 1.0 is assumed.
	 */
	//	sys.eta = 1.0;
	/*
	 * Shear rate
	 */
	sys.shear_rate = 1.0;
	/*
	 * Temperature
	 */
	sys.kb_T=1.0;
	/*
	 * Simulation terminate at this value
	 */
	shear_strain = 100.;
	/*
	 * dt: the time step to integrate the equation of motion.
	 * We need to give a good criterion to give.
	 */
	sys.dt = 1e-5 / sys.shear_rate;
	/* 
	 * dt_mid: the intermediate time step for the mid-point 
	 * algortithm. dt/dt_mid = dt_ratio
	 * Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 * ASD code from Brady has dt_ratio=150
	 */
	sys.dt_ratio = 100;
	sys.dt_mid = sys.dt/sys.dt_ratio;
	/*
	 * The time steps finishing simulation.
	 */
	ts_max = (int)(shear_strain / sys.dt);
	/*
	 * Range of lubrication force
	 */
	sys.lub_max = 2.5;
	sys.sq_lub_max = sys.lub_max*sys.lub_max; // square of lubrication cutoff length.
	/*
	 * Contact force parameters
	 *
	 */
	sys.kn = 100; // normal spring constant
	sys.kt = 100; // tangential spring constant
	/*
	 * Particles are spined by the background vorticity.
	 * Small friction coeffient may not stop the sliding between surfaces
	 * of spinning particles, when normal force is small.
	 * We should also estimate the effect of lubrication torque.
	 *
	 */
	sys.mu_static = 1; // static friction coeffient
	sys.mu_dynamic = 0.8; // dynamic friction coeffient
	/*
	 * This is a threshold velocity to swich from dynamic friction to
	 * static friction. But this is a temporal provísional.
	 * There is no reference to give this value.
	 */
	sys.dynamic_friction_critical_velocity = 0.01;
	/*
	 * Visualization
	 */
	/*
	 * For yaplot output data,
	 * rotation of disk (2D) is visualized by cross.
	 */
	if (sys.dimension == 2)
		sys.draw_rotation_2d = true;
	else
		sys.draw_rotation_2d = false;
	/*
	 * snapshot for yaplot data.
	 */
	interval_snapshot = 100;
	/*
	 * The bond width indicates the force strength.
	 */
	yap_force_factor = 0.05;
	/*
	 * The middle height of the simulation box is set to the flow zero level.
	 */
	origin_zero_flow = true;
	/*
	 * Set simulation name and name of output files.
	 */
	sys.prepareSimulationName();
	string yap_filename = "yap_" + sys.simu_name + ".yap";
	string vel_filename = "rheo_" + sys.simu_name + ".dat";
	string vpy_filename = "vpy_" + sys.simu_name + ".dat";
	fout_yap.open(yap_filename.c_str());
	fout_rheo.open(vel_filename.c_str());
	fout_vpy.open(vpy_filename.c_str());
}

void
Simulation::importInitialPositionFile(){
	fstream file_import;
	file_import.open( filename_import_positions.c_str());
	vec3d pos;
	while (true){
		file_import >> pos.x >> pos.y >> pos.z;
		if (file_import.eof())
			break;
		initial_positions.push_back(pos);
	}
	file_import.close();
}

/*
 * Main simulation
 */
void
Simulation::SimulationMain(int argc, const char * argv[]){
	SetParameters(argc, argv);
	importInitialPositionFile();
	sys.n = (int)initial_positions.size();
	cerr << "N = " << sys.n  << endl;
	sys.prepareSimulation();
	for (int i=0; i < sys.n; i++){
		sys.position[i] = initial_positions[i];
		sys.angle[i] = 0;
	}

	sys.checkNewInteraction();
	//	int count = 0;
	//	while(count++<3){
	double time=0;
	while(true){
		sys.timeEvolution(interval_snapshot);
		outputRheologyData();
		output_yap();
		time+=(double)interval_snapshot;
		output_vpython(time);
	}

}

/*
 *
 */
void
Simulation::outputRheologyData(){
	sys.calcStress();
	/*
	 * Output the sum of the normal forces.
	 */
	fout_rheo << sys.dt * sys.ts << ' ';// 1
	fout_rheo << sys.mean_lub_stress[2] + sys.mean_contact_stress[2]  << ' ' ; //2
	fout_rheo << sys.mean_lub_stress[0] << ' ' ; //3
	fout_rheo << sys.mean_lub_stress[1] << ' ' ; //4
	fout_rheo << sys.mean_lub_stress[2] << ' ' ; //5
	fout_rheo << sys.mean_lub_stress[3] << ' ' ; //6
	fout_rheo << sys.mean_lub_stress[4] << ' ' ; //7
	fout_rheo << sys.mean_contact_stress[0] << ' ' ; //8
	fout_rheo << sys.mean_contact_stress[1] << ' ' ; //9
	fout_rheo << sys.mean_contact_stress[2] << ' ' ; //10
	fout_rheo << sys.mean_contact_stress[3] << ' ' ; //11
	fout_rheo << sys.mean_contact_stress[4] << ' ' ; //12
	fout_rheo << endl;
}

vec3d
Simulation::shiftUpCoordinate(double x, double y, double z){
	if (origin_zero_flow){
		z += sys.lz2;
		if (z > sys.lz2){
			x += - sys.shear_disp;
		if ( x < - sys.lx2)
			x += sys.lx;
			z -=  sys.lz;
		}
	}
	return vec3d(x,y,z);
}

void
Simulation::drawLine(char type , vec3d pos, vec3d vec, ofstream &fout){
	fout << type << ' ';
	fout << pos.x << ' '<< pos.y << ' '<< pos.z << ' ';
	fout << pos.x + vec.x << ' '<< pos.y + vec.y << ' '<< pos.z + vec.z << endl;
}

void
Simulation::drawLine2(char type , vec3d pos1, vec3d pos2, ofstream &fout){
	vec3d seg = pos2 - pos1;
	fout << type << ' ';
	if (seg.z > sys.lz2){
		pos2.z -= sys.lz;
		pos2.x -= sys.shear_disp;
		seg = pos2 - pos1;
	} else if (seg.z < -sys.lz2){
		pos2.z += sys.lz;
		pos2.x += sys.shear_disp;
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

void
Simulation::drawLine(double x0, double y0, double z0,
			  double x1, double y1, double z1,
			  ofstream &fout){
	fout << 'l' << ' ';
	fout << x0 << ' ' << y0 << ' ' << z0 << ' ';
	fout << x1 << ' ' << y1 << ' ' << z1 << endl;
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
	fout_yap << "r 1\n";
	fout_yap << "@ " << color_white << endl;
	vec3d pos;
	for (int i=0; i < sys.n; i++){
		pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
								sys.position[i].y - sys.ly2,
								sys.position[i].z - sys.lz2);
		fout_yap << "c " << pos.x << ' ' << pos.y << ' ' << pos.z << endl;
	}

	/* Layer 4: Orientation of particle (2D simulation)
	 */
	
	if (sys.draw_rotation_2d){
		fout_yap << "y 5\n";
		fout_yap << "@ " << color_white << endl;
		for (int i=0; i < sys.n; i++){
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
		for (int k=0; k < sys.num_interaction; k++){
			if ( sys.interaction[k].contact){
				if (sys.interaction[k].static_friction)
					fout_yap << "@ " << color_green << endl;
				else
					fout_yap << "@ " << color_orange << endl;
				
				int i = sys.interaction[k].particle_num[0];
				pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
										sys.position[i].y - sys.ly2,
										sys.position[i].z - sys.lz2);
				fout_yap << "r " << yap_force_factor*sys.interaction[k].f_tangent.norm()  << endl;
				drawLine('s', pos, -sys.interaction[k].nr_vec, fout_yap);
				int j = sys.interaction[k].particle_num[1];
				pos = shiftUpCoordinate(sys.position[j].x - sys.lx2,
										sys.position[j].y - sys.ly2,
										sys.position[j].z - sys.lz2);
				drawLine('s', pos, sys.interaction[k].nr_vec, fout_yap);
			}
		}
	}
	/* Layer 3: Normal
	 * Lubrication + contact force
	 */
	fout_yap << "y 3\n";
	fout_yap << "@ " << color_yellow << endl;
	for (int k=0; k < sys.num_interaction; k++){
		if ( sys.interaction[k].active ){
			int i = sys.interaction[k].particle_num[0];
			int j = sys.interaction[k].particle_num[1];
			fout_yap << "r " << yap_force_factor*sys.interaction[k].valNormalForce() << endl;
			pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
									sys.position[i].y - sys.ly2,
									sys.position[i].z - sys.lz2);
			drawLine('s', pos, -sys.interaction[k].nr_vec, fout_yap);
			pos = shiftUpCoordinate(sys.position[j].x - sys.lx2,
									sys.position[j].y - sys.ly2,
									sys.position[j].z - sys.lz2);
			drawLine('s', pos, sys.interaction[k].nr_vec, fout_yap);
		}
	}

	/* Layer 3: Normal
	 * Lubrication + contact force
	 */
//	fout_yap << "y 4\n";
//	fout_yap << "@ " << color_yellow << endl;
//	for (int i=0; i < sys.n; i++){
//		for (int j = i+1; j < sys.n; j++){
//			double f_ij = 0;
//			if (sys.lub == true){
//				f_ij += sys.lubricationForceFactor(i, j);
//			}
//			if (contact_pair[i][j] != -1){
//				f_ij += -fc[contact_pair[i][j]].f_normal;
//			}
//			if (f_ij > 0){
//				fout_yap << "r " << yap_force_factor*f_ij << endl;
//				vec3d pos1 = shiftUpCoordinate(sys.position[i].x - sys.lx2,
//											   sys.position[i].y - sys.ly2,
//											   sys.position[i].z - sys.lz2);
//				vec3d pos2 = shiftUpCoordinate(sys.position[j].x - sys.lx2,
//											   sys.position[j].y - sys.ly2,
//											   sys.position[j].z - sys.lz2);
//				
//				drawLine2('s', pos1, pos2, fout_yap);
//			}
//		}
//	}

		
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



/* Output data for vpython visualization.
 *
 */
void
Simulation::output_vpython(double time){
	vec3d pos;
	fout_vpy << "time: " << time << endl;
	for (int i=0; i < sys.n; i++){
		pos = shiftUpCoordinate(sys.position[i].x - sys.lx2,
								sys.position[i].y - sys.ly2,
								sys.position[i].z - sys.lz2);
		fout_vpy << i << ' ' << pos.x << ' ' << pos.y << ' ' << pos.z << endl;
	}
	fout_vpy << endl;
}
