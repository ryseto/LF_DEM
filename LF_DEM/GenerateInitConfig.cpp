//
//  GenerateInitConfig.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/6/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include <stdlib.h> // necessary for Linux
#include "GenerateInitConfig.h"
#include "Simulation.h"

#define RANDOM ( rand_gen.rand() ) // RNG uniform [0,1]

using namespace std;

int
GenerateInitConfig::generate(int argc, const char * argv[]){
	setParameters(argc, argv);
	sys.np(np);
	sys.allocateRessources();
	sys.lx(lx);
	sys.ly(ly);
	sys.lz(lz);
	sys.setSystemVolume();
	sys.volume_fraction = volume_fraction;
	sys.lub_max = 2.5;
	sys.in_predictor = false;
	putRandom();
	sys.setupSystemForGenerateInit();
	grad = new vec3d [np];
	prev_grad = new vec3d [np];
	step_size = 10;
	gradientDescent();
	step_size /= 4.;
	gradientDescent();
	step_size /= 4.;
	gradientDescent();
	step_size /= 4.;
	gradientDescent();
	step_size /= 4.;
	gradientDescent();
	step_size /= 4.;
	gradientDescent();
	step_size /= 4.;
	gradientDescent();
	step_size /= 4.;
	gradientDescent();
	
	int count = 0;
	double energy = 0;
	ofstream fout;
	fout.open("energy_decay.dat");
	double energy_previous = 0;
	double diff_energy;
	do {
		energy = zeroTMonteCarloSweep();
		cerr << energy << endl;
		if (count % 100 == 0) {
			diff_energy = energy-energy_previous;
			fout << count << ' ' << energy << ' ' <<  diff_energy << endl;
			energy_previous = energy;
		}
		count ++;
	} while (abs(diff_energy) > 0.01);

	//sys.setupSystem(initial_positions, radii);
	//	outputDataHeader(fout_particle);
	//	outputConfigurationData();
	
	vector<vec3d> positions;
	vector<double> radii;
	positions.resize(sys.np());
	radii.resize(sys.np());
	for (int i=0; i<sys.np(); i++){
		positions[i] = sys.position[i];
		radii[i] = sys.radius[i];
	}
	Simulation simulation;
	simulation.RelaxationZeroShear(positions, radii, lx, ly, lz);
	//	solveOverlap();
	outputPositionData(positions,
					   radii);
	
	delete [] grad;
	delete [] prev_grad;
	
	return 0;
}

void
GenerateInitConfig::outputPositionData(vector<vec3d> &positions,
									   vector<double> &radii){
	ofstream fout;
	ostringstream ss_posdatafilename;
	ss_posdatafilename << "D" << dimension;
	ss_posdatafilename << "N" << np;
	ss_posdatafilename << "VF" << volume_fraction;
	if (number_ratio == 1) {
		ss_posdatafilename << "Mono";
	} else {
		ss_posdatafilename << "Poly" << a2 << "_" << number_ratio ;
	}

	if (dimension == 2) {
		if (lx_lz == 1) {
			ss_posdatafilename << "Square"; // square
		} else {
			ss_posdatafilename << "L" << (int)(10*lx_lz) << "_" << 10;
		}
	} else {
		if (lx_lz == 1 && ly_lz == 1) {
			ss_posdatafilename << "Cubic"; //
		} else {
			ss_posdatafilename << "L" << (int)(10*lx_lz) << "_" << (int)(10*ly_lz) << "_" << 10;
		}
	}
	ss_posdatafilename << "_" << rand_seed << ".dat";
	cerr << ss_posdatafilename.str() << endl;
	fout.open(ss_posdatafilename.str().c_str());
	fout << "# np1 np2 vf lx ly lz" << endl;
	fout << "# " << np1 << ' ' << np2 << ' ' << volume_fraction << ' ';
	fout << lx << ' ' << ly << ' ' << lz << endl;
	for (int i = 0; i < np ; i++){
		fout << positions[i].x << ' ';
		fout << positions[i].y << ' ';
		fout << positions[i].z << ' ';
		fout << radii[i] << endl;
	}
	fout.close();
}


double
GenerateInitConfig::computeGradient(){
	
	int i, j;
	double r, rcont;
	vec3d nr_vec;
	
	double amp, amp2;
	double energy=0.;
	
	for(i=0; i<np; i++){
		for(int u=0;u<sys.dimension;u++){
			grad[i].reset();
		}
	}
	
	
	for (int k=0; k<sys.num_interaction; k++){
		if(sys.interaction[k].contact){
			i = sys.interaction[k].par_num[0];
			j = sys.interaction[k].par_num[1];
			r = sys.interaction[k].r();
			rcont = sys.interaction[k].ro;
			nr_vec = sys.interaction[k].nr_vec;
			
			amp=(1./rcont-1./r); // negative
			amp2=4.*amp/rcont;
			
			grad[i] -= r*nr_vec*amp2;
			grad[j] += r*nr_vec*amp2;
			
			energy += 2*r*amp*amp;
		}
	}
	
	return energy;
	
}



void
GenerateInitConfig::moveAlongGradient(vec3d *g, int dir){
	double grad_norm;
	double gradient_power = 0.9;
	vec3d step;
	
	grad_norm = 0;
	for (int i=0; i<np;i++) {
		for(int u=0; u<sys.dimension; u++) {
			grad_norm += g[i].sq_norm();
		}
	}
	
	if (grad_norm != 0) {
		double rescale = pow(grad_norm, gradient_power);
		for(int i=0; i<np;i++){
			step = -dir*g[i]*step_size/rescale;
			sys.displacement(i, step);
		}
		sys.checkNewInteraction();
		sys.updateInteractions();
	}
	
}

void
GenerateInitConfig::storeGradient(){
	for (int i=0;i<np;i++) {
		prev_grad[i] = grad[i];
	}
}

double
GenerateInitConfig::gradientDescent(){
	double old_running_energy;
	double running_energy;
	double relative_en;

	long long int steps = 0;

	cerr << endl << " Gradient Descent..." << endl;
	
	storeGradient();
	
	running_energy=computeGradient();
	
	cerr << "  Starting Energy "<< running_energy/np << endl;
	
	do {
		old_running_energy = running_energy;
		
		moveAlongGradient(grad, 1);
		storeGradient();
		running_energy=computeGradient();
		
		relative_en=(old_running_energy-running_energy)/(old_running_energy+running_energy);
		
		if (steps%100 == 0) {
			cerr << "    Steps = " << steps << " :::   Energy : "<< running_energy/np << endl;
		}
		
		steps++;
		
	} while(relative_en > 1e-6);
	
	
	if(relative_en < 0) {
		cerr << "    Steps = " << steps;
		cerr << " :::   Last Step Upwards. Old Energy : " << old_running_energy/np;
		cerr << " New Energy : " << running_energy/np;
		cerr << " Relative : " << relative_en << endl;
		cerr << "      Reverting last step..." << endl;
		moveAlongGradient(prev_grad, -1);
		return old_running_energy;
	}

	if (relative_en > 0
		&& relative_en < 1e-6) {
		cerr << "    Steps = " << steps;
		cerr << " :::   Stuck: too slow (Relative energy difference : " << relative_en  << endl;
		cerr << "  Ending Energy "<< running_energy << endl<< endl;
		return running_energy;
	}
	return running_energy;
}

void
GenerateInitConfig::putRandom(){
	rand_gen.seed(rand_seed);
	for (int i=0; i < np; i++) {
		sys.position[i].x = lx*RANDOM;
		sys.position[i].z = lz*RANDOM;
		if (dimension == 2) {
			sys.position[i].y = ly2;
		} else {
			sys.position[i].y = ly*RANDOM;
		}
		if (i < np1) {
			sys.radius[i] = a1;
		} else {
			sys.radius[i] = a2;
		}
	}
	sys.setRadiusMax(a2);
}

void
GenerateInitConfig::updateInteractions(int i){
	vector <Interaction*> inter_list;
	for (set<Interaction*>::iterator it = sys.interaction_list[i].begin();
		 it != sys.interaction_list[i].end(); it ++) {
		inter_list.push_back(*it);
	}
	for (unsigned int k=0; k<inter_list.size(); k++) {
		bool desactivated;
		inter_list[k]->updateState(desactivated);
		if (desactivated) {
			sys.deactivated_interaction.push(inter_list[k]->label);
		}
	}
}

int
GenerateInitConfig::overlapNumber(int i){
	int overlaps = 0;
	for (set<Interaction*>::iterator it = sys.interaction_list[i].begin();
		 it != sys.interaction_list[i].end(); it ++) {
		if ((*it)->r() < (*it)->ro) {
			overlaps++;
		}
	}
	return overlaps;
}

double
GenerateInitConfig::particleEnergy(int i){
	double energy = 0;
	for (set<Interaction*>::iterator it = sys.interaction_list[i].begin();
		 it != sys.interaction_list[i].end(); it ++) {
		double r = (*it)->r();
		double rcont = (*it)->ro;
		if (r < rcont) {
			double amp = (r/rcont-1); // negative
			energy += 2*amp*amp;
		}
	}
	return energy;
}

double
GenerateInitConfig::zeroTMonteCarloSweep(){
	int steps = 0;
	int init_overlaps = 0;
	double init_energy = 0;
	for(int i=0; i<np; i++) {
	 	init_overlaps += overlapNumber(i);
	}
	for(int i=0; i<np; i++) {
	 	init_energy += particleEnergy(i);
	}
	while (steps < np) {
		int moved_part = (int)(RANDOM*np);
		//		int overlap_pre_move = overlapNumber(moved_part);
		double energy_pre_move = particleEnergy(moved_part);
		vec3d trial_move;
		if (sys.dimension == 3) {
			trial_move = randUniformSphere(0.002);
		} else {
			trial_move = randUniformCircle(0.002);
		}
		trial_move *= RANDOM;
		sys.displacement(moved_part, trial_move);
		updateInteractions(moved_part);
		
		//		int overlap_post_move = overlapNumber(moved_part);
		double energy_post_move = particleEnergy(moved_part);
		
		//		if( overlap_pre_move <= overlap_post_move ){
		if(energy_pre_move < energy_post_move) {
			sys.displacement(moved_part, -trial_move);
			updateInteractions(moved_part);
		}
		steps ++;
	}
	
	sys.checkNewInteraction();
	sys.updateInteractions();
	
	int final_overlaps = 0;
	double final_energy = 0;
	for(int i=0; i<np; i++) {
		final_overlaps += overlapNumber(i);
	}
	for(int i=0; i<np; i++) {
	 	final_energy += particleEnergy(i);
	}
	
	cerr << " MC sweep : init energy " << init_energy/np << " final energy " << final_energy/np;
	cerr << " init overlaps " << init_overlaps << " final overlaps " << final_overlaps << endl;
	
	return final_energy;
}

vec3d
GenerateInitConfig::randUniformSphere(double r){
	double z = 2*RANDOM - 1;
	double phi = 2*M_PI*RANDOM;
	double sin_theta = sqrt(1-z*z);
	return vec3d(r*sin_theta*cos(phi), r*sin_theta*sin(phi), r*z);
}

vec3d
GenerateInitConfig::randUniformCircle(double r){
	double phi = 2*M_PI*RANDOM;
	return vec3d(r*cos(phi), 0, r*sin(phi));
}

template<typename T>
T readStdinDefault(T default_value,	string message){
	string input;
	T value;
	cerr << message << "[" << default_value << "]: ";
	getline(cin, input);
	if ( !input.empty() ) {
		istringstream stream( input );
		stream >> value;
	} else {
		value = default_value;
	}
	cerr << value << endl;
	return value;
}

void
GenerateInitConfig::setParameters(int argc, const char * argv[]){
	/*
	 *  Read parameters from standard input
	 *
	 */
	np = readStdinDefault(200, "number of particle");
	dimension = readStdinDefault(3, "dimension (2 or 3)");
	sys.dimension = dimension;
	if (dimension == 2){
		volume_fraction =  readStdinDefault(0.7, "volume_fraction");
	} else {
		volume_fraction =  readStdinDefault(0.5, "volume_fraction");
	}
	lx_lz = readStdinDefault(1.0 , "Lx/Lz [1]: ");
	if (dimension == 3){
		ly_lz = 1;
		ly_lz = readStdinDefault(1.0 , "Ly/Lz [1]: ");
	}
	char m_p_disperse = readStdinDefault('m' , "(m)onodisperse or (p)olydisperse");
	number_ratio = 1.0; // mono
	a1 = 1.0;
	if ( m_p_disperse == 'p'){
		cerr << "a1 = 1.0" << endl;
		do {
			a2 = readStdinDefault(1.4 , "a2 (a2>a1)");
		} while (a2 < a1);
		do{
			number_ratio = readStdinDefault(0.5, "n1/(n1+n2)");
		} while ( number_ratio < 0 || number_ratio > 1);
	}
	rand_seed = readStdinDefault(1, "random seed");
	/*
	 *  Calculate parameters
	 */
	double np1_tmp = np*number_ratio;
	if (np1_tmp - (int)np1_tmp <= 0.5) {
		np1 = (int)np1_tmp;
	} else {
		np1 = (int)np1_tmp+1;
	}
	np2 = np - np1;
	double pvolume1, pvolume2;
	if (dimension == 2) {
		pvolume1 = M_PI*np1;
		pvolume2 = M_PI*a2*a2*np2;
	} else {
		pvolume1 = (4.0/3)*M_PI*np1;
		pvolume2 = (4.0/3)*M_PI*a2*a2*a2*np2;
	}
	double pvolume = pvolume1 + pvolume2;
	if (dimension == 2) {
		lz = sqrt(pvolume / (lx_lz*volume_fraction ));
		lx = lz*lx_lz;
		ly = 0.0;
	} else {
		lz = pow( pvolume / (lx_lz*ly_lz*volume_fraction), 1.0/3);
		lx = lz*lx_lz;
		ly = lz*ly_lz;
	}
	lx2 = lx/2;
	ly2 = ly/2;
	lz2 = lz/2;
	cerr << "np = " << np1+np2 << endl;
	cerr << "np1 : np2 " << np1  << ":" << np2 << endl;
	cerr << "vf = " << volume_fraction << endl;
	cerr << "box =" << lx << ' ' << ly << ' ' << lz << endl;
	
	sys.np(np);
//	if (np2 > 0) {
//		sys.poly = true;
//	} else {
//		sys.poly = false;
//	}
	cerr << "np = " << np << endl;
	if (ly == 0) {
		sys.dimension = 2;
	} else {
		sys.dimension = 3;
	}
	cerr << "dimension = " << sys.dimension << endl;

}

void
GenerateInitConfig::setSystemParameters(){
	double viscosity_solvent = 0.001;
	double radius_of_particle = 1e-6;
	double cf_amp_dl0 = 1;
	double unit_of_length = radius_of_particle; // = radius of smaller particle (a0)
	double unit_of_velocity = sys.shear_rate*unit_of_length;
	double unit_of_force = 6*M_PI*viscosity_solvent*radius_of_particle*unit_of_velocity;
	double cf_amp = cf_amp_dl0*6*M_PI*viscosity_solvent*radius_of_particle*radius_of_particle;
	sys.cf_amp_dl = cf_amp/unit_of_force;
	cerr << "unit_of_force = " << unit_of_force << endl;
	cerr << "cf_amp_dl = " << sys.cf_amp_dl << endl;
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
	sys.kn = 5000;
	sys.overlap_target = 0.03;
	sys.disp_tan_target = 0.03;
	/*
	 * Colloidal force parameter
	 * Short range repulsion is assumed.
	 * cf_amp_dl0: cf_amp_dl at shearrate = 1
	 */
	cf_amp_dl0 = 1; // dimensionless
	sys.cf_range_dl = 0.1; // dimensionless
	/*
	 * mu_static: static friction coeffient
	 * mu_dynamic: dynamic friction coeffient
	 */
	sys.mu_static = 10;
	/*
	 * Output interval
	 */
	sys.strain_interval_output = 0.05;
	/*
	 *  Data output
	 */
	/*
	 * The middle height of the simulation box is set to the flow zero level.
	 */
}


