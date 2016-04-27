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
#ifndef USE_DSFMT
#define RANDOM ( rand_gen.rand() ) // RNG uniform [0,1]
#endif
#ifdef USE_DSFMT
#define RANDOM ( dsfmt_genrand_close_open(&rand_gen) ) // RNG uniform [0,1]
#endif
using namespace std;

int GenerateInitConfig::generate(int rand_seed_, int config_type)
{
	if (config_type == 10) {
		cerr << "generate magnetic configuration" <<endl;
		magnetic_config = true;
	} else if (config_type == 2) {
		cerr << "generate circular wide gap " <<endl;
		circulargap_config = true;
	} else if (config_type == 4) {
		cerr << "generate winding walls" <<endl;
		winding_wall_config = true;
	} else if (config_type == 3) {
        cerr << "generate flat walls" <<endl;
        parallel_wall_config = true;
    }
	setParameters();
	rand_seed = rand_seed_;
	if (circulargap_config) {
		/* Note:
		 * Wall particles are placed at
		 *   r = radius_in - a
		 *   r = radius_out + a;
		 * Mobile partilces can be in radius_in < r < radius_out
		 */
        np_wall1 = ((cg_radius_in-radius_wall_particle)*2*M_PI)/(2*radius_wall_particle);
		np_wall2 = ((cg_radius_out+radius_wall_particle)*2*M_PI)/(2*radius_wall_particle);
		cerr << "np_in = " << np_wall1 << endl;
		cerr << "np_out = " << np_wall2 << endl;
		np_movable = np;
		np_fix = np_wall1+np_wall2;
		np += np_fix;
		cerr << "np " << np << endl;
	} else if (parallel_wall_config) {
		/* Note:
		 * Wall particles are placed at
		 *   z = z_bot - a
		 *   z = z_top + a;
		 * Mobile partilces can be in radius_in < r < radius_out
		 */
        np_wall1 = lx/(2*radius_wall_particle);
		np_wall2 = lx/(2*radius_wall_particle);
        np_movable = np;
		np_fix = np_wall1+np_wall2;
		np += np_fix;
	} else if (winding_wall_config) {
		np_wall1 = (cg_radius_out+cg_radius_in)*M_PI/2/1.5+1;
		np_wall2 = (cg_radius_out+cg_radius_in)*M_PI/2/1.5+1;
		np_movable = np;
		np_fix = np_wall1+np_wall2;
		np += np_fix;
	} else {
		np_movable = np;
	}
	sys.set_np(np);
    sys.set_np_mobile(np_movable);
	sys.friction = false;
	sys.repulsiveforce = false;
	sys.p.interaction_range = 2.5;
	sys.p.lub_max_gap = 0.5;
	sys.allocateRessourcesPreConfiguration();

	sys.setBoxSize(lx, ly, lz);
	sys.setSystemVolume();
	sys.in_predictor = false;
	sys.p.integration_method = 0;
	putRandom();
	double inflate_ratio = 1.03;
	bool relax = true;
	for (int i=0; i<np; i++) {
		sys.radius[i] *= inflate_ratio;
	}
	sys.setInteractions_GenerateInitConfig();
	if (relax) {
		grad = new vec3d [np];
		prev_grad = new vec3d [np];
		step_size = 5;
		gradientDescent();
		step_size /= 4.;
		gradientDescent();
		step_size /= 4.;
		gradientDescent();
	}
	// step_size /= 4.;
	// gradientDescent();
	// step_size /= 4.;
	// gradientDescent();
	// step_size /= 4.;
	// gradientDescent();
	// step_size /= 4.;
	// gradientDescent();
	// step_size /= 4.;
	// gradientDescent();
	int count = 0;
	double energy = 0;
	ofstream fout;
	fout.open("energy_decay.dat");
	double energy_previous = 0;
	double diff_energy = 99999;
	int cnt = 0;

	if (relax) {
		do {
			energy = zeroTMonteCarloSweep();
			cerr << energy << endl;
			diff_energy = energy-energy_previous;
			energy_previous = energy;
			if (count%100 == 0) {
				fout << count << ' ' << energy << ' ' <<  diff_energy << endl;
			}
			count ++;
			cerr << "diff = " << abs(diff_energy) << endl;
			if (abs(diff_energy) < 1e-10) {
				cnt++;
			} else {
				cnt = 0;
			}
		} while (cnt < 3);
	}
	//deflate
	for (int i=0; i<np_movable; i++) {
		if (i < np1) {
			sys.radius[i] = a1;
		} else {
			sys.radius[i] = a2;
		}
	}
	for (int i=np_movable; i<np; i++) {
		sys.radius[i] = radius_wall_particle;
	}
//	position.resize(np);
//	radius.resize(np);
//	for (int i=0; i<sys.get_np(); i++) {
//		position[i] = sys.position[i];
//		radius[i] = sys.radius[i];
//	}
	outputPositionData();
	delete [] grad;
	delete [] prev_grad;
	return 0;
}

void GenerateInitConfig::outputPositionData()
{
	ofstream fout;
	ofstream fout_yap;
	ostringstream ss_posdatafilename;
	if (sys.twodimension) {
		ss_posdatafilename << "D2";
	} else {
		ss_posdatafilename << "D3";
	}
	ss_posdatafilename << "N" << np_movable;
	ss_posdatafilename << "VF" << volume_fraction;
	if (disperse_type == 'm') {
		ss_posdatafilename << "Mono";
	} else if (disperse_type == 'b') {
		if (vf_ratio > 0) {
			ss_posdatafilename << "Bidi" << a2 << "_" << vf_ratio;
		} else {
			ss_posdatafilename << "Bidi" << a2 << "_nr" << -vf_ratio;
		}
	} else {
		cerr << "disperse_type is wrong." << endl;
		exit(1);
	}
    if (circulargap_config) {
        ss_posdatafilename << "cylinders" << cg_ratio_radii; // square
    } else if (parallel_wall_config) {
		ss_posdatafilename << "shearwalls"; // square
	} else if (winding_wall_config) {
		ss_posdatafilename << "windingwalls"; // square
	} else {
		if (sys.twodimension) {
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
    }
	vector<double> magnetic_susceptibility;
	if (magnetic_config) {
		for (int i=0; i<np; i++) {
			if (i < np/2) {
				magnetic_susceptibility.push_back(1);
			} else {
				double size_factor = sys.radius[i]*sys.radius[i]*sys.radius[i];
				magnetic_susceptibility.push_back(-1*size_factor);
			}
		}
	}
	
	ss_posdatafilename << "_" << rand_seed << "_.dat";
	cerr << ss_posdatafilename.str() << endl;
	fout.open(ss_posdatafilename.str().c_str());
	ss_posdatafilename << ".yap";
	fout_yap.open(ss_posdatafilename.str().c_str());
    if (circulargap_config || winding_wall_config) {
		fout << "# np1 np2 vf lx ly lz np_wall1 np_wall2 radius_in radius_out" << endl;
		fout << std::setprecision(15);
        fout << "# " << np1 << ' ' << np2 << ' ' << volume_fraction << ' ';
        fout << lx << ' ' << ly << ' ' << lz << ' ';
        fout << np_wall1 << ' ' << np_wall2 << ' ';
        fout << cg_radius_in << ' ' << cg_radius_out << endl;
    } else if (parallel_wall_config) {
		fout << "# np1 np2 vf lx ly lz np_wall1 np_wall2 z_bot z_top" << endl;
		fout << std::setprecision(15);
        fout << "# " << np1 << ' ' << np2 << ' ' << volume_fraction << ' ';
        fout << lx << ' ' << ly << ' ' << lz << ' ';
        fout << np_wall1 << ' ' << np_wall1 << ' ';
        fout << z_bot << ' ' << z_top  << endl;
	} else {
        fout << "# np1 np2 vf lx ly lz vf1 vf2 disp" << endl;
		fout << std::setprecision(15);
        fout << "# " << np1 << ' ' << np2 << ' ' << volume_fraction << ' ';
        fout << lx << ' ' << ly << ' ' << lz << ' ';
        fout << volume_fraction1 << ' ' << volume_fraction2 << ' ' << 0 << endl;
    }
	
	for (int i = 0; i<np; i++) {
		fout << std::setprecision(15);
		fout << sys.position[i].x << ' ';
		fout << sys.position[i].y << ' ';
		fout << sys.position[i].z << ' ';
		fout << sys.radius[i] << ' ';
		if (magnetic_config) {
			fout << 0 << ' ';
			fout << 0 << ' ';
			fout << 0 << ' ';
			fout << magnetic_susceptibility[i];
		}
		fout << endl;
		fout_yap << "r ";
		fout_yap << sys.radius[i] << endl;
		fout_yap << "c ";
		fout_yap << sys.position[i].x << ' ';
		fout_yap << sys.position[i].y << ' ';
		fout_yap << sys.position[i].z << endl;
//		fout_yap << "t ";
//		fout_yap << position[i].x << ' ';
//		fout_yap << position[i].y << ' ';
//		fout_yap << position[i].z << ' ';
//		fout_yap << i << endl;
	}
	fout_yap << "@ 4 \n";
	fout_yap << "y 4 \n";
	for (int k=0; k<sys.nb_interaction; k++) {
		unsigned int i, j;
		sys.interaction[k].get_par_num(i, j);
		vec3d d_pos = sys.position[i]-sys.position[j];
		if (d_pos.norm() < 10){
			fout_yap << "l ";
			fout_yap << sys.position[i].x << ' ';
			fout_yap << sys.position[i].y << ' ';
			fout_yap << sys.position[i].z << ' ';
			fout_yap << sys.position[j].x << ' ';
			fout_yap << sys.position[j].y << ' ';
			fout_yap << sys.position[j].z << endl;
		}
	}
	fout.close();
}

double GenerateInitConfig::computeGradient()
{
	for(int i=0; i<np; i++) {
		grad[i].reset();
	}
	unsigned int i, j;
	double r, rcont;
	double amp, amp2;
	double energy = 0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_contact()) {
			sys.interaction[k].get_par_num(i, j);
			r = sys.interaction[k].r;
            rcont = sys.interaction[k].ro;
			const vec3d& nr_vec = sys.interaction[k].nvec;
			amp = (1/rcont-1/r); // negative
			amp2 = 4*amp/rcont;
			grad[i] -= r*nr_vec*amp2;
			grad[j] += r*nr_vec*amp2;
			energy += 2*r*amp*amp;
		}
	}
	return energy;
}

void GenerateInitConfig::moveAlongGradient(vec3d* g, int dir)
{
	double grad_norm;
	double gradient_power = 0.5;
	vec3d step;
	grad_norm = 0;
	for (int i=0; i<np_movable; i++) {
		grad_norm += g[i].sq_norm();
	}
	if (grad_norm != 0) {
		double rescale = pow(grad_norm, gradient_power);
		for (int i=0; i<np_movable; i++) {
			step = -dir*g[i]*step_size/rescale;
			sys.displacement(i, step);
		}
		sys.checkNewInteraction();
		sys.updateInteractions();
	}
}

void GenerateInitConfig::storeGradient()
{
	for (int i=0; i<np; i++) {
		prev_grad[i] = grad[i];
	}
}

double GenerateInitConfig::gradientDescent()
{
	double old_running_energy;
	double running_energy;
	double relative_en;
	long long int steps = 0;
	cerr << endl << " Gradient Descent..." << endl;
	storeGradient();
	running_energy = computeGradient();
	cerr << "  Starting Energy " << running_energy/np << endl;
	do {
		old_running_energy = running_energy;
		moveAlongGradient(grad, 1);
		storeGradient();
		running_energy = computeGradient();
		relative_en = (old_running_energy-running_energy)/(old_running_energy+running_energy);
		if (steps%100 == 0) {
			cerr << "    Steps = " << steps << " :::   Energy : " << running_energy/np << endl;
		}
		steps++;
		events.clear();
	} while(relative_en > 1e-6);
	if (relative_en < 0) {
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

void GenerateInitConfig::putRandom()
{
    sys.allocatePositionRadius();
#ifndef USE_DSFMT
    rand_gen.seed(rand_seed);
#endif
#ifdef USE_DSFMT
    dsfmt_init_gen_rand(&rand_gen, rand_seed) ; // hash of time and clock trick from MersenneTwister code v1.0 by Richard J. Wagner
#endif
	if (circulargap_config) {
		vec3d r_center(lx/2, 0, lz/2);
        int i = 0;
        while (i < np_movable) {
            vec3d pos(lx*RANDOM, 0, lz*RANDOM);
            double r = (pos-r_center).norm();
            if (r > cg_radius_in+2*radius_wall_particle && r < cg_radius_out-2*radius_wall_particle) {
                sys.position[i] = pos;
                if (i < np1) {
                    sys.radius[i] = a1;
                } else {
                    sys.radius[i] = a2;
                }
                i++;
            }
        }
        for (i=0; i<np_wall1; i++){
            double t = i*(2*M_PI/np_wall1);
            vec3d pos = r_center + (cg_radius_in-radius_wall_particle)*vec3d(cos(t), 0, sin(t));
            sys.position[i+np_movable] = pos;
            sys.radius[i+np_movable] = radius_wall_particle;
        }
        for (i=0; i<np_wall2; i++){
            double t = i*(2*M_PI/np_wall2);
            vec3d pos = r_center + (cg_radius_out+radius_wall_particle)*vec3d(cos(t), 0, sin(t));
            sys.position[i+np_movable+np_wall1] = pos;
            sys.radius[i+np_movable+np_wall1] = radius_wall_particle;
        }
        cerr << np_wall1 << ' ' << np_wall2 << endl;
    } else if (parallel_wall_config) {
        int i = 0;
        while (i < np_movable) {
            vec3d pos(lx*RANDOM, 0, lz*RANDOM);
			double a;
			if (i < np1) {
				a = a1;
			} else {
				a = a2;
			}
            if (pos.z > z_bot+a && pos.z < z_top-a) {
				pos.cerr();
				sys.position[i] = pos;
				sys.radius[i] = a;
				i++;
            }
        }
        double delta_x = lx/np_wall1;
        for (i=0; i<np_wall1; i++){
            vec3d pos(1+delta_x*i, 0, z_bot-radius_wall_particle);
            sys.position[i+np_movable] = pos;
            sys.radius[i+np_movable] = radius_wall_particle;
        }
        for (i=0; i<np_wall2; i++){
            vec3d pos(1+delta_x*i, 0, z_top+radius_wall_particle);
            sys.position[i+np_movable+np_wall1] = pos;
            sys.radius[i+np_movable+np_wall1] = radius_wall_particle;
        }
        cerr << "*" << endl;
        cerr << np_wall1 << ' ' << np_wall2 << endl;
	} else if (winding_wall_config) {
		int i = 0;
		vec3d r_center1(0, 0, 0);
		vec3d r_center2(lx/2, 0, lx/2);
		vec3d r_center3(lx, 0, 0);
		while (i < np_movable) {
			vec3d pos(lx*RANDOM, 0, lz*RANDOM);
			double r;
			if (pos.z > pos.x) {
				r = (pos-r_center1).norm();
			} else if (pos.z < -pos.x + lx) {
				r = (pos-r_center2).norm();
			} else {
				r = (pos-r_center3).norm();
			}
			if (r > cg_radius_in+2*radius_wall_particle && r < cg_radius_out-2*radius_wall_particle) {
				sys.position[i] = pos;
				if (i < np1) {
					sys.radius[i] = a1;
				} else {
					sys.radius[i] = a2;
				}
				i++;
			}
		}

		double dl = (cg_radius_in+cg_radius_out)*M_PI/2/np_wall1;
		double l;
		
		l = 0;
		while (l + dl < cg_radius_in*M_PI/2) {
			double theta = l/cg_radius_in;
			double theta0 = 3*M_PI/4;
			vec3d u_vec(cos(theta0-theta), 0, sin(theta0-theta));
			sys.position[i] = r_center3+cg_radius_in*u_vec;
			if (sys.position[i].x > lx) {
				sys.position[i].x -= lx;
			}
			sys.radius[i] = radius_wall_particle;
			l += dl;
			i++;
		}
		//l -= dl;
		while (l + dl<= (cg_radius_out+cg_radius_in)*M_PI/2) {
			
			double l1 = cg_radius_in*M_PI/2;
			double theta = (l-l1)/cg_radius_out;
			double theta0 = 5*M_PI/4;
			vec3d u_vec(cos(theta0+theta), 0, sin(theta0+theta));
			sys.position[i] = r_center2+cg_radius_out*u_vec;
			if (sys.position[i].x > lx) {
				sys.position[i].x -= lx;
			}
			sys.radius[i] = radius_wall_particle;
			l += dl;
			i++;
		}

		
		l = 0;
		while (l + dl < cg_radius_out*M_PI/2) {

			double theta = l/cg_radius_out;
			double theta0 = 3*M_PI/4;
			vec3d u_vec(cos(theta0-theta), 0, sin(theta0-theta));
			sys.position[i] = r_center3+cg_radius_out*u_vec;
			if (sys.position[i].x > lx) {
				sys.position[i].x -= lx;
			}
			sys.radius[i] = radius_wall_particle;
			l += dl;
			i++;
		}
		//l -= dl;
		while (l + dl <= (cg_radius_out+cg_radius_in)*M_PI/2) {

			double l1 = cg_radius_out*M_PI/2;
			double theta = (l-l1)/cg_radius_in;
			double theta0 = 5*M_PI/4;
			vec3d u_vec(cos(theta0+theta), 0, sin(theta0+theta));
			sys.position[i] = r_center2+cg_radius_in*u_vec;
			if (sys.position[i].x > lx) {
				sys.position[i].x -= lx;
			}
			sys.radius[i] = radius_wall_particle;
			l += dl;
			i++;
		}
		cerr << np_wall1 << " " << np_wall2<< endl;
		cerr << np << endl;
		cerr << sys.position.size() << endl;
	} else {
        for (int i=0; i<np_movable; i++) {
            sys.position[i].x = lx*RANDOM;
            sys.position[i].z = lz*RANDOM;
            if (sys.twodimension) {
                sys.position[i].y = ly_half;
            } else {
                sys.position[i].y = ly*RANDOM;
            }
            if (i < np1) {
				sys.radius[i] = a1;
            } else {
                sys.radius[i] = a2;
            }
        }
	}
}

void GenerateInitConfig::updateInteractions(int i)
{
	vector <Interaction*> inter_list;
	for (auto&& inter : sys.interaction_list[i]){
		inter_list.push_back(inter);
	}

	for (auto& il : inter_list) {
		if (il->is_active()) {
			bool desactivated = false;
			il->updateState(desactivated);
			if (desactivated) {
				sys.deactivated_interaction.push(il->get_label());
			}
		}
	}
}

int GenerateInitConfig::overlapNumber(int i)
{
	int overlaps = 0;
	for (auto&& inter : sys.interaction_list[i]){
		if (inter->is_overlap()) {
			overlaps++;
		}
	}
	return overlaps;
}

double GenerateInitConfig::particleEnergy(int i)
{
	double energy = 0;
	for (auto&& inter : sys.interaction_list[i]){
		if (inter->is_overlap()) {
			double amp = inter->get_a_reduced()*(1/inter->ro-1/inter->r); // negative
			energy += inter->r*amp*amp;
		}
	}
	return energy;
}

double GenerateInitConfig::zeroTMonteCarloSweep()
{
	int steps = 0;
	int init_overlaps = 0;
	double init_energy = 0;
	for(int i=0; i<np_movable; i++) {
	 	init_overlaps += overlapNumber(i);
	}
	for(int i=0; i<np_movable; i++) {
	 	init_energy += particleEnergy(i);
	}
	double dx = 0.04;
	while (steps < np_movable) {
		int moved_part = (int)(RANDOM*np_movable);
		//		int overlap_pre_move = overlapNumber(moved_part);
		double energy_pre_move = particleEnergy(moved_part);
		vec3d trial_move;
		if (sys.twodimension) {
			trial_move = randUniformCircle(dx);
		} else {
			trial_move = randUniformSphere(dx);
		}
		trial_move *= RANDOM;
		sys.displacement(moved_part, trial_move);
		updateInteractions(moved_part);
		//		int overlap_post_move = overlapNumber(moved_part);
		double energy_post_move = particleEnergy(moved_part);
		//		if( overlap_pre_move <= overlap_post_move ){
		if (energy_pre_move < energy_post_move) {
			sys.displacement(moved_part, -trial_move);
			updateInteractions(moved_part);
		}
		//sys.checkNewInteraction();
		steps ++;
		events.clear();
	}
	sys.boxset.update();
	sys.checkNewInteraction();
	sys.updateInteractions();
	int final_overlaps = 0;
	double final_energy = 0;
	for(int i=0; i<np_movable; i++) {
		final_overlaps += overlapNumber(i);
	}
	for(int i=0; i<np_movable; i++) {
	 	final_energy += particleEnergy(i);
	}
	cerr << " MC sweep : init energy " << init_energy/np << " final energy " << final_energy/np;
	cerr << " init overlaps " << init_overlaps << " final overlaps " << final_overlaps << endl;
	return final_energy/np;
}

vec3d GenerateInitConfig::randUniformSphere(double r)
{
	double z = 2*RANDOM-1;
	double phi = 2*M_PI*RANDOM;
	double sin_theta = sqrt(1-z*z);
	return vec3d(r*sin_theta*cos(phi), r*sin_theta*sin(phi), r*z);
}

vec3d GenerateInitConfig::randUniformCircle(double r)
{
	double phi = 2*M_PI*RANDOM;
	return vec3d(r*cos(phi), 0, r*sin(phi));
}

template<typename T> T readStdinDefault(T default_value, string message)
{
	string input;
	T value;
	cerr << message << "[" << default_value << "]: ";
	getline(cin, input);
	if (!input.empty()) {
		istringstream stream( input );
		stream >> value;
	} else {
		value = default_value;
	}
	cerr << value << endl;
	return value;
}

void GenerateInitConfig::setParameters()
{
	/*
	 *  Read parameters from standard input
	 *
	 */
	np = readStdinDefault(500, "number of particle");
    if (circulargap_config || parallel_wall_config) {
        sys.twodimension = true;
    } else {
        int dimension = readStdinDefault(2, "dimension (2 or 3)");
        if (dimension == 2) {
            sys.twodimension = true;
        } else {
			sys.twodimension = false;
		}
    }
	if (sys.twodimension) {
		volume_fraction = readStdinDefault(0.78, "volume_fraction");
	} else {
		volume_fraction = readStdinDefault(0.5, "volume_fraction");
	}
    if (circulargap_config || parallel_wall_config) {
		lx_lz = 1.0;
	} else if (winding_wall_config) {
		lx_lz = 0.5;
    } else {
        lx_lz = readStdinDefault(1.0 , "Lx/Lz [1]: "); // default value needs to be float number.
        if (!sys.twodimension) {
            ly_lz = readStdinDefault(1.0 , "Ly/Lz [1]: "); // default value needs to be float number.
        }
    }
	disperse_type = readStdinDefault('b' , "(m)onodisperse or (b)idisperse");
	a1 = 1;
	a2 = 1;
	volume_fraction1 = volume_fraction; // mono
	if (disperse_type == 'b') {
		cerr << "a1 = 1.0" << endl;
		do {
			a2 = readStdinDefault(1.4, "a2 (a2>a1)");
		} while (a2 < a1);
		do {
			vf_ratio = readStdinDefault(0.5, "volume fraction ratio of smaller particle (If negative, number ratio)");
		} while (vf_ratio < -1 || vf_ratio > 1);
	} else {
		vf_ratio = 1;
	}
	
	//rand_seed = readStdinDefault(1, "random seed");
	/*
	 *  Calculate parameters
	 */
	
	double total_volume;
	double pvolume1, pvolume2;
	if (sys.twodimension) {
		pvolume1 = M_PI*a1*a1;
		pvolume2 = M_PI*a2*a2;
	} else {
		pvolume1 = (4.0/3)*M_PI*a1*a1*a1;
		pvolume2 = (4.0/3)*M_PI*a2*a2*a2;
	}
	
	if (vf_ratio > 0) {
		volume_fraction1 = volume_fraction*vf_ratio;
		if (disperse_type == 'b') {
			volume_fraction2 = volume_fraction-volume_fraction1;
		} else {
			volume_fraction2 = 0;
		}
		if (np > 0) {
			total_volume = np/(volume_fraction1/pvolume1+volume_fraction2/pvolume2);
			double np1_tmp = volume_fraction1*total_volume/pvolume1;
			if (np1_tmp-(int)np1_tmp <= 0.5) {
				np1 = (int)np1_tmp;
			} else {
				np1 = (int)np1_tmp+1;
			}
			np2 = np-np1;
			double pvolume = np1*pvolume1+np2*pvolume2;
			if (sys.twodimension) {
				lz = sqrt(pvolume/(lx_lz*volume_fraction));
				lx = lz*lx_lz;
				ly = 0;
			} else {
				lz = pow(pvolume/(lx_lz*ly_lz*volume_fraction), 1.0/3);
				lx = lz*lx_lz;
				ly = lz*ly_lz;
			}
		} else {
			lz = readStdinDefault(10, "lz");
			lx = lz*lx_lz;
			ly = lz*ly_lz;
			double pvolume1_ = lx*ly*lz*volume_fraction1;
			double pvolume2_ = lx*ly*lz*volume_fraction2;
			np1 = (int)(pvolume1_/pvolume1+0.5);
			np2 = (int)(pvolume2_/pvolume2+0.5);
			np = np1+np2;
		}
	} else {
		double nf_ratio = -vf_ratio;
		np1 = nf_ratio*np;
		np2 = np-np1;
		double pvolume = np1*pvolume1+np2*pvolume2;
		if (sys.twodimension) {
			lz = sqrt(pvolume/(lx_lz*volume_fraction));
			lx = lz*lx_lz;
			ly = 0;
		} else {
			lz = pow(pvolume/(lx_lz*ly_lz*volume_fraction), 1.0/3);
			lx = lz*lx_lz;
			ly = lz*ly_lz;
		}
	}

    if (circulargap_config) {
		double area_particle = np1*pvolume1+np2*pvolume2;
		double area_gap = area_particle/volume_fraction;
		cerr << "area_particle = " << area_particle << endl;
		cerr << "area_gap = " << area_gap << endl;
		cg_ratio_radii = readStdinDefault(0.5, "radius ratio (R_in/R_out)");
		cg_radius_out = sqrt(area_gap/(M_PI*(1-cg_ratio_radii*cg_ratio_radii)));
		cg_radius_in = cg_ratio_radii*cg_radius_out;
		cerr << cg_radius_in << endl;
		cerr << cg_radius_out << endl;
		radius_wall_particle = readStdinDefault(1.0, "wall particle size");
		lz = 2*cg_radius_out+5;
		lx = lz*lx_lz;
		ly = 0;
	} else if (winding_wall_config) {
		double area_particle = np1*pvolume1+np2*pvolume2;
		//double area_gap = area_particle/volume_fraction;
		double rout = readStdinDefault(1, "r_out (>sqrt{2}/2)");
		double width = rout - sqrt(2)/2;
		double a = 1;
		//lz = (-sqrt(2)*a+sqrt(2*a*a + 4*sqrt(2)*width*area_particle/(M_PI*volume_fraction)))/(2*sqrt(2)*width);
		lz = (a+sqrt(a*a+4*width*area_particle/(M_PI*volume_fraction*sqrt(2))))/(2*width);
		lx = 2*lz;
		cg_radius_out = lz*(sqrt(2)/2+width);
		cg_radius_in = lz*(sqrt(2)/2-width);
		lz += 2;
		cerr << "lz = " << lz << endl;
		cerr << "cg_radius_out = " << cg_radius_out << endl;
		cerr << "cg_radius_in = " << cg_radius_in << endl;
		double area1 = M_PI*(cg_radius_out-a)*(cg_radius_out-a)/2;
		double area2 = M_PI*(cg_radius_in+a)*(cg_radius_in+a)/2;
		cerr << " area fraction = " << area_particle/(area1-area2);
    } else if (parallel_wall_config) {
        lz += 10;
		radius_wall_particle = readStdinDefault(1.0, "wall particle size");

        z_bot = 5;
        z_top = lz-5;
    }
	lx_half = lx/2;
	ly_half = ly/2;
	lz_half = lz/2;
	cerr << "np = " << np1+np2 << endl;
	cerr << "np1 : np2 " << np1 << ":" << np2 << endl;
	cerr << "vf1 = " << volume_fraction << endl;
	cerr << "vf2 = " << volume_fraction2 << endl;
	cerr << "box =" << lx << ' ' << ly << ' ' << lz << endl;
	cerr << "radius_wall_particle =" << radius_wall_particle << endl;
}
