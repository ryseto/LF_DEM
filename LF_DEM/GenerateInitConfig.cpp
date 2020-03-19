//
//  GenerateInitConfig.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/6/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include <tuple>
#include "GenerateInitConfig.h"
#include "SystemHelperFunctions.h"
#include "Configuration.h"
#include "ParameterSetFactory.h"

#ifndef USE_DSFMT
#define RANDOM ( rand_gen.rand() ) // RNG uniform [0,1]
#endif
#ifdef USE_DSFMT
#define RANDOM ( dsfmt_genrand_close_open(&rand_gen) ) // RNG uniform [0,1]
#endif
using namespace std;

GenerateInitConfig::GenerateInitConfig():
z_top(-1),
z_bot(-1),
a1(1),
a2(1),
circulargap_config(false),
parallel_wall_config(false),
winding_wall_config(false),
bottom_wall_config(false),
filter_mesh_config(false),
cg_radius_in(-1),
cg_radius_out(-1),
cg_ratio_radii(-1),
np_wall1(0),
np_wall2(0),
np_fix(0),
np_movable(0),
radius_wall_particle(-1),
nb_pin(-1),
symmetry_check(false)
{
	cerr << "GenerateInitConfig" << endl;
}

template<typename T>
void GenerateInitConfig::baseSetup(T &conf, bool is2d, double inflate_ratio)
{
	std::tie(conf.position, conf.radius) = putRandom(is2d);
	for (int i=0; i<np; i++) {
		conf.radius[i] *= inflate_ratio;
	}
	if (is2d) {
		conf.angle.resize(conf.position.size(), 0);
	}
	conf.lx = lx;
	conf.ly = ly;
	conf.lz = lz;
	conf.lees_edwards_disp = {0, 0, 0};
	conf.volume_or_area_fraction = volume_fraction;
}

void GenerateInitConfig::generateBasic(int rand_seed_, double volume_frac, unsigned N, bool bidimensional)
{
	Simulation simu;
	setParametersBasic(simu, volume_frac, N, bidimensional);
	rand_seed = rand_seed_;
	cerr << "rand_seed = " << rand_seed_ << endl;
	auto &sys = simu.getSys();
	sys.imposed_flow = std::make_shared<Geometry::ImposedDeformation>();
	double contact_ratio = 0.05;
	double min_gap = -0.01;
	double inflate_ratio = 1-min_gap;
	struct base_shear_configuration c;
	np_movable = np;
	baseSetup(c, sys.twodimension, inflate_ratio);
	sys.setupConfiguration(c, Parameters::ControlVariable::rate);
	sys.interaction->checkNewInteractions();
	sys.interaction->updateInteractions();
	auto contact_nb = countNumberOfContact(sys);
	int cnt_iteration = 0;
	while(((double)contact_nb.first)/sys.get_np() > contact_ratio && evaluateMinGap(sys) < min_gap) {
		sys.timeEvolution(sys.get_time()+2, -1);
		std::cout << "." << cnt_iteration << std::flush;
		contact_nb = countNumberOfContact(sys);
		std::cerr << (double)contact_nb.first << std::endl;
		if (max_iteration > 0
			&& cnt_iteration++ > max_iteration) {
			std::cout << "max iteration " << std::flush;
			break;
		}
	}
	
	for (int i=0; i<np_movable; i++) {
		if (i < np1) {
			sys.conf->radius[i] = a1;
		} else {
			sys.conf->radius[i] = a2;
		}
	}
	for (int i=np_movable; i<np; i++) {
		sys.conf->radius[i] = radius_wall_particle;
	}
	outputPositionData(sys);
}

int GenerateInitConfig::generate(int rand_seed_, double volume_frac_gen_, double cluster_phi_,
								 int config_type)
{
	if (config_type == 2) {
		cerr << "generate circular wide gap " <<endl;
		circulargap_config = true;
	} else if (config_type == 4) {
		cerr << "generate winding walls" <<endl;
		winding_wall_config = true;
	} else if (config_type == 3) {
		cerr << "generate flat walls" <<endl;
		parallel_wall_config = true;
	} else if (config_type == 5) {
		cerr << "generate a flat bottom" <<endl;
		bottom_wall_config = true;
	} else if (config_type == 6) {
		cerr << "filter mesh at right side" <<endl;
		filter_mesh_config = true;
	}
	Simulation simu;
	setParameters(simu, volume_frac_gen_);
	if (rand_seed_ == -999) {
		rand_seed_ = 137;
		symmetry_check = true;
	}
	rand_seed = rand_seed_;
	cerr << "rand_seed = " << rand_seed_ << endl;
	cluster_phi = cluster_phi_;
	auto &sys = simu.getSys();
	sys.imposed_flow = std::make_shared<Geometry::ImposedDeformation>();
	double contact_ratio = 0.05;
	double min_gap = -0.01;
	double inflate_ratio = 1-min_gap;
	
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
		struct circular_couette_configuration c;
		baseSetup(c, sys.twodimension, inflate_ratio);
		c.np_wall1 = np_wall1;
		c.np_wall2 = np_wall2;
		c.radius_in = cg_radius_in;
		c.radius_out = cg_radius_out;
		cerr << "np " << np << endl;
		sys.setupConfiguration(c, Parameters::ControlVariable::rate);
	} else if (parallel_wall_config || bottom_wall_config) {
		/* Note:
		 * Wall particles are placed at
		 *   z = z_bot - a
		 *   z = z_top + a;
		 * Mobile partilces can be in z_bot < r < z_top
		 */
		np_wall1 = lx/(2*radius_wall_particle);
		np_wall2 = lx/(2*radius_wall_particle);
		np_movable = np;
		np_fix = np_wall1 + np_wall2 + 2*nb_pin;
		np += np_fix;
		struct fixed_velo_configuration c;
		baseSetup(c, sys.twodimension, inflate_ratio);
		c.np_wall1 = np_wall1;
		c.np_wall2 = np_wall2;
		c.z_bot = z_bot;
		c.z_top = z_top;
		sys.setupConfiguration(c, Parameters::ControlVariable::rate);
	} else if (filter_mesh_config) {
		np_wall1 =nb_pin;
		np_wall2 = 0;
		np_movable = np;
		np_fix = nb_pin;
		np += nb_pin;
		struct fixed_velo_configuration c;
		baseSetup(c, sys.twodimension, inflate_ratio);
		c.np_wall1 = nb_pin;
		c.np_wall2 = 0;
		sys.setupConfiguration(c, Parameters::ControlVariable::rate);
	} else if (winding_wall_config) {
		np_wall1 = (cg_radius_out+cg_radius_in)*M_PI/2/1.5+1;
		np_wall2 = (cg_radius_out+cg_radius_in)*M_PI/2/1.5+1;
		np_movable = np;
		np_fix = np_wall1+np_wall2;
		np += np_fix;
		struct circular_couette_configuration c; //@@@
		baseSetup(c, sys.twodimension, inflate_ratio);
		c.np_wall1 = np_wall1;
		c.np_wall2 = np_wall2;
		c.radius_in = cg_radius_in;
		c.radius_out = cg_radius_out;
		sys.setupConfiguration(c, Parameters::ControlVariable::rate);
	} else {
		struct base_shear_configuration c;
		np_movable = np;
		baseSetup(c, sys.twodimension, inflate_ratio);
		sys.setupConfiguration(c, Parameters::ControlVariable::rate);
	}
	
	sys.interaction->checkNewInteractions();
	sys.interaction->updateInteractions();
	auto contact_nb = countNumberOfContact(sys);
	int cnt_iteration = 0;
	while(((double)contact_nb.first)/sys.get_np() > contact_ratio && evaluateMinGap(sys) < min_gap) {
		sys.timeEvolution(sys.get_time()+2, -1);
		std::cout << "." << cnt_iteration << std::flush;
		contact_nb = countNumberOfContact(sys);
		std::cerr << (double)contact_nb.first << std::endl;
		if (max_iteration > 0
			&& cnt_iteration++ > max_iteration) {
			std::cout << "max iteration " << std::flush;
			break;
		}
	}
	
	///@@@@ Symmetry check
	if (symmetry_check) {
		for (int i=0; i<np_movable/2; i++) {
			sys.conf->position[i+np_movable/2].x = sys.conf->position[i].x;
			sys.conf->position[i+np_movable/2].z = lz-sys.conf->position[i].z;
		}
	}
	for (int i=0; i<np_movable; i++) {
		if (i < np1) {
			sys.conf->radius[i] = a1;
		} else {
			sys.conf->radius[i] = a2;
		}
	}
	for (int i=np_movable; i<np; i++) {
		sys.conf->radius[i] = radius_wall_particle;
	}
	if (bottom_wall_config) {
		np -= np_wall2;
		np_wall2 = 0;
	}
	outputPositionData(sys);
	return 0;
}

void GenerateInitConfig::outputPositionData(const System &sys)
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
		ss_posdatafilename << "walls"; // square
	} else if (bottom_wall_config) {
		ss_posdatafilename << "bottom"; // square
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
	} else if (parallel_wall_config || bottom_wall_config) {
		fout << "# np1 np2 vf lx ly lz np_wall1 np_wall2 z_bot z_top" << endl;
		fout << std::setprecision(15);
		fout << "# " << np1 << ' ' << np2 << ' ' << volume_fraction << ' ';
		fout << lx << ' ' << ly << ' ' << lz << ' ';
		fout << np_wall1 << ' ' << np_wall2 << ' ';
		fout << z_bot << ' ' << z_top  << endl;
	} else if (filter_mesh_config) {
		fout << "# np1 np2 vf lx ly lz np_wall1 np_wall2" << endl;
		fout << std::setprecision(15);
		fout << "# " << np1 << ' ' << np2 << ' ' << volume_fraction << ' ';
		fout << lx << ' ' << ly << ' ' << lz << ' ';
		fout << np_wall1 << ' ' << np_wall2 << endl;
	} else {
		fout << "# np1 np2 vf lx ly lz vf1 vf2 dispx dispy" << endl;
		fout << std::setprecision(15);
		fout << "# " << np1 << ' ' << np2 << ' ' << volume_fraction << ' ';
		fout << lx << ' ' << ly << ' ' << lz << ' ';
		fout << volume_fraction1 << ' ' << volume_fraction2 << ' ';
		fout << 0 << ' ' << 0 << endl; 
	}
	
	for (int i = 0; i<np; i++) {
		fout << std::setprecision(15);
		fout << sys.conf->position[i].x << ' ';
		fout << sys.conf->position[i].y << ' ';
		fout << sys.conf->position[i].z << ' ';
		fout << sys.conf->radius[i] << ' ';
		fout << endl;
		fout_yap << "r ";
		fout_yap << sys.conf->radius[i] << endl;
		fout_yap << "c ";
		fout_yap << sys.conf->position[i].x << ' ';
		fout_yap << sys.conf->position[i].y << ' ';
		fout_yap << sys.conf->position[i].z << endl;
		//		fout_yap << "t ";
		//		fout_yap << position[i].x << ' ';
		//		fout_yap << position[i].y << ' ';
		//		fout_yap << position[i].z << ' ';
		//		fout_yap << i << endl;
	}
	fout_yap << "@ 4 \n";
	fout_yap << "y 4 \n";
	for (const auto &inter: *(sys.interaction)) {
		unsigned int i, j;
		std::tie(i, j) = inter->get_par_num();
		vec3d d_pos = sys.conf->position[i]-sys.conf->position[j];
		if (d_pos.norm() < 10){
			fout_yap << "l ";
			fout_yap << sys.conf->position[i].x << ' ';
			fout_yap << sys.conf->position[i].y << ' ';
			fout_yap << sys.conf->position[i].z << ' ';
			fout_yap << sys.conf->position[j].x << ' ';
			fout_yap << sys.conf->position[j].y << ' ';
			fout_yap << sys.conf->position[j].z << endl;
		}
	}
	fout.close();
}

std::pair<std::vector<vec3d>, std::vector<double>> GenerateInitConfig::putRandom(bool twodimension)
{
	std::vector<vec3d> position (np);
	std::vector<double> radius (np);
	
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
				position[i] = pos;
				if (i < np1) {
					radius[i] = a1;
				} else {
					radius[i] = a2;
				}
				i++;
			}
		}
		for (i=0; i<np_wall1; i++){
			double t = i*(2*M_PI/np_wall1);
			vec3d pos = r_center+(cg_radius_in-radius_wall_particle)*vec3d(cos(t), 0, sin(t));
			position[i+np_movable] = pos;
			radius[i+np_movable] = radius_wall_particle;
		}
		for (i=0; i<np_wall2; i++){
			double t = i*(2*M_PI/np_wall2);
			vec3d pos = r_center + (cg_radius_out+radius_wall_particle)*vec3d(cos(t), 0, sin(t));
			position[i+np_movable+np_wall1] = pos;
			radius[i+np_movable+np_wall1] = radius_wall_particle;
		}
		cerr << np_wall1 << ' ' << np_wall2 << endl;
	} else if (parallel_wall_config || bottom_wall_config) {
		int i = 0;
		double shift_up = 2;
		double cluster_radius;
		if (bottom_wall_config) {
			shift_up = 5;
		}
		int num_cl = 0;
		if (cluster_phi != 0) {
			if (cluster_phi > 0) {
				num_cl = 1;
			} else {
				num_cl = 2;
			}
			cluster_radius = sqrt((np1*a1*a1 + np2*a2*a2)/(num_cl* abs(cluster_phi)));
			cerr << "cluster_radius = " << cluster_radius << endl;
		}
		while (i < np_movable) {
			double a = (i < np1 ? a1 : a2);
			vec3d pos = {lx*RANDOM, 0, lz*RANDOM};
			if (num_cl == 0) {
				if (pos.z > z_bot+shift_up && pos.z < z_top-shift_up) {
					position[i] = pos;
					radius[i] = a;
					i++;
				}
			} else {
				if (num_cl ==  1) {
					vec3d pos_center(lx_half, 0, lz_half);
					if ( (pos_center - pos).norm()< cluster_radius) {
						position[i] = pos;
						radius[i] = a;
						i++;
					}
				} else if (num_cl == 2) {
					double z1 = (lz-cluster_radius)/3;
					double z2 = (2*lz+cluster_radius)/3;
					double z_center = (i < np_movable/2 ? z1 : z2);
					vec3d pos_center(lx_half, 0, z_center);
					if ( (pos_center - pos).norm()< cluster_radius) {
						position[i] = pos;
						radius[i] = a;
						i++;
					}
				}
			}
		}
		double delta_x = lx/np_wall1;
		for (i=0; i<np_wall1; i++){
			double x = radius_wall_particle+delta_x*i;
			double z = z_bot;
			vec3d pos(x, 0, z);
			position[i+np_movable] = pos;
			radius[i+np_movable] = radius_wall_particle;
		}
		for (i=0; i<np_wall1; i++){
			double x = radius_wall_particle+delta_x*i;
			double z = z_top;
			vec3d pos(x, 0, z);
			position[i+np_movable+np_wall1] = pos;
			radius[i+np_movable+np_wall1] = radius_wall_particle;
		}
		delta_x = lx/nb_pin;
		for (i=0; i<nb_pin; i++){
			double x = delta_x*i + delta_x/2;
			double z = z_bot+2*radius_wall_particle;
			vec3d pos(x, 0, z);
			position[2*i+np_movable+np_wall1] = pos;
			radius[2*i+np_movable+np_wall1] = radius_wall_particle;
			z = z_top-2*radius_wall_particle;
			pos.set(x, 0, z);
			position[2*i+1+np_movable+np_wall1] = pos;
			radius[2*i+1+np_movable+np_wall1] = radius_wall_particle;
		}
//		for (i=0; i<np_wall2; i++){
//			vec3d pos(radius_wall_particle+delta_x*i, 0, z_top);
//			if (wall_pin_interval > 0 && i%wall_pin_interval == wall_pin_interval-1) {
//				pos -= del;
//			}
//			position[i+np_movable+np_wall1] = pos;
//			radius[i+np_movable+np_wall1] = radius_wall_particle;
//		}
		cerr << np_wall1 << ' ' << np_wall2 << endl;
		
	} else if (filter_mesh_config) {
		int i = 0;
		while (i < np_movable) {
			double a;
			if (i < np1) {
				a = a1;
			} else {
				a = a2;
			}
			vec3d pos(lx*RANDOM, 0, lz*RANDOM);
			position[i] = pos;
			radius[i] = a;
			i++;
		}
		double delta_z = lz/nb_pin;
		for (i=0; i<np_wall1; i++) {
			double z = delta_z*i + delta_z/2;
			double x = lx-radius_wall_particle;
			vec3d pos(x, 0, z);
			position[i+np_movable] = pos;
			radius[i+np_movable] = radius_wall_particle;
		}
		//		for (i=0; i<np_wall2; i++){
		//			vec3d pos(radius_wall_particle+delta_x*i, 0, z_top);
		//			if (wall_pin_interval > 0
		//				pos -= del;
		//			}
		//			position[i+np_movable+np_wall1] = pos;
		//			radius[i+np_movable+np_wall1] = radius_wall_particle;
		//		}
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
				position[i] = pos;
				if (i < np1) {
					radius[i] = a1;
				} else {
					radius[i] = a2;
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
			position[i] = r_center3+cg_radius_in*u_vec;
			if (position[i].x > lx) {
				position[i].x -= lx;
			}
			radius[i] = radius_wall_particle;
			l += dl;
			i++;
		}
		//l -= dl;
		while (l + dl<= (cg_radius_out+cg_radius_in)*M_PI/2) {
			
			double l1 = cg_radius_in*M_PI/2;
			double theta = (l-l1)/cg_radius_out;
			double theta0 = 5*M_PI/4;
			vec3d u_vec(cos(theta0+theta), 0, sin(theta0+theta));
			position[i] = r_center2+cg_radius_out*u_vec;
			if (position[i].x > lx) {
				position[i].x -= lx;
			}
			radius[i] = radius_wall_particle;
			l += dl;
			i++;
		}
		
		l = 0;
		while (l + dl < cg_radius_out*M_PI/2) {
			double theta = l/cg_radius_out;
			double theta0 = 3*M_PI/4;
			vec3d u_vec(cos(theta0-theta), 0, sin(theta0-theta));
			position[i] = r_center3+cg_radius_out*u_vec;
			if (position[i].x > lx) {
				position[i].x -= lx;
			}
			radius[i] = radius_wall_particle;
			l += dl;
			i++;
		}
		//l -= dl;
		while (l + dl <= (cg_radius_out+cg_radius_in)*M_PI/2) {
			double l1 = cg_radius_out*M_PI/2;
			double theta = (l-l1)/cg_radius_in;
			double theta0 = 5*M_PI/4;
			vec3d u_vec(cos(theta0+theta), 0, sin(theta0+theta));
			position[i] = r_center2+cg_radius_in*u_vec;
			if (position[i].x > lx) {
				position[i].x -= lx;
			}
			radius[i] = radius_wall_particle;
			l += dl;
			i++;
		}
		cerr << np_wall1 << " " << np_wall2<< endl;
		cerr << np << endl;
		cerr << position.size() << endl;
	} else {
		for (int i=0; i<np_movable; i++) {
			position[i].x = lx*RANDOM;
			position[i].z = lz*RANDOM;
			if (twodimension) {
				position[i].y = ly_half;
			} else {
				position[i].y = ly*RANDOM;
			}
			if (i < np1) {
				radius[i] = a1;
			} else {
				radius[i] = a2;
			}
		}
	}
	return make_pair(position, radius);
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

void GenerateInitConfig::setParameters(Simulation &simu, double volume_frac_init)
{
	/*
	 *  Read parameters from standard input
	 *
	 */
	Parameters::ParameterSetFactory PFactory(Dimensional::Unit::hydro);
	simu.sys.p = PFactory.getParameterSet();
	
	auto &sys = simu.getSys();
	simu.sys.p.integration_method = 0;
	simu.sys.p.disp_max = 5e-3;
	simu.sys.p.lub.model = "none";
	simu.sys.p.contact.friction_model = Interactions::FrictionModel::frictionless;
	simu.sys.p.contact.kn = 1;
	simu.sys.p.contact.relaxation_time_tan = 1e-4;
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
	
	if (volume_frac_init == 0) {
		if (sys.twodimension) {
			volume_fraction = readStdinDefault(0.78, "volume_fraction");
		} else {
			volume_fraction = readStdinDefault(0.5, "volume_fraction");
		}
	} else {
		cerr << "volume_fraction is set to " << volume_frac_init << endl;
		volume_fraction = volume_frac_init;
	}
	if (circulargap_config) {
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
	} else if (parallel_wall_config || bottom_wall_config ) {
		radius_wall_particle = readStdinDefault(1.0, "wall particle size");
		nb_pin = readStdinDefault(0, "number of wall pin");
		double delta = 1e-6;
		z_bot = 0 + delta;
		z_top = lz - delta;
	} else if (filter_mesh_config) {
		radius_wall_particle = readStdinDefault(1.0, "wall particle size");
		nb_pin = readStdinDefault(0, "number of wall pin");
	}
	lx_half = lx/2;
	ly_half = ly/2;
	lz_half = lz/2;
	
	max_iteration = readStdinDefault(-1, "max iteration number");
	
	cerr << "np = " << np1+np2 << endl;
	cerr << "np1 : np2 " << np1 << ":" << np2 << endl;
	cerr << "vf1 = " << volume_fraction << endl;
	cerr << "vf2 = " << volume_fraction2 << endl;
	cerr << "box =" << lx << ' ' << ly << ' ' << lz << endl;
	cerr << "radius_wall_particle =" << radius_wall_particle << endl;
}


void GenerateInitConfig::setParametersBasic(Simulation &simu, double volume_frac, unsigned N, bool bidimensional)
{
	/*
	 *  Read parameters from standard input
	 *
	 */
	Parameters::ParameterSetFactory PFactory(Dimensional::Unit::hydro);
	simu.sys.p = PFactory.getParameterSet();
	
	auto &sys = simu.getSys();
	simu.sys.p.integration_method = 0;
	simu.sys.p.disp_max = 5e-3;
	simu.sys.p.lub.model = "none";
	simu.sys.p.contact.friction_model = Interactions::FrictionModel::frictionless;
	simu.sys.p.contact.kn = 1;
	simu.sys.p.contact.relaxation_time_tan = 1e-4;
	np = N;
	sys.twodimension = bidimensional;
	
	volume_fraction = volume_frac;

	lx_lz = 1;
	ly_lz = 1;

	disperse_type = 'b';
	a1 = 1;
	a2 = 1.4;
	vf_ratio = 0.5;
	volume_fraction1 = volume_fraction; // mono
	
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
	
	volume_fraction1 = volume_fraction*vf_ratio;
	volume_fraction2 = volume_fraction-volume_fraction1;
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
	}
	lx_half = lx/2;
	ly_half = ly/2;
	lz_half = lz/2;
	
	max_iteration = -1;
	
	cerr << "np = " << np1+np2 << endl;
	cerr << "np1 : np2 " << np1 << ":" << np2 << endl;
	cerr << "vf1 = " << volume_fraction << endl;
	cerr << "vf2 = " << volume_fraction2 << endl;
	cerr << "box =" << lx << ' ' << ly << ' ' << lz << endl;
}
