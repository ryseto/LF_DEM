//
//  GenerateInitConfig.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/6/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "GenerateInitConfig.h"

int GenerateInitConfig::generate(int argc, const char * argv[]){
	setParameters(argc, argv);
	position.resize(np);
	radius.resize(np);
	putRandom();
	solveOverlap();
	outputPositionData();
	return 0;
}

void GenerateInitConfig::outputPositionData(){
	ofstream fout;
	ostringstream ss_posdatafilename;
	ss_posdatafilename << "D" << dimension;
	ss_posdatafilename << "N" << np;
	ss_posdatafilename << "VF" << volume_fraction;
	if (number_ratio == 1){
		ss_posdatafilename << "Mono";
	} else {
		ss_posdatafilename << "Poly" << a2 << "_" << number_ratio ;
	}
	
	if (dimension == 2){
		if (lx_lz == 1){
			ss_posdatafilename << "Square"; // square
		} else {
			ss_posdatafilename << "L" << (int)(10*lx_lz) << "_" << 10;
		}
	} else {
		if (lx_lz == 1 && ly_lz == 1){
			ss_posdatafilename << "Qubic"; //
		} else {
			ss_posdatafilename << "L" << (int)(10*lx_lz) << "_" << (int)(10*ly_lz) << "_" << 10;
		}
	}
	ss_posdatafilename << "_" << random_seed << ".dat";
	cerr << ss_posdatafilename.str() << endl;
	fout.open(ss_posdatafilename.str().c_str());
	fout << "# np1 np2 vf lx ly lz" << endl;
	fout << "# " << np1 << ' ' << np2 << ' ' << volume_fraction << ' ';
	fout << lx << ' ' << ly << ' ' << lz << endl;
	for (int i = 0; i < np ; i++){
		fout << position[i].x << ' ';
		fout << position[i].y << ' ';
		fout << position[i].z << ' ';
		fout << radius[i] << endl;
	}
	fout.close();
}


void GenerateInitConfig::solveOverlap(){
	int cc = 0;
	vector<int> previous_overlap;
	previous_overlap.resize(np);
	double dd = 0.01;
	vec3d delta_translation;
	while (true){
		int i = lrand48() % np;
		if (dimension == 2){
			double rand_angle = 2*M_PI*drand48();
			delta_translation.set(dd*cos(rand_angle), 0, dd*sin(rand_angle));
		} else {
			delta_translation = randUniformSphere(dd);
		}
		position[i] += delta_translation;
		position[i].periodicBoundaryBox(lx, ly, lz);
		int overlap = -1;
		if (is_contact(i, previous_overlap[i] )){
			/* First check is the overlap last time is dissolved?
			 */
			overlap = previous_overlap[i];
		} else {
			/* If previous overlap is resolved.
			 * search other overlap.
			 */
			for (int j = 0; j < np ; j++){
				if (j != i){
					if (is_contact(i, j)){
						previous_overlap[i] = j;
						overlap = j;
						break;
					}
				}
			}
		}
		if (overlap == -1){
			if (cc > 10000){
				if (cc % 1000 == 0 ){
					if ( checkOverlap() == false){
						break;
					}
				}
			}
			cc ++;
		} else {
			position[i] += dd*dr;
			position[i].periodicBoundaryBox(lx, ly, lz);
			position[overlap] -= dd*dr;
			position[overlap].periodicBoundaryBox(lx, ly, lz);
		}
	}
}

void GenerateInitConfig::putRandom(){
	srand48(random_seed);
	for (int i=0; i < np; i++){
		position[i].x = lx*drand48();
		position[i].z = lz*drand48();
		if (dimension == 2){
			position[i].y = ly2;
		} else {
			position[i].y = ly*drand48();
		}
		if (i < np1){
			radius[i] = a1;
		} else {
			radius[i] = a2;
		}
	}
}

bool GenerateInitConfig::checkOverlap(){
	static int i_previous = 0;
	static int j_previous = 1;
	if ( is_contact(i_previous, j_previous ) ){
		return true;
	}
	for (int i = 0; i < np ; i++){
		for (int j = i+1; j < np ; j++){
			if ( is_contact(i,j) ){
				i_previous = i ;
				j_previous = j ;
				return true;
			}
		}
	}
	return false;
}

bool GenerateInitConfig::is_contact(int i, int j){
	double contact_distance = radius[i] + radius[j];
	double sq_contact_distance = contact_distance*contact_distance;
	dr.x = position[i].x - position[j].x;
	dr.z = position[i].z - position[j].z;
	if (dr.z > contact_distance ){
		dr.z -= lz;
	} else if (dr.z < -contact_distance){
		dr.z += lz;
	}
	if (abs(dr.z) < contact_distance){
		if(dr.x > contact_distance){
			dr.x -= lx;
		} else if(dr.x < - contact_distance){
			dr.x += lx;
		}
		if (abs(dr.x) < contact_distance){
			if (dimension == 3){
				dr.y = position[i].y - position[j].y;
				if (dr.y > contact_distance ){
					dr.y -= ly;
				} else if (dr.y < -contact_distance){
					dr.y += ly;
				}
				if (abs(dr.y) < contact_distance){
					if (dr.sq_norm() < sq_contact_distance){
						return true;
					} 
				}
			} else {
				if (dr.sq_norm_xz() < sq_contact_distance){
					return true;
				}
			}
		}
	}
	return false;
}


vec3d GenerateInitConfig::randUniformSphere(double r){
	double z = 2*drand48() - 1.0;
	double phi = 2*M_PI*drand48();
	double sin_theta = sqrt(1.0-z*z);
	return vec3d( r*sin_theta*cos(phi), r*sin_theta*sin(phi), r*z);
}


template<typename T> T readStdinDefault(T default_value,
										string message){
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

void GenerateInitConfig::setParameters(int argc, const char * argv[]){
	/*
	 *  Read parameters from standard input
	 *
	 */
	np = readStdinDefault(200, "number of particle");
	dimension = readStdinDefault(3, "dimension (2 or 3)");
	volume_fraction =  readStdinDefault(0.5, "volume_fraction");
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
	random_seed = readStdinDefault(1, "random seed");
	/*
	 *  Calculate parameters
	 */
	double np1_tmp = np*number_ratio;
	if (np1_tmp - (int)np1_tmp <= 0.5){
		np1 = (int)np1_tmp;
	} else {
		np1 = (int)np1_tmp+1;
	}
	np2 = np - np1;
	double pvolume1, pvolume2;
	if (dimension == 2){
		pvolume1 = M_PI*np1;
		pvolume2 = M_PI*a2*a2*np2;
	} else {
		pvolume1 = (4.0/3)*M_PI*np1;
		pvolume2 = (4.0/3)*M_PI*a2*a2*a2*np2;
	}
	double pvolume = pvolume1 + pvolume2;
	if (dimension == 2){
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
}

