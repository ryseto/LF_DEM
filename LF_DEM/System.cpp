//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "System.h"
#include <sstream>
#include <Accelerate/Accelerate.h>


System::~System(){
	delete [] position;
	delete [] angle;
	delete [] velocity;
	delete [] ang_velocity;
	delete [] force;
	delete [] torque;
	delete [] work;
	delete [] ipiv;
};


void System::init(){
	lx2 = 0.5*lx;
	ly2 = 0.5*ly;
	lz2 = 0.5*lz;
	ostringstream ss_simu_name;
	ss_simu_name << "D" << dimension << "L" << lz << "vf" << volume_fraction <<  "ms" << mu_static << "md" << mu_dynamic << "lub" << lubcore ;
	
	simu_name = ss_simu_name.str();
	cerr << simu_name << endl;
	sq_critical_velocity = dynamic_friction_critical_velocity * dynamic_friction_critical_velocity;
	
	vel_difference = shear_rate*lz;
}

/* Set number of particles.
 * Allocate vectors for the state.
 */
void System::setNumberParticle(int num_particle){
	n = num_particle;
	position = new vec3d [num_particle];
	n3 = 3*n;
	i_position = new int * [num_particle];
	for (int i = 0; i<num_particle; i++ ){
		i_position[i] = new int [3];
	}
	angle = new double [num_particle];
	velocity = new vec3d [num_particle];
	ang_velocity = new vec3d [num_particle];
	force = new vec3d [num_particle];
	torque = new vec3d [num_particle];
	res = new double [9*num_particle*num_particle];
	mov = new double [9*num_particle*num_particle];
	
	/* for dgesv_ or dsysv_
	 */
	b_vector = new double [n3];
	x_vector = new double [n3];
	lwork =n3*4;
	work = new double [n3*4];
	ipiv= new int [n3];
	UPLO = 'U';
	nrhs= 1;

	lda = n3;
	ldb = n3;
}

void System::forceReset(){
	for (int i=0; i < n; i++){
		force[i].reset();
	}
}
void System::torqueReset(){
	if (friction){
		for (int i=0; i < n; i++){
			torque[i].reset();
		}
	}
}

inline vec3d randUniformSphere(double r){
	double z = 2*drand48() - 1.0;
	double phi = 2*M_PI*drand48();
	double sin_theta = sqrt(1.0-z*z);
	return vec3d( r*sin_theta*cos(phi),  r*sin_theta*sin(phi), r*z);
}

bool System::nooverlap(){
	bool no_overlap = true;
	static int i_previous = 1;
	static int j_previous = 5;
	if ( checkContact(i_previous, j_previous) < 4){
		cerr << "." ;
		return false;
	}
		
	for (int i = 0; i < n ; i++){
		for (int j = i+1; j < n ; j++){
			if ( checkContact(i, j) < 4){
				no_overlap = false;
				i_previous = i ;
				j_previous = j ;
				
				break;
			}
		}
	}
	return no_overlap;
	
}


/*
 * To prepare an initial configuration.
 */
void System::setRandomPosition(){
	vec3d trial_pos;
	for (int i=0; i < n; i++){
		trial_pos.x = lx*drand48();
		trial_pos.z = lz*drand48();
		if (dimension == 2){
			trial_pos.y = ly2;
		} else {
			trial_pos.y = ly*drand48();
		}
		position[i] = trial_pos;
		angle[i] = 0;
	}

	int cc = 0;
	vector<int> previous_overlap;
	previous_overlap.resize(n);
	double dd=0.01;
	while (true){
		int i = lrand48() % n;
		if (dimension == 2){
			double rand_angle = 2 * M_PI * drand48();
			vec3d shift(dd*cos(rand_angle), 0, dd*sin(rand_angle));
			displacement(i, shift.x ,0,shift.z );
		} else {
			vec3d shift = randUniformSphere(dd);
			displacement(i, shift.x ,shift.y ,shift.z );
		}
		int overlap = -1;
		
		if (checkContact(i, previous_overlap[i]) < 4){
			overlap =  previous_overlap[i];
		} else {
			for (int j = 0; j < n ; j++){
				if (j != i){
					if ( checkContact(i, j) < 4){
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
					if (nooverlap()){
						break;
					}
				}
			}
			cc ++;
		} else {
			displacement(i, dd*dx, dd*dy, dd*dz);
			displacement(overlap, -dd*dx, -dd*dy, -dd*dz);
		}
	}
}


void System::updateVelocity(){
	vec3d U_inf(0, 0, 0);
	for (int i=0; i < n; i++){
		U_inf.x = shear_rate*position[i].z;
		velocity[i] = (1.0/eta)*force[i] + U_inf;
	}
	if(friction){
		vec3d O_inf(0, 0.5*shear_rate, 0);
		for (int i=0; i < n; i++){
			ang_velocity[i] = (1.33333/eta)*torque[i] + O_inf;
		}
	}
}

void resmatrix(double *res, double *nvec, int ii, int jj, double alpha, int n3){
	int ii3 = 3*ii;
	int jj3 = 3*jj;
	int jj3_1 = jj3+1;
	int jj3_2 = jj3+2;
	double alpha_n1n0 = alpha*nvec[1]*nvec[0];
	double alpha_n2n1 = alpha*nvec[2]*nvec[1];
	double alpha_n0n2 = alpha*nvec[0]*nvec[2];

	res[ n3*ii3   + jj3   ]   += alpha*nvec[0]*nvec[0]; // 00
	res[ n3*ii3   + jj3_1 ]   += alpha_n1n0; // 10
	res[ n3*ii3   + jj3_2 ]   += alpha_n0n2; // 20
	
	res[ n3*(ii3+1) + jj3   ] += alpha_n1n0; // 01
	res[ n3*(ii3+1) + jj3_1 ] += alpha*nvec[1]*nvec[1]; //11
	res[ n3*(ii3+1) + jj3_2 ] += alpha_n2n1; // 21
	res[ n3*(ii3+2) + jj3   ] += alpha_n0n2; // 02
	res[ n3*(ii3+2) + jj3_1 ] += alpha_n2n1; // 12
	res[ n3*(ii3+2) + jj3_2 ] += alpha*nvec[2]*nvec[2]; //22
}


double System::lubricationForceFactor(int i, int j){
	double r_sq = sq_distance(i,j);
	if( r_sq < sq_lub_max){
		double r = sqrt(r_sq);
		double h = r - lubcore;
		vec3d nv(dx/r, dy/r, dz/r);
		vec3d rel_vel = velocity[j] - velocity[i];
		if (position[i].z -position[j].z > lz2){
			rel_vel.x += vel_difference;
		} else if (position[i].z -position[j].z < -lz2){
			rel_vel.x -= vel_difference;
		}
		double alpha = 1.0/(4*h);
		double force = abs(alpha*dot( rel_vel , nv));
		return force;
	} else {
		return 0;
	}
}


void System::updateVelocityLubrication(){
	double tmp[3];
	double nvec[3];
	for (int k = 0; k < n3*n3; k++){
		res[k]=0;
	}
	for (int k = 0;k < n3; k++){
		b_vector[k] = 0;
	}
	
	for (int i = 0 ; i < n; i ++){
		int i3 = 3*i;
		res[n3*(i3  ) + i3  ] = 1;
		res[n3*(i3+1) + i3+1] = 1;
		res[n3*(i3+2) + i3+2] = 1;
	}

	double r, r_sq;
	if (lub){
		for (int i = 0 ; i < n - 1; i ++){
			for (int j = i+1 ; j < n; j ++){
				r_sq = sq_distance(i,j);
				if( r_sq < sq_lub_max){
					r = sqrt(r_sq);
					nvec[0] = dx/r;
					nvec[1] = dy/r;
					nvec[2] = dz/r;
					double h = r - lubcore;
					double alpha = - 1/(4*h);
					if ( h > 0){
						// (i, j) (k,l) --> res[ n3*(3*i+l) + 3*j+k ]
						resmatrix(res, nvec, i, i, -alpha, n3);
						resmatrix(res, nvec, i, j, +alpha, n3);
						resmatrix(res, nvec, j, j, -alpha, n3);
						resmatrix(res, nvec, j, i, +alpha, n3);
						double tmp1 = alpha*shear_rate*dz*nvec[0];
						tmp[0] = tmp1*nvec[0];
						tmp[1] = tmp1*nvec[1];
						tmp[2] = tmp1*nvec[2];
						b_vector[3*i]   += tmp[0];
						b_vector[3*i+1] += tmp[1];
						b_vector[3*i+2] += tmp[2];
						b_vector[3*j]   -= tmp[0];
						b_vector[3*j+1] -= tmp[1];
						b_vector[3*j+2] -= tmp[2];
					}
				}
			}
		}
	}
	/* F = R (V - V_inf)
	 *
	 * (V - V_inf) = M F
	 * A.x = b_vector
	 * b_vector[n] : r-h-s vector
	 * atimes (int n, static double *x, double *b, void *param) :
	 *        calc matrix-vector product A.x = b.
	 */
	for (int i = 0; i < n; i++){
		int i3 = 3*i;
		b_vector[i3] += force[i].x;
		b_vector[i3+1] += force[i].y;
		b_vector[i3+2] += force[i].z;
	}
	// LU
	dgesv_(&n3, &nrhs, res, &lda, ipiv, b_vector, &ldb, &info);
	//	dsysv_(&UPLO, &n3, &nrhs, res, &lda, ipiv, b_vector, &ldb, work, &lwork, &info);
	for (int i = 0; i < n; i++){
		int i3 = 3*i;
		velocity[i].x = b_vector[i3] + shear_rate*position[i].z;
		velocity[i].y = b_vector[i3+1];
		velocity[i].z = b_vector[i3+2];
	}
	if(friction){
		double delta_omega =  0.5*shear_rate;
		for (int i=0; i < n; i++){
			ang_velocity[i] = (1.33333/eta)*torque[i];
			ang_velocity[i].y += delta_omega;
		}
	}
}


void System::displacement(int i, const double &dx_, const double &dy_, const double &dz_){
	position[i].x += dx_;
	position[i].y += dy_;
	position[i].z += dz_;
	
	if (position[i].z > lz ){
		position[i].z -= lz;
		position[i].x -= x_shift;
	} else if ( position[i].z < 0 ){
		position[i].z += lz;
		position[i].x += x_shift;
	}
	if ( position[i].x > lx ){
		position[i].x -= lx;
	} else if (position[i].x < 0 ){
		position[i].x += lx;
	}
	if (dimension == 3){
		if ( position[i].y > ly ){
			position[i].y -= ly;
		} else if (position[i].y < 0 ){
			position[i].y += ly;
		}
	}
	
	
}


void System::deltaTimeEvolution(){
	x_shift += vel_difference*dt;
	if (x_shift > lx){
		x_shift -= lx;
	}
	for (int i=0; i < n; i++){
		
		displacement(i, velocity[i].x*dt,velocity[i].y*dt,velocity[i].z*dt);
	}

	if (friction){
		if (dimension == 2){
			for (int i=0; i < n; i++){
				angle[i] += ang_velocity[i].y*dt;
			}
		}
	}
	
}



/* Distance between particle i and particle j
 */
double System::distance(int i, int j){
	return sqrt(sq_distance(i,j));
}

/* Square norm of vector (dx, dy, dz)
 */
double System::sq_norm(){
	if (dz > lz2 ){
		dz -= lz;
		dx -= x_shift;
	} else if (dz < -lz2){
		dz += lz;
		dx += x_shift;
	}
	
	while(dx > lx2){
		dx -= lx;
	}
	while(dx < - lx2){
		dx += lx;
	}
	if (dimension == 3){
		if (dy > ly2 ){
			dy -= ly;
		} else if (dy < -ly2){
			dy += ly;
		}
		return dx*dx + dy*dy + dz*dz;
	} else {
		return dx*dx + dz*dz;
		
	}
}

/* Square distance between particle i and particle j
 */
double System::sq_distance(int i, int j){
	dx = position[i].x - position[j].x;
	dy = position[i].y - position[j].y;
	dz = position[i].z - position[j].z;
	return sq_norm();
}

/* 
 * Square distance between position pos and particle i
 */
double System::sq_distance(vec3d &pos , int i){
	dx = position[i].x - pos.x;
	dy = position[i].y - pos.y;
	dz = position[i].z - pos.z;
	return sq_norm();
}

double System::checkContact(int i, int j){
	dx = position[i].x - position[j].x;
	dz = position[i].z - position[j].z;
	if (dz > lz2 ){
		dz -= lz;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
		dx -= x_shift;
	} else if (dz < -lz2){
		dz += lz;
		dx += x_shift;
	}
	if (abs(dz) < 2){
		while(dx > lx2){
			dx -= lx;
		}
		while(dx < - lx2){
			dx += lx;
		}
		if (abs(dx) < 2){
			if (dimension == 3){
				dy = position[i].y - position[j].y;
				if (dy > ly2 ){
					dy -= ly;
				} else if (dy < -ly2){
					dy += ly;
				}
				if (abs(dy) < 2){
					return dx*dx + dy*dy + dz*dz;
				}
			} else {
				return dx*dx + dz*dz;
			}
		}
	}
	return 100;
}

