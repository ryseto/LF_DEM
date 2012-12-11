//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
#include <sstream>

System::~System(){
	delete [] position;
	delete [] angle;
	delete [] velocity;
	delete [] ang_velocity;
	delete [] force;
	delete [] torque;
#ifdef CHOLMOD
	delete [] diag_values;
	delete [] off_diag_values;
	delete [] ploc;
	cholmod_free_dense(&v, &c);
	cholmod_free_dense(&rhs_b, &c);
	cholmod_free_sparse(&sparse_res, &c);
	cholmod_finish(&c);
#else
	delete [] work;
	delete [] ipiv;
#endif
};

void System::prepareSimulationName(){
	ostringstream ss_simu_name;
	if (dimension == 2){
		ss_simu_name << "D" << dimension << "L" << lx << "_" <<lz ;
	} else {
		ss_simu_name << "D" << dimension << "L" << lx << "_" << ly << "_" << lz ;
	}
	if (friction == true){
		ss_simu_name << "vf" << volume_fraction ;
		ss_simu_name << "fs" << mu_static << "fd" << mu_dynamic;
	} else {
		ss_simu_name << "vf" << volume_fraction ;
	}
	if (lub == true){
		ss_simu_name << "lub" << lubcore ;
	}
	simu_name = ss_simu_name.str();
	cerr << simu_name << endl;
	
}

/* Set number of particles.
 * Allocate vectors for the state.
 */
void System::prepareSimulation(){
	ts = 0;
	lx2 = 0.5*lx;
	ly2 = 0.5*ly;
	lz2 = 0.5*lz;
	shear_disp = 0;
	vel_difference = shear_rate*lz;
	sq_critical_velocity = dynamic_friction_critical_velocity * dynamic_friction_critical_velocity;
	position = new vec3d [n];
	n3 = 3*n;
	angle = new double [n];
	velocity = new vec3d [n];
	ang_velocity = new vec3d [n];
	force = new vec3d [n];
	torque = new vec3d [n];
	stress = new double* [n];
	for (int i=0; i < n; i++){
		stress[i] = new double [5];
	}
	
	double O_inf_y = 0.5*shear_rate/2.0;
	for (int i=0; i < n; i++){
		ang_velocity[i].set(0, O_inf_y, 0);
		torque[i].reset();
	}
	
	fb = new BrownianForce(this);
	maxnum_interactionpair = (int)(12*n);
	
	num_interaction = 0;
	interaction = new Interaction [maxnum_interactionpair];
	for (int k=0; k < maxnum_interactionpair ; k++){
		interaction[k].init(this);
	}
	
	initInteractionPair();
	
#ifdef CHOLMOD
	cholmod_start (&c) ;
	stype = -1; // 1 is symmetric, stored upper triangular (UT), -1 is LT
	sorted = 0;		/* TRUE if columns sorted, FALSE otherwise*/
	packed = 1;		/* TRUE if matrix packed, FALSE otherwise */
	xtype = CHOLMOD_REAL;
	diag_values = new double [6*n];
	off_diag_values = new vector <double> [3];
	ploc = new int [n+1];
	fb->init();
#else
	/* for dgesv_ or dsysv_
	 */
	res = new double [9*n*n];
	rhs_b = new double [n3];
	lwork =n3*4;
	work = new double [n3*4];
	ipiv= new int [n3];
	UPLO = 'U';
	nrhs= 1;
	lda = n3;
	ldb = n3;
#endif
}

void System::initInteractionPair(){
	interaction_pair = new int * [n];
	for (int i=0; i < n; i++){
		interaction_pair[i] = new int [n];
	}
	for (int i=0; i < n-1; i++){
		for (int j=i+1; j < n; j++){
			interaction_pair[i][j] = -1;
		}
	}
}

void System::timeEvolution(int time_step){
	int ts_next = ts + time_step;
	while (ts < ts_next){
		checkNewInteraction();
		for (int k = 0; k < num_interaction; k++){
			
			interaction[k].normalElement();
		}
		forceReset();
		torqueReset();
		calcContactForces();
		if (lub){
			// Lubrication dynamics
			updateVelocityLubrication();
		} else {
			// Free-draining approximation
			updateVelocity();
		}
		
		deltaTimeEvolution();
		if (friction){
			updateContactForceConfig();
		}
		//		sys.checkInteractionEnd();
		
		ts ++;
	}
}

/* Check the distance between separating particles.
 * i < j
 *
 * A patch-up prescription to aboid
 * contact_pair[i][j] < 0 indicates separating particles to be checked.
 * contact_pair[i][j] = -1, the particles are near contact. So every time step, distance should be checked.a
 * contact_pair[i][j] < -1, the particles have some distance.
 */
void System::checkNewInteraction(){
	//	for (int k = 0; k < num_interaction; k++){
	//		int i = interaction[k].particle_num[0];
	//		int j = interaction[k].particle_num[1];
	//		double sq_distance = sq_distanceToCheckContact(i, j);
	//		interaction[k].r = sqrt(sq_distance);
	//		//		cerr << interaction[k].r << endl;
	//	}
	for (int k = 0; k < num_interaction; k++){
		if ( interaction[k].active ){
			if( interaction[k].r > lub_max){
				interaction[k].active = false;
				interaction_pair[interaction[k].particle_num[0]][interaction[k].particle_num[1]] = -1;
				deactivated_interaction.push(k);
			}else if( interaction[k].contact == true ){
				if ( interaction[k].r > 2){
					interaction[k].contact = false;
				}
			} else {
				if ( interaction[k].r < 2){
					// activate new contact
					interaction[k].newContact();
				}
			}
		}
	}
	for (int i=0; i < n-1; i++){
		for (int j=i+1; j < n; j++){
			if ( interaction_pair[i][j] == -1){
				double sq_distance = sq_distanceToCheckContact(i, j);
				if ( sq_distance < sq_lub_max){
					int interaction_new;
					if (deactivated_interaction.empty()){
						// add an interaction object.
						interaction_new = num_interaction;
						num_interaction ++;
					} else {
						// fill a deactivated interaction object.
						interaction_new = deactivated_interaction.front();
						deactivated_interaction.pop();
					}
					interaction[interaction_new].create(i, j);
					interaction_pair[i][j] = interaction_new;
				}
			}
		}
	}
}

void System::calcContactForces(){
	if (friction){
		for (int k = 0; k < num_interaction; k++){
			interaction[k].calcContactInteraction();
		}
	} else {
		for (int k=0; k < num_interaction; k++){
			interaction[k].calcContactInteractionNoFriction();
		}
	}
}



void System::updateContactForceConfig(){
	
	for (int k = 0; k < num_interaction; k++){
		interaction[k].incrementContactTangentialDisplacement();
	}
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

void System::stressReset(){
	for (int i=0; i < n; i++){
		for (int j=0; j < 5; j++){
			stress[i][j]=0;
		}
	}
}

/*
 * Free-draining approximation
 */
void System::updateVelocity(){
	vec3d U_inf(0, 0, 0);
	for (int i=0; i < n; i++){
		U_inf.x = shear_rate*position[i].z;
		velocity[i] = (1.0/eta)*force[i] + U_inf;
	}
	if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < n; i++){
			ang_velocity[i] = (1.33333/eta)*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
	}
}


#ifdef CHOLMOD
//off-diagonal terms
void System::appendToColumn(double *nvec, int jj, double alpha){
	int jj3   = 3*jj;
	int jj3_1 = jj3+1;
	int jj3_2 = jj3+2;
	double alpha_n0 = alpha*nvec[0];
	double alpha_n1 = alpha*nvec[1];
	double alpha_n2 = alpha*nvec[2];
	double alpha_n1n0 = alpha_n0*nvec[1];
	double alpha_n2n1 = alpha_n1*nvec[2];
	double alpha_n0n2 = alpha_n2*nvec[0];
	
	rows.push_back(jj3);
	rows.push_back(jj3_1);
	rows.push_back(jj3_2);
	
	off_diag_values[0].push_back(alpha_n0*nvec[0]); // 00
	off_diag_values[0].push_back(alpha_n1n0); // 10
	off_diag_values[0].push_back(alpha_n0n2); // 20
	off_diag_values[1].push_back(alpha_n1n0); // 01
	off_diag_values[1].push_back(alpha_n1*nvec[1]); //11
	off_diag_values[1].push_back(alpha_n2n1); // 21
	off_diag_values[2].push_back(alpha_n0n2); // 02
	off_diag_values[2].push_back(alpha_n2n1); // 12
	off_diag_values[2].push_back(alpha_n2*nvec[2]); //22
}

// diagonal terms
void System::addToDiag(double *nvec, int ii, double alpha){
	int ii6 = 6*ii;
	
	double alpha_n0 = alpha*nvec[0];
	double alpha_n1 = alpha*nvec[1];
	double alpha_n2 = alpha*nvec[2];
	double alpha_n1n0 = alpha_n0*nvec[1];
	double alpha_n2n1 = alpha_n1*nvec[2];
	double alpha_n0n2 = alpha_n2*nvec[0];
	
	diag_values[ii6]   += alpha_n0*nvec[0]; // 00
	diag_values[ii6+1] += alpha_n1n0; // 10
	diag_values[ii6+2] += alpha_n0n2; // 20
	
	diag_values[ii6+3] += alpha_n1*nvec[1]; //11
	diag_values[ii6+4] += alpha_n2n1; // 21
	
	diag_values[ii6+5] += alpha_n2*nvec[2]; //22
	
}

void System::fillSparseResmatrix(){
	// fill
	for(int j = 0; j < n; j++){
		int j3 = 3*j;
		int j6 = 6*j;
		
		((int*)sparse_res->p)[j3  ] = j6   + 3*ploc[j];
		((int*)sparse_res->p)[j3+1] = j6+3 + 2*ploc[j] +   ploc[j+1];
		((int*)sparse_res->p)[j3+2] = j6+5 +   ploc[j] + 2*ploc[j+1];
		
		int pj3   = ((int*)sparse_res->p)[j3];
		int pj3_1 = ((int*)sparse_res->p)[j3+1];
		int pj3_2 = ((int*)sparse_res->p)[j3+2];
		
		// diagonal blocks row indices
		((int*)sparse_res->i)[ pj3 ]       = j3;
		((int*)sparse_res->i)[ pj3 + 1 ]   = j3+1;
		((int*)sparse_res->i)[ pj3 + 2 ]   = j3+2;
		
		((int*)sparse_res->i)[ pj3_1 ]     = j3+1;
		((int*)sparse_res->i)[ pj3_1 + 1 ] = j3+2;
		
		((int*)sparse_res->i)[ pj3_2 ]     = j3+2;
		
		// diagonal blocks row values
		((double*)sparse_res->x)[ pj3 ]       = diag_values[j6];
		((double*)sparse_res->x)[ pj3 + 1 ]   = diag_values[j6+1];
		((double*)sparse_res->x)[ pj3 + 2 ]   = diag_values[j6+2];
		
		((double*)sparse_res->x)[ pj3_1 ]     = diag_values[j6+3];
		((double*)sparse_res->x)[ pj3_1 + 1 ] = diag_values[j6+4];
		
		((double*)sparse_res->x)[ pj3_2 ]     = diag_values[j6+5];
		
		//    cout << j3+2 <<" " << diag_values[j6+5]<< " " << ((int*)sparse_res->p)[j3] << " " << ((int*)sparse_res->p)[j3+1] << " " << ((int*)sparse_res->p)[j3+2] <<endl;
		// off-diagonal blocks row indices and values
		for(int k = ploc[j]; k < ploc[j+1]; k++){
			int u = k - ploc[j];
			((int*)sparse_res->i)[ pj3   + u + 3 ] = rows[k];
			((int*)sparse_res->i)[ pj3_1 + u + 2 ] = rows[k];
			((int*)sparse_res->i)[ pj3_2 + u + 1 ] = rows[k];
			
			((double*)sparse_res->x)[ pj3   + u + 3 ] = off_diag_values[0][k];
			((double*)sparse_res->x)[ pj3_1 + u + 2 ] = off_diag_values[1][k];
			((double*)sparse_res->x)[ pj3_2 + u + 1 ] = off_diag_values[2][k];
		}
	}
	((int*)sparse_res->p)[n3] = ((int*)sparse_res->p)[n3-1] + 1;
}

#else
void fillResmatrix(double *res, double *nvec, int ii, int jj, double alpha, int n3){
	int ii3     = 3*ii;
	int n3ii3   = n3*ii3;
	int n3ii3_1 = n3*(ii3+1);
	int n3ii3_2 = n3*(ii3+2);
	int jj3   = 3*jj;
	int jj3_1 = jj3+1;
	int jj3_2 = jj3+2;
	double alpha_n0 = alpha*nvec[0];
	double alpha_n1 = alpha*nvec[1];
	double alpha_n2 = alpha*nvec[2];
	double alpha_n1n0 = alpha_n0*nvec[1];
	double alpha_n2n1 = alpha_n1*nvec[2];
	double alpha_n0n2 = alpha_n2*nvec[0];
	
	res[n3ii3   + jj3  ] += alpha_n0*nvec[0]; // 00
	res[n3ii3   + jj3_1] += alpha_n1n0; // 10
	res[n3ii3   + jj3_2] += alpha_n0n2; // 20
	res[n3ii3_1 + jj3  ] += alpha_n1n0; // 01
	res[n3ii3_1 + jj3_1] += alpha_n1*nvec[1]; //11
	res[n3ii3_1 + jj3_2] += alpha_n2n1; // 21
	res[n3ii3_2 + jj3  ] += alpha_n0n2; // 02
	res[n3ii3_2 + jj3_1] += alpha_n2n1; // 12
	res[n3ii3_2 + jj3_2] += alpha_n2*nvec[2]; //22
}
#endif


#ifdef CHOLMOD
/*
 * fills resistance matrix and part of the rhs force coming from lubrication
 *
 */
void System::buildLubricationTerms(){
	for (int k = 0; k < 6*n; k++){
		diag_values[k] = 0.;
	}
	rows.clear();
	off_diag_values[0].clear();
	off_diag_values[1].clear();
	off_diag_values[2].clear();
	for (int i = 0; i < n; i ++){
		int i6=6*i;
		diag_values[i6  ] = 1.;
		diag_values[i6+3] = 1.;
		diag_values[i6+5] = 1.;
	}
	
	for (int i = 0; i < n - 1; i ++){
		ploc[i] = (unsigned int)rows.size();
		for (int j = i+1 ; j < n; j ++){
			double h = 0;
			double nvec[3];
			int k = interaction_pair[i][j];
			if( k != -1){
				nvec[0] = interaction[k].nr_vec.x;
				nvec[1] = interaction[k].nr_vec.y;
				nvec[2] = interaction[k].nr_vec.z;
				h = interaction[k].r - lubcore;
				if(h > 0){
					double alpha = - 1/(4*h);
					// (i, j) (k,l) --> res[ n3*(3*i+l) + 3*j+k ]
					addToDiag(nvec, i, -alpha);
					addToDiag(nvec, j, -alpha);
					appendToColumn(nvec, j, +alpha);
					double alpha_gd_dz_n0 = alpha*shear_rate*interaction[k].r_vec.z*nvec[0];
					
					
					double alpha_gd_dz_n0_n[] = { \
						alpha_gd_dz_n0*nvec[0],
						alpha_gd_dz_n0*nvec[1],
						alpha_gd_dz_n0*nvec[2]};
					
					((double*)rhs_b->x)[3*i  ] += alpha_gd_dz_n0_n[0];
					((double*)rhs_b->x)[3*i+1] += alpha_gd_dz_n0_n[1];
					((double*)rhs_b->x)[3*i+2] += alpha_gd_dz_n0_n[2];
					((double*)rhs_b->x)[3*j  ] -= alpha_gd_dz_n0_n[0];
					((double*)rhs_b->x)[3*j+1] -= alpha_gd_dz_n0_n[1];
					((double*)rhs_b->x)[3*j+2] -= alpha_gd_dz_n0_n[2];
				} else {
					cerr << i << ' ' << j << ' ' << endl;
					cerr << "h<0 : " << h <<   endl;
					//				cerr << "r = " << interaction[k].r << endl;
					position[i].cerr();
					position[j].cerr();
					//				interaction[k].nr_vec.cerr();
					
					exit(1);
				}
			}
		}
	}
	ploc[n-1] = (unsigned int)rows.size();
	ploc[n] = (unsigned int)rows.size();
	
}
#endif

#ifdef CHOLMOD
void System::buildBrownianTerms(){
	// add Brownian force
	fb->add_to(rhs_b);
}
// #else
// void System::buildBrownianTerms(){
// 	// add Brownian force
// 	//	fb->generate(rhs_b); // right now not working, as it relies on Cholesky factor
// }
#endif

#ifdef CHOLMOD
void System::buildContactTerms(){
	// add contact force
	for (int i = 0; i < n; i++){
		int i3 = 3*i;
		((double*)rhs_b->x)[i3] += force[i].x;
		((double*)rhs_b->x)[i3+1] += force[i].y;
		((double*)rhs_b->x)[i3+2] += force[i].z;
	}
}
// #else
// void System::buildContactTerms(){
// 	// add contact force
// 	for (int i = 0; i < n; i++){
// 		int i3 = 3*i;
// 		rhs_b[i3] += force[i].x;
// 		rhs_b[i3+1] += force[i].y;
// 		rhs_b[i3+2] += force[i].z;
// 	}
// }
#endif

#ifdef CHOLMOD
void System::updateVelocityLubrication(){
	rhs_b = cholmod_zeros(n3, 1, xtype, &c);
	buildLubricationTerms();
	// allocate
	int nzmax;  // non-zero values
	nzmax = 6*n; // diagonal blocks
	for(int s=0; s<3; s++){
		nzmax += off_diag_values[s].size();  // off-diagonal
	}
	sparse_res = cholmod_allocate_sparse(n3, n3, nzmax, sorted, packed, stype,xtype, &c);
	fillSparseResmatrix();
	
	L = cholmod_analyze (sparse_res, &c);
	cholmod_factorize (sparse_res, L, &c);
	
	buildContactTerms();
	
	if (brownian){
		buildBrownianTerms();
	}
	
	v = cholmod_solve (CHOLMOD_A, L, rhs_b, &c) ;
	for (int i = 0; i < n; i++){
		int i3 = 3*i;
		velocity[i].x = ((double*)v->x)[i3] + shear_rate*position[i].z;
		velocity[i].y = ((double*)v->x)[i3+1];
		velocity[i].z = ((double*)v->x)[i3+2];
		//velocity[i].cerr();
	}
	if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < n; i++){
			ang_velocity[i] = (1.33333/eta)*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
	}
	
	cholmod_free_sparse(&sparse_res, &c);
	cholmod_free_factor(&L, &c);
	cholmod_free_dense(&rhs_b, &c);
	cholmod_free_dense(&v, &c);
}

#else

/*
 * By using lapack
 */
void System::updateVelocityLubrication(){
	for (int k = 0;k < n3; k++){
		rhs_b[k] = 0.;
	}
	for (int k = 0; k < n3*n3; k++){
		res[k] = 0.;
	}
	for (int i = 0 ; i < n; i ++){
		int i3 = 3*i;
		res[n3*(i3  ) + i3  ] = 1.;
		res[n3*(i3+1) + i3+1] = 1.;
		res[n3*(i3+2) + i3+2] = 1.;
	}
	if (lub){
		for (int i = 0 ; i < n - 1; i ++){
			for (int j = i+1 ; j < n; j ++){
				double r_sq = sq_distance(i,j);
				double r;
				if( r_sq < sq_lub_max){
					r = sqrt(r_sq);
					double nvec[] = {dx/r, dy/r, dz/r};
					double h = r - lubcore;
					double alpha = - 1/(4*h);
					if (h > 0){
						// (i, j) (k,l) --> res[ n3*(3*i+l) + 3*j+k ]
						fillResmatrix(res, nvec, i, i, -alpha, n3);
						fillResmatrix(res, nvec, i, j, +alpha, n3);
						fillResmatrix(res, nvec, j, j, -alpha, n3);
						fillResmatrix(res, nvec, j, i, +alpha, n3);
						double alpha_gd_dz_n0 = alpha*shear_rate*dz*nvec[0];
						double alpha_gd_dz_n0_n[] = { \
							alpha_gd_dz_n0*nvec[0],
							alpha_gd_dz_n0*nvec[1],
							alpha_gd_dz_n0*nvec[2]};
						rhs_b[3*i  ] += alpha_gd_dz_n0_n[0];
						rhs_b[3*i+1] += alpha_gd_dz_n0_n[1];
						rhs_b[3*i+2] += alpha_gd_dz_n0_n[2];
						rhs_b[3*j  ] -= alpha_gd_dz_n0_n[0];
						rhs_b[3*j+1] -= alpha_gd_dz_n0_n[1];
						rhs_b[3*j+2] -= alpha_gd_dz_n0_n[2];
					}
				}
			}
		}
	}
	for (int i = 0; i < n; i++){
		int i3 = 3*i;
		rhs_b[i3  ] += force[i].x;
		rhs_b[i3+1] += force[i].y;
		rhs_b[i3+2] += force[i].z;
	}
	
	dsysv_(&UPLO, &n3, &nrhs, res, &lda, ipiv, rhs_b, &ldb, work, &lwork, &info);
	for (int i = 0; i < n; i++){
		int i3 = 3*i;
		velocity[i].x = rhs_b[i3] + shear_rate*position[i].z;
		velocity[i].y = rhs_b[i3+1];
		velocity[i].z = rhs_b[i3+2];
	}
	if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < n; i++){
			ang_velocity[i] = 1.33333*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
	}
}
#endif

void System::displacement(int i, const double &dx_, const double &dy_, const double &dz_){
	position[i].x += dx_;
	position[i].y += dy_;
	position[i].z += dz_;
	
	if (position[i].z > lz ){
		position[i].z -= lz;
		position[i].x -= shear_disp;
	} else if ( position[i].z < 0 ){
		position[i].z += lz;
		position[i].x += shear_disp;
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
	shear_disp += vel_difference*dt;
	if (shear_disp > lx){
		shear_disp -= lx;
	}
	for (int i=0; i < n; i++){
		displacement(i, velocity[i].x*dt, velocity[i].y*dt, velocity[i].z*dt);
	}
	
	if (draw_rotation_2d){
		for (int i=0; i < n; i++){
			angle[i] += ang_velocity[i].y*dt;
		}
	}
}

/*
 * Distance between particle i and particle j
 */
double System::distance(int i, int j){
	return sqrt(sq_distance(i,j));
}

/*
 * Square norm of vector (dx, dy, dz)
 */

/*
 * Square distance between particle i and particle j
 */
double System::sq_distance(int i, int j){
	double dx = position[i].x - position[j].x;
	double dy = position[i].y - position[j].y;
	double dz = position[i].z - position[j].z;
	if (dz > lz2 ){
		dz -= lz;
		dx -= shear_disp;
	} else if (dz < -lz2){
		dz += lz;
		dx += shear_disp;
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

/*
 * Distance less than 'lub_max' can be calculated.
 * Otherwize it returns 1000.
 */
double System::sq_distanceToCheckContact(int i, int j){
	
	
	double dx = position[i].x - position[j].x;
	double dz = position[i].z - position[j].z;
	if (dz > lz2 ){
		dz -= lz;
		dx -= shear_disp;
	} else if (dz < -lz2){
		dz += lz;
		dx += shear_disp;
	}
	if (abs(dz) < lub_max){
		while(dx > lx2){
			dx -= lx;
		}
		while(dx < -lx2){
			dx += lx;
		}
		if (abs(dx) < lub_max){
			if (dimension == 3){
				double dy = position[i].y - position[j].y;
				if (dy > ly2 ){
					dy -= ly;
				} else if (dy < -ly2){
					dy += ly;
				}
				if (abs(dy) < lub_max){
					return dx*dx + dy*dy + dz*dz;
				}
			} else {
				return dx*dx + dz*dz;
			}
		}
	}
	return 1000;
}


double System::lubricationForceFactor(int i, int j){
	//double r_sq = sq_distance(i,j);
	int k = interaction_pair[i][j];
	if(interaction[k].active){
		double h = interaction[k].r  - lubcore;
		if ( h > 0){
			vec3d nv = interaction[k].nr_vec;
			vec3d rel_vel = velocity[j] - velocity[i];
			double zi_zj = position[i].z - position[j].z;
			if (zi_zj > lz2){
				rel_vel.x += vel_difference;
			} else if (zi_zj < -lz2){
				rel_vel.x -= vel_difference;
			}
			double alpha = 1.0/(4*h);
			return abs(alpha*dot(rel_vel , nv));
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}

void System::lubricationStress(int i, int j){
	//	double r_sq = sq_distance(i, j);
	int k = interaction_pair[i][j];
	if(interaction[k].active){
		
		//		double r = sqrt(r_sq);
		//		double h = r - lubcore;
		double h = interaction[k].r  - lubcore;
		//		vec3d nv(dx/r, dy/r, dz/r);
		vec3d nv =interaction[k].nr_vec;
		
		vec3d rel_vel = velocity[j] - velocity[i];
		double zi_zj = position[i].z - position[j].z;
		if (zi_zj > lz2){
			rel_vel.x += vel_difference;
		} else if (zi_zj < -lz2){
			rel_vel.x -= vel_difference;
		}
		double alpha;
		if (h > 0){
			alpha = 1.0/(4*h);
			vec3d f_lub_ij = -alpha*dot(rel_vel , nv)*nv;
			
			double Sxx = 2*(f_lub_ij.x * nv.x);
			double Sxy = f_lub_ij.x * nv.y + f_lub_ij.y * nv.x ;
			double Sxz = f_lub_ij.x * nv.z + f_lub_ij.z * nv.x ;
			double Syz = f_lub_ij.y * nv.z + f_lub_ij.z * nv.y ;
			double Syy = 2*(f_lub_ij.y * nv.y);
			stress[i][0] += Sxx;
			stress[j][0] += Sxx;
			
			stress[i][1] += Sxy;
			stress[j][1] += Sxy;
			
			stress[i][2] += Sxz;
			stress[j][2] += Sxz;
			
			stress[i][3] += Syz;
			stress[j][3] += Syz;
			
			stress[i][4] += Syy;
			stress[j][4] += Syy;
		}
		
		
	}
}


void System::calcStressAverage(){
	double totalStress[5] = {0,0,0,0,0};
	for (int i=0; i< n ; i++){
		for (int k=0; k < 5; k++){
			totalStress[k] += stress[i][k];
		}
	}
	for (int k=0; k < 5; k++){
		mean_stress[k] = totalStress[k]/n;
	}
}










