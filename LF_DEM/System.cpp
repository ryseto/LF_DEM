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
	if (!position)
		delete [] position;
	if (!radius)
		delete [] radius;
	if (!angle)
		delete [] angle;
	if (!velocity)
		delete [] velocity;
	if (!relative_velocity)
		delete [] relative_velocity;
	if (!ang_velocity)
		delete [] ang_velocity;
	if (!force)
		delete [] force;
	if (!torque)
		delete [] torque;
	if (!interaction_list)
		delete [] interaction_list;
	if (!interaction_partners)
		delete [] interaction_partners;
	
#ifdef CHOLMOD
	if (!diag_values)
		delete [] diag_values;
	if (!off_diag_values )
		delete [] off_diag_values;
	if (!ploc)
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

void
System::allocateRessources(){
	
	position = new vec3d [n];
	radius = new double [n];
	angle = new double [n];
	velocity = new vec3d [n];
	relative_velocity = new vec3d [n];
	ang_velocity = new vec3d [n];
	force = new vec3d [n];
	torque = new vec3d [n];
	contactstress = new double* [n];
	
	for (int i=0; i < n; i++){
		position[i].x=0.;
		position[i].y=0.;
		position[i].z=0.;
		radius[i]=0.;
		velocity[i].x=0.;
		velocity[i].y=0.;
		velocity[i].z=0.;
		relative_velocity[i].x=0.;
		relative_velocity[i].y=0.;
		relative_velocity[i].z=0.;
	}
	
	for (int i=0; i < n; i++){
		contactstress[i] = new double [5];
	}
	
	lubstress = new double* [n];
	for (int i=0; i < n; i++){
		lubstress[i] = new double [5];
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
	interaction_list = new set <Interaction*> [n];
	interaction_partners = new set <int> [n];
	
	
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

void
System::initializeBoxing(){// need to know radii first
	
	double max_radius=0.;
	for (int i=0; i < n; i++){
		if(radius[i]>max_radius){
			max_radius=radius[i];
		}
	}
	
	boxset = new BoxSet(lub_max*max_radius, this);
	for (int i=0; i < n; i++){
		boxset->box(i);
	}
}

void
System::timeEvolution(int time_step){
	int ts_next = ts + time_step;
	checkNewInteraction();
	while (ts < ts_next){
		forceReset();
		torqueReset();
		calcContactForces();
		if (lubrication){
			// Lubrication dynamics
			if(brownian){
				updateVelocityLubricationBrownian();
			}
			else{
				updateVelocityLubrication();
			}
		} else {
			// Free-draining approximation
			updateVelocity();
		}
		deltaTimeEvolution();
		ts ++;
	}
}

void
System::checkNewInteraction(){
	vector<int>::iterator it;
	vector<int>::iterator it_beg;
	vector<int>::iterator it_end;
	vec3d pos_diff;
	int zshift;
	double sq_dist;
	for (int i=0; i < n-1; i++){
		
		it_beg = boxset->neighborhood_begin(i);
		it_end = boxset->neighborhood_end(i);
		
		for (it = it_beg; it != it_end; it++){
			int j=*it;
			if(j>i){
				if ( interaction_partners[i].find(j) == interaction_partners[i].end() ){
					
					// this is done in 3 steps because we need each information for Interaction creation
					pos_diff = position[j] - position[i];
					periodize_diff(&pos_diff, &zshift);
					sq_dist = pos_diff.sq_norm();
					
					double sq_dist_lim = sq_lub_max * 0.25 * ( radius[i] + radius[j] ) * ( radius[i] + radius[j] );
					if ( sq_dist < sq_dist_lim){
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
						// new interaction
						interaction[interaction_new].activate(i, j, +pos_diff, sqrt(sq_dist), zshift);
					}
				}
			}
		}
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
void
System::updateInteractions(){
	
	for (int k = 0; k < num_interaction; k++){
		int switch_off = interaction[k].update();
		if(switch_off)
			deactivated_interaction.push(k);
	}
	
}


void
System::forceReset(){
	for (int i=0; i < n; i++){
		force[i].reset();
	}
}

void
System::torqueReset(){
	for (int i=0; i < n; i++){
		torque[i].reset();
	}
}

void
System::stressReset(){
	for (int i=0; i < n; i++){
		for (int j=0; j < 5; j++){
			lubstress[i][j]=0;
			contactstress[i][j]=0;
		}
	}
}

/*
 * Free-draining approximation
 */
void
System::updateVelocity(){
	vec3d U_inf(0, 0, 0);
	for (int i=0; i < n; i++){
		U_inf.x = shear_rate*position[i].z;
		relative_velocity[i] = (1.0/eta)*force[i];
		velocity[i] = relative_velocity[i] + U_inf;
		
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
void
System::appendToColumn(double *nvec, int jj, double alpha){
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
void
System::addToDiag(double *nvec, int ii, double alpha){
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

void
System::fillSparseResmatrix(){
	
	allocateSparseResmatrix();
	
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
void
fillResmatrix(double *res, double *nvec, int ii, int jj, double alpha, int n3){
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
void
System::addStokesDrag(){
	for (int i = 0; i < n; i ++){
		int i6=6*i;
		double d_value = bgf_factor*radius[i];
		diag_values[i6  ] = d_value;
		diag_values[i6+3] = d_value;
		diag_values[i6+5] = d_value;
	}
}



// We solve A*(U-Uinf) = Gtilde*Einf
// This method computes elements of matrix A and vector Gtilde*Einf
void
System::buildLubricationTerms(){
	for (int k = 0; k < 6*n; k++){
		diag_values[k] = 0.;
	}
	rows.clear();
	off_diag_values[0].clear();
	off_diag_values[1].clear();
	off_diag_values[2].clear();
	
	addStokesDrag();
	
	double XAii, XAjj, XAij, XAji;
	//	double XGii, XGjj, XGij, XGji;
	
	double nvec[3];
	double GEi[3];
	double GEj[3];
	
	set<Interaction*>::iterator it;
	int j;
	Interaction *inter;
	for (int i = 0; i < n - 1; i ++){
		ploc[i] = (unsigned int)rows.size();
		
		for (it = interaction_list[i].begin() ; it != interaction_list[i].end(); it ++){
			inter=*it;
		 	j=inter->partner(i);
			if(j>i){
				
				nvec[0] = inter->nr_vec.x; // nvec defined from i to j
				nvec[1] = inter->nr_vec.y;
				nvec[2] = inter->nr_vec.z;
				
				inter->XA(XAii, XAij, XAji, XAjj);
				
				// (i, j) (k,l) --> res[ n3*(3*i+l) + 3*j+k ]
				addToDiag(nvec, i, inter->a0 * XAii);
				addToDiag(nvec, j, inter->a1 * XAjj);
				appendToColumn(nvec, j, 0.5 * inter->ro * XAji);
				
				
				inter->GE(GEi, GEj);  // G*E_\infty term
				
				for(int u=0; u<3; u++){
					((double*)rhs_b->x)[ 3*i + u ] += GEi[ u ];
					((double*)rhs_b->x)[ 3*j + u ] += GEj[ u ];
				}
			}
		}
	}
	
	ploc[n-1] = (unsigned int)rows.size();
	ploc[n] = (unsigned int)rows.size();
	
}

void
System::buildBrownianTerms(){
	// add Brownian force
	fb->add_to(rhs_b);
}
// #else
// void System::buildBrownianTerms(){
// 	// add Brownian force
// 	//	fb->generate(rhs_b); // right now not working, as it relies on Cholesky factor
// }

void
System::buildContactTerms(){
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

void
System::allocateSparseResmatrix(){
	// allocate
	int nzmax;  // non-zero values
	nzmax = 6*n; // diagonal blocks
	for(int s=0; s<3; s++){
		nzmax += off_diag_values[s].size();  // off-diagonal
	}
	sparse_res = cholmod_allocate_sparse(n3, n3, nzmax, sorted, packed, stype,xtype, &c);
	
}



void
System::print_res(){ // testing
	cout << " Diag " << endl;
	for(int i=0;i<n;i++){
		int ii6=6*i;
		cout << i << " " << diag_values[ii6] << " " << diag_values[ii6+1]<< " " << diag_values[ii6+2]<< " " << diag_values[ii6+3]<< " " << diag_values[ii6+4]<< " " << diag_values[ii6+5] << endl;
	}
	
	cout << endl<< " OffDiag " << endl;
	for(int i=0;i<off_diag_values[0].size();i++){
		cout << off_diag_values[0][i] << " " << off_diag_values[1][i] << " " << off_diag_values[2][i]<< endl;
	}
	
	
}

void
System::updateVelocityLubrication(){
	rhs_b = cholmod_zeros(n3, 1, xtype, &c);
	buildLubricationTerms();
	
	fillSparseResmatrix();
	
	L = cholmod_analyze (sparse_res, &c);
	cholmod_factorize (sparse_res, L, &c);
	if(c.status){
		// Cholesky decomposition has failed: usually because matrix is incorrectly found to be positive-definite
		// It is very often enough to force another preconditioner to solve the problem.
		cerr << " factorization failed. forcing simplicial algorithm... " << endl;
		c.supernodal = CHOLMOD_SIMPLICIAL;
		L = cholmod_analyze (sparse_res, &c);
		cholmod_factorize (sparse_res, L, &c) ;
		cerr << " factorization status " << c.status << " final_ll ( 0 is LDL, 1 is LL ) " <<  c.final_ll <<endl;
		c.supernodal = CHOLMOD_SUPERNODAL;
	}
	
	buildContactTerms();
	
	v = cholmod_solve (CHOLMOD_A, L, rhs_b, &c) ;
	
	/********** testing *
	 double m1 [2] = {-1,0};
	 double p1 [2] = {1,0};
	 cholmod_dense *r = cholmod_copy_dense(rhs_b, &c);
	 cholmod_sdmult(sparse_res, 0, m1, p1, v, r, &c);
	 cout << " cholmod residu " << cholmod_norm_dense(r,0, &c) << endl;
	 cholmod_free_dense(&r, &c);
	 * end testing *************/
	
	/* TEST IMPLEMENTATION
	 * SDFF : Stokes drag force factor:
	 * SDFF = 1.0 : full drag forces from the undisturbed background flow.
	 * SDFF = 0.0 : no drag force from the undisturbed background flow.
	 */
	for (int i = 0; i < n; i++){
		int i3 = 3*i;
		relative_velocity[i].x = ((double*)v->x)[i3];
		relative_velocity[i].y = ((double*)v->x)[i3+1];
		relative_velocity[i].z = ((double*)v->x)[i3+2];
		
		velocity[i].x = relative_velocity[i].x + shear_rate*position[i].z;
		velocity[i].y = relative_velocity[i].y;
		velocity[i].z = relative_velocity[i].z;
	}
	
	if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < n; i++){
			ang_velocity[i] = 1.33333*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
	}
	
	cholmod_free_sparse(&sparse_res, &c);
	cholmod_free_factor(&L, &c);
	cholmod_free_dense(&rhs_b, &c);
	cholmod_free_dense(&v, &c);
}

void System::updateVelocityLubricationBrownian(){
	
	rhs_b = cholmod_zeros(n3, 1, xtype, &c);
	buildLubricationTerms();
	
	fillSparseResmatrix();
	
	L = cholmod_analyze (sparse_res, &c);
	cholmod_factorize (sparse_res, L, &c);
	if(c.status){
		// Cholesky decomposition has failed: usually because matrix is incorrectly found to be positive-definite
		// It is very often enough to force another preconditioner to solve the problem.
		cerr << " factorization failed. forcing simplicial algorithm... " << endl;
		c.supernodal = CHOLMOD_SIMPLICIAL;
		L = cholmod_analyze (sparse_res, &c);
		cholmod_factorize (sparse_res, L, &c) ;
		cerr << " factorization status " << c.status << " final_ll ( 0 is LDL, 1 is LL ) " <<  c.final_ll <<endl;
		c.supernodal = CHOLMOD_SUPERNODAL;
	}
	
	buildContactTerms();
	
	v_nonBrownian = cholmod_solve (CHOLMOD_A, L, rhs_b, &c) ;
	// now the Brownian part of the velocity:
	// mid-point algortithm a la Banchio & Brady
	brownian_force = fb->generate();
	v_Brownian_init = cholmod_solve (CHOLMOD_A, L, brownian_force, &c) ;
	
	// move particles to intermediate point
	for (int i=0; i < n; i++){
		int i3 = 3*i;
		displacement(i, ((double*)v_Brownian_init->x)[i3]*dt_mid, ((double*)v_Brownian_init->x)[i3+1]*dt_mid, ((double*)v_Brownian_init->x)[i3+2]*dt_mid);
	}
	
	updateInteractions();
	
	// rebuild new R_FU
	cholmod_free_factor(&L, &c);
	cholmod_free_sparse(&sparse_res, &c);
	buildLubricationTerms();
	fillSparseResmatrix();
	L = cholmod_analyze (sparse_res, &c);
	cholmod_factorize (sparse_res, L, &c);
	if(c.status){
		// Cholesky decomposition has failed: usually because matrix is incorrectly found to be positive-definite
		// It is very often enough to force another preconditioner to solve the problem.
		cerr << " factorization failed. forcing simplicial algorithm... " << endl;
		c.supernodal = CHOLMOD_SIMPLICIAL;
		L = cholmod_analyze (sparse_res, &c);
		cholmod_factorize (sparse_res, L, &c) ;
		cerr << " factorization status " << c.status << " final_ll ( 0 is LDL, 1 is LL ) " <<  c.final_ll <<endl;
		c.supernodal = CHOLMOD_SUPERNODAL;
	}
	
	// get the intermediate brownian velocity
	v_Brownian_mid = cholmod_solve (CHOLMOD_A, L, brownian_force, &c) ;
	
	/* testing
	 for (int i=0; i < n; i++){
	 int i3 = 3*i;
	 cout << ((double*)v_Brownian_init->x)[i3] << " " << ((double*)v_Brownian_init->x)[i3+1] << " " << ((double*)v_Brownian_init->x)[i3+2] << " "  <<0.5*dt_ratio*(((double*)v_Brownian_mid->x)[i3] - ((double*)v_Brownian_init->x)[i3] ) << " " << 0.5*dt_ratio*(((double*)v_Brownian_mid->x)[i3+1] - ((double*)v_Brownian_init->x)[i3+1] ) << " " << 0.5*dt_ratio*(((double*)v_Brownian_mid->x)[i3+2] - ((double*)v_Brownian_init->x)[i3+2] ) << " "  << ((double*)v_nonBrownian->x)[i3] << " " << ((double*)v_nonBrownian->x)[i3+1] << " " << ((double*)v_nonBrownian->x)[i3+2] <<endl;
	 }
	 getchar();
	 */
	
	// move particles back to initial point
	for (int i=0; i < n; i++){
		int i3 = 3*i;
		displacement(i, -((double*)v_Brownian_init->x)[i3]*dt_mid, -((double*)v_Brownian_init->x)[i3+1]*dt_mid, -((double*)v_Brownian_init->x)[i3+2]*dt_mid);
	}
	
	updateInteractions();
	
	// update total velocity
	// first term is hydrodynamic + contact velocities
	// second term is Brownian velocities
	// third term is Brownian drift
	// fourth term for vx is the shear rate
	for (int i = 0; i < n; i++){
		int i3 = 3*i;
		relative_velocity[i].x = ((double*)v_nonBrownian->x)[i3] + ((double*)v_Brownian_init->x)[i3] + 0.5*dt_ratio*(((double*)v_Brownian_mid->x)[i3] - ((double*)v_Brownian_init->x)[i3] );
		relative_velocity[i].y = ((double*)v_nonBrownian->x)[i3+1] + ((double*)v_Brownian_init->x)[i3+1] + 0.5*dt_ratio*(((double*)v_Brownian_mid->x)[i3+1] - ((double*)v_Brownian_init->x)[i3+1] );
		relative_velocity[i].z = ((double*)v_nonBrownian->x)[i3+2] + ((double*)v_Brownian_init->x)[i3+2] + 0.5*dt_ratio*(((double*)v_Brownian_mid->x)[i3+2] - ((double*)v_Brownian_init->x)[i3+2] );
		
		velocity[i].x = relative_velocity[i].x + shear_rate*position[i].z;
		velocity[i].y = relative_velocity[i].y;
		velocity[i].z = relative_velocity[i].z;
	}
	
	if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < n; i++){
			ang_velocity[i] = 1.33333*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
	}
	
	cholmod_free_sparse(&sparse_res, &c);
	cholmod_free_factor(&L, &c);
	cholmod_free_dense(&rhs_b, &c);
	cholmod_free_dense(&v_nonBrownian, &c);
	cholmod_free_dense(&v_Brownian_init, &c);
	cholmod_free_dense(&v_Brownian_mid, &c);
}


#else

/*
 * By using lapack
 */
void
System::updateVelocityLubrication(){
	for (int k = 0;k < n3; k++){
		rhs_b[k] = 0.;
	}
	for (int k = 0; k < n3*n3; k++){
		res[k] = 0.;
	}

	for (int i = 0 ; i < n; i ++){
		int i3 = 3*i;
		double d_value = bgf_factor*radius[i];
		res[n3*(i3  ) + i3  ] = d_value;
		res[n3*(i3+1) + i3+1] = d_value;
		res[n3*(i3+2) + i3+2] = d_value;
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
	
	/********** testing *
	 double * rhs_b_cpy = new double [n3];
	 for (int i = 0; i < n3; i++){
	 rhs_b_cpy[i]= rhs_b[i];
	 }
	 * end testing ***********/
	dsysv_(&UPLO, &n3, &nrhs, res, &lda, ipiv, rhs_b, &ldb, work, &lwork, &info);
	
	/********** testing *
	 int inc=1;
	 double beta=-1.;
	 double alpha=1.;
	 enum CBLAS_ORDER order=CblasRowMajor;
	 enum CBLAS_UPLO ul=CblasUpper;
	 cblas_dsymv(order, ul, n3, alpha, res, lda, rhs_b, inc, beta, rhs_b_cpy,inc);
	 double infty_norm=0.;
	 for (int i = 0; i < n3; i++){
	 cout << rhs_b_cpy[i] << endl;
	 if(fabs(rhs_b_cpy[i])>infty_norm)
	 infty_norm=fabs(rhs_b_cpy[i]);
	 }
	 
	 cout << " lapack residu " << infty_norm << endl;
	 getchar();
	 delete [] rhs_b_cpy;
	 * end testing ***********/
	
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

void
System::displacement(int i, const double &dx_, const double &dy_, const double &dz_){
	position[i].x += dx_;
	position[i].y += dy_;
	position[i].z += dz_;
	
	periodize(&(position[i]));
	boxset->box(i);
}


// [0,l]
void
System::periodize(vec3d *pos){
	if (pos->z > lz() ){
		pos->z -= lz();
		pos->x -= shear_disp;
	} else if ( pos->z < 0 ){
		pos->z += lz();
		pos->x += shear_disp;
	}
	while ( pos->x > lx() ){
		pos->x -= lx();
	}
	while (pos->x < 0 ){
		pos->x += lx();
	}
	if (dimension == 3){
		if ( pos->y > ly() ){
			pos->y -= ly();
		} else if (pos->y < 0 ){
			pos->y += ly();
		}
	}
}

// [-l/2,l/2]
void
System::periodize_diff(vec3d *pos_diff){
	if (pos_diff->z > lz2() ){
		pos_diff->z -= lz();
		pos_diff->x -= shear_disp;
	} else if ( pos_diff->z < -lz2() ){
		pos_diff->z += lz();
		pos_diff->x += shear_disp;
	}
	while ( pos_diff->x > lx2() ){
		pos_diff->x -= lx();
	}
	while (pos_diff->x < -lx2() ){
		pos_diff->x += lx();
	}
	if (dimension == 3){
		if ( pos_diff->y > ly2() ){
			pos_diff->y -= ly();
		} else if (pos_diff->y < -ly2() ){
			pos_diff->y += ly();
		}
	}
}
// periodize + give z_shift= number of boundaries crossed in z-direction
void
System::periodize_diff(vec3d *pos_diff, int *zshift){
	if (pos_diff->z > lz2() ){
		pos_diff->z -= lz();
		pos_diff->x -= shear_disp;
		(*zshift)=-1;
	} else if ( pos_diff->z < -lz2() ){
		pos_diff->z += lz();
		pos_diff->x += shear_disp;
		(*zshift)=+1;
	}
	else{
		(*zshift)=0;
	}
	while ( pos_diff->x > lx2() ){
		pos_diff->x -= lx();
	}
	while (pos_diff->x < -lx2() ){
		pos_diff->x += lx();
	}
	if (dimension == 3){
		if ( pos_diff->y > ly2() ){
			pos_diff->y -= ly();
		} else if (pos_diff->y < -ly2() ){
			pos_diff->y += ly();
		}
	}
}



void
System::deltaTimeEvolution(){
	
	// evolve PBC
	shear_disp += vel_difference*dt;
	if (shear_disp > lx()){
		shear_disp -= lx();
	}
	
	// move particles
	for (int i=0; i < n; i++){
		displacement(i, velocity[i].x*dt, velocity[i].y*dt, velocity[i].z*dt);
	}
	if (draw_rotation_2d){
		for (int i=0; i < n; i++){
			angle[i] += ang_velocity[i].y*dt;
		}
	}
	// update boxing system
	boxset->update();
	
	checkNewInteraction();
	updateInteractions();
	
}

/*
 * Distance between particle i and particle j
 */
double
System::distance(int i, int j){
	return sqrt(sq_distance(i,j));
}

/*
 * Square distance between particle i and particle j
 */
double
System::sq_distance(int i, int j){
	vec3d pos_diff = position[j] - position[i];
	
	periodize_diff(&pos_diff);
	
	if (dimension == 3){
		return pos_diff.sq_norm();
	} else {
	 	return pos_diff.sq_norm_xz();
	}
}

void
System::calcStress(){
	stressReset();
	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			interaction[k].addLubricationStress();
			if (interaction[k].contact){
				interaction[k].addContactStress();
			}
		}
	}
	
	double total_lub_stress[5] = {0,0,0,0,0};
	double total_contact_stress[5] = {0,0,0,0,0};
	for (int i=0; i < n ; i++){
		for (int k=0; k < 5; k++){
			total_lub_stress[k] += lubstress[i][k];
			total_contact_stress[k] += contactstress[i][k];
		}
	}
	/*
	 *  mean_hydro_stress
	 *  The term 5.0/9 is the one-body part
	 *
	 */
	for (int k=0; k < 5; k++){
		mean_hydro_stress[k] = total_lub_stress[k] / n + (5.0/9)*bgf_factor;
		mean_contact_stress[k] = total_contact_stress[k] / n;
	}
}

void
System::calcContactForces(){
	if (friction){
		for (int k=0; k < num_interaction; k++){
			interaction[k].calcContactInteraction();
		}
	} else {
		for (int k=0; k < num_interaction; k++){
			interaction[k].calcContactInteractionNoFriction();
		}
	}
	
}
