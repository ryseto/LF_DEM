//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
#include <sstream>
#include <BelosCGIteration.hpp>

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

	if (!total_force)
		delete [] total_force;
	if (!lubrication_force)
		delete [] lubrication_force;
	if (!contact_force)
		delete [] contact_force;

	if (!brownian_force)
	  delete [] brownian_force;

	if (!torque)
		delete [] torque;

	for (int i=0; i < np; i++){
		delete [] lubstress[i];
	}
	delete [] lubstress;

	for (int i=0; i < np; i++)
		delete [] contactstress[i];
	delete [] contactstress;

	if (brownian){
		for (int i=0; i < np; i++)
			delete [] brownianstress[i];
		delete [] brownianstress;
	}

	if (!interaction_list)
		delete [] interaction_list;
	if (!interaction_partners)
		delete [] interaction_partners;
	if (!fb){
		delete [] fb;
	}
	
	if (!diag_values)
		delete [] diag_values;
	if (!off_diag_values )
		delete [] off_diag_values;
	if (!ploc)
		delete [] ploc;
#ifdef CHOLMOD
	cholmod_free_dense(&v, &c);
	cholmod_free_dense(&total_rhs, &c);
	//	cholmod_free_dense(&nonbrownian_rhs, &c); 
	if (brownian){
		cholmod_free_dense(&brownian_rhs, &c);
	}
	//cholmod_free_dense(&contact_rhs, &c);
	//cholmod_free_dense(&lubrication_rhs, &c);
	cholmod_free_sparse(&sparse_res, &c);
	cholmod_finish(&c);
#endif
};


void
System::allocateRessources(){
	position = new vec3d [np];
	radius = new double [np];
	angle = new double [np];
	velocity = new vec3d [np];
	relative_velocity = new vec3d [np];
	ang_velocity = new vec3d [np];
	total_force = new vec3d [np];
	lubrication_force = new vec3d [np];
	contact_force = new vec3d [np];
	brownian_force = new vec3d [np];
	torque = new vec3d [np];
	lub_force = new vec3d [np];
	contact_number.resize(np);
	contactstress = new double* [np];
	for (int i=0; i < np; i++){
		contactstress[i] = new double [5];
	}
	lubstress = new double* [np];
	for (int i=0; i < np; i++){
		lubstress[i] = new double [5];
	}
	brownianstress = new double* [np];
	for (int i=0; i < np; i++){
		brownianstress[i] = new double [5];
	}
	fb = new BrownianForce(this);

	maxnum_interactionpair = (int)(12*np);
	interaction = new Interaction [maxnum_interactionpair];
	interaction_list = new set <Interaction*> [np];
	interaction_partners = new set <int> [np];
	
	dof = 3;
	linalg_size = dof*np;

#ifdef TRILINOS
	int numlhs = 1;
	int numrhs = 1;
	Map = rcp(new Epetra_Map(linalg_size, 0, Comm));
	v = rcp( new VEC(*Map, numlhs) );
	lubrication_rhs = rcp( new VEC(*Map, numrhs) );
	//	sparse_res = rcp( new MAT(linalg_size) );
	stokes_equation = rcp( new Belos::LinearProblem < SCAL, VEC, MAT > ( sparse_res, v, lubrication_rhs ) ) ;
#endif

	diag_values = new double [6*np];
	off_diag_values = new vector <double> [3];
	ploc = new int [np+1];

	
}

void
System::setupSystem(const vector<vec3d> &initial_positions,
					const vector <double> &radii){
	allocateRessources();

	for (int i=0; i < np; i++){
		position[i] = initial_positions[i];
		radius[i] = radii[i];
		angle[i] = 0;
	}
	
	for (int k=0; k < maxnum_interactionpair ; k++){
		interaction[k].init(this);
	}
	for (int i=0; i < np; i++){
		velocity[i].x=0.;
		velocity[i].y=0.;
		velocity[i].z=0.;
		relative_velocity[i].x=0.;
		relative_velocity[i].y=0.;
		relative_velocity[i].z=0.;
	}
	double O_inf_y = 0.5*shear_rate/2.0;
	for (int i=0; i < np; i++){
		ang_velocity[i].set(0, O_inf_y, 0);
		torque[i].reset();
	}
	initializeBoxing();
	checkNewInteraction();

	dt = dt * 1.0/radius_max;
	shear_strain = 0;
	num_interaction = 0;
	
	if (kb_T == 0){
		brownian = false;
	} else {
		brownian = true;
		fb->init();
	}
	sq_critical_velocity = \
	dynamic_friction_critical_velocity*dynamic_friction_critical_velocity;
	sq_lub_max = lub_max*lub_max; // square of lubrication cutoff length.
	ts = 0;
	shear_disp = 0;
	vel_difference = shear_rate*_lz;
	/*
	 * dt_mid: the intermediate time step for the mid-point
	 * algortithm. dt/dt_mid = dt_ratio
	 * Banchio/Brady (J Chem Phys) gives dt_ratio=100
	 * ASD code from Brady has dt_ratio=150
	 */
	dt_mid = dt/dt_ratio;


	// linear algebra
#ifdef CHOLMOD
	cholmod_start (&c) ;
	stype = -1; // 1 is symmetric, stored upper triangular (UT), -1 is LT
	sorted = 0;		/* TRUE if columns sorted, FALSE otherwise*/
	packed = 1;		/* TRUE if matrix packed, FALSE otherwise */
	xtype = CHOLMOD_REAL;
	L=NULL;
#endif
#ifdef TRILINOS
	RCP<ParameterList> solverParams = parameterList();
	int blocksize = 40;
	int maxiters = 400;
	double tol = 1.e-8;
	solverParams->set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
	solverParams->set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
	solverParams->set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    solverParams->set( "Verbosity", Belos::Errors + Belos::Warnings );
	solver = factory.create ("CG", solverParams);
#endif
	
}

void
System::initializeBoxing(){// need to know radii first
	
	double max_radius=0.;
	for (int i=0; i < np; i++){
		if(radius[i]>max_radius){
			max_radius=radius[i];
		}
	}

	boxset = new BoxSet(lub_max*max_radius, this);
	for (int i=0; i < np; i++){
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
			} else{
				updateVelocityLubrication();
			}
		} else {
			// Free-draining approximation
			updateVelocity();
		}
		if ( ts == ts_next - 1){
			/* Evaluation of the state of suspension.
			 * The following evaluations are valid
			 * when the positions and velocities are consistent.
			 * For that, we need to evaluate them before updating the positions.
			 */
			calcStress();
			analyzeState();

		}
		deltaTimeEvolution();
		ts ++;
		shear_strain += shear_rate*dt;
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
	for (int i=0; i < np-1; i++){
		it_beg = boxset->neighborhood_begin(i);
		it_end = boxset->neighborhood_end(i);
		for (it = it_beg; it != it_end; it++){
			int j=*it;
			if(j>i){
				if ( interaction_partners[i].find(j) == interaction_partners[i].end() ){
					// this is done in 3 steps because we need each information for Interaction creation
					pos_diff = position[j] - position[i];
					periodize_diff(pos_diff, zshift);
					sq_dist = pos_diff.sq_norm();
					double ri_rj_2 = 0.5*(radius[i] + radius[j]);
					double sq_dist_lim = sq_lub_max * ri_rj_2 * ri_rj_2;
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
		bool switch_off = interaction[k].update();
		if(switch_off)
			deactivated_interaction.push(k);
	}
}

void
System::forceReset(){
	for (int i=0; i < np; i++){
		total_force[i].reset();
		lubrication_force[i].reset();
		contact_force[i].reset();
		brownian_force[i].reset();
	}
}

void
System::torqueReset(){
	for (int i=0; i < np; i++){
		torque[i].reset();
	}
}

void
System::stressReset(){
	for (int i=0; i < np; i++){
		for (int j=0; j < 5; j++){
			lubstress[i][j]=0;
			contactstress[i][j]=0;
			brownianstress[i][j]=0;
		}
	}
}

/*
 * Free-draining approximation
 */
void
System::updateVelocity(){
	vec3d U_inf(0, 0, 0);
	for (int i=0; i < np; i++){
		U_inf.x = shear_rate*position[i].z;
		relative_velocity[i] = (1.0/eta)*total_force[i];
		velocity[i] = relative_velocity[i] + U_inf;
	}
	if(friction){
		double O_inf_y = 0.5*shear_rate;
		for (int i=0; i < np; i++){
			ang_velocity[i] = (1.33333/eta)*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
	}
}

//off-diagonal terms

void
System::appendToColumn(const vec3d &nvec, int jj, double alpha){
	int jj3   = 3*jj;
	int jj3_1 = jj3+1;
	int jj3_2 = jj3+2;
	double alpha_n0 = alpha*nvec.x;
	double alpha_n1 = alpha*nvec.y;
	double alpha_n2 = alpha*nvec.z;
	double alpha_n1n0 = alpha_n0*nvec.y;
	double alpha_n2n1 = alpha_n1*nvec.z;
	double alpha_n0n2 = alpha_n2*nvec.x;
	
	rows.push_back(jj3);
	rows.push_back(jj3_1);
	rows.push_back(jj3_2);
	
	off_diag_values[0].push_back(alpha_n0*nvec.x); // 00
	off_diag_values[0].push_back(alpha_n1n0); // 10
	off_diag_values[0].push_back(alpha_n0n2); // 20
	off_diag_values[1].push_back(alpha_n1n0); // 01
	off_diag_values[1].push_back(alpha_n1*nvec.y); //11
	off_diag_values[1].push_back(alpha_n2n1); // 21
	off_diag_values[2].push_back(alpha_n0n2); // 02
	off_diag_values[2].push_back(alpha_n2n1); // 12
	off_diag_values[2].push_back(alpha_n2*nvec.z); //22
}

// diagonal terms

void
System::addToDiag(const vec3d &nvec, int ii, double alpha){
	int ii6 = 6*ii;
	
	double alpha_n0 = alpha*nvec.x;
	double alpha_n1 = alpha*nvec.y;
	double alpha_n2 = alpha*nvec.z;
	double alpha_n1n0 = alpha_n0*nvec.y;
	double alpha_n2n1 = alpha_n1*nvec.z;
	double alpha_n0n2 = alpha_n2*nvec.x;
	
	diag_values[ii6]   += alpha_n0*nvec.x; // 00
	diag_values[ii6+1] += alpha_n1n0; // 10
	diag_values[ii6+2] += alpha_n0n2; // 20
	
	diag_values[ii6+3] += alpha_n1*nvec.y; //11
	diag_values[ii6+4] += alpha_n2n1; // 21
	
	diag_values[ii6+5] += alpha_n2*nvec.z; //22
}


void
System::addStokesDrag(){
	for (int i = 0; i < np; i ++){
		int i6=6*i;
		double d_value = bgf_factor*radius[i];
		diag_values[i6  ] = d_value;
		diag_values[i6+3] = d_value;
		diag_values[i6+5] = d_value;
	}
}

#ifdef CHOLMOD
/*************** Cholmod Matrix Filling *************
Cholmod matrices we are using are defined in column major order (index j is column index)

Cholmod matrices are defined as follows:
- all values are stored in array x ( size nzmax )
- locations of values are encoded in array p ( size np ):
  values corresponding to column j are x[ p[j] ]  to x[ p[j+1] - 1 ]
- corresponding rows are stored in array i ( size nzmax ):
  rows corresponding to column j are i[ p[j] ]  to i[ p[j+1] - 1 ]

Hence: 
with p[j] < a < p[j+1]-1  
           . . . . j . . . . . .
        .|         .            |
        .|         .            |
        .|         .            |  
     i[a]| . . . .x[a]          |
        .|                      |
        .|                      |


*****************************************************/
void
System::fillSparseResmatrix(){
	

	allocateSparseResmatrix();
	
	// fill
	for(int j = 0; j < np; j++){
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
	((int*)sparse_res->p)[np3] = ((int*)sparse_res->p)[np3-1] + 1;
}
#endif

#ifdef TRILINOS
/*************** Epetra_CrsMatrix Filling *************
Epetra_CrsMatrix we are using are defined in row major order.

Epetra_CrsMatrix elements are not accessed directly for filling.
Instead we use user friendly methods, that take one row at a time.

*****************************************************/

void
System::fillSparseResmatrix(){
	
	allocateSparseResmatrix();
	
	int * diag_indices_j3 = new int [3];
	int * diag_indices_j3_1 = new int [2];
	int * diag_indices_j3_2;

	for(int j = 0; j < np; j++){
		int j3 = 3*j;
		int j6 = 6*j;
		

	    /******* diagonal-block values ******/

		diag_indices_j3[0] = j3;
		diag_indices_j3[1] = j3+1;
		diag_indices_j3[2] = j3+2;
		diag_indices_j3_1[0] = j3+1;
		diag_indices_j3_1[1] = j3+2;
		*diag_indices_j3_2 = j3+2;


		sparse_res->InsertMyValues(j3   , 3 , &(diag_values[j6])   , diag_indices_j3  );
		sparse_res->InsertMyValues(j3+1 , 2 , &(diag_values[j6+3]) , diag_indices_j3_1);
		sparse_res->InsertMyValues(j3+1 , 1 , &(diag_values[j6+5]) , diag_indices_j3_2); 


	    /******* off-diagonal-block values ******/
		int nb_entries = ploc[j+1] - ploc[j];
		

		sparse_res->InsertMyValues(j3   , nb_entries , &(off_diag_values[ 0 ][ ploc[j] ]) , &(rows[ ploc[j] ]) );
		sparse_res->InsertMyValues(j3+1 , nb_entries , &(off_diag_values[ 1 ][ ploc[j] ]) , &(rows[ ploc[j] ]) );
		sparse_res->InsertMyValues(j3+1 , nb_entries , &(off_diag_values[ 2 ][ ploc[j] ]) , &(rows[ ploc[j] ]) ); 
	}
	sparse_res->FillComplete();

}
#endif


// We solve A*(U-Uinf) = Gtilde*Einf ( in Jeffrey's notations )
// This method computes elements of matrix A and vector Gtilde*Einf
void
System::buildLubricationTerms(){
	cerr << " Now building new matrix " << endl;
	for (int k = 0; k < 6*np; k++){
		diag_values[k] = 0.;
	}
	rows.clear();
	off_diag_values[0].clear();
	off_diag_values[1].clear();
	off_diag_values[2].clear();
	
	addStokesDrag();
	
	double XAii, XAjj, XAij, XAji;
	//	double XGii, XGjj, XGij, XGji;
	double GEi[3];
	double GEj[3];
	
	set<Interaction*>::iterator it;
	int j;
	Interaction *inter;
	for (int i = 0; i < np - 1; i ++){
		ploc[i] = (unsigned int)rows.size();
		for (it = interaction_list[i].begin() ; it != interaction_list[i].end(); it ++){
			inter=*it;
		 	j=inter->partner(i);
			if(j>i){
				inter->XA(XAii, XAij, XAji, XAjj);
				// (i, j) (k,l) --> res[ n3*(3*i+l) + 3*j+k ]
				addToDiag(inter->nr_vec, i, inter->a0 * XAii);
				addToDiag(inter->nr_vec, j, inter->a1 * XAjj);
				appendToColumn(inter->nr_vec, j, 0.5 * inter->ro * XAji);
				inter->GE(GEi, GEj);  // G*E_\infty term
				for(int u=0; u<3; u++){
					// ((double*)lubrication_rhs->x)[ 3*i + u ] += GEi[ u ];
					// ((double*)lubrication_rhs->x)[ 3*j + u ] += GEj[ u ];
#ifdef CHOLMOD
					((double*)total_rhs->x)[ 3*i + u ] += GEi[ u ];
					((double*)total_rhs->x)[ 3*j + u ] += GEj[ u ];
#endif
				}
			}
		}
	}
	
	ploc[np-1] = (unsigned int)rows.size();
	ploc[np] = (unsigned int)rows.size();
}

#ifdef CHOLMOD
void
System::updateResistanceMatrix(){
		cerr << " Now updating matrix " << endl;
	size_t nrow = np3;
	size_t ncol = 1;
	size_t nzmax_od = 6;
	size_t nzmax_d = 3;

	int sort = 1;
	int pack = 1;
	int styp = 0;
	cholmod_sparse *update_vector_off_diag = cholmod_allocate_sparse ( nrow, ncol, nzmax_od, sort, pack, styp, xtype, &c); 
	cholmod_sparse *update_vector_diag = cholmod_allocate_sparse ( nrow, ncol, nzmax_d, sort, pack, styp, xtype, &c); ; // sparse packed vector. it has to be sorted. 
	//	for(int i=0; i<nrow; i++){
	((int*)update_vector_off_diag->p)[0] = 0;
	((int*)update_vector_off_diag->p)[1] = nzmax_od;
	((int*)update_vector_diag->p)[0] = 0;
	((int*)update_vector_diag->p)[1] = nzmax_d;


  	double XAii, XAjj, XAij, XAji;
  	double prev_XAii, prev_XAjj, prev_XAij, prev_XAji;
	int update;  // true for update, false for downdate
	double off_diag_diff;
	double sqrt_off_diag_diff;
	double diag_diff_comp_i, diag_diff_comp_j;
	double sqrt_diag_diff_comp_i, sqrt_diag_diff_comp_j;

	double nvec[3];

	for (int k = 0; k < num_interaction; k++){
		if( interaction[k].active || interaction[k].just_off() ){
			
			int i = interaction[k].particle_num[0];
			int j = interaction[k].particle_num[1];
			int i3 = 3*i;
			int j3 = 3*j; 

			nvec[0] = interaction[k].nr_vec.x;
			nvec[1] = interaction[k].nr_vec.y;
			nvec[2] = interaction[k].nr_vec.z;


			// DOWNDATE PART
			interaction[k].prev_XA(prev_XAii, prev_XAij, prev_XAji, prev_XAjj);
			// UPDATE PART
			interaction[k].XA(XAii, XAij, XAji, XAjj);


			off_diag_diff = 0.5 * interaction[k].ro * ( XAji - prev_XAji );

			if ( off_diag_diff > 0. ){
				update = 1;
				sqrt_off_diag_diff = sqrt( off_diag_diff );
			}
			else{
				update = 0;
				sqrt_off_diag_diff = sqrt( - off_diag_diff );
			}
			for(int u=0; u<3; u++) {
				((int*)update_vector_off_diag->i)[u] = i3+u;
				((double*)update_vector_off_diag->x)[u] = sqrt_off_diag_diff * nvec[u] ;
				((int*)update_vector_off_diag->i)[3+u] = j3+u;
				((double*)update_vector_off_diag->x)[3+u] = sqrt_off_diag_diff * nvec[u] ;
			}
			
			cholmod_updown(update, update_vector_off_diag, L, &c);

			// Now compensate the terms introduced on the diagonal, and add the diagonal update

			if(update){
				diag_diff_comp_i = interaction[k].a0 * ( XAii - prev_XAii ) - off_diag_diff ;
				diag_diff_comp_j = interaction[k].a1 * ( XAjj - prev_XAjj ) - off_diag_diff ;
			}
			else{
				diag_diff_comp_i = interaction[k].a0 * ( XAii - prev_XAii ) + off_diag_diff ;
				diag_diff_comp_j = interaction[k].a1 * ( XAjj - prev_XAjj ) + off_diag_diff ;
			}
			
			// first i
			if ( diag_diff_comp_i > 0. ){
				update = 1;
				sqrt_diag_diff_comp_i = sqrt( diag_diff_comp_i );
			}
			else{
				update = 0;
				sqrt_diag_diff_comp_i = sqrt( - diag_diff_comp_i );
			}


			for(int u=0; u<3; u++) {
				((int*)update_vector_diag->i)[u] = i3+u;
				((double*)update_vector_diag->x)[u] = sqrt_diag_diff_comp_i * nvec[u] ;
			}
			
			cholmod_updown(update, update_vector_diag, L, &c);

			// second j
			if ( diag_diff_comp_i > 0. ){
				update = 1;
				sqrt_diag_diff_comp_i = sqrt( diag_diff_comp_i );
			}
			else{
				update = 0;
				sqrt_diag_diff_comp_i = sqrt( - diag_diff_comp_i );
			}


			for(int u=0; u<3; u++) {
				((int*)update_vector_diag->i)[u] = j3+u;
				((double*)update_vector_diag->x)[u] = sqrt_diag_diff_comp_j * nvec[u] ;
			}

			cholmod_updown(update, update_vector_diag, L, &c);

		}
	}
}

#endif

// void
// System::calcLubricationForce(){
// 	for (int k = 0; k < 6*np; k++){
// 		diag_values[k] = 0.;
// 	}
// 	off_diag_values[0].clear();
// 	off_diag_values[1].clear();
// 	off_diag_values[2].clear();

// 	double GEi[3];
// 	double GEj[3];
// 	double XAii, XAjj, XAij, XAji;

// 	set<Interaction*>::iterator it;
// 	int j;
// 	Interaction *inter;
// 	for (int i = 0; i < np - 1; i ++){
// 		ploc[i] = (unsigned int)rows.size();
// 		for (it = interaction_list[i].begin() ; it != interaction_list[i].end(); it ++){
// 			inter=*it;
// 		 	j=inter->partner(i);
// 			if(j>i){
// 				inter->XA(XAii, XAij, XAji, XAjj);
// 				// (i, j) (k,l) --> res[ n3*(3*i+l) + 3*j+k ]
// 				addToDiag(inter->nr_vec, i, inter->a0 * XAjj);
// 				addToDiag(inter->nr_vec, j, inter->a1 * XAjj);
// 				appendToColumn(inter->nr_vec, j, 0.5 * inter->ro * XAji);
// 				inter->GE(GEi, GEj);  // G*E_\infty term
// 				for(int u=0; u<3; u++){
// 					// ((double*)lubrication_rhs->x)[ 3*i + u ] += GEi[ u ];
// 					// ((double*)lubrication_rhs->x)[ 3*j + u ] += GEj[ u ];
// 					((double*)total_rhs->x)[ 3*i + u ] += GEi[ u ];
// 					((double*)total_rhs->x)[ 3*j + u ] += GEj[ u ];
// 				}
// 			}
// 		}
// 	}
// 	for (int i=0; i < np; i++)
// 		lub_force[i].reset();
	
// 	for (int i=0; i < np; i++){
// 		for (int j=0; j < np; i++){
			
			
			
// 		}
// 	}
	
	
// }

#ifdef CHOLMOD
void
System::buildContactTerms(){
	// add contact force
	for (int i = 0; i < np; i++){
		int i3 = 3*i;
		// ((double*)contact_rhs->x)[i3] += contact_force[i].x;
		// ((double*)contact_rhs->x)[i3+1] += contact_force[i].y;
		// ((double*)contact_rhs->x)[i3+2] += contact_force[i].z;
		((double*)total_rhs->x)[i3] += contact_force[i].x;
		((double*)total_rhs->x)[i3+1] += contact_force[i].y;
		((double*)total_rhs->x)[i3+2] += contact_force[i].z;
	}
}
#endif
#ifdef CHOLMOD
void
System::allocateSparseResmatrix(){
	// allocate
	int nzmax;  // non-zero values
	nzmax = 6*np; // diagonal blocks
	for(int s=0; s<3; s++){
		nzmax += off_diag_values[s].size();  // off-diagonal
	}
	sparse_res = cholmod_allocate_sparse(np3, np3, nzmax, sorted, packed, stype,xtype, &c);
}
#endif
#ifdef TRILINOS
void
System::allocateSparseResmatrix(){
	// allocate
	int * nz = new int [linalg_size];  // non-zero values
	
	for(int i=0;i<linalg_size;i++){
	  nz[i] = dof; // diagonal
	  nz[i] += ploc[i+1]-1 - ploc[i];  // off-diagonal
	}
	sparse_res = rcp( new Epetra_CrsMatrix(Copy, *Map, nz));
	delete [] nz;
}
#endif


void
System::print_res(){ // testing
	cout << " Diag " << endl;
	for(int i=0;i<np;i++){
		int ii6=6*i;
		cout << i << " " << diag_values[ii6] << " " << diag_values[ii6+1]<< " " << diag_values[ii6+2]<< " " << diag_values[ii6+3]<< " " << diag_values[ii6+4]<< " " << diag_values[ii6+5] << endl;
	}
	
	cout << endl<< " OffDiag " << endl;
	for(int i=0;i<off_diag_values[0].size();i++){
		cout << off_diag_values[0][i] << " " << off_diag_values[1][i] << " " << off_diag_values[2][i]<< endl;
	}
}

#ifdef CHOLMOD
void
System::updateVelocityLubrication(){
	// contact_rhs = cholmod_zeros(n3, 1, xtype, &c);
	// lubrication_rhs = cholmod_zeros(n3, 1, xtype, &c);
	total_rhs = cholmod_zeros(np3, 1, xtype, &c);

	if( L == NULL ) {
		buildLubricationTerms();
		fillSparseResmatrix();
		factorizeResistanceMatrix();
	}
	else{
	  updateResistanceMatrix();
	}

	buildContactTerms();

	v = cholmod_solve (CHOLMOD_A, L, total_rhs, &c) ;
	// v_lub = cholmod_solve (CHOLMOD_A, L, lubrication_rhs, &c) ;
	// v_cont = cholmod_solve (CHOLMOD_A, L, contact_rhs, &c) ;
	
	/* TEST IMPLEMENTATION
	 * SDFF : Stokes drag force factor:
	 * SDFF = 1.0 : full drag forces from the undisturbed background flow.
	 * SDFF = 0.0 : no drag force from the undisturbed background flow.
	 */
	for (int i = 0; i < np; i++){
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
		for (int i=0; i < np; i++){
			ang_velocity[i] = 1.33333*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
	}
	
	cholmod_free_sparse(&sparse_res, &c);
	//	cholmod_free_factor(&L, &c);
	cholmod_free_dense(&total_rhs, &c);
	// cholmod_free_dense(&lubrication_rhs, &c);
	// cholmod_free_dense(&contact_rhs, &c);
	cholmod_free_dense(&v, &c);
}
#endif

#ifdef CHOLMOD
void
System::factorizeResistanceMatrix(){
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

}
#endif

#ifdef CHOLMOD
void System::updateVelocityLubricationBrownian(){
	total_rhs = cholmod_zeros(np3, 1, xtype, &c);
	buildLubricationTerms();
	
	fillSparseResmatrix();

	factorizeResistanceMatrix();

	buildContactTerms();
	
	v_nonBrownian = cholmod_solve (CHOLMOD_A, L, total_rhs, &c) ;

	// now the Brownian part of the velocity:
	// mid-point algortithm (see Melrose & Ball), modified (intermediate tstep) a la Banchio & Brady

	brownian_rhs = fb->generate();

	v_Brownian_init = cholmod_solve (CHOLMOD_A, L, brownian_rhs, &c) ;
	
	// move particles to intermediate point
	for (int i=0; i < np; i++){
		int i3 = 3*i;
		displacement(i, ((double*)v_Brownian_init->x)[i3]*dt_mid, ((double*)v_Brownian_init->x)[i3+1]*dt_mid, ((double*)v_Brownian_init->x)[i3+2]*dt_mid);
	}
	
	updateInteractions();
	
	// rebuild new R_FU
	cholmod_free_factor(&L, &c);
	cholmod_free_sparse(&sparse_res, &c);
	buildLubricationTerms();
	fillSparseResmatrix();

	factorizeResistanceMatrix();
	
	// get the intermediate brownian velocity
	v_Brownian_mid = cholmod_solve (CHOLMOD_A, L, brownian_rhs, &c) ;
	
	// move particles back to initial point, and update interactions
	//
	// Note that, although it looks like a complete reversal of the initial move (as it should be), 
	// the final state we obtain can be slightly different than the initial one, as the 1st move's update of the interaction
	// might switch off some of them. The 2nd move's update is not able to switch them back on.
	for (int i=0; i < np; i++){
		int i3 = 3*i;
		displacement(i, -((double*)v_Brownian_init->x)[i3]*dt_mid, -((double*)v_Brownian_init->x)[i3+1]*dt_mid, -((double*)v_Brownian_init->x)[i3+2]*dt_mid);
	}
	updateInteractions();

	// update total velocity
	// first term is hydrodynamic + contact velocities
	// second term is Brownian velocities
	// third term is Brownian drift
	// fourth term for vx is the shear rate
	for (int i = 0; i < np; i++){
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
		for (int i=0; i < np; i++){
			ang_velocity[i] = 1.33333*torque[i];
			ang_velocity[i].y += O_inf_y;
		}
	}
	
	cholmod_free_sparse(&sparse_res, &c);
	cholmod_free_factor(&L, &c);
	cholmod_free_dense(&total_rhs, &c);
	cholmod_free_dense(&v_nonBrownian, &c);
	cholmod_free_dense(&v_Brownian_init, &c);
	cholmod_free_dense(&v_Brownian_mid, &c);
}
#endif

#ifdef CHOLMOD
void System::computeBrownianStress(){
	double stresslet_i[5];
	double stresslet_j[5];
	for (int k=0; k < 5; k ++){
		stresslet_i[k] = 0.;
		stresslet_j[k] = 0.;
	}
	double vi [3];
	double vj [3];

	total_rhs = cholmod_zeros(np3, 1, xtype, &c);
	buildLubricationTerms();
	
	fillSparseResmatrix();

	factorizeResistanceMatrix();

	// now the Brownian part of the velocity:
	// mid-point algortithm (see Melrose & Ball), modified (intermediate tstep) a la Banchio & Brady

	brownian_rhs = fb->generate();

	v_Brownian_init = cholmod_solve (CHOLMOD_A, L, brownian_rhs, &c) ;
	
	for (int k = 0; k < num_interaction; k++){
		for (int u=0; u < 5; u ++){
			stresslet_i[u] = 0.;
			stresslet_j[u] = 0.;
		}

		int i = interaction[k].particle_num[0];
		int j = interaction[k].particle_num[1];
		int i3 = 3*i;
		int j3 = 3*j;
		vi[0] = ((double*)v_Brownian_init->x)[i3  ];
		vi[1] = ((double*)v_Brownian_init->x)[i3+1];
		vi[2] = ((double*)v_Brownian_init->x)[i3+2];
		
		vj[0] = ((double*)v_Brownian_init->x)[j3  ];
		vj[1] = ((double*)v_Brownian_init->x)[j3+1];
		vj[2] = ((double*)v_Brownian_init->x)[j3+2];
		interaction[k].pairStresslet(vi, vj, stresslet_i, stresslet_j);

		for (int u=0; u < 5; u++){
			brownianstress[i][u] += stresslet_i[u];
			brownianstress[j][u] += stresslet_j[u];
		}
	}
	

	// move particles to intermediate point
	for (int i=0; i < np; i++){
		int i3 = 3*i;
		displacement(i, ((double*)v_Brownian_init->x)[i3]*dt_mid, ((double*)v_Brownian_init->x)[i3+1]*dt_mid, ((double*)v_Brownian_init->x)[i3+2]*dt_mid);
	}
	
	updateInteractions();
	
	// rebuild new R_FU
	cholmod_free_factor(&L, &c);
	cholmod_free_sparse(&sparse_res, &c);
	buildLubricationTerms();
	fillSparseResmatrix();

	factorizeResistanceMatrix();
	
	// get the intermediate brownian velocity
	v_Brownian_mid = cholmod_solve (CHOLMOD_A, L, brownian_rhs, &c) ;
	
	for (int k = 0; k < num_interaction; k++){
		for (int u=0; u < 5; u ++){
			stresslet_i[u] = 0.;
			stresslet_j[u] = 0.;
		}

		int i = interaction[k].particle_num[0];
		int j = interaction[k].particle_num[1];
		int i3 = 3*i;
		int j3 = 3*j;
		vi[0] = ((double*)v_Brownian_mid->x)[i3  ];
		vi[1] = ((double*)v_Brownian_mid->x)[i3+1];
		vi[2] = ((double*)v_Brownian_mid->x)[i3+2];
		
		vj[0] = ((double*)v_Brownian_mid->x)[j3  ];
		vj[1] = ((double*)v_Brownian_mid->x)[j3+1];
		vj[2] = ((double*)v_Brownian_mid->x)[j3+2];
		interaction[k].pairStresslet(vi, vj, stresslet_i, stresslet_j);

		for (int u=0; u < 5; u++){
			brownianstress[i][u] -= stresslet_i[u];
			brownianstress[j][u] -= stresslet_j[u];
		}
	}

	for (int i=0; i < np; i++){
		for (int u=0; u < 5; u++){
			brownianstress[i][u] *= 0.5*dt_ratio;
		}
	}
	
	// move particles back to initial point, and update interactions
	//
	// Note that, although it looks like a complete reversal of the initial move (as it should be), 
	// the final state we obtain can be slightly different than the initial one, as the 1st move's update of the interaction
	// might switch off some of them. The 2nd move's update is not able to switch them back on.
	for (int i=0; i < np; i++){
		int i3 = 3*i;
		displacement(i, -((double*)v_Brownian_init->x)[i3]*dt_mid, -((double*)v_Brownian_init->x)[i3+1]*dt_mid, -((double*)v_Brownian_init->x)[i3+2]*dt_mid);
	}
	updateInteractions();

	cholmod_free_sparse(&sparse_res, &c);
	cholmod_free_factor(&L, &c);
	cholmod_free_dense(&total_rhs, &c);
	cholmod_free_dense(&v_nonBrownian, &c);
	cholmod_free_dense(&v_Brownian_init, &c);
	cholmod_free_dense(&v_Brownian_mid, &c);
}
#endif

void
System::displacement(int i, const double &dx_, const double &dy_, const double &dz_){
	position[i].x += dx_;
	position[i].y += dy_;
	position[i].z += dz_;
	periodize( position[i]);
	boxset->box(i);
}


// [0,l]
void
System::periodize(vec3d &pos){
	if (pos.z > lz() ){
		pos.z -= lz();
		pos.x -= shear_disp;
	} else if ( pos.z < 0 ){
		pos.z += lz();
		pos.x += shear_disp;
	}
	while ( pos.x > lx() ){
		pos.x -= lx();
	}
	while (pos.x < 0 ){
		pos.x += lx();
	}
	if (dimension == 3){
		if ( pos.y > ly() ){
			pos.y -= ly();
		} else if (pos.y < 0 ){
			pos.y += ly();
		}
	}
}

// [-l/2,l/2]
void
System::periodize_diff(vec3d &pos_diff){
	if (pos_diff.z > lz2() ){
		pos_diff.z -= lz();
		pos_diff.x -= shear_disp;
	} else if ( pos_diff.z < -lz2() ){
		pos_diff.z += lz();
		pos_diff.x += shear_disp;
	}
	while ( pos_diff.x > lx2() ){
		pos_diff.x -= lx();
	}
	while (pos_diff.x < -lx2() ){
		pos_diff.x += lx();
	}
	if (dimension == 3){
		if ( pos_diff.y > ly2() ){
			pos_diff.y -= ly();
		} else if (pos_diff.y < -ly2() ){
			pos_diff.y += ly();
		}
	}
}

// periodize + give z_shift= number of boundaries crossed in z-direction
void
System::periodize_diff(vec3d &pos_diff, int &zshift){
	if (pos_diff.z > lz2() ){
		pos_diff.z -= lz();
		pos_diff.x -= shear_disp;
		zshift = -1;
	} else if ( pos_diff.z < -lz2() ){
		pos_diff.z += lz();
		pos_diff.x += shear_disp;
		zshift = +1;
	} else{
		zshift = 0;
	}
	while ( pos_diff.x > lx2() ){
		pos_diff.x -= lx();
	}
	while (pos_diff.x < -lx2() ){
		pos_diff.x += lx();
	}
	if (dimension == 3){
		if ( pos_diff.y > ly2() ){
			pos_diff.y -= ly();
		} else if (pos_diff.y < -ly2() ){
			pos_diff.y += ly();
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
	for (int i=0; i < np; i++){
		displacement(i, velocity[i].x*dt, velocity[i].y*dt, velocity[i].z*dt);
	}
	if (draw_rotation_2d){
		for (int i=0; i < np; i++){
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
	periodize_diff(pos_diff);
	if (dimension == 3){
		return pos_diff.sq_norm();
	} else {
	 	return pos_diff.sq_norm_xz();
	}
}

void
System::calcStress(){
	stressReset();

	if(brownian)
	  computeBrownianStress();

	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			interaction[k].addLubricationStress();
			if (interaction[k].contact){
				interaction[k].addContactStress();
			}
			interaction[k].evaluateLubricationForce();
		}
	}
	for(int k=0; k < 10 ; k++){
		cnt_contact_number[k] = 0;
	}
	for (int i =0; i < np; i++){
		cnt_contact_number[ contact_number[i] ] ++;
	}
	
	for (int k=0; k < 5; k++){
		total_lub_stress[k] = 0;
		total_contact_stress[k] = 0;
		total_brownian_stress[k] = 0;
	}
	
	
	for (int i=0; i < np; i++){
		for (int k=0; k < 5; k++){
			total_lub_stress[k] += lubstress[i][k];
			total_contact_stress[k] += contactstress[i][k];
			total_brownian_stress[k] += brownianstress[i][k];
		}
	}
	
	/*
	 * The term 5.0/9 is the one-body part
	 *
	 */
	total_stress_bgf = 0;
	for (int i=0; i < np; i++){
		double a = radius[i];
		total_stress_bgf += (5.0/9)*bgf_factor*a*a*a;
	}

}

void
System::analyzeState(){
	gap_min = lz();
	double sum_overlap = 0;
	int cnt_overlap = 0;
	max_age = 0;
	double sum_age = 0;
	int cnt_age = 0;
	
	for (int i=0; i < np; i++){
		contact_number[i] = 0;
	}
	
	total_contact = 0;
	// for analysis
	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			
			if (interaction[k].gap() < 0){
				sum_overlap +=interaction[k].gap();
				cnt_overlap ++;
			}
			
			if (interaction[k].gap() < gap_min){
				gap_min = interaction[k].gap();
			}
			
			if (interaction[k].age() > 0 ){
				if (max_age < interaction[k].age()){
					max_age = interaction[k].age();
				}
				sum_age = interaction[k].age();
				cnt_age ++;
			}
			interaction[k].recordTrajectory();
			
			if (interaction[k].near){
				total_contact ++;
				contact_number[interaction[k].particle_num[0]] ++;
				contact_number[interaction[k].particle_num[1]] ++;
			}
		}
	}
	ave_overlap = sum_overlap / cnt_overlap;
	ave_age = sum_age / cnt_age;

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

void
System::setSystemVolume(){
	if (dimension == 2){
		system_volume = _lx*_lz*2*radius_max;
	} else {
		system_volume = _lx*_ly*_lz;
	}
}







