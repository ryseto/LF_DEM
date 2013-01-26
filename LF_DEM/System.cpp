//
//  System.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 11/14/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include "System.h"
#include <sstream>
#ifdef TRILINOS
#include <BelosCGIteration.hpp>
#endif
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
	
	if (!v_lub_cont)
		delete [] v_lub_cont;

#ifdef CHOLMOD
	if (!diag_values)
		delete [] diag_values;
	if (!off_diag_values )
		delete [] off_diag_values;
	if (!ploc)
		delete [] ploc;


	cholmod_free_dense(&chol_v_lub_cont, &chol_c);
	cholmod_free_dense(&chol_rhs_lub_cont, &chol_c);
	//	cholmod_free_dense(&chol_nonbrownian_rhs, &chol_c); 
	if (brownian){
		cholmod_free_dense(&chol_brownian_rhs, &chol_c);
	}
	//cholmod_free_dense(&chol_contact_rhs, &chol_c);
	//cholmod_free_dense(&chol_lubrication_rhs, &chol_c);
	cholmod_free_sparse(&chol_rfu_matrix, &chol_c);
	cholmod_finish(&chol_c);
#endif

#ifdef TRILINOS
	for (int i=0; i < linalg_size; i++)
		delete [] columns[i];
	delete [] columns;

	delete [] columns_nb;

	for (int i=0; i < linalg_size; i++)
		delete [] values[i];
	delete [] values;	
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

	int maxnum_interactionpair_per_particle = 12;
	maxnum_interactionpair = (int)(maxnum_interactionpair_per_particle*np);
	interaction = new Interaction [maxnum_interactionpair];
	interaction_list = new set <Interaction*> [np];
	interaction_partners = new set <int> [np];
	
	dof = 3;
	linalg_size = dof*np;
	
	
#ifdef TRILINOS
	columns_max_nb = dof*maxnum_interactionpair_per_particle;
	int numlhs = 1;
	int numrhs = 1;
	Map = rcp(new Epetra_Map(linalg_size, 0, Comm));
	tril_v_lub_cont = rcp( new VEC(*Map, numlhs) );
	tril_rhs_lub_cont = rcp( new VEC(*Map, numrhs) );
	//	tril_rfu_matrix = rcp( new MAT(linalg_size) );
	tril_rfu_matrix = rcp( new Epetra_CrsMatrix(Copy, *Map, maxnum_interactionpair_per_particle*dof + dof ));
	tril_stokes_equation = rcp( new Belos::LinearProblem < SCAL, VEC, MAT > ( tril_rfu_matrix, tril_v_lub_cont, tril_rhs_lub_cont ) ) ;

	columns = new int* [linalg_size];
	for (int i=0; i < linalg_size; i++){
		columns[i] = new int [columns_max_nb];
		for(int j=0; j < columns_max_nb; j++){
		  columns[i][j] = -1;
		}
	}
	values = new double* [linalg_size];
	for (int i=0; i < linalg_size; i++){
		values[i] = new double [columns_max_nb];
		for(int j=0; j < columns_max_nb; j++){
		  values[i][j] = 0.;
		}
	}
	columns_nb = new int [linalg_size];
	for (int i=0; i < linalg_size; i++)
		columns_nb[i] = 0;


#endif
#ifdef CHOLMOD
	diag_values = new double [6*np];
	off_diag_values = new vector <double> [3];
	ploc = new int [np+1];
	v_lub_cont = new double [linalg_size];
#endif
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
	cholmod_start (&chol_c) ;
	stype = -1; // 1 is symmetric, stored upper triangular (UT), -1 is LT
	sorted = 0;		/* TRUE if columns sorted, FALSE otherwise*/
	packed = 1;		/* TRUE if matrix packed, FALSE otherwise */
	xtype = CHOLMOD_REAL;
	chol_L=NULL;
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
	tril_solver = tril_factory.create ("CG", solverParams);
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
#ifdef CHOLMOD
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
#endif

//off-diagonal terms
#ifdef TRILINOS
void
System::appendToRow(const vec3d &nvec, int ii, int jj, double alpha){
	int ii3   = 3*ii;
	int ii3_1 = ii3+1;
	int ii3_2 = ii3+2;
	int jj3   = 3*jj;
	int jj3_1 = jj3+1;
	int jj3_2 = jj3+2;

	double alpha_n0 = alpha*nvec.x;
	double alpha_n1 = alpha*nvec.y;
	double alpha_n2 = alpha*nvec.z;
	double alpha_n1n0 = alpha_n0*nvec.y;
	double alpha_n2n1 = alpha_n1*nvec.z;
	double alpha_n0n2 = alpha_n2*nvec.x;
	

	// declare ii and jj new columns, and update column nb
	int last_col_nb = columns_nb[ii3];

	columns[ii3  ][last_col_nb  ] = jj3  ;
	columns[ii3  ][last_col_nb+1] = jj3_1;
	columns[ii3  ][last_col_nb+2] = jj3_2;
	columns[ii3_1][last_col_nb  ] = jj3  ;
	columns[ii3_1][last_col_nb+1] = jj3_1;
	columns[ii3_1][last_col_nb+2] = jj3_2;
	columns[ii3_2][last_col_nb  ] = jj3  ;
	columns[ii3_2][last_col_nb+1] = jj3_1;
	columns[ii3_2][last_col_nb+2] = jj3_2;

	columns[jj3  ][last_col_nb  ] = ii3  ;
	columns[jj3  ][last_col_nb+1] = ii3_1;
	columns[jj3  ][last_col_nb+2] = ii3_2;
	columns[jj3_1][last_col_nb  ] = ii3  ;
	columns[jj3_1][last_col_nb+1] = ii3_1;
	columns[jj3_1][last_col_nb+2] = ii3_2;
	columns[jj3_2][last_col_nb  ] = ii3  ;
	columns[jj3_2][last_col_nb+1] = ii3_1;
	columns[jj3_2][last_col_nb+2] = ii3_2;

	columns_nb[ii3] += 3;
	columns_nb[ii3_1] += 3;
	columns_nb[ii3_2] += 3;
	columns_nb[jj3] += 3;
	columns_nb[jj3_1] += 3;
	columns_nb[jj3_2] += 3;

	// set values	
	values[ii3  ][last_col_nb  ] = alpha_n0*nvec.x; // 00
	values[ii3  ][last_col_nb+1] = alpha_n1n0;      // 01
	values[ii3  ][last_col_nb+2] = alpha_n0n2;      // 02
	values[ii3_1][last_col_nb  ] = alpha_n1n0;      // 10
	values[ii3_1][last_col_nb+1] = alpha_n1*nvec.y; // 11
	values[ii3_1][last_col_nb+2] = alpha_n2n1;      // 12
	values[ii3_2][last_col_nb  ] = alpha_n0n2;      // 20
	values[ii3_2][last_col_nb+1] = alpha_n2n1;      // 21
	values[ii3_2][last_col_nb+2] = alpha_n2*nvec.z;      // 22

	values[jj3  ][last_col_nb  ] = alpha_n0*nvec.x; // 00
	values[jj3  ][last_col_nb+1] = alpha_n1n0;      // 01
	values[jj3  ][last_col_nb+2] = alpha_n0n2;      // 02
	values[jj3_1][last_col_nb  ] = alpha_n1n0;      // 10
	values[jj3_1][last_col_nb+1] = alpha_n1*nvec.y; // 11
	values[jj3_1][last_col_nb+2] = alpha_n2n1;      // 12
	values[jj3_2][last_col_nb  ] = alpha_n0n2;      // 20
	values[jj3_2][last_col_nb+1] = alpha_n2n1;      // 21
	values[jj3_2][last_col_nb+2] = alpha_n2*nvec.z;      // 22

}
#endif



// diagonal terms

void
System::addToDiag(const vec3d &nvec, int ii, double alpha){

	
	double alpha_n0 = alpha*nvec.x;
	double alpha_n1 = alpha*nvec.y;
	double alpha_n2 = alpha*nvec.z;
	double alpha_n1n0 = alpha_n0*nvec.y;
	double alpha_n2n1 = alpha_n1*nvec.z;
	double alpha_n0n2 = alpha_n2*nvec.x;
	
#ifdef CHOLMOD
	int ii6 = 6*ii;
	diag_values[ii6]   += alpha_n0*nvec.x; // 00
	diag_values[ii6+1] += alpha_n1n0; // 10
	diag_values[ii6+2] += alpha_n0n2; // 20
	
	diag_values[ii6+3] += alpha_n1*nvec.y; //11
	diag_values[ii6+4] += alpha_n2n1; // 21
	
	diag_values[ii6+5] += alpha_n2*nvec.z; //22
#endif
#ifdef TRILINOS
	int iidof = dof*ii;
	values[iidof  ][0] += alpha_n0*nvec.x; // 00
	values[iidof  ][1] += alpha_n1n0; // 01
	values[iidof  ][2] += alpha_n0n2; // 02

	values[iidof+1][0] += alpha_n1n0; // 10
	values[iidof+1][1] += alpha_n1*nvec.y; // 11
	values[iidof+1][2] += alpha_n2n1; // 12

	values[iidof+2][0] += alpha_n0n2; // 20
	values[iidof+2][1] += alpha_n2n1; // 21
	values[iidof+2][2] += alpha_n2*nvec.z; // 20
#endif
}

#ifdef CHOLMOD
void
System::addStokesDrag(){
	for (int i = 0; i < np; i ++){
		int i6=6*i;
		double d_value = bgf_factor*radius[i];
		diag_values[i6  ] += d_value;
		diag_values[i6+3] += d_value;
		diag_values[i6+5] += d_value;
	}
}
#endif

#ifdef TRILINOS
void
System::addStokesDrag(){
	for (int i = 0; i < np; i ++){
		int idof=dof*i;
		double d_value = bgf_factor*radius[i];
		for (int j = 0; j < dof; j ++){
			values[idof+j][j] = d_value;
		}
	}
}
#endif


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
		
		((int*)chol_rfu_matrix->p)[j3  ] = j6   + 3*ploc[j];
		((int*)chol_rfu_matrix->p)[j3+1] = j6+3 + 2*ploc[j] +   ploc[j+1];
		((int*)chol_rfu_matrix->p)[j3+2] = j6+5 +   ploc[j] + 2*ploc[j+1];
		
		int pj3   = ((int*)chol_rfu_matrix->p)[j3];
		int pj3_1 = ((int*)chol_rfu_matrix->p)[j3+1];
		int pj3_2 = ((int*)chol_rfu_matrix->p)[j3+2];
		
		// diagonal blocks row indices
		((int*)chol_rfu_matrix->i)[ pj3 ]       = j3;
		((int*)chol_rfu_matrix->i)[ pj3 + 1 ]   = j3+1;
		((int*)chol_rfu_matrix->i)[ pj3 + 2 ]   = j3+2;
		
		((int*)chol_rfu_matrix->i)[ pj3_1 ]     = j3+1;
		((int*)chol_rfu_matrix->i)[ pj3_1 + 1 ] = j3+2;
		
		((int*)chol_rfu_matrix->i)[ pj3_2 ]     = j3+2;
		
		// diagonal blocks row values
		((double*)chol_rfu_matrix->x)[ pj3 ]       = diag_values[j6];
		((double*)chol_rfu_matrix->x)[ pj3 + 1 ]   = diag_values[j6+1];
		((double*)chol_rfu_matrix->x)[ pj3 + 2 ]   = diag_values[j6+2];
		
		((double*)chol_rfu_matrix->x)[ pj3_1 ]     = diag_values[j6+3];
		((double*)chol_rfu_matrix->x)[ pj3_1 + 1 ] = diag_values[j6+4];
		
		((double*)chol_rfu_matrix->x)[ pj3_2 ]     = diag_values[j6+5];
		
		//    cout << j3+2 <<" " << diag_values[j6+5]<< " " << ((int*)chol_rfu_matrix->p)[j3] << " " << ((int*)chol_rfu_matrix->p)[j3+1] << " " << ((int*)chol_rfu_matrix->p)[j3+2] <<endl;
		// off-diagonal blocks row indices and values
		for(int k = ploc[j]; k < ploc[j+1]; k++){
			int u = k - ploc[j];
			((int*)chol_rfu_matrix->i)[ pj3   + u + 3 ] = rows[k];
			((int*)chol_rfu_matrix->i)[ pj3_1 + u + 2 ] = rows[k];
			((int*)chol_rfu_matrix->i)[ pj3_2 + u + 1 ] = rows[k];
			
			((double*)chol_rfu_matrix->x)[ pj3   + u + 3 ] = off_diag_values[0][k];
			((double*)chol_rfu_matrix->x)[ pj3_1 + u + 2 ] = off_diag_values[1][k];
			((double*)chol_rfu_matrix->x)[ pj3_2 + u + 1 ] = off_diag_values[2][k];
		}
	}
	((int*)chol_rfu_matrix->p)[np3] = ((int*)chol_rfu_matrix->p)[np3-1] + 1;
}
#endif

#ifdef TRILINOS
/*************** Epetra_CrsMatrix Filling *************
Epetra_CrsMatrix we are using are defined in row major order.

Epetra_CrsMatrix must be stored completely.

Epetra_CrsMatrix elements are not accessed directly for filling.
Instead we use user friendly methods, that take one row at a time.

*****************************************************/

void
System::fillSparseResmatrix(){
	
	allocateSparseResmatrix();
	
	int * diag_indices_j3 = new int [3];
	int * diag_indices_j3_1 = new int [2];
	int diag_indices_j3_2;

	for(int i = 0; i < linalg_size; i++){
		tril_rfu_matrix->InsertGlobalValues(i, columns_nb[i] , values[i], columns[i]);
	}
	tril_rfu_matrix->FillComplete();

}
#endif


// We solve A*(U-Uinf) = Gtilde*Einf ( in Jeffrey's notations )
// This method computes elements of matrix A and vector Gtilde*Einf
#ifdef TRILINOS
void
System::buildLubricationTerms(){
	for (int i = 0; i < linalg_size; i++){
		for (int j = 0; j < columns_max_nb; j++){
			columns[i][j] = -1;
			values[i][j] = 0.;
		}
	}

	// declare the diagonal blocks
	for (int i = 0; i < np; i++){
		int idof = dof*i;
		for (int j = 0; j < dof; j++){
			columns[idof  ][j] = idof+j;
			columns[idof+1][j] = idof+j;
			columns[idof+2][j] = idof+j;
		}
		columns_nb[idof  ] = dof;
		columns_nb[idof+1] = dof;
		columns_nb[idof+2] = dof;
	}

	addStokesDrag();
	
	double XAii, XAjj, XAij, XAji;
	//	double XGii, XGjj, XGij, XGji;
	double GEi[3];
	double GEj[3];
	
	set<Interaction*>::iterator it;
	int j;
	Interaction *inter;
	for (int i = 0; i < np - 1; i ++){
		for (it = interaction_list[i].begin() ; it != interaction_list[i].end(); it ++){
			inter=*it;
		 	j=inter->partner(i);
			if(j>i){
				inter->XA(XAii, XAij, XAji, XAjj);
				// (i, j) (k,l) --> res[ n3*(3*i+l) + 3*j+k ]
				addToDiag(inter->nr_vec, i, inter->a0 * XAii);
				addToDiag(inter->nr_vec, j, inter->a1 * XAjj);
				appendToRow(inter->nr_vec, i, j, 0.5 * inter->ro * XAji);
				inter->GE(GEi, GEj);  // G*E_\infty term
				for(int u=0; u<3; u++){
					tril_rhs_lub_cont->SumIntoMyValue( 3*i + u, 0, GEi[ u ]);
					tril_rhs_lub_cont->SumIntoMyValue( 3*j + u, 0, GEj[ u ]);
				}
			}
		}
	}

}
#endif

// We solve A*(U-Uinf) = Gtilde*Einf ( in Jeffrey's notations )
// This method computes elements of matrix A and vector Gtilde*Einf
#ifdef CHOLMOD
void
System::buildLubricationTerms(){
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
					((double*)chol_rhs_lub_cont->x)[ 3*i + u ] += GEi[ u ];
					((double*)chol_rhs_lub_cont->x)[ 3*j + u ] += GEj[ u ];
				}
			}
		}
	}


	ploc[np-1] = (unsigned int)rows.size();
	ploc[np] = (unsigned int)rows.size();
}


void
System::updateResistanceMatrix(){
	size_t nrow = np3;
	size_t ncol = 1;
	size_t nzmax_od = 6;
	size_t nzmax_d = 3;

	int sort = 1;
	int pack = 1;
	int styp = 0;
	cholmod_sparse *update_vector_off_diag = cholmod_allocate_sparse ( nrow, ncol, nzmax_od, sort, pack, styp, xtype, &chol_c); 
	cholmod_sparse *update_vector_diag = cholmod_allocate_sparse ( nrow, ncol, nzmax_d, sort, pack, styp, xtype, &chol_c); ; // sparse packed vector. it has to be sorted. 

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
			
			cholmod_updown(update, update_vector_off_diag, chol_L, &chol_c);

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
			
			cholmod_updown(update, update_vector_diag, chol_L, &chol_c);

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

			cholmod_updown(update, update_vector_diag, chol_L, &chol_c);

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

void
System::buildContactTerms(){
	// add contact force
	for (int i = 0; i < np; i++){
		int i3 = 3*i;
#ifdef CHOLMOD
		((double*)chol_rhs_lub_cont->x)[i3] += contact_force[i].x;
		((double*)chol_rhs_lub_cont->x)[i3+1] += contact_force[i].y;
		((double*)chol_rhs_lub_cont->x)[i3+2] += contact_force[i].z;
#endif
#ifdef TRILINOS
		tril_rhs_lub_cont->SumIntoMyValue( 3*i    , 0, contact_force[i].x);
		tril_rhs_lub_cont->SumIntoMyValue( 3*i + 1, 0, contact_force[i].y);
		tril_rhs_lub_cont->SumIntoMyValue( 3*i + 2, 0, contact_force[i].z);
#endif
	}
}

#ifdef CHOLMOD
void
System::allocateSparseResmatrix(){
	// allocate
	int nzmax;  // non-zero values
	nzmax = 6*np; // diagonal blocks
	for(int s=0; s<3; s++){
		nzmax += off_diag_values[s].size();  // off-diagonal
	}
	chol_rfu_matrix = cholmod_allocate_sparse(np3, np3, nzmax, sorted, packed, stype,xtype, &chol_c);
}
#endif
#ifdef TRILINOS
void
System::allocateSparseResmatrix(){
	// // allocate
	// int * nz = new int [linalg_size];  // non-zero values
	
	// for(int i=0;i<linalg_size;i++){
	//   nz[i] = dof; // diagonal
	//   nz[i] += ploc[i+1]-1 - ploc[i];  // off-diagonal
	// }
	// tril_rfu_matrix = rcp( new Epetra_CrsMatrix(Copy, *Map, nz));
	// delete [] nz;
}
#endif

#ifdef CHOLMOD
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
#endif

void
System::updateVelocityLubrication(){
  
#ifdef CHOLMOD
	chol_rhs_lub_cont = cholmod_zeros(np3, 1, xtype, &chol_c);
#endif
#ifdef TRILINOS
	tril_rhs_lub_cont->PutScalar(0.);
	tril_rfu_matrix->PutScalar(0.);
#endif

	buildLubricationTerms();
	fillSparseResmatrix();
	

	/***** if you want to switch to update procedure for CHOLMOD (SLOW!),
		   replace above 3 lines by the following
	if( chol_L == NULL ) {
		buildLubricationTerms();
		fillSparseResmatrix();
		factorizeResistanceMatrix();
	}
	else{
	  updateResistanceMatrix();
	}
	*****/ 

	buildContactTerms();

#ifdef CHOLMOD
	factorizeResistanceMatrix();
	chol_v_lub_cont = cholmod_solve (CHOLMOD_A, chol_L, chol_rhs_lub_cont, &chol_c) ;
	v_lub_cont = (double*)chol_v_lub_cont->x;
#endif
#ifdef TRILINOS
	bool set_success = tril_stokes_equation->setProblem();
	if(!set_success){
		cerr << "ERROR:  Belos::LinearProblem failed to set up correctly" << endl;
		exit(1);
	}
	tril_solver->setProblem (tril_stokes_equation);
	Belos::ReturnType ret = tril_solver->solve();
	if( ret != Belos::Converged )
		cerr << " Warning: Belos::Solver did not converge" << endl;
	int iter_steps = tril_solver->getNumIters();
	cout << " iterations " << iter_steps << endl;
	tril_v_lub_cont->ExtractCopy(&v_lub_cont);
#endif
	
	/* TEST IMPLEMENTATION
	 * SDFF : Stokes drag force factor:
	 * SDFF = 1.0 : full drag forces from the undisturbed background flow.
	 * SDFF = 0.0 : no drag force from the undisturbed background flow.
	 */
	for (int i = 0; i < np; i++){
		int i3 = 3*i;
		relative_velocity[i].x = v_lub_cont[i3];
		relative_velocity[i].y = v_lub_cont[i3+1];
		relative_velocity[i].z = v_lub_cont[i3+2];
		
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
#ifdef CHOLMOD
	cholmod_free_factor(&chol_L, &chol_c);
	cholmod_free_sparse(&chol_rfu_matrix, &chol_c);
	cholmod_free_dense(&chol_rhs_lub_cont, &chol_c);
	cholmod_free_dense(&chol_v_lub_cont, &chol_c);
#endif
	
}

#ifdef CHOLMOD
void
System::factorizeResistanceMatrix(){
	chol_L = cholmod_analyze (chol_rfu_matrix, &chol_c);
	cholmod_factorize (chol_rfu_matrix, chol_L, &chol_c);
	if(chol_c.status){
		// Cholesky decomposition has failed: usually because matrix is incorrectly found to be positive-definite
		// It is very often enough to force another preconditioner to solve the problem.
		cerr << " factorization failed. forcing simplicial algorithm... " << endl;
		chol_c.supernodal = CHOLMOD_SIMPLICIAL;
		chol_L = cholmod_analyze (chol_rfu_matrix, &chol_c);
		cholmod_factorize (chol_rfu_matrix, chol_L, &chol_c) ;
		cerr << " factorization status " << chol_c.status << " final_ll ( 0 is LDL, 1 is LL ) " <<  chol_c.final_ll <<endl;
		chol_c.supernodal = CHOLMOD_SUPERNODAL;
	}

}
#endif


void System::updateVelocityLubricationBrownian(){
#ifdef CHOLMOD
	chol_rhs_lub_cont = cholmod_zeros(np3, 1, xtype, &chol_c);
	buildLubricationTerms();
	
	fillSparseResmatrix();

	factorizeResistanceMatrix();

	buildContactTerms();
	
	chol_v_nonBrownian = cholmod_solve (CHOLMOD_A, chol_L, chol_rhs_lub_cont, &chol_c) ;

	// now the Brownian part of the velocity:
	// mid-point algortithm (see Melrose & Ball), modified (intermediate tstep) a la Banchio & Brady

	chol_brownian_rhs = fb->generate();

	chol_v_Brownian_init = cholmod_solve (CHOLMOD_A, chol_L, chol_brownian_rhs, &chol_c) ;
	
	// move particles to intermediate point
	for (int i=0; i < np; i++){
		int i3 = 3*i;
		displacement(i, ((double*)chol_v_Brownian_init->x)[i3]*dt_mid, ((double*)chol_v_Brownian_init->x)[i3+1]*dt_mid, ((double*)chol_v_Brownian_init->x)[i3+2]*dt_mid);
	}
	
	updateInteractions();
	
	// rebuild new R_FU
	cholmod_free_factor(&chol_L, &chol_c);
	cholmod_free_sparse(&chol_rfu_matrix, &chol_c);
	buildLubricationTerms();
	fillSparseResmatrix();

	factorizeResistanceMatrix();
	
	// get the intermediate brownian velocity
	chol_v_Brownian_mid = cholmod_solve (CHOLMOD_A, chol_L, chol_brownian_rhs, &chol_c) ;
	
	// move particles back to initial point, and update interactions
	//
	// Note that, although it looks like a complete reversal of the initial move (as it should be), 
	// the final state we obtain can be slightly different than the initial one, as the 1st move's update of the interaction
	// might switch off some of them. The 2nd move's update is not able to switch them back on.
	for (int i=0; i < np; i++){
		int i3 = 3*i;
		displacement(i, -((double*)chol_v_Brownian_init->x)[i3]*dt_mid, -((double*)chol_v_Brownian_init->x)[i3+1]*dt_mid, -((double*)chol_v_Brownian_init->x)[i3+2]*dt_mid);
	}
	updateInteractions();

	// update total velocity
	// first term is hydrodynamic + contact velocities
	// second term is Brownian velocities
	// third term is Brownian drift
	// fourth term for vx is the shear rate
	for (int i = 0; i < np; i++){
		int i3 = 3*i;
		relative_velocity[i].x = ((double*)chol_v_nonBrownian->x)[i3] + ((double*)chol_v_Brownian_init->x)[i3] + 0.5*dt_ratio*(((double*)chol_v_Brownian_mid->x)[i3] - ((double*)chol_v_Brownian_init->x)[i3] );
		relative_velocity[i].y = ((double*)chol_v_nonBrownian->x)[i3+1] + ((double*)chol_v_Brownian_init->x)[i3+1] + 0.5*dt_ratio*(((double*)chol_v_Brownian_mid->x)[i3+1] - ((double*)chol_v_Brownian_init->x)[i3+1] );
		relative_velocity[i].z = ((double*)chol_v_nonBrownian->x)[i3+2] + ((double*)chol_v_Brownian_init->x)[i3+2] + 0.5*dt_ratio*(((double*)chol_v_Brownian_mid->x)[i3+2] - ((double*)chol_v_Brownian_init->x)[i3+2] );
		
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
	
	cholmod_free_sparse(&chol_rfu_matrix, &chol_c);
	cholmod_free_factor(&chol_L, &chol_c);
	cholmod_free_dense(&chol_rhs_lub_cont, &chol_c);
	cholmod_free_dense(&chol_v_nonBrownian, &chol_c);
	cholmod_free_dense(&chol_v_Brownian_init, &chol_c);
	cholmod_free_dense(&chol_v_Brownian_mid, &chol_c);
#endif
}



void System::computeBrownianStress(){
#ifdef CHOLMOD
	double stresslet_i[5];
	double stresslet_j[5];
	for (int k=0; k < 5; k ++){
		stresslet_i[k] = 0.;
		stresslet_j[k] = 0.;
	}
	vec3d vi;
	vec3d vj;

	chol_rhs_lub_cont = cholmod_zeros(np3, 1, xtype, &chol_c);
	buildLubricationTerms();
	
	fillSparseResmatrix();

	factorizeResistanceMatrix();

	// now the Brownian part of the velocity:
	// mid-point algortithm (see Melrose & Ball), modified (intermediate tstep) a la Banchio & Brady

	chol_brownian_rhs = fb->generate();

	chol_v_Brownian_init = cholmod_solve (CHOLMOD_A, chol_L, chol_brownian_rhs, &chol_c) ;
	
	for (int k = 0; k < num_interaction; k++){
		for (int u=0; u < 5; u ++){
			stresslet_i[u] = 0.;
			stresslet_j[u] = 0.;
		}

		int i = interaction[k].particle_num[0];
		int j = interaction[k].particle_num[1];
		int i3 = 3*i;
		int j3 = 3*j;

		vi.x = ((double*)chol_v_Brownian_init->x)[i3  ];
		vi.y = ((double*)chol_v_Brownian_init->x)[i3+1];
		vi.z = ((double*)chol_v_Brownian_init->x)[i3+2];
		
		vj.x = ((double*)chol_v_Brownian_init->x)[j3  ];
		vj.y = ((double*)chol_v_Brownian_init->x)[j3+1];
		vj.z = ((double*)chol_v_Brownian_init->x)[j3+2];

		interaction[k].pairStresslet(vi, vj, stresslet_i, stresslet_j);

		for (int u=0; u < 5; u++){
			brownianstress[i][u] += stresslet_i[u];
			brownianstress[j][u] += stresslet_j[u];
		}
	}
	

	// move particles to intermediate point
	for (int i=0; i < np; i++){
		int i3 = 3*i;
		displacement(i, ((double*)chol_v_Brownian_init->x)[i3]*dt_mid, ((double*)chol_v_Brownian_init->x)[i3+1]*dt_mid, ((double*)chol_v_Brownian_init->x)[i3+2]*dt_mid);
	}
	
	updateInteractions();
	
	// rebuild new R_FU
	cholmod_free_factor(&chol_L, &chol_c);
	cholmod_free_sparse(&chol_rfu_matrix, &chol_c);
	buildLubricationTerms();
	fillSparseResmatrix();

	factorizeResistanceMatrix();
	
	// get the intermediate brownian velocity
	chol_v_Brownian_mid = cholmod_solve (CHOLMOD_A, chol_L, chol_brownian_rhs, &chol_c) ;
	
	for (int k = 0; k < num_interaction; k++){
		for (int u=0; u < 5; u ++){
			stresslet_i[u] = 0.;
			stresslet_j[u] = 0.;
		}

		int i = interaction[k].particle_num[0];
		int j = interaction[k].particle_num[1];
		int i3 = 3*i;
		int j3 = 3*j;

		vi.x = ((double*)chol_v_Brownian_mid->x)[i3  ];
		vi.y = ((double*)chol_v_Brownian_mid->x)[i3+1];
		vi.z = ((double*)chol_v_Brownian_mid->x)[i3+2];
		
		vj.x = ((double*)chol_v_Brownian_mid->x)[j3  ];
		vj.y = ((double*)chol_v_Brownian_mid->x)[j3+1];
		vj.z = ((double*)chol_v_Brownian_mid->x)[j3+2];

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
		displacement(i, -((double*)chol_v_Brownian_init->x)[i3]*dt_mid, -((double*)chol_v_Brownian_init->x)[i3+1]*dt_mid, -((double*)chol_v_Brownian_init->x)[i3+2]*dt_mid);
	}
	updateInteractions();

	cholmod_free_sparse(&chol_rfu_matrix, &chol_c);
	cholmod_free_factor(&chol_L, &chol_c);
	cholmod_free_dense(&chol_rhs_lub_cont, &chol_c);
	cholmod_free_dense(&chol_v_nonBrownian, &chol_c);
	cholmod_free_dense(&chol_v_Brownian_init, &chol_c);
	cholmod_free_dense(&chol_v_Brownian_mid, &chol_c);
#endif
}


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







