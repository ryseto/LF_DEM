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
	if (!relative_velocity_lub_cont)
		delete [] relative_velocity_lub_cont;
	if (!relative_velocity_brownian)
		delete [] relative_velocity_brownian;

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

	if (!interaction_list)
		delete [] interaction_list;
	if (!interaction_partners)
		delete [] interaction_partners;
	if (!fb){
		delete [] fb;
	}
	
	if (!v_lub_cont)
		delete [] v_lub_cont;

	if(brownian){
	    if (!v_Brownian_init)
		delete [] v_Brownian_init;

	    if (!v_Brownian_mid)
		delete [] v_Brownian_mid;
	}

	delete stokes_solver;
};


void
System::allocateRessources(){
	position = new vec3d [np];
	radius = new double [np];
	angle = new double [np];
	velocity = new vec3d [np];
	relative_velocity = new vec3d [np];
	relative_velocity_lub_cont = new vec3d [np];
	relative_velocity_brownian = new vec3d [np];
	ang_velocity = new vec3d [np];
	total_force = new vec3d [np];
	lubrication_force = new vec3d [np];
	contact_force = new vec3d [np];
	brownian_force = new vec3d [np];
	torque = new vec3d [np];
	lub_force = new vec3d [np];
	nearing_number.resize(np);
	lubstress.resize(np);
	contactstress.resize(np);
	brownianstress.resize(np);


	int maxnum_interactionpair_per_particle = 15;
	maxnum_interactionpair = (int)(maxnum_interactionpair_per_particle*np);
	interaction = new Interaction [maxnum_interactionpair];
	interaction_list = new set <Interaction*> [np];
	interaction_partners = new set <int> [np];
	
	lubstress2.resize(np);
	
	dof = 3;
	linalg_size = dof*np;
	
	v_lub_cont = new double [linalg_size];
	if(brownian){
	    v_Brownian_init = new double [linalg_size];
	    v_Brownian_mid = new double [linalg_size];
	}
	stokes_solver = new StokesSolver(np, brownian);
	fb = new BrownianForce(this); 
}

void
System::setupSystem(const vector<vec3d> &initial_positions,
					const vector <double> &radii){

	if (kb_T == 0){
		brownian = false;
	} else {
		brownian = true;
	}

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
		relative_velocity_lub_cont[i].x=0.;
		relative_velocity_lub_cont[i].y=0.;
		relative_velocity_lub_cont[i].z=0.;
		relative_velocity_brownian[i].x=0.;
		relative_velocity_brownian[i].y=0.;
		relative_velocity_brownian[i].z=0.;
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

	stokes_solver->initialize();
	
	// initialize the brownian force after the solver, as it assumes
	// the cholmod_common of the solver is already initialized
	if(brownian){
		fb->init();
	}

	brownianstress_calc_nb = 0;
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
	boxset->update();
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
					// distance is done in 3 steps because we need each information for Interaction creation
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
		for (int u=0; u < 5; u++){
			lubstress[i].elm[u]=0;
			contactstress[i].elm[u]=0;
			lubstress2[i].elm[u]=0;
		}
	}
}
void
System::stressBrownianReset(){
	for (int i=0; i < np; i++){
		for (int u=0; u < 5; u++){
			brownianstress[i].elm[u]=0;
		}
	}
	brownianstress_calc_nb = 0;
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

void
System::addStokesDrag(){
    for (int i = 0; i < np; i ++){
	double d_value = bgf_factor*radius[i];
	stokes_solver->addToDiag_RFU(i, d_value);    
    }
}


// We solve A*(U-Uinf) = Gtilde*Einf ( in Jeffrey's notations )
// This method computes:
//  - elements of matrix A 
//  - vector Gtilde*Einf if rhs is true (default behavior)
void
System::buildLubricationTerms(bool rhs=true){
    
    double XAii, XAjj, XAij, XAji;

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
		  
		  stokes_solver->addToDiagBlock_RFU(inter->nr_vec, i, inter->a0 * XAii);
		  stokes_solver->addToDiagBlock_RFU(inter->nr_vec, j, inter->a1 * XAjj);
		  stokes_solver->appendToOffDiagBlock_RFU(inter->nr_vec, i, j, 0.5 * inter->ro * XAji);
		  
		  if(rhs){
		    inter->GE(GEi, GEj);  // G*E_\infty term
		    for(int u=0; u<3; u++){
			  stokes_solver->addToRHS( 3*i + u, GEi[u] );
			  stokes_solver->addToRHS( 3*j + u, GEj[u] );
		    }
		  }
		  
	    }
	  }
	  stokes_solver->doneBlocks(i);
    }
}

void
System::buildContactTerms(){
    // add contact force
    for (int i = 0; i < np; i++){
	int i3 = 3*i;
	stokes_solver->addToRHS(i3  , contact_force[i].x);
	stokes_solver->addToRHS(i3+1, contact_force[i].y);
	stokes_solver->addToRHS(i3+2, contact_force[i].z);
    }
}



void
System::addToDiagBlock_RFU(cholmod_sparse *chol_rfu_submatrix, const vec3d &nvec, int ii, double alpha){
    double alpha_n0 = alpha*nvec.x;
    double alpha_n1 = alpha*nvec.y;
    double alpha_n2 = alpha*nvec.z;
    double alpha_n1n0 = alpha_n0*nvec.y;
    double alpha_n2n1 = alpha_n1*nvec.z;
    double alpha_n0n2 = alpha_n2*nvec.x;
    
	int pii3_0;
	int pii3_1;
	int pii3_2;

	if(ii==1){
	  pii3_0 = 15;
	  pii3_1 = 18;
	  pii3_2 = 20;
	}
	if(ii==0){
	  pii3_0 = 0;
	  pii3_1 = 6;
	  pii3_2 = 11;
	}



	((double*)chol_rfu_submatrix->x)[pii3_0  ] += alpha_n0*nvec.x; 
	((double*)chol_rfu_submatrix->x)[pii3_0+1] += alpha_n1n0; // 10
	((double*)chol_rfu_submatrix->x)[pii3_0+2] += alpha_n0n2; // 20
	
	((double*)chol_rfu_submatrix->x)[pii3_1  ] += alpha_n1*nvec.y; //11
	((double*)chol_rfu_submatrix->x)[pii3_1+1] += alpha_n2n1; // 21
	
  	((double*)chol_rfu_submatrix->x)[pii3_2  ]  += alpha_n2*nvec.z; //22

}


void
System::addToOffDiagBlock_RFU(cholmod_sparse *chol_rfu_submatrix, const vec3d &nvec, double alpha){


    double alpha_n0 = alpha*nvec.x;
    double alpha_n1 = alpha*nvec.y;
    double alpha_n2 = alpha*nvec.z;
    double alpha_n1n0 = alpha_n0*nvec.y;
    double alpha_n2n1 = alpha_n1*nvec.z;
    double alpha_n0n2 = alpha_n2*nvec.x;
    

	int p0 = 3;
	int p1 = 8;
	int p2 = 12;

    
	((double*)chol_rfu_submatrix->x)[p0  ] = alpha_n0*nvec.x; // 00
    ((double*)chol_rfu_submatrix->x)[p0+1] = alpha_n1n0; // 10
    ((double*)chol_rfu_submatrix->x)[p0+2] = alpha_n0n2; // 20
	((double*)chol_rfu_submatrix->x)[p1  ] = alpha_n1n0; // 01
	((double*)chol_rfu_submatrix->x)[p1+1] = alpha_n1*nvec.y; //11
	((double*)chol_rfu_submatrix->x)[p1+2] = alpha_n2n1; // 21
	((double*)chol_rfu_submatrix->x)[p2  ] = alpha_n0n2; // 02
	((double*)chol_rfu_submatrix->x)[p2+1] = alpha_n2n1; // 12
	((double*)chol_rfu_submatrix->x)[p2+2] = alpha_n2*nvec.z; //22
}

void 
System::buildBrownianTerms(double *force){

    double XAii, XAjj, XAij, XAji;

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

		  //		  cout << endl<<" part " << i << " " << j << endl;
		  cholmod_sparse *chol_rfu_submatrix;
		  cholmod_dense *chol_pair_force;
		  cholmod_dense *rand_vec;
		  cholmod_factor *chol_L ;
		  cholmod_factor *L_copy ;
		  cholmod_sparse *L_sparse ;
		  cholmod_common chol_c = stokes_solver->chol_c;
		  
		  //		  cholmod_start (&chol_c) ;
		  
		  int sorted = 1;
		  int packed = 1;
		  int stype = -1;
		  chol_rfu_submatrix = cholmod_allocate_sparse(6, 6, 21, sorted, packed, stype,CHOLMOD_REAL, &chol_c);
		  rand_vec = cholmod_zeros(6, 1, CHOLMOD_REAL, &chol_c);
		  chol_pair_force = cholmod_zeros(6, 1, CHOLMOD_REAL, &chol_c);
		  
		  //		  cout << "allocated "<< endl;
		  
		  ((int*)chol_rfu_submatrix->p)[0]=0;
		  ((int*)chol_rfu_submatrix->p)[1]=6;
		  ((int*)chol_rfu_submatrix->p)[2]=11;
		  ((int*)chol_rfu_submatrix->p)[3]=15;
		  ((int*)chol_rfu_submatrix->p)[4]=18;
		  ((int*)chol_rfu_submatrix->p)[5]=20;
		  
		  int loc=0;
		  for(int a=0; a<6; a++){
			for(int b=a; b<6; b++){
			  ((int*)chol_rfu_submatrix->i)[loc]=b;
			  ((double*)chol_rfu_submatrix->x)[loc]=0;
			  loc++;
			}
		  }

		  // Stokes drag
		  double d_value = bgf_factor*radius[i];
		  ((double*)chol_rfu_submatrix->x)[0]=d_value;
		  ((double*)chol_rfu_submatrix->x)[6]=d_value;
		  ((double*)chol_rfu_submatrix->x)[11]=d_value;
		  d_value = bgf_factor*radius[j];
		  ((double*)chol_rfu_submatrix->x)[15]=d_value;
		  ((double*)chol_rfu_submatrix->x)[18]=d_value;
		  ((double*)chol_rfu_submatrix->x)[20]=d_value;



		  inter->XA(XAii, XAij, XAji, XAjj);
		  
		  addToDiagBlock_RFU(chol_rfu_submatrix, inter->nr_vec, 0, inter->a0 * XAii);
		  addToDiagBlock_RFU(chol_rfu_submatrix, inter->nr_vec, 1, inter->a1 * XAjj);
		  addToOffDiagBlock_RFU(chol_rfu_submatrix, inter->nr_vec, 0.5 * inter->ro * XAji);

		  // cout << endl<< " chol subrfu " << endl;
		  // for(int a = 0; a < 6; a++){
		  // 	for(int k =((int*)chol_rfu_submatrix->p)[a] ; k < ((int*)chol_rfu_submatrix->p)[a+1]; k++){
		  // 	  cout << a << " " << ((int*)chol_rfu_submatrix->i)[k] << " " << ((double*)chol_rfu_submatrix->x)[k] << endl;
		  // 	}
		  // }

		  double ** L_dense;
		  L_dense = new double* [6];
		  for(int a=0; a<6; a++){
			L_dense[a] = new double [6];
			for(int b=a; b<6; b++){
			  L_dense[a][b]=0.;
			}
		  }

		  double aux, ad;
		  //		  cout << endl;
		  for(int a=0; a<6; a++){
			int pa = ((int*)chol_rfu_submatrix->p)[a];
			for(int b=a; b<6; b++){
			  aux = ((double*)chol_rfu_submatrix->x)[ pa + b-a ]; 
			  //			  cout << a << " " << b << " " << aux  << endl;
			  for(int c = a-1;c>=0; c--){
				aux = aux - L_dense[a][c]*L_dense[b][c];
				//				cout << a << " " << b << " " <<c << " " << aux  << endl;
			  }
			  if (a == b){
				if (aux <= 0){
				  if ( aux > 1e-10 ){
					aux = 0.;
				  }
				  cout << "choldc failed" << endl;
				  cout <<  aux << " " << a << endl;
				  exit(1);
				}
				ad = sqrt(aux);
				L_dense[a][a] = ad;
			  }
			  else{
				L_dense[b][a] = aux/ad;
			  }
			  //			  cout << " L " << b << " " << a <<  " " << L_dense[b][a] << endl;
			}
		  }


		  // L_dense test
		  // cout << " chol test sys " << endl;
		  // for(int a=0; a<6; a++){
		  // 	int pa = ((int*)chol_rfu_submatrix->p)[a];
		  // 	for(int b=a; b<6; b++){
		  // 	  cout << a << " " << b << " " << ((double*)chol_rfu_submatrix->x)[ pa + b-a ] << " "; 
		  // 	  aux =0.;
		  // 	  for(int c=0; c<=a; c++){
		  // 		aux += L_dense[a][c]*L_dense[b][c];
		  // 	  }
		  // 	  cout << aux << endl;
		  // 	}
		  // }
		  // getchar();
		  //		  chol_L = cholmod_analyze (chol_rfu_submatrix, &chol_c);
		  //		  cout << "analyzed "<< endl;
		  //		  cholmod_factorize (chol_rfu_submatrix, chol_L, &chol_c);
		  //		  cout << "factorizeed "<< endl;
		  //		  L_copy = cholmod_copy_factor(chol_L, &chol_c); // sadly it seems we have to make a copy. Is there a way to avoid this?
		  //		  L_sparse = cholmod_factor_to_sparse(L_copy, &chol_c);
		  //		  cout << "copied "<< endl;
		  
		  for(int a=0; a<6; a++){
			//			((double*)rand_vec->x)[a] = r_gen.randNorm(0., 1.);
			((double*)rand_vec->x)[a] = 1.;
		  }
		  

		  //		  cout << " force sys " << endl;
		  double sqrt_kbT2_dt = sqrt(2*kb_T/dt);
		  for(int a=0; a<3; a++){
			for(int b=0; b<=a; b++){
			  force[3*i+a] += sqrt_kbT2_dt*L_dense[a][b]*((double*)rand_vec->x)[b];
			  //			  cout << a << " " << b << " "<< L_dense[a][b] << " "<< force[3*i+a]<< endl;		
			}
			//		  cout << sqrt_kbT2_dt*((double*)rand_vec->x)[a] << endl;
			//		  			  cout <<  force[3*i+a] << endl;
		  }
		  for(int a=3; a<6; a++){
			for(int b=0; b<=a; b++){
			  force[3*j+a-3] += sqrt_kbT2_dt*L_dense[a][b]*((double*)rand_vec->x)[b];
			  //			  cout << a << " " << b << " "<< L_dense[a][b] << " "<< force[3*j+a-3] << endl;		
			}
			//		  cout << sqrt_kbT2_dt*((double*)rand_vec->x)[a] << endl;
			//		  cout <<  force[3*j+a-3] << endl;
		  }

		  for(int a=0; a<6; a++)
			delete [] L_dense[a];
		  delete [] L_dense;


		  // double sqrt_temp2_dt [2] = {sqrt(2*kb_T/dt),0};
		  // double zero [2] = {0,0};
		  // cholmod_sdmult(L_sparse, 0, sqrt_temp2_dt, zero, rand_vec, chol_pair_force, &chol_c);
		  // cout << "multiplied "<< endl;

		  // for(int a=0; a<3; a++){
		  // 	force[3*i+a] += ((double*)chol_pair_force->x)[a  ];
		  // 	force[3*j+a] += ((double*)chol_pair_force->x)[a+3];
		  // }

		  
		  
		  cholmod_free_dense(&chol_pair_force, &chol_c);
		  cholmod_free_dense(&rand_vec, &chol_c);
		  //		  cholmod_free_factor(&chol_L, &chol_c);
		  //		  cholmod_free_factor(&L_copy, &chol_c);
		  cholmod_free_sparse(&chol_rfu_submatrix, &chol_c);
		  //		  cholmod_free_sparse(&L_sparse, &chol_c);
		  
	    }
	  }
    }
	

}

void
System::updateVelocityLubrication(){

    stokes_solver->resetRHS();
    stokes_solver->prepareNewBuild_RFU("direct");
//	stokes_solver->prepareNewBuild_RFU("iterative");

    addStokesDrag();
    buildLubricationTerms();

    stokes_solver->complete_RFU();

    buildContactTerms();
	
    stokes_solver->solve(v_lub_cont);
	
	/* TEST IMPLEMENTATION
	 * SDFF : Stokes drag force factor:
	 * SDFF = 1.0 : full drag forces from the undisturbed background flow.
	 * SDFF = 0.0 : no drag force from the undisturbed background flow.
	 */
    for (int i = 0; i < np; i++){
	int i3 = 3*i;
	relative_velocity_lub_cont[i].x = v_lub_cont[i3];
	relative_velocity_lub_cont[i].y = v_lub_cont[i3+1];
	relative_velocity_lub_cont[i].z = v_lub_cont[i3+2];
	
	velocity[i].x = relative_velocity_lub_cont[i].x + shear_rate*position[i].z;
	velocity[i].y = relative_velocity_lub_cont[i].y;
	velocity[i].z = relative_velocity_lub_cont[i].z;
    }
	
    if(friction){
	double O_inf_y = 0.5*shear_rate;
	for (int i=0; i < np; i++){
	    ang_velocity[i] = 1.33333*torque[i];
	    ang_velocity[i].y += O_inf_y;
	}
    }
    
    stokes_solver->solvingIsDone();
}



// On the way, we compute the Brownian Stress, as this needs a lot of averaging to converge
void System::updateVelocityLubricationBrownian(){

    stresslet stresslet_i_init;
    stresslet stresslet_j_init;
    stresslet stresslet_i_mid;
    stresslet stresslet_j_mid;
    stresslet *step_stresslet = new stresslet [np];

    for (int i=0; i < np; i++){
	  for (int u=0; u < 5; u++){
	    step_stresslet[i].elm[u] = 0.;
	  }
    }

	
    vec3d vi;
    vec3d vj;


    stokes_solver->resetRHS();
    stokes_solver->prepareNewBuild_RFU("direct");

    addStokesDrag();
    buildLubricationTerms();

    stokes_solver->complete_RFU();

    buildContactTerms();
	
    stokes_solver->solve(v_lub_cont);


    // now the Brownian part of the velocity:
    // mid-point algortithm (see Melrose & Ball), modified (intermediate tstep) a la Banchio & Brady
    // 
    // we do not call solvingIsDone() before new solve(), because 
    // R_FU has not changed, so same factorization is safely used


 	 // double *bforce = new double [np3];
	 // for(int i=0; i<np3; i++){
	 //   bforce[i] = 0.;
	 // }
	 // buildBrownianTerms(bforce);

	 // fb->generate();

	 // stokes_solver->setRHS( bforce );

	 // delete [] bforce;
	stokes_solver->setRHS( fb->generate() );

	stokes_solver->solve( v_Brownian_init );

	// Brownian Stress: term R_SU * v_Brownian_init
    for (int k = 0; k < num_interaction; k++){
	  for (int u=0; u < 5; u ++){
 		stresslet_i_init.elm[u] = 0.;
		stresslet_j_init.elm[u] = 0.;
	  }
	  
	  int i = interaction[k].particle_num[0];
	  int j = interaction[k].particle_num[1];
	  int i3 = 3*i;
	  int j3 = 3*j;
	  
	  vi.x = v_Brownian_init[i3  ];
	  vi.y = v_Brownian_init[i3+1];
	  vi.z = v_Brownian_init[i3+2];
	  
	  vj.x = v_Brownian_init[j3  ];
	  vj.y = v_Brownian_init[j3+1];
	  vj.z = v_Brownian_init[j3+2];
	  
	  interaction[k].pairVelocityStresslet(vi, vj, stresslet_i_init, stresslet_j_init);
	  
	  for (int u=0; u < 5; u++){
	    step_stresslet[i].elm[u] += stresslet_i_init.elm[u];
	    step_stresslet[j].elm[u] += stresslet_j_init.elm[u];
	  }
    }


    // move particles to intermediate point
    for (int i=0; i < np; i++){
	  int i3 = 3*i;
	  displacement(i, v_Brownian_init[i3]*dt_mid, v_Brownian_init[i3+1]*dt_mid, v_Brownian_init[i3+2]*dt_mid);
    }
    updateInteractions();
	
    // build new Resistance matrix after move
	stokes_solver->solvingIsDone();
	//	stokes_solver->convertDirectToIterative();
	//  stokes_solver->prepareNewBuild_RFU("iterative");
    stokes_solver->prepareNewBuild_RFU("direct");
    addStokesDrag();
    buildLubricationTerms(false); // false: don't modify rhs
    stokes_solver->complete_RFU();

    // get the intermediate brownian velocity
    stokes_solver->solve(v_Brownian_mid);
    stokes_solver->solvingIsDone();


	// Brownian Stress: term  -R_SU_mid * v_Brownian_mid
    for (int k = 0; k < num_interaction; k++){
	  for (int u=0; u < 5; u ++){
		stresslet_i_mid.elm[u] = 0.;
		stresslet_j_mid.elm[u] = 0.;
	  }
	  
	  int i = interaction[k].particle_num[0];
	  int j = interaction[k].particle_num[1];
	  int i3 = 3*i;
	  int j3 = 3*j;
	  
	  vi.x = v_Brownian_mid[i3  ];
	  vi.y = v_Brownian_mid[i3+1];
	  vi.z = v_Brownian_mid[i3+2];
	  
	  vj.x = v_Brownian_mid[j3  ];
	  vj.y = v_Brownian_mid[j3+1];
	  vj.z = v_Brownian_mid[j3+2];
	  
	  interaction[k].pairVelocityStresslet(vi, vj, stresslet_i_mid, stresslet_j_mid);
	  
	  for (int u=0; u < 5; u++){
	    step_stresslet[i].elm[u] -= stresslet_i_mid.elm[u];
	    step_stresslet[j].elm[u] -= stresslet_j_mid.elm[u];
	  }
    }
	// Brownian Stress: multiply all by  0.5*dt/dt_mid
    for (int i=0; i < np; i++){
	  for (int u=0; u < 5; u++){
	    step_stresslet[i].elm[u] *= 0.5*dt_ratio;
		brownianstress[i].elm[u] += step_stresslet[i].elm[u];
	  }
    }
	brownianstress_calc_nb ++;

    // move particles back to initial point, and update interactions
    //
    // Note that, although it looks like a complete reversal of the initial move (as it should be), 
    // the final state we obtain can be slightly different than the initial one, as the 1st move's update of the interaction
    // might switch off some of them. The 2nd move's update is not able to switch them back on.
    for (int i=0; i < np; i++){
	  int i3 = 3*i;
	  displacement(i, -v_Brownian_init[i3]*dt_mid, -v_Brownian_init[i3+1]*dt_mid, -v_Brownian_init[i3+2]*dt_mid);
    }
    updateInteractions();

    // update total velocity
    // first term is hydrodynamic + contact velocities
    // second term is Brownian velocities
    // third term is Brownian drift
    // fourth term for vx is the shear rate
    for (int i = 0; i < np; i++){
	  int i3 = 3*i;
	  relative_velocity_lub_cont[i].x = v_lub_cont[i3];
	  relative_velocity_lub_cont[i].y = v_lub_cont[i3+1];
	  relative_velocity_lub_cont[i].z = v_lub_cont[i3+2];
	  
	  relative_velocity_brownian[i].x = v_Brownian_init[i3] + 0.5*dt_ratio*( v_Brownian_mid[i3] - v_Brownian_init[i3] );
	  relative_velocity_brownian[i].y = v_Brownian_init[i3+1] + 0.5*dt_ratio*( v_Brownian_mid[i3+1] - v_Brownian_init[i3+1] ); 
	  relative_velocity_brownian[i].z = v_Brownian_init[i3+2] + 0.5*dt_ratio*( v_Brownian_mid[i3+2] - v_Brownian_init[i3+2] );

	  velocity[i].x = relative_velocity_lub_cont[i].x + relative_velocity_brownian[i].x + shear_rate*position[i].z;
	  velocity[i].y = relative_velocity_lub_cont[i].y + relative_velocity_brownian[i].y;
	  velocity[i].z = relative_velocity_lub_cont[i].z + relative_velocity_brownian[i].z;
    }
	
    if(friction){
	double O_inf_y = 0.5*shear_rate;
	for (int i=0; i < np; i++){
	    ang_velocity[i] = 1.33333*torque[i];
	    ang_velocity[i].y += O_inf_y;
	}
    }



 	// vb_avg = 0;
	// vb_avg_nb = 0;
	
	// for(int i=0; i < np; i++){
	//   vb_avg += relative_velocity_lub_cont[i].sq_norm();
	//   vb_avg_nb++;
	// }
	// cout << " vlub_avg " << sqrt(vb_avg)/vb_avg_nb << " ";




	// vb_avg = 0;
	// vb_avg_nb = 0;
	
	// for(int i=0; i < np; i++){
	//   vb_avg += relative_velocity_brownian[i].sq_norm();
	//   vb_avg_nb++;
	// }
	// cout << " vb_avg " << sqrt(vb_avg)/vb_avg_nb << " ";

	delete [] step_stresslet;

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
	
	ksi_avg = 1.;
	ksi_avg_nb = 0;
	checkNewInteraction();
	updateInteractions();
	//	cout << " ksi_min " << ksi_avg << endl;

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

	for (int k = 0; k < num_interaction; k++){
		if (interaction[k].active){
			interaction[k].evaluateLubricationForce();
			interaction[k].addLubricationStress();
			if (interaction[k].contact){
				interaction[k].addContactStress();
			}
		}
	}
	for(int u=0; u < 10 ; u++){
		cnt_nearing_number[u] = 0;
	}
	for (int i =0; i < np; i++){
		cnt_nearing_number[ nearing_number[i] ] ++;
	}
	
	for (int u=0; u < 5; u++){
		total_lub_stress[u] = 0;
		total_contact_stress[u] = 0;
		total_brownian_stress[u] = 0;
	}
	
	
	for (int i=0; i < np; i++){
		for (int u=0; u < 5; u++){
			total_lub_stress[u] += lubstress[i].elm[u];
			total_contact_stress[u] += contactstress[i].elm[u];
			total_brownian_stress[u] += brownianstress[i].elm[u]/brownianstress_calc_nb;
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

	stressBrownianReset();
}

void
System::analyzeState(){
	gap_min = lz();
	double sum_overlap = 0;
	int cnt_overlap = 0;
	max_nearing_time = 0;
	
	for (int i=0; i < np; i++){
		nearing_number[i] = 0;
	}
	
	num_nearing = 0;
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

			interaction[k].recordTrajectory();
			
			if (interaction[k].near){
				num_nearing ++;
				nearing_number[interaction[k].particle_num[0]] ++;
				nearing_number[interaction[k].particle_num[1]] ++;
				if ( max_nearing_time < interaction[k].nearingTime()){
					max_nearing_time = interaction[k].nearingTime();
				}
			}
		}
	}
	
	double sum_nearing_time = 0;
	for(int k = 0; k < nearing_time_record.size(); k++){
		sum_nearing_time += nearing_time_record[k];
	}
	ave_overlap = sum_overlap / cnt_overlap;
	if (nearing_time_record.size() > 0){
		ave_nearing_time = sum_nearing_time / nearing_time_record.size();
	} else {
		ave_nearing_time = 0;
	}
	nearing_time_record.clear();
	
	
	
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

