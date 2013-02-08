#include "BrownianForce.h"

#define GRANDOM ( r_gen.randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.

using namespace std;

// BrownianForce::BrownianForce(double kb_T, cholmod_dense *brownian_force_vector, cholmod_common *c_, System *sys_){
// 	sys = sys_;
// 	force = *brownian_force_vector;
// 	c = c_;// }

BrownianForce::BrownianForce(System *sys_){
	sys = sys_;
	kb_T = sys->kb_T;
	kb_T2= 2*kb_T;

	forces = new double [3*sys->numpart()];
	pair_resistance_matrix = new double* [6];
	L_factor = new double* [6];
	for(int i=0;i<6;i++){
	  pair_resistance_matrix[i] = new double [6];
	  L_factor[i] = new double [6];
	}
}

BrownianForce::~BrownianForce(){

  delete [] forces;

  for(int i=0;i<6;i++){
	delete [] pair_resistance_matrix[i];
	delete [] L_factor[i];
  }
  delete [] pair_resistance_matrix;
  delete [] L_factor;
}

void
BrownianForce::init(){

}


double*
BrownianForce::generate(){
    double XAii, XAjj, XAij, XAji;

    set<Interaction*>::iterator it;
    int j;
    Interaction *inter;
	int np = sys->numpart();
	int np3 = 3*np; 
	double rand_vec6 [6];
	
    for (int i = 0; i < np3; i ++){
	  forces[i] = 0.;
	}
    for (int i = 0; i < np - 1; i ++){

	  for (it = (sys->interaction_list[i]).begin() ; it != (sys->interaction_list[i]).end(); it++){
	    inter=*it;
	    j=inter->partner(i);

	    if(j>i){
		  for(int a=0; a<6; a++){
			for(int b=0; b<6; b++){
			  L_factor[a][b]=0.;
			  pair_resistance_matrix[a][b]=0.;
			}
		  }

		  // Stokes drag
		  double d_value = sys->bgf_factor*sys->radius[i];
		  for(int a=0; a<3; a++){
			pair_resistance_matrix[a][a]=d_value;
		  }
		  d_value = sys->bgf_factor*sys->radius[j];
		  for(int a=3; a<6; a++){
			pair_resistance_matrix[a][a]=d_value;
		  }

		  inter->XA(XAii, XAij, XAji, XAjj);
		  
		  addToDiagBlock_RFU(inter->nr_vec, 0, inter->a0 * XAii);
		  addToDiagBlock_RFU(inter->nr_vec, 1, inter->a1 * XAjj);
		  addToOffDiagBlock_RFU(inter->nr_vec, 0.5 * inter->ro * XAji);

		  factorize();
		  
		  double sqrt_kbT2_dt = sqrt(kb_T2/sys->dt);

		  for(int a=0; a<6; a++){
		    rand_vec6[a] = sqrt_kbT2_dt * r_gen.randNorm(0., 1.);
		  }
		  

		  for(int a=0; a<3; a++){
			for(int b=0; b<=a; b++){
			  forces[3*i+a] += L_factor[a][b]*rand_vec6[b];
			}
		  }
		  for(int a=3; a<6; a++){
			for(int b=0; b<=a; b++){
			  forces[3*j+a-3] += L_factor[a][b]*rand_vec6[b];
			}
		  }

		}
	  }
    }
	// for(int i=0;i<np3;i++){
	//   cout << forces[i] << endl;
	// }
	return forces;
}



void
BrownianForce::addToDiagBlock_RFU(const vec3d &nvec, int ii, double alpha){
    double alpha_n0 = alpha*nvec.x;
    double alpha_n1 = alpha*nvec.y;
    double alpha_n2 = alpha*nvec.z;
    double alpha_n1n0 = alpha_n0*nvec.y;
    double alpha_n2n1 = alpha_n1*nvec.z;
    double alpha_n0n2 = alpha_n2*nvec.x;
    
	int ii3   = 3*ii;
	int ii3_1 = ii3 + 1;
	int ii3_2 = ii3 + 2;

	pair_resistance_matrix[ii3  ][ii3  ] += alpha_n0*nvec.x; 
	pair_resistance_matrix[ii3_1][ii3  ] += alpha_n1n0; // 10
	pair_resistance_matrix[ii3_2][ii3  ] += alpha_n0n2; // 20
	
	pair_resistance_matrix[ii3_1][ii3_1] += alpha_n1*nvec.y; //11
	pair_resistance_matrix[ii3_2][ii3_1] += alpha_n2n1; // 21
	
  	pair_resistance_matrix[ii3_2][ii3_2]  += alpha_n2*nvec.z; //22

}


void
BrownianForce::addToOffDiagBlock_RFU(const vec3d &nvec, double alpha){


    double alpha_n0 = alpha*nvec.x;
    double alpha_n1 = alpha*nvec.y;
    double alpha_n2 = alpha*nvec.z;
    double alpha_n1n0 = alpha_n0*nvec.y;
    double alpha_n2n1 = alpha_n1*nvec.z;
    double alpha_n0n2 = alpha_n2*nvec.x;
    

	pair_resistance_matrix[3][0] = alpha_n0*nvec.x; // 00
    pair_resistance_matrix[4][0] = alpha_n1n0; // 10
    pair_resistance_matrix[5][0] = alpha_n0n2; // 20
	pair_resistance_matrix[3][1] = alpha_n1n0; // 01
	pair_resistance_matrix[4][1] = alpha_n1*nvec.y; //11
	pair_resistance_matrix[5][1] = alpha_n2n1; // 21
	pair_resistance_matrix[3][2] = alpha_n0n2; // 02
	pair_resistance_matrix[4][2] = alpha_n2n1; // 12
	pair_resistance_matrix[5][2] = alpha_n2*nvec.z; //22
}


void
BrownianForce::factorize(){
  double fact, sqrt_fact;
  for(int a=0; a<6; a++){
	for(int b=a; b<6; b++){
	  fact = pair_resistance_matrix[b][a]; 
	  for(int c=0; c<a; c++){
		fact = fact - L_factor[a][c]*L_factor[b][c];
	  }
	  if (a == b){
		if (fact <= 0.){
		  cout << " Error : BrownianForce::CholeskyDecomposition() failed." << endl;
		  exit(1);
		}
		sqrt_fact = sqrt(fact);
		L_factor[a][a] = sqrt_fact;
	  }
	  else{
		L_factor[b][a] = fact/sqrt_fact;
	  }
	}
  }
  
}
