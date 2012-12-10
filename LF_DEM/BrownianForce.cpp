#include "BrownianForce.h"


#define RANDOM ( r_gen.rand() ) // RNG uniform in [0,1]
#define GRANDOM ( r_gen.randNorm(0., 1.) ) // RNG gaussian with mean 0. and variance 1.

using namespace std;

// BrownianForce::BrownianForce(double kb_T, cholmod_dense *brownian_force_vector, cholmod_common *c_, System *sys_){
// 	sys = sys_;
// 	force = *brownian_force_vector;
// 	c = c_;
// }

BrownianForce::BrownianForce(System *sys_){
	sys = sys_;
	kb_T = sys->kb_T;
	kb_T2= 2*kb_T;
}
BrownianForce::~BrownianForce(){
	cholmod_free_dense(&forces, c);
	cholmod_free_dense(&rand_vec, c);
}

void
BrownianForce::init(){
	c=&(sys->c);
	rand_vec = cholmod_zeros(3*sys->numpart(), 1, CHOLMOD_REAL, c);
	forces = cholmod_zeros(3*sys->numpart(), 1, CHOLMOD_REAL, c);
}

void
BrownianForce::add_to(cholmod_dense *force_vec){
	generate();
  
	int n3=3*sys->numpart();
	for(int i=0; i<n3;i++){
	  ((double*)force_vec->x)[i] += ((double*)forces->x)[i];
	}
}

void
BrownianForce::generate(){

	cholmod_factor* L_copy = cholmod_copy_factor(sys->L, c); // sadly it seems we have to make a copy. Is there a way to avoid this?
	L_sparse = cholmod_factor_to_sparse(L_copy, c);
	
	double temperature [2] = {kb_T2,0};
	double zero [2] = {0,0};
	cholmod_sdmult(L_sparse, 0, temperature, zero, random_vector(), forces, c);

	// cout << endl<< endl << " brownian force "<< endl<< endl;
	// cout << " kbt " << kb_T2 << endl;
	// for (int i = 0; i < 20; i++){
	// 	cout << ((double*)forces->x)[i] << endl;
	// }

	cholmod_free_sparse(&L_sparse, c);
	cholmod_free_factor(&L_copy, c);

}


void
BrownianForce::generate(double* f){
	generate();
	int n3 = 3*sys->numpart();
	for(int i=0; i<n3; i++){
		f[i] = ((double*)forces->x)[i];
	}
}

cholmod_dense*
BrownianForce::random_vector(){
	int n3 = 3*sys->numpart();
	for(int i=0; i<n3; i++){
		((double*)rand_vec->x)[i] = GRANDOM;
	}

	return rand_vec;
}
