#include "fstream"
#include "StokesSolver.h"
using namespace std;
#define DELETE(x) if(x){delete [] x; x = NULL;}

/******************************************************
 *                                                     *
 *                   Public Methods                    *
 *                                                     *
 ******************************************************/

StokesSolver::~StokesSolver()
{
	if (!odbrows_table) {
		delete [] odbrows_table;
	}
  if (!odbrows_table_mf) {
		delete [] odbrows_table_mf;
	}
  if (!odbrows_table_ff) {
		delete [] odbrows_table_ff;
	}
	if (!chol_solution) {
		cholmod_free_dense(&chol_solution, &chol_c);
	}
	if (!chol_rhs) {
		cholmod_free_dense(&chol_rhs, &chol_c);
	}
//	if (brownian) {
//		cholmod_free_dense(&chol_brownian_rhs, &chol_c);
//	}
	if (!chol_res_matrix) {
		cholmod_free_sparse(&chol_res_matrix, &chol_c);
	}
	if (chol_init) {
		cholmod_finish(&chol_c);
	}
}

void StokesSolver::init(int n)
{
	np = n;
	mobile_particle_nb = np;
	//brownian = is_brownian;
	// initializing values that can be changed later
	chol_init = false;
	//	FTcoupling = false;
}

void StokesSolver::initialize()
{
	// CHOLMOD parameters
	stype = -1; // 1 is symmetric, stored upper triangular (UT), -1 is LT
	sorted = 0;		/* TRUE if columns sorted, FALSE otherwise*/
	packed = 1;		/* TRUE if matrix packed, FALSE otherwise */
	xtype = CHOLMOD_REAL;
	// resistance matrix characteristics (see header for matrix description)
	dblocks_size = 18*mobile_particle_nb;
	allocateRessources();
	chol_L_to_be_freed = false;
}

/************* Matrix filling methods **********************/

// Diagonal Blocks Terms, FT/UW version
void StokesSolver::addToDiagBlock(const vec3d& nvec, int ii,
								  double scaledXA, double scaledYA,
								  double scaledYB, double scaledYC)
{

  if(ii > mobile_particle_nb){
    // FF matrix
    addToDBlock(dblocks_ff[ii-mobile_particle_nb], nvec, scaledXA, scaledYA, scaledYB, scaledYC);
  }
  else{
    // MM matrix
    addToDBlock(dblocks[ii], nvec, scaledXA, scaledYA, scaledYB, scaledYC);
  }
}


// Diagonal Blocks Terms, FT/UW version
void StokesSolver::addToDBlock(struct DBlock &b, const vec3d& nvec,
								  double scaledXA, double scaledYA,
								  double scaledYB, double scaledYC)
{
	double n0n0 = nvec.x*nvec.x;
	double n0n1 = nvec.x*nvec.y;
	double n0n2 = nvec.x*nvec.z;
	double n1n1 = nvec.y*nvec.y;
	double n1n2 = nvec.y*nvec.z;
	double n2n2 = nvec.z*nvec.z;
	double one_n0n0 = 1-n0n0;
	double one_n1n1 = 1-n1n1;
	double one_n2n2 = 1-n2n2;

	// (*,0)
	b.col0[0] += scaledXA*n0n0 + scaledYA*one_n0n0; // 00 element of the dblock
	b.col0[1] += (scaledXA-scaledYA)*n0n1; // 10
	b.col0[2] += (scaledXA-scaledYA)*n0n2; // 20
	b.col0[3] += 0;                    // 30
	b.col0[4] += -scaledYB*nvec.z;     // 40
	b.col0[5] += +scaledYB*nvec.y;      // 50
	// (*,1)
	b.col1[0] += scaledXA*n1n1 + scaledYA*one_n1n1;        // 11
	b.col1[1] += (scaledXA-scaledYA)*n1n2;        // 21
	b.col1[2] += 0;                    // 41
	b.col1[3] += -scaledYB*nvec.x;     // 51
	// (*,2)
	b.col2[0] += scaledXA*n2n2 + scaledYA*one_n2n2;        // 22
	b.col2[1] += 0;                    // 32
	// (*,3)
	b.col3[0] += scaledYC*one_n0n0;    // 33
	b.col3[1] += -scaledYC*n0n1;       // 43
	b.col3[2] += -scaledYC*n0n2;       // 53
	// (*,4)
	b.col4[0] += scaledYC*one_n1n1;    // 44
	b.col4[1] += -scaledYC*n1n2;       // 54
	// (*,5)
	b.col5[0] += scaledYC*one_n2n2;    // 55
}


// Off-Diagonal Blocks Terms, FT/UW version
void StokesSolver::setOffDiagBlock(const vec3d& nvec, int jj,
								   double scaledXA,
								   double scaledYA, double scaledYB,
								   double scaledYBtilde, double scaledYC)
{
  // FF Matrix
  if(mobile_matrix_done){
    odbrows_ff.push_back(6*(jj-mobile_particle_nb));
    odblocks_ff[odbrows_ff.size()-1] = buildODBlock(nvec, scaledXA, scaledYA, scaledYB, scaledYBtilde, scaledYC);
  }
  else{
    // MM matrix
    if(jj < mobile_particle_nb){
      odbrows.push_back(6*jj);
      odblocks[odbrows.size()-1] = buildODBlock(nvec, scaledXA, scaledYA, scaledYB, scaledYBtilde, scaledYC);
    }
    // MF matrix
    else{
      odbrows_mf.push_back(6*(jj-mobile_particle_nb));
      odblocks_mf[odbrows_mf.size()-1] = buildODBlock(nvec, scaledXA, scaledYA, scaledYB, scaledYBtilde, scaledYC);
    }
  }
	return;
}


/*************** Cholmod Matrix Filling *************
 Cholmod matrices we are using are defined in column major order (index j is column index)

 Cholmod matrices are defined as follows:
 - all values are stored in array x ( size nzmax )
 - locations of values are encoded in array p ( size np ):
 values corresponding to column j are x[ p[j] ]  to x[ p[j+1] - 1 ]
 - corresponding rows are stored in array i ( size nzmax ):
 rows corresponding to column j are i[ p[j] ]  to i[ p[j+1] - 1 ]

 Hence:
 with p[j]-1 < a < p[j+1]
        . . . . j . . . . . .
     .|         .            |
     .|         .            |
     .|         .            |
 i[a] | . . . .x[a]          |
     .|                      |
     .|                      |


 *****************************************************/

void StokesSolver::completeResistanceMatrix()
{
	// this function is commented, but you are strongly advised to read
	// the description of storage in the header file first :)

    // declare the last 2 values of odbrows_table
    int size = mobile_particle_nb; // temporary
    odbrows_table[size-1] = (unsigned int)odbrows.size();
    odbrows_table[size] = (unsigned int)odbrows.size();
    allocateResistanceMatrix();
    // fill
	for (int j=0; j<size; j++) {
		// associated with particle j are 6 columns in the matrix:
		// { 6j, ... , 6j+5 }
		int j6 = 6*j;
		int j21 = 21*j;
		// the number of non-zero elements before column 6j is:
		// - 21*j from j diagonal blocks
		// - 36*odbrows_table[j] from odbrows_table[j] off-diagonal blocks
		//
		// the number of non-zero elements before column 6j+1 is:
		// - number of non-zero before column 6j + number of non-zero in column 6*j
		// (in 6j: 6 elements in diagonal block, plus 6*(odbrows_table[j+1]-odbrows_table[j])
		//
		// for 6j+2 --> 6j+5: same idea
		int od_nzero_nb = 6*(odbrows_table[j+1]-odbrows_table[j]);
		((int*)chol_res_matrix->p)[j6  ] = j21   + 36*odbrows_table[j];
		((int*)chol_res_matrix->p)[j6+1] = ((int*)chol_res_matrix->p)[j6] + 6 + od_nzero_nb;
		((int*)chol_res_matrix->p)[j6+2] = ((int*)chol_res_matrix->p)[j6+1] + 5 + od_nzero_nb;
		((int*)chol_res_matrix->p)[j6+3] = ((int*)chol_res_matrix->p)[j6+2] + 4 + od_nzero_nb;
		((int*)chol_res_matrix->p)[j6+4] = ((int*)chol_res_matrix->p)[j6+3] + 3 + od_nzero_nb;
		((int*)chol_res_matrix->p)[j6+5] = ((int*)chol_res_matrix->p)[j6+4] + 2 + od_nzero_nb;
		int pj6   = ((int*)chol_res_matrix->p)[j6];
		int pj6_1 = ((int*)chol_res_matrix->p)[j6+1];
		int pj6_2 = ((int*)chol_res_matrix->p)[j6+2];
		int pj6_3 = ((int*)chol_res_matrix->p)[j6+3];
		int pj6_4 = ((int*)chol_res_matrix->p)[j6+4];
		int pj6_5 = ((int*)chol_res_matrix->p)[j6+5];
		// diagonal block row indices (21)
		((int*)chol_res_matrix->i)[pj6  ] = j6;   // column j6
		((int*)chol_res_matrix->i)[pj6+1] = j6+1;
		((int*)chol_res_matrix->i)[pj6+2] = j6+2;
		((int*)chol_res_matrix->i)[pj6+3] = j6+3;
		((int*)chol_res_matrix->i)[pj6+4] = j6+4;
		((int*)chol_res_matrix->i)[pj6+5] = j6+5;
		((int*)chol_res_matrix->i)[pj6_1  ] = j6+1;    // column j6+1
		((int*)chol_res_matrix->i)[pj6_1+1] = j6+2;
		((int*)chol_res_matrix->i)[pj6_1+2] = j6+3;
		((int*)chol_res_matrix->i)[pj6_1+3] = j6+4;
		((int*)chol_res_matrix->i)[pj6_1+4] = j6+5;
		((int*)chol_res_matrix->i)[pj6_2  ] = j6+2;    // column j6+2
		((int*)chol_res_matrix->i)[pj6_2+1] = j6+3;
		((int*)chol_res_matrix->i)[pj6_2+2] = j6+4;
		((int*)chol_res_matrix->i)[pj6_2+3] = j6+5;
		((int*)chol_res_matrix->i)[pj6_3  ] = j6+3;    // column j6+3
		((int*)chol_res_matrix->i)[pj6_3+1] = j6+4;
		((int*)chol_res_matrix->i)[pj6_3+2] = j6+5;
		((int*)chol_res_matrix->i)[pj6_4  ] = j6+4;    // column j6+4
		((int*)chol_res_matrix->i)[pj6_4+1] = j6+5;
		((int*)chol_res_matrix->i)[pj6_5  ] = j6+5;    // column j6+5
		// diagonal blocks row values (21)
		((double*)chol_res_matrix->x)[pj6  ] = dblocks[j].col0[0];   // column j6
		((double*)chol_res_matrix->x)[pj6+1] = dblocks[j].col0[1];
		((double*)chol_res_matrix->x)[pj6+2] = dblocks[j].col0[2];
		((double*)chol_res_matrix->x)[pj6+3] = dblocks[j].col0[3];
		((double*)chol_res_matrix->x)[pj6+4] = dblocks[j].col0[4];
		((double*)chol_res_matrix->x)[pj6+5] = dblocks[j].col0[5];
		((double*)chol_res_matrix->x)[pj6_1  ] = dblocks[j].col1[0];   // column j6+1
		((double*)chol_res_matrix->x)[pj6_1+1] = dblocks[j].col1[1];
		((double*)chol_res_matrix->x)[pj6_1+2] = -dblocks[j].col0[4];  // anti-symmetry
		((double*)chol_res_matrix->x)[pj6_1+3] = dblocks[j].col1[2];
		((double*)chol_res_matrix->x)[pj6_1+4] = dblocks[j].col1[3];
		((double*)chol_res_matrix->x)[pj6_2  ] = dblocks[j].col2[0];   // column j6+2
		((double*)chol_res_matrix->x)[pj6_2+1] = -dblocks[j].col0[5];   // anti-symmetry
		((double*)chol_res_matrix->x)[pj6_2+2] = -dblocks[j].col1[3];   // anti-symmetry
		((double*)chol_res_matrix->x)[pj6_2+3] = dblocks[j].col2[1];
		((double*)chol_res_matrix->x)[pj6_3  ] = dblocks[j].col3[0];   // column j6+3
		((double*)chol_res_matrix->x)[pj6_3+1] = dblocks[j].col3[1];
		((double*)chol_res_matrix->x)[pj6_3+2] = dblocks[j].col3[2];
		((double*)chol_res_matrix->x)[pj6_4  ] = dblocks[j].col4[0];   // column j6+4
		((double*)chol_res_matrix->x)[pj6_4+1] = dblocks[j].col4[1];
		((double*)chol_res_matrix->x)[pj6_5  ] = dblocks[j].col5[0];   // column j6+5
		/*****  2  : off-diagonal blocks row indices and values ***********/
		// 36 non-zero elements per block
		for (int k = odbrows_table[j]; k<odbrows_table[j+1]; k++) {
			int u = 6*(k-odbrows_table[j]);
			// we are filling the "k-odbFrows_table[j]"th off-diag block of the column.
			// For column j6, for exemple, the indices of the non-zero values are:
			// pj6 for all non-zero elements before column j6,
			// + 6 for the diagonal block of column j6
			// + u (=6*(k-odbFrows_table[j])) for the off-diag blocks of j6
			// + index inside the current block
			for (int s=0; s<6; s++) {
				((int*)chol_res_matrix->i)[pj6+6 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_1+5 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_2+4 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_3+3 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_4+2 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_5+1 + u +s] = odbrows[k]+s;
			}

			((double*)chol_res_matrix->x)[pj6+6 + u   ]   = odblocks[k].col0[0]; // A   // column j6
			((double*)chol_res_matrix->x)[pj6+6 + u +1]   = odblocks[k].col0[1];
			((double*)chol_res_matrix->x)[pj6+6 + u +2]   = odblocks[k].col0[2];
			((double*)chol_res_matrix->x)[pj6+6 + u +3]   = odblocks[k].col0[3]; // B
			((double*)chol_res_matrix->x)[pj6+6 + u +4]   = odblocks[k].col0[4];
			((double*)chol_res_matrix->x)[pj6+6 + u +5]   = odblocks[k].col0[5];
			//
			((double*)chol_res_matrix->x)[pj6_1+5 + u   ] = odblocks[k].col0[1]; // A  // column j6+1
			((double*)chol_res_matrix->x)[pj6_1+5 + u +1] = odblocks[k].col1[0];
			((double*)chol_res_matrix->x)[pj6_1+5 + u +2] = odblocks[k].col1[1];
			((double*)chol_res_matrix->x)[pj6_1+5 + u +3] = -odblocks[k].col0[4]; // B
			((double*)chol_res_matrix->x)[pj6_1+5 + u +4] = odblocks[k].col1[2];
			((double*)chol_res_matrix->x)[pj6_1+5 + u +5] = odblocks[k].col1[3];
			//
			((double*)chol_res_matrix->x)[pj6_2+4 + u   ] = odblocks[k].col0[2]; // A // column j6+2
			((double*)chol_res_matrix->x)[pj6_2+4 + u +1] = odblocks[k].col1[1];
			((double*)chol_res_matrix->x)[pj6_2+4 + u +2] = odblocks[k].col2[0];
			((double*)chol_res_matrix->x)[pj6_2+4 + u +3] = -odblocks[k].col0[5]; // B
			((double*)chol_res_matrix->x)[pj6_2+4 + u +4] = -odblocks[k].col1[3];
			((double*)chol_res_matrix->x)[pj6_2+4 + u +5] = odblocks[k].col2[1];
			//
			((double*)chol_res_matrix->x)[pj6_3+3 + u   ] = odblocks[k].col3[0]; // Btilde   // column j6+3
			((double*)chol_res_matrix->x)[pj6_3+3 + u +1] = odblocks[k].col3[1];
			((double*)chol_res_matrix->x)[pj6_3+3 + u +2] = odblocks[k].col3[2];
			((double*)chol_res_matrix->x)[pj6_3+3 + u +3] = odblocks[k].col3[3]; // C
			((double*)chol_res_matrix->x)[pj6_3+3 + u +4] = odblocks[k].col3[4];
			((double*)chol_res_matrix->x)[pj6_3+3 + u +5] = odblocks[k].col3[5];

			((double*)chol_res_matrix->x)[pj6_4+2 + u   ] = -odblocks[k].col3[1]; // Btilde // column j6+4
			((double*)chol_res_matrix->x)[pj6_4+2 + u +1] = odblocks[k].col4[0];
			((double*)chol_res_matrix->x)[pj6_4+2 + u +2] = odblocks[k].col4[1];
			((double*)chol_res_matrix->x)[pj6_4+2 + u +3] = odblocks[k].col3[4]; // C
			((double*)chol_res_matrix->x)[pj6_4+2 + u +4] = odblocks[k].col4[2];
			((double*)chol_res_matrix->x)[pj6_4+2 + u +5] = odblocks[k].col4[3];
			//
			((double*)chol_res_matrix->x)[pj6_5+1 + u   ] = -odblocks[k].col3[2]; // Btilde // column j6+5
			((double*)chol_res_matrix->x)[pj6_5+1 + u +1] = -odblocks[k].col4[1];
			((double*)chol_res_matrix->x)[pj6_5+1 + u +2] = odblocks[k].col5[0];
			((double*)chol_res_matrix->x)[pj6_5+1 + u +3] = odblocks[k].col3[5]; // C
			((double*)chol_res_matrix->x)[pj6_5+1 + u +4] = odblocks[k].col4[3];
			((double*)chol_res_matrix->x)[pj6_5+1 + u +5] = odblocks[k].col5[1];
		}
	}
  ((int*)chol_res_matrix->p)[6*size] = ((int*)chol_res_matrix->p)[6*size-1]+1;

	// printResistanceMatrix(cout, "dense");getchar();
	factorizeResistanceMatrix();
}



void StokesSolver::resetResistanceMatrix(int nb_of_interactions,
										 double *reset_resmat_dblocks)
{
	//setSolverType(solver_type);
	odblocks_nb = nb_of_interactions;

	for(auto &b: dblocks){
		resetDBlock(b, reset_resmat_dblocks);
	}
	odbrows.clear();
  odblocks.resize(odblocks_nb);
  for(auto &b: odblocks){
    resetODBlock(b);
  }
	odbrows_table[0] = 0;
	mobile_matrix_done = false;

// for the mixed problem
	int odblocks_mf_nb = 0;
  odbrows_mf.clear();
  odblocks_mf.resize(odblocks_mf_nb);
  for(auto &b: odblocks_mf){
    resetODBlock(b);
  }
	odbrows_table_mf[0] = 0;

	// for the mixed problem
	for(auto &b: dblocks_ff){
		resetDBlock(b, reset_resmat_dblocks);
	}
	int odblocks_ff_nb = 0;
  odbrows_ff.clear();
  odblocks_ff.resize(odblocks_ff_nb);
  for(auto &b: odblocks_ff){
    resetODBlock(b);
  }
	odbrows_table_ff[0] = 0;
}

void StokesSolver::resetRHS()
{
	int size = chol_rhs->nrow;
	for (int i=0; i<size; i++){
		((double*)chol_rhs->x)[i] = 0;
	}
}

void StokesSolver::resetRHStorque()
{
	int size = chol_rhs->nrow;
	for (int i=3; i<size; i+=6){
		((double*)chol_rhs->x)[i] = 0;
		((double*)chol_rhs->x)[i+1] = 0;
		((double*)chol_rhs->x)[i+2] = 0;
	}
}

void StokesSolver::addToRHSForce(int i, double* force_i)
{
	int i6 = 6*i;
	for (int u=0; u<3; u++) {
		((double*)chol_rhs->x)[i6+u] += force_i[u];
	}
}

void StokesSolver::addToRHSForce(int i, const vec3d& force_i)
{
	int i6 = 6*i;
	((double*)chol_rhs->x)[i6] += force_i.x;
	((double*)chol_rhs->x)[i6+1] += force_i.y;
	((double*)chol_rhs->x)[i6+2] += force_i.z;
}

void StokesSolver::addToRHSTorque(int i, double* torque_i)
{
	int i6_3 = 6*i+3;
	for (int u=0; u<3; u++) {
		((double*)chol_rhs->x)[i6_3+u] += torque_i[u];
	}
}

void StokesSolver::addToRHSTorque(int i, const vec3d& torque_i)
{
	int i6_3 = 6*i+3;
	((double*)chol_rhs->x)[i6_3] += torque_i.x;
	((double*)chol_rhs->x)[i6_3+1] += torque_i.y;
	((double*)chol_rhs->x)[i6_3+2] += torque_i.z;
}

void StokesSolver::addToRHS(double* rhs)
{
	int size = chol_rhs->nrow;
	for (int i=0; i<size; i++) {
		((double*)chol_rhs->x)[i] += rhs[i];
	}
}

void StokesSolver::setRHS(double* rhs)
{
	int size = chol_rhs->nrow;
	for (int i=0; i<size; i++) {
		((double*)chol_rhs->x)[i] = rhs[i];
	}
}

void StokesSolver::setRHSForce(int i, const vec3d& force_i)
{
	int i6 = 6*i;
	((double*)chol_rhs->x)[i6] = force_i.x;
	((double*)chol_rhs->x)[i6+1] = force_i.y;
	((double*)chol_rhs->x)[i6+2] = force_i.z;
}

void StokesSolver::setRHSTorque(int i, const vec3d& torque_i)
{
	int i6_3 = 6*i+3;
	((double*)chol_rhs->x)[i6_3] = torque_i.x;
	((double*)chol_rhs->x)[i6_3+1] = torque_i.y;
	((double*)chol_rhs->x)[i6_3+2] = torque_i.z;
}

void StokesSolver::getRHS(double* rhs)
{
	int size = chol_rhs->nrow;
	for (int i=0; i<size; i++) {
		rhs[i] = ((double*)chol_rhs->x)[i];
	}
}


// Computes X = L*RHS
void StokesSolver::compute_LTRHS(double* X)
{
	/*
	 Cholmod gives a factorizationof a permutated resistance
	 matrix Lc*Lc^T = P*RFU*P^T

	 That means P*L = Lc, with  L*L^T = RFU

	 So for a rhs Y:
	 X = L*Y = P^T*Lc*Y
		*/
	if (!chol_L->is_ll) {
		cerr << " The factorization is LDL^T. compute_LTRHS(double* X) only works for LL^T factorization." << endl;
	}
	double alpha[] = {1, 0};
	double beta[] = {0, 0};
	int transpose = 0;
	cholmod_factor* chol_L_copy = cholmod_copy_factor(chol_L, &chol_c);
	cholmod_sparse* chol_L_sparse = cholmod_factor_to_sparse(chol_L_copy, &chol_c);
	cholmod_sdmult(chol_L_sparse, transpose, alpha, beta, chol_rhs, chol_Psolution, &chol_c); // chol_Psolution = Lc*Y
	chol_solution = cholmod_solve(CHOLMOD_Pt, chol_L, chol_Psolution, &chol_c); // chol_solution = P^T*chol_Psolution

	int size = chol_solution->nrow;
	for (int i=0; i<size; i++) {
		X[i] = ((double*)chol_solution->x)[i];
	}
	cholmod_free_sparse(&chol_L_sparse, &chol_c);
	cholmod_free_dense(&chol_solution, &chol_c);
	cholmod_free_factor(&chol_L_copy, &chol_c);
}

// Finds solutions to L^T X = RHS
void StokesSolver::solve_LT(double* X)
{
	chol_PTsolution = cholmod_solve(CHOLMOD_Lt, chol_L, chol_rhs, &chol_c) ;
	chol_solution = cholmod_solve(CHOLMOD_Pt, chol_L, chol_PTsolution, &chol_c) ;
	int size = chol_solution->nrow;
	for (int i=0; i<size; i++) {
		X[i] = ((double*)chol_solution->x)[i];
	}
	cholmod_free_dense(&chol_solution, &chol_c);
	cholmod_free_dense(&chol_PTsolution, &chol_c);
}

void StokesSolver::solve_LT(vec3d* X, vec3d* ang_X)
{
	chol_PTsolution = cholmod_solve(CHOLMOD_Lt, chol_L, chol_rhs, &chol_c) ;
	chol_solution = cholmod_solve(CHOLMOD_Pt, chol_L, chol_PTsolution, &chol_c) ;
	int size = chol_solution->nrow/6;
	for (int i=0; i<size; i++) {
		int i6 =6*i;
		X[i].x = ((double*)chol_solution->x)[i6];
		X[i].y = ((double*)chol_solution->x)[i6+1];
		X[i].z = ((double*)chol_solution->x)[i6+2];
		ang_X[i].x = ((double*)chol_solution->x)[i6+3];
		ang_X[i].y = ((double*)chol_solution->x)[i6+4];
		ang_X[i].z = ((double*)chol_solution->x)[i6+5];
	}
	cholmod_free_dense(&chol_solution, &chol_c);
	cholmod_free_dense(&chol_PTsolution, &chol_c);
}

void StokesSolver::solve(vec3d* velocity, vec3d* ang_velocity)
{
	chol_solution = cholmod_solve(CHOLMOD_A, chol_L, chol_rhs, &chol_c);
	int size = chol_solution->nrow/6;
	for (int i=0; i<size; i++) {
		int i6 =6*i;
		velocity[i].x = ((double*)chol_solution->x)[i6];
		velocity[i].y = ((double*)chol_solution->x)[i6+1];
		velocity[i].z = ((double*)chol_solution->x)[i6+2];
		ang_velocity[i].x = ((double*)chol_solution->x)[i6+3];
		ang_velocity[i].y = ((double*)chol_solution->x)[i6+4];
		ang_velocity[i].z = ((double*)chol_solution->x)[i6+5];
	}
	cholmod_free_dense(&chol_solution, &chol_c);
}

void StokesSolver::solve(double* velocity)
{
	chol_solution = cholmod_solve(CHOLMOD_A, chol_L, chol_rhs, &chol_c);
	int size = chol_solution->nrow;
	for (int i=0; i<size; i++) {
		velocity[i] = ((double*)chol_solution->x)[i];
	}
	cholmod_free_dense(&chol_solution, &chol_c);
}

// testing function, don't use it in production code, very slow and unclean
void StokesSolver::multiplySolutionByResMat(double* vec)
{
	chol_solution = cholmod_solve(CHOLMOD_A, chol_L, chol_rhs, &chol_c);
	cholmod_dense* r;
   	r = cholmod_copy_dense(chol_rhs, &chol_c);
	double one[] = {1, 0};
	double zero[] = {0, 0};
	cholmod_sdmult(chol_res_matrix, 0, one, zero, chol_solution, r, &chol_c);
	int size = r->nrow;
	for (int i=0; i<size; i++) {
			vec[i] = ((double*)r->x)[i];
	}
	cholmod_free_dense(&r, &chol_c);
	cholmod_free_dense(&chol_solution, &chol_c);
}


void StokesSolver::solvingIsDone()
{
	cholmod_free_factor(&chol_L, &chol_c);
	cholmod_free_sparse(&chol_res_matrix, &chol_c);
	//	cholmod_free_dense(&chol_rhs, &chol_c);
}

/******************************************************
 *                                                     *
 *                  Private Methods                    *
 *                                                     *
 ******************************************************/

void StokesSolver::allocateRessources()
{
		dblocks.resize(mobile_particle_nb);
		odbrows_table = new int [mobile_particle_nb+1];
		odbrows_table_mf = new int [mobile_particle_nb+1];
		odbrows_table_ff = new int [np-mobile_particle_nb+1];
		dblocks_ff.resize(np-mobile_particle_nb);

		cholmod_start(&chol_c);
		chol_init = true;
		int size = 6*mobile_particle_nb;
		chol_rhs = cholmod_allocate_dense(size, 1, size, xtype, &chol_c);
		chol_Psolution = cholmod_allocate_dense(size, 1, size, xtype, &chol_c); // used for Brownian motion
		for (int i=0; i<size; i++) {
			((double*)chol_rhs->x)[i] = 0;
		}
		chol_L = NULL;
}

void StokesSolver::allocateResistanceMatrix()
{
	// allocate
	int nzmax; // non-zero values
	int size = 6*dblocks.size();

	nzmax = 21*dblocks.size(); // diagonal blocks
	nzmax += 36*odblocks_nb;  // off-diagonal
	chol_res_matrix = cholmod_allocate_sparse(size, size, nzmax, sorted, packed, stype, xtype, &chol_c);
}

void StokesSolver::doneBlocks(int i)
{
  if (i==mobile_particle_nb) {
    mobile_matrix_done = true;
  }
  if (mobile_matrix_done) {
    i -= mobile_matrix_done;
    odbrows_table_ff[i+1] = (unsigned int)odbrows_ff.size();
  } else {
    odbrows_table[i+1] = (unsigned int)odbrows.size();
    odbrows_table_mf[i+1] = (unsigned int)odbrows_mf.size();
  }

}

// odblocks fillings, for FT/UW version
struct ODBlock StokesSolver::buildODBlock(const vec3d &nvec,
							 double scaledXA,
							 double scaledYA, double scaledYB,
							 double scaledYBtilde, double scaledYC)
{
	double n0n0 = nvec.x*nvec.x;
	double n0n1 = nvec.x*nvec.y;
	double n0n2 = nvec.x*nvec.z;
	double n1n1 = nvec.y*nvec.y;
	double n1n2 = nvec.y*nvec.z;
	double n2n2 = nvec.z*nvec.z;
	double one_n0n0 = 1-n0n0;
	double one_n1n1 = 1-n1n1;
	double one_n2n2 = 1-n2n2;
  struct ODBlock block;

	block.col0[0] = scaledXA*n0n0+scaledYA*one_n0n0; // column 0
	block.col0[1] = (scaledXA-scaledYA)*n0n1;
	block.col0[2] = (scaledXA-scaledYA)*n0n2;
	block.col0[3] = 0;
	block.col0[4] = -scaledYB*nvec.z;
	block.col0[5] = scaledYB*nvec.y;
	block.col1[0] = scaledXA*n1n1+scaledYA*one_n1n1; // column 1
	block.col1[1] = (scaledXA-scaledYA)*n1n2;
	block.col1[2] = 0;
	block.col1[3] = -scaledYB*nvec.x;
	block.col2[0] = scaledXA*n2n2+scaledYA*one_n2n2; // column 2
	block.col2[1] = 0;
	block.col3[0] = 0; // column 3
	block.col3[1] = scaledYBtilde*nvec.z;
	block.col3[2] = -scaledYBtilde*nvec.y;
	block.col3[3] = scaledYC*one_n0n0;
	block.col3[4] = -scaledYC*n0n1;
	block.col3[5] = -scaledYC*n0n2;
	block.col4[0] = 0; // column 4
	block.col4[1] = scaledYBtilde*nvec.x;
	block.col4[2] = scaledYC*one_n1n1;
	block.col4[3] = -scaledYC*n1n2;
	block.col5[0] = 0; // column 5
	block.col5[1] = scaledYC*one_n2n2;
  return block;
}


void StokesSolver::factorizeResistanceMatrix()
{
	/*reference code */
	chol_c.supernodal = CHOLMOD_SUPERNODAL;
	chol_L = cholmod_analyze(chol_res_matrix, &chol_c);
	cholmod_factorize(chol_res_matrix, chol_L, &chol_c);
    if (chol_c.status) {
		// Cholesky decomposition has failed: usually because matrix is incorrectly found to be positive-definite
		// It is very often enough to force another preconditioner to solve the problem.
		cerr << " factorization failed. forcing simplicial algorithm... " << endl;
		chol_c.supernodal = CHOLMOD_SIMPLICIAL;
		chol_L = cholmod_analyze(chol_res_matrix, &chol_c);
		cholmod_factorize(chol_res_matrix, chol_L, &chol_c);
		cerr << " factorization status " << chol_c.status << " final_ll ( 0 is LDL, 1 is LL ) " << chol_c.final_ll <<endl;
		chol_c.supernodal = CHOLMOD_SUPERNODAL;
    }
}

#ifdef CHOLMOD_EXTRA
/*
 setSpInvPreconditioner() :
 A sparse inverse (left-)preconditioner.
 cholmod-extra routine by Jaakko Luttinen

 */
void StokesSolver::setSpInvPreconditioner()
{
	cholmod_sparse* sparse_inv = cholmod_spinv(chol_L, &chol_c);
	cholmod_free_sparse(&sparse_inv, &chol_c);
}
#endif

// testing
void StokesSolver::printResistanceMatrix(ostream& out, string sparse_or_dense)
{
	if (sparse_or_dense == "sparse") {
		//		out << endl<< " chol res " << endl;
		int size = chol_res_matrix->nrow;
		for (int i = 0; i<size; i++) {
			for (int k =((int*)chol_res_matrix->p)[i]; k<((int*)chol_res_matrix->p)[i+1]; k++) {
				out << i << " " << ((int*)chol_res_matrix->i)[k] << " " << ((double*)chol_res_matrix->x)[k] << endl;
			}
		}
	}
	if (sparse_or_dense == "dense") {
		cholmod_dense* dense_res = cholmod_sparse_to_dense(chol_res_matrix, &chol_c);
		int size = chol_res_matrix->nrow;
		for (int i = 0; i<size; i++) {
			// if(i==0){
			// 	for (int j = 0; j < size/6; j++) {
			// 		out << j << "\t \t \t \t \t \t" ;
			// 	}
			// 	out << endl;
			// }
			for (int j = 0; j<size; j++) {
				out << setprecision(3) << ((double*)dense_res->x)[i+j*size] << "\t" ;
			}
			out << endl;
		}
		out << endl;
		cholmod_free_dense(&dense_res, &chol_c);
	}
}

void StokesSolver::printFactor(ostream &out)
{
	cholmod_factor* chol_L_copy = cholmod_copy_factor(chol_L, &chol_c);
	cholmod_sparse* chol_L_sparse = cholmod_transpose(cholmod_factor_to_sparse(chol_L_copy, &chol_c), 1, &chol_c);
	cholmod_dense* chol_L_dense = cholmod_sparse_to_dense(chol_L_sparse, &chol_c);
	//	cholmod_sparse* chol_PTL_sparse = cholmod_dense_to_sparse(cholmod_solve(CHOLMOD_Pt, chol_L, chol_L_dense, &chol_c), 1, &chol_c) ; // chol_solution = P^T*chol_Psolution
	//	cholmod_dense* chol_LT_dense = cholmod_sparse_to_dense(cholmod_transpose(chol_L_sparse, 1, &chol_c), &chol_c);
	int transpose = 1;
	double alpha[] = {1, 0};
	double beta[] = {0, 0};
	cholmod_dense *dense_res = cholmod_sparse_to_dense(chol_res_matrix, &chol_c);
	cholmod_sdmult(chol_L_sparse, transpose, alpha, beta, chol_L_dense, dense_res, &chol_c);
	int size = chol_res_matrix->nrow;

	// out << " Cholesky factor" << endl;
	// for (int i = 0; i < size; i++) {
	// 		// if(i==0){
	// 		// 	for (int j = 0; j < size/6; j++) {
	// 		// 		out << j << "\t \t \t \t \t \t" ;
	// 		// 	}
	// 		// 	out << endl;
	// 		// }
	// 	for (int j = 0; j < size; j++) {
	// 		out <<setprecision(3) <<  ((double*)chol_L_dense->x)[i+j*size] << "\t" ;
	// 	}
	// 	out << endl;
	// }
	// out << endl;
	//		out << " Cholesky squared  " << endl;
	for (int i = 0; i<size; i++) {
		// if(i==0){
		// 	for (int j = 0; j < size/6; j++) {
		// 		out << j << "\t \t \t \t \t \t" ;
		// 	}
		// 	out << endl;
		// }
		for (int j=0; j<size; j++) {
			out <<setprecision(3) << ((double*)dense_res->x)[i+j*size] << "\t" ;
		}
		out << endl;
	}
	out << endl;
	cholmod_free_sparse(&chol_L_sparse, &chol_c);
	//		cholmod_free_dense(&chol_LT_dense, &chol_c);
	//		cholmod_free_sparse(&chol_PTL_sparse, &chol_c);
	cholmod_free_factor(&chol_L_copy, &chol_c);
	cholmod_free_dense(&dense_res, &chol_c);
}

// testing
void StokesSolver::printRHS()
{
	int size = chol_rhs->nrow;

	for (int i=0; i<size; i++) {
		cout << i << " (part " << " " << (i-i%6)/6 << " )  " << ((double*)chol_rhs->x)[i] << endl;
	}
}
