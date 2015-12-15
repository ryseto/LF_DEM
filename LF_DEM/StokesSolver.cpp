#include <fstream>
#include "StokesSolver.h"

using namespace std;

/******************************************************
 *                                                     *
 *                   Public Methods                    *
 *                                                     *
 ******************************************************/

StokesSolver::~StokesSolver()
{
	if (!chol_solution) {
		cholmod_free_dense(&chol_solution, &chol_c);
	}
	if (!chol_rhs) {
		cholmod_free_dense(&chol_rhs, &chol_c);
	}
	if (!chol_res_matrix) {
		cholmod_free_sparse(&chol_res_matrix, &chol_c);
	}
	cholmod_finish(&chol_c);
}

void StokesSolver::init(int n)
{
	np = n;
	mobile_particle_nb = np;

	// CHOLMOD parameters
	stype = -1; // 1 is symmetric, stored upper triangular (UT), -1 is LT
	sorted = 0;		/* TRUE if columns sorted, FALSE otherwise*/
	packed = 1;		/* TRUE if matrix packed, FALSE otherwise */
	xtype = CHOLMOD_REAL;
	// resistance matrix characteristics (see header for matrix description)
	dblocks_size = 18*mobile_particle_nb;
	allocateRessources();
	chol_L_to_be_freed = false;

	odb_layout.resize(6);
	odb_layout[0] = vector<int>({0,1,2,4,5});
	odb_layout[1] = vector<int>({0,1,2,3,5});
	odb_layout[2] = vector<int>({0,1,2,3,4});
	odb_layout[3] = vector<int>({1,2,3,4,5});
	odb_layout[4] = vector<int>({0,2,3,4,5});
	odb_layout[5] = vector<int>({0,1,3,4,5});
	db_layout.resize(6);
	db_layout[0] = vector<int>({0,1,2,4,5});
	db_layout[1] = vector<int>({1,2,3,5});
	db_layout[2] = vector<int>({2,3,4});
	db_layout[3] = vector<int>({3,4,5});
	db_layout[4] = vector<int>({4,5});
	db_layout[5] = vector<int>({5});

	dblocks_cntnonzero.resize(6);
	for (int i=0; i<6; i++) {
		dblocks_cntnonzero[i] = db_layout[i].size();
	}
}

/************* Matrix filling methods **********************/

// Diagonal Blocks Terms, FT/UW version
void StokesSolver::addToDiagBlock(const vec3d& nvec, int ii,
								  double scaledXA, double scaledYA,
								  double scaledYB, double scaledYC)
{

	if (ii > mobile_particle_nb) {
		// FF matrix
		addToDBlock(dblocks_ff[ii-mobile_particle_nb], nvec, scaledXA, scaledYA, scaledYB, scaledYC);
	} else {
		// MM matrix
		addToDBlock(dblocks[ii], nvec, scaledXA, scaledYA, scaledYB, scaledYC);
	}
}

// Diagonal Blocks Terms, FT/UW version
void StokesSolver::addToDBlock(struct DBlock& b, const vec3d& nvec,
							   double scaledXA, double scaledYA,double scaledYB, double scaledYC)
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
	b.col0[3] += -scaledYB*nvec.z;     // 40
	b.col0[4] += +scaledYB*nvec.y;      // 50
	// (*,1)
	b.col1[0] += scaledXA*n1n1 + scaledYA*one_n1n1;        // 11
	b.col1[1] += (scaledXA-scaledYA)*n1n2;        // 21
	b.col1[2] += -scaledYB*nvec.x;     // 51
	// (*,2)
	b.col2[0] += scaledXA*n2n2 + scaledYA*one_n2n2;        // 22
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
	if (mobile_matrix_done) {
		// FF Matrix
		odbrows_ff.push_back(6*(jj-mobile_particle_nb));
		odblocks_ff[odbrows_ff.size()-1] = buildODBlock(nvec, scaledXA, scaledYA, scaledYB, scaledYBtilde, scaledYC);
	} else {
		if (jj < mobile_particle_nb) {
			// MM matrix
			odbrows.push_back(6*jj);
			odblocks[odbrows.size()-1] = buildODBlock(nvec, scaledXA, scaledYA, scaledYB, scaledYBtilde, scaledYC);
		} else {
			// MF matrix
			odbrows_mf.push_back(6*(jj-mobile_particle_nb));
			odblocks_mf[odbrows_mf.size()-1] = buildODBlock(nvec, scaledXA, scaledYA, scaledYB, scaledYBtilde, scaledYC);
		}
	}
	return;
}

void StokesSolver::insertDBlockValues(double *matrix_x,
										 const vector<int>& index_chol_ix,
										 const struct DBlock& b)
{
	int pcol0 = index_chol_ix[0];
	int pcol1 = index_chol_ix[1];
	int pcol2 = index_chol_ix[2];
	int pcol3 = index_chol_ix[3];
	int pcol4 = index_chol_ix[4];
	int pcol5 = index_chol_ix[5];
	// diagonal blocks row values (21)
	matrix_x[pcol0  ] = b.col0[0];   // column j6
	matrix_x[pcol0+1] = b.col0[1];
	matrix_x[pcol0+2] = b.col0[2];
	matrix_x[pcol0+3] = b.col0[3];
	matrix_x[pcol0+4] = b.col0[4];
	matrix_x[pcol1  ] = b.col1[0];   // column j6+1
	matrix_x[pcol1+1] = b.col1[1];
	matrix_x[pcol1+2] = -b.col0[3];  // anti-symmetry
	matrix_x[pcol1+3] = b.col1[2];
	matrix_x[pcol2  ] = b.col2[0];   // column j6+2
	matrix_x[pcol2+1] = -b.col0[4];   // anti-symmetry
	matrix_x[pcol2+2] = -b.col1[2];   // anti-symmetry
	matrix_x[pcol3  ] = b.col3[0];   // column j6+3
	matrix_x[pcol3+1] = b.col3[1];
	matrix_x[pcol3+2] = b.col3[2];
	matrix_x[pcol4  ] = b.col4[0];   // column j6+4
	matrix_x[pcol4+1] = b.col4[1];
	matrix_x[pcol5  ] = b.col5[0];   // column j6+5
}

void StokesSolver::insertODBlockValues(double *matrix_x,
										  const vector<int>& index_chol_ix,
										  const struct ODBlock& b)
{
	int start_index = index_chol_ix[0];
	matrix_x[start_index  ]   = b.col0[0]; // A   // column j6
	matrix_x[start_index +1]   = b.col0[1];
	matrix_x[start_index +2]   = b.col0[2];
	matrix_x[start_index +3]   = b.col0[3]; // B
	matrix_x[start_index +4]   = b.col0[4];
	//
	start_index = index_chol_ix[1];
	matrix_x[start_index   ] = b.col0[1]; // A  // column j6+1
	matrix_x[start_index +1] = b.col1[0];
	matrix_x[start_index +2] = b.col1[1];
	matrix_x[start_index +3] = -b.col0[3]; // B
	matrix_x[start_index +4] = b.col1[2];
	//
	start_index = index_chol_ix[2];
	matrix_x[start_index   ] = b.col0[2]; // A // column j6+2
	matrix_x[start_index +1] = b.col1[1];
	matrix_x[start_index +2] = b.col2[0];
	matrix_x[start_index +3] = -b.col0[4]; // B
	matrix_x[start_index +4] = -b.col1[2];
	//
	start_index = index_chol_ix[3];
	matrix_x[start_index   ] = b.col3[0];// Btilde   // column j6+3
	matrix_x[start_index +1] = b.col3[1];
	matrix_x[start_index +2] = b.col3[2]; // C
	matrix_x[start_index +3] = b.col3[3];
	matrix_x[start_index +4] = b.col3[4];
	//
	start_index = index_chol_ix[4];
	matrix_x[start_index   ] = -b.col3[0]; // Btilde // column j6+4
	matrix_x[start_index +1] = b.col4[0];
	matrix_x[start_index +2] = b.col3[3]; // C
	matrix_x[start_index +3] = b.col4[1];
	matrix_x[start_index +4] = b.col4[2];
	//
	start_index = index_chol_ix[5];
	matrix_x[start_index   ] = -b.col3[1]; // Btilde // column j6+5
	matrix_x[start_index +1] = -b.col4[0];
	matrix_x[start_index +2] = b.col3[4]; // C
	matrix_x[start_index +3] = b.col4[2];
	matrix_x[start_index +4] = b.col5[0];
}

void StokesSolver::insertBlockColumnIndices(int *matrix_p, const vector<int> &pvalues){
	/**
			Insert the starting indices (pvalues) for the 6 columns corresponding to a column of blocks in a cholmod_sparse::p array.
			You must to gives a pointer to cholmod_sparse::p[first_column] as matrix_p, *not* the bare cholmod_sparse::p
	*/
	for (int col=0; col<6; col++) {
		matrix_p[col] = pvalues[col];
	}
}


void StokesSolver::insertDBlockRows(int *matrix_i, const vector<int> &index_values, int top_row_nb){
	/**
			Insert row numbers for the diagonal block with top row top_row_nb in a cholmod_sparse::i array.
			You need to provide the location of the 6 columns of the block in the cholmod_sparse::i array
			through a vector of indices index_values.
	*/
	for (int col=0; col<6; col++) {
		int slim = dblocks_cntnonzero[col];
		int index_start = index_values[col];
		for (int s=0; s<slim; s++) {
			matrix_i[index_start+s] = top_row_nb + db_layout[col][s];
		}
	}
}

void StokesSolver::insertODBlockRows(int *matrix_i, const vector<int> &index_values, int top_row_nb)
{
	/**
			Insert row numbers for an off-diagonal block with top row top_row_nb in a cholmod_sparse::i array.
			You need to provide the location of the 6 columns of the block in the cholmod_sparse::i array
			through a vector of indices index_values.
	*/
	for (int col=0; col<6; col++) {
		int index_start = index_values[col];
		vector<int> &layout = odb_layout[col];
		for (int s=0; s<5; s++) {
			matrix_i[index_start+s] = top_row_nb + layout[s];
		}
	}
}

void StokesSolver::insertDBlock(cholmod_sparse *matrix, const vector<int> &index_chol_ix, int top_row_nb, const struct DBlock &diagblock)
{
	/**
			Insert the DBlock diagblock in a cholmod_sparse matrix.
			You must provide the location of the 6 columns of the block in the cholmod_sparse::i array
			through a vector of indices index_chol_ix, and the nb of the topmost row (equivalently leftmost column) of the block.
	*/
	insertBlockColumnIndices((int*)matrix->p+top_row_nb, index_chol_ix);
	insertDBlockRows((int*)matrix->i, index_chol_ix, top_row_nb); // left_col_nb is also the top row nb on a diag block
	insertDBlockValues((double*)matrix->x, index_chol_ix, diagblock);

}
void StokesSolver::insertODBlock(cholmod_sparse *matrix, const vector<int> &index_chol_ix, int top_row_nb, const struct ODBlock &offdiagblock)
{
	/**
			Insert the ODBlock offdiagblock in a cholmod_sparse matrix.
			You must provide the location of the 6 columns of the block in the cholmod_sparse::i array
			through a vector of indices index_chol_ix, and the nb of of the topmost row of the block.
	*/
	insertODBlockRows((int*)matrix->i, index_chol_ix, top_row_nb);
	insertODBlockValues((double*)matrix->x, index_chol_ix, offdiagblock);

}


void StokesSolver::completeResistanceMatrix(){
	allocateResistanceMatrix();

	completeResistanceMatrix_MobileMobile();
	factorizeResistanceMatrix();

	// completeResistanceMatrix_MobileFixed();
	// completeResistanceMatrix_FixedFixed();
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

void StokesSolver::completeResistanceMatrix_MobileMobile()
{
	// this function is commented, but you are strongly advised to read
	// the description of storage in the header file first :)

	int size = mobile_particle_nb;

// the vector index_chol_ix tracks the indices in the i and x arrays of the cholmod matrix for the 6 columns
	vector<int> index_chol_ix;
	index_chol_ix.resize(6);


	for (int j=0; j<size; j++) {
		/******* Initialize index_chol_ix for this column of blocks **************/
		// associated with particle j are 6 columns in the matrix:
		// { 6j, ... , 6j+5 }
		// the number of non-zero elements before column 6j is:
		// - 18*j from j diagonal blocks
		// - 30*odbrows_table[j] from odbrows_table[j] off-diagonal blocks
		//
		// the number of non-zero elements before column 6j+1 is:
		// - number of non-zero before column 6j + number of non-zero in column 6*j
		// (in 6j: dblocks_cntnonzero[0] elements in diagonal block, plus 5*(odbrows_table[j+1]-odbrows_table[j])
		//
		// for 6j+2 --> 6j+5: same idea

		int j6 = 6*j;
		int od_nzero_nb = 5*(odbrows_table[j+1]-odbrows_table[j]);
		index_chol_ix[0] = 18*j+30*odbrows_table[j];
		for (int col=1; col<6; col++) { // set the starting indices for cols 0-5
			index_chol_ix[col] = index_chol_ix[col-1]+dblocks_cntnonzero[col-1]+od_nzero_nb; // nb before previous + elements in previous
		}

		/********* 1: Insert the diagonal blocks elements *********/
		insertDBlock(chol_res_matrix, index_chol_ix, j6, dblocks[j]);
		for (int col=0; col<6; col++) {
			index_chol_ix[col] +=  dblocks_cntnonzero[col];
		}
		/********  2  : off-diagonal blocks blocks elements ***********/
		for (int k = odbrows_table[j]; k<odbrows_table[j+1]; k++) {
			insertODBlock(chol_res_matrix, index_chol_ix, odbrows[k], odblocks[k]);
			for (int col=0; col<6; col++) {
				index_chol_ix[col] += 5;// 5 non-zero elements per columns in odblocks
			}
		}
	}
	((int*)chol_res_matrix->p)[6*size] = ((int*)chol_res_matrix->p)[6*size-1]+1;
}




void StokesSolver::resetResistanceMatrix(int nb_of_interactions,
										 const vector<struct DBlock> &reset_resmat_dblocks)
{
	odblocks_nb = nb_of_interactions;

	for (unsigned int i=0; i<dblocks.size(); i++) {
		dblocks[i] = reset_resmat_dblocks[i];
	}
	odbrows.clear();
	odblocks.resize(odblocks_nb);
	for (auto &b: odblocks) {
		resetODBlock(b);
	}
	odbrows_table[0] = 0;
	mobile_matrix_done = false;

	// for the mixed problem
	int odblocks_mf_nb = 0;
	odbrows_mf.clear();
	odblocks_mf.resize(odblocks_mf_nb);
	for (auto &b: odblocks_mf) {
		resetODBlock(b);
	}
	odbrows_table_mf[0] = 0;

	// for the mixed problem
	for (unsigned int i=0; i<dblocks_ff.size(); i++) {
		dblocks[i] = reset_resmat_dblocks[i+mobile_particle_nb];
	}
	int odblocks_ff_nb = 0;
	odbrows_ff.clear();
	odblocks_ff.resize(odblocks_ff_nb);
	for (auto &b: odblocks_ff) {
		resetODBlock(b);
	}
	odbrows_table_ff[0] = 0;
}

void StokesSolver::resetRHS()
{
	int size = chol_rhs->nrow;
	for (int i=0; i<size; i++) {
		((double*)chol_rhs->x)[i] = 0;
	}
}

void StokesSolver::resetRHStorque()
{
	int size = chol_rhs->nrow;
	for (int i=3; i<size; i+=6) {
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

// // Finds solutions to L^T X = RHS
// void StokesSolver::solve_LT(double* X)
// {
// 	chol_PTsolution = cholmod_solve(CHOLMOD_Lt, chol_L, chol_rhs, &chol_c) ;
// 	chol_solution = cholmod_solve(CHOLMOD_Pt, chol_L, chol_PTsolution, &chol_c) ;
// 	int size = chol_solution->nrow;
// 	for (int i=0; i<size; i++) {
// 		X[i] = ((double*)chol_solution->x)[i];
// 	}
// 	cholmod_free_dense(&chol_solution, &chol_c);
// 	cholmod_free_dense(&chol_PTsolution, &chol_c);
// }
//
// void StokesSolver::solve_LT(vec3d* X, vec3d* ang_X)
// {
// 	chol_PTsolution = cholmod_solve(CHOLMOD_Lt, chol_L, chol_rhs, &chol_c) ;
// 	chol_solution = cholmod_solve(CHOLMOD_Pt, chol_L, chol_PTsolution, &chol_c) ;
// 	int size = chol_solution->nrow/6;
// 	for (int i=0; i<size; i++) {
// 		int i6 =6*i;
// 		X[i].x = ((double*)chol_solution->x)[i6];
// 		X[i].y = ((double*)chol_solution->x)[i6+1];
// 		X[i].z = ((double*)chol_solution->x)[i6+2];
// 		ang_X[i].x = ((double*)chol_solution->x)[i6+3];
// 		ang_X[i].y = ((double*)chol_solution->x)[i6+4];
// 		ang_X[i].z = ((double*)chol_solution->x)[i6+5];
// 	}
// 	cholmod_free_dense(&chol_solution, &chol_c);
// 	cholmod_free_dense(&chol_PTsolution, &chol_c);
// }

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
	odbrows_table.resize(mobile_particle_nb+1);
	odbrows_table_mf.resize(mobile_particle_nb+1);
	odbrows_table_ff.resize(np-mobile_particle_nb+1);
	dblocks_ff.resize(np-mobile_particle_nb);
	cholmod_start(&chol_c);
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
	nzmax = 18*dblocks.size(); // diagonal blocks
	nzmax += 30*odblocks_nb;  // off-diagonal
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
	block.col0[3] = -scaledYB*nvec.z;
	block.col0[4] = scaledYB*nvec.y;
	block.col1[0] = scaledXA*n1n1+scaledYA*one_n1n1; // column 1
	block.col1[1] = (scaledXA-scaledYA)*n1n2;
	block.col1[2] = -scaledYB*nvec.x;
	block.col2[0] = scaledXA*n2n2+scaledYA*one_n2n2; // column 2
	block.col3[0] = scaledYBtilde*nvec.z; // column 3
	block.col3[1] = -scaledYBtilde*nvec.y;
	block.col3[2] = scaledYC*one_n0n0;
	block.col3[3] = -scaledYC*n0n1;
	block.col3[4] = -scaledYC*n0n2;
	block.col4[0] = scaledYBtilde*nvec.x; // column 4
	block.col4[1] = scaledYC*one_n1n1;
	block.col4[2] = -scaledYC*n1n2;
	block.col5[0] = scaledYC*one_n2n2; // column 5
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

void StokesSolver::printFactor(ostream& out)
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
