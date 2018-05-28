#include <fstream>
#include "StokesSolver.h"
#ifdef USE_GPU
#include <cuda_runtime_api.h>
#endif
#ifndef CUDA_DEVICE
#define CUDA_DEVICE 0
#endif

using namespace std;

/******************************************************
 *                                                     *
 *                   Public Methods                    *
 *                                                     *
 ******************************************************/

StokesSolver::StokesSolver():
chol_L(nullptr),
chol_rhs(nullptr),
chol_res_matrix(nullptr),
chol_res_matrix_mf(nullptr),
chol_res_matrix_ff(nullptr),
chol_solution(nullptr),
chol_solveE_workspace(nullptr),
chol_solveY_workspace(nullptr),
chol_vel_mob(nullptr),
chol_vel_fix(nullptr),
chol_force_mob(nullptr),
chol_force_fix(nullptr),
chol_Psolution(nullptr)
{
	CHOL_FUNC(start) (&chol_c);
}

StokesSolver::~StokesSolver()
{
	if (chol_solution) {
		CHOL_FUNC(free_dense) (&chol_solution, &chol_c);
	}
	if (chol_solveE_workspace) {
		CHOL_FUNC(free_dense) (&chol_solveE_workspace, &chol_c);
	}
	if (chol_solveY_workspace) {
		CHOL_FUNC(free_dense) (&chol_solveY_workspace, &chol_c);
	}
	if (chol_rhs) {
		CHOL_FUNC(free_dense) (&chol_rhs, &chol_c);
	}
	if (chol_vel_mob) {
		CHOL_FUNC(free_dense) (&chol_vel_mob, &chol_c);
	}
	if (chol_vel_fix) {
		CHOL_FUNC(free_dense) (&chol_vel_fix, &chol_c);
	}
	if (chol_force_mob) {
		CHOL_FUNC(free_dense) (&chol_force_mob, &chol_c);
	}
	if (chol_force_fix) {
		CHOL_FUNC(free_dense) (&chol_force_fix, &chol_c);
	}
	if (chol_res_matrix) {
		CHOL_FUNC(free_sparse) (&chol_res_matrix, &chol_c);
	}
	if (chol_res_matrix_mf) {
		CHOL_FUNC(free_sparse) (&chol_res_matrix_mf, &chol_c);
	}
	if (chol_res_matrix_ff) {
		CHOL_FUNC(free_sparse) (&chol_res_matrix_ff, &chol_c);
	}
	if (chol_L) {
		CHOL_FUNC(free_factor) (&chol_L, &chol_c);
	}
	CHOL_FUNC(finish) (&chol_c);
}

void StokesSolver::init(chol_int np_total, chol_int np_mobile)
{
	np = np_total;
	mobile_particle_nb = np_mobile;
	fixed_particle_nb = np_total-mobile_particle_nb;
	// resistance matrix characteristics (see header for matrix description)
	dblocks_size = 18*mobile_particle_nb;
	to_be_factorized = true;
	allocateRessources();

	odb_layout.resize(6);
	odb_layout[0] = {0,1,2,4,5};
	odb_layout[1] = {0,1,2,3,5};
	odb_layout[2] = {0,1,2,3,4};
	odb_layout[3] = {1,2,3,4,5};
	odb_layout[4] = {0,2,3,4,5};
	odb_layout[5] = {0,1,3,4,5};
	db_layout.resize(6);
	db_layout[0] = {0,1,2,4,5};
	db_layout[1] = {1,2,3,5};
	db_layout[2] = {2,3,4};
	db_layout[3] = {3,4,5};
	db_layout[4] = {4,5};
	db_layout[5] = {5};

	dblocks_cntnonzero.resize(6);
	for (unsigned i=0; i<6; i++) {
		dblocks_cntnonzero[i] = (int)db_layout[i].size();
	}
}

/************* Matrix filling methods **********************/

void StokesSolver::addToDiagBlock(chol_int ii, const struct DBlock &b)
{
	if (ii >= mobile_particle_nb) {
		// FF matrix
		dblocks_ff[ii-mobile_particle_nb] += b;
	} else {
		// MM matrix
		dblocks[ii] += b;
	}
}

void StokesSolver::addResistanceBlocks(chol_int i,
                                       chol_int j,
                                       const std::pair<struct DBlock, struct DBlock> &DiagBlocks_i_and_j,
                                       const struct ODBlock& ODBlock_ij)
{
	addToDiagBlock(i, DiagBlocks_i_and_j.first);
	addToDiagBlock(j, DiagBlocks_i_and_j.second);
	setOffDiagBlock(j, ODBlock_ij);
}

// Off-Diagonal Blocks Terms, FT/UW version
void StokesSolver::setOffDiagBlock(chol_int jj, const struct ODBlock& b)
{
	if (mobile_matrix_done) {
		// FF Matrix
		odbrows_ff.push_back(6*(jj-mobile_particle_nb));
		auto i = odbrows_ff.size()-1;
		odblocks_ff[i] = b;
	} else {
		if (jj < mobile_particle_nb) {
			// MM matrix
			odbrows.push_back(6*jj);
			auto i = odbrows.size()-1;
			odblocks[i] = b;
		} else {
			// MF matrix
			odbrows_mf.push_back(6*(jj-mobile_particle_nb));
			auto i = odbrows_mf.size()-1;
			odblocks_mf[i] = b;
		}
	}
}

void StokesSolver::insertDBlockValues(double *matrix_x,
									  const vector<chol_int>& index_chol_ix,
									  const struct DBlock& b)
{
	auto pcol0 = index_chol_ix[0];
	auto pcol1 = index_chol_ix[1];
	auto pcol2 = index_chol_ix[2];
	auto pcol3 = index_chol_ix[3];
	auto pcol4 = index_chol_ix[4];
	auto pcol5 = index_chol_ix[5];

	// diagonal blocks row values (21)
	matrix_x[pcol0  ] =  b.col0[0];   // column j6
	matrix_x[pcol0+1] =  b.col0[1];
	matrix_x[pcol0+2] =  b.col0[2];
	matrix_x[pcol0+3] =  b.col0[3];
	matrix_x[pcol0+4] =  b.col0[4];
	matrix_x[pcol1  ] =  b.col1[0];   // column j6+1
	matrix_x[pcol1+1] =  b.col1[1];
	matrix_x[pcol1+2] = -b.col0[3];  // anti-symmetry
	matrix_x[pcol1+3] =  b.col1[2];
	matrix_x[pcol2  ] =  b.col2[0];   // column j6+2
	matrix_x[pcol2+1] = -b.col0[4];   // anti-symmetry
	matrix_x[pcol2+2] = -b.col1[2];   // anti-symmetry
	matrix_x[pcol3  ] =  b.col3[0];   // column j6+3
	matrix_x[pcol3+1] =  b.col3[1];
	matrix_x[pcol3+2] =  b.col3[2];
	matrix_x[pcol4  ] =  b.col4[0];   // column j6+4
	matrix_x[pcol4+1] =  b.col4[1];
	matrix_x[pcol5  ] =  b.col5[0];   // column j6+5
}

void StokesSolver::insertODBlockValues(double *matrix_x,
									   const vector<chol_int>& index_chol_ix,
									   const struct ODBlock& b)
{
	auto start_index = index_chol_ix[0];
	matrix_x[start_index   ] =  b.col0[0]; // A   // column j6
	matrix_x[start_index +1] =  b.col0[1];
	matrix_x[start_index +2] =  b.col0[2];
	matrix_x[start_index +3] =  b.col0[3]; // B
	matrix_x[start_index +4] =  b.col0[4];
	//
	start_index = index_chol_ix[1];
	matrix_x[start_index   ] =  b.col0[1]; // A  // column j6+1
	matrix_x[start_index +1] =  b.col1[0];
	matrix_x[start_index +2] =  b.col1[1];
	matrix_x[start_index +3] = -b.col0[3]; // B
	matrix_x[start_index +4] =  b.col1[2];
	//
	start_index = index_chol_ix[2];
	matrix_x[start_index   ] =  b.col0[2]; // A // column j6+2
	matrix_x[start_index +1] =  b.col1[1];
	matrix_x[start_index +2] =  b.col2[0];
	matrix_x[start_index +3] = -b.col0[4]; // B
	matrix_x[start_index +4] = -b.col1[2];
	//
	start_index = index_chol_ix[3];
	matrix_x[start_index   ] =  b.col3[0];// Btilde   // column j6+3
	matrix_x[start_index +1] =  b.col3[1];
	matrix_x[start_index +2] =  b.col3[2]; // C
	matrix_x[start_index +3] =  b.col3[3];
	matrix_x[start_index +4] =  b.col3[4];
	//
	start_index = index_chol_ix[4];
	matrix_x[start_index   ] = -b.col3[0]; // Btilde // column j6+4
	matrix_x[start_index +1] =  b.col4[0];
	matrix_x[start_index +2] =  b.col3[3]; // C
	matrix_x[start_index +3] =  b.col4[1];
	matrix_x[start_index +4] =  b.col4[2];
	//
	start_index = index_chol_ix[5];
	matrix_x[start_index   ] = -b.col3[1]; // Btilde // column j6+5
	matrix_x[start_index +1] = -b.col4[0];
	matrix_x[start_index +2] =  b.col3[4]; // C
	matrix_x[start_index +3] =  b.col4[2];
	matrix_x[start_index +4] =  b.col5[0];
}

void StokesSolver::insertBlockColumnIndices(chol_int *matrix_p, const vector<chol_int>& pvalues)
{
	/**
	 Insert the starting indices (pvalues) for the 6 columns corresponding to a column of blocks in a cholmod_sparse::p array.
	 You must give a pointer to cholmod_sparse::p[first_column] as matrix_p, *not* the bare cholmod_sparse::p
	 */
	for (unsigned col=0; col<6; col++) {
		matrix_p[col] = pvalues[col];
	}
}

void StokesSolver::insertDBlockRows(chol_int *matrix_i, const vector<chol_int>& index_values, chol_int top_row_nb)
{
	/**
	 Insert row numbers for the diagonal block with top row top_row_nb in a cholmod_sparse::i array.
	 You need to provide the location of the 6 columns of the block in the cholmod_sparse::i array
	 through a vector of indices index_values.
	 */
	for (unsigned col=0; col<6; col++) {
		auto slim = dblocks_cntnonzero[col];
		auto index_start = index_values[col];
		for (decltype(slim) s=0; s<slim; s++) {
			matrix_i[index_start+s] = top_row_nb + db_layout[col][s];
		}
	}
}

void StokesSolver::insertODBlockRows(chol_int *matrix_i, const vector<chol_int>& index_values, chol_int top_row_nb)
{
	/**
	 Insert row numbers for an off-diagonal block with top row top_row_nb in a cholmod_sparse::i array.
	 You need to provide the location of the 6 columns of the block in the cholmod_sparse::i array
	 through a vector of indices index_values.
	 */
	for (unsigned col=0; col<6; col++) {
		auto index_start = index_values[col];
		auto &layout = odb_layout[col];
		for (unsigned s=0; s<5; s++) {
			matrix_i[index_start+s] = top_row_nb + layout[s];
		}
	}
}

void StokesSolver::insertDBlock(cholmod_sparse *matrix, const vector<chol_int>& index_chol_ix,
								chol_int top_row_nb, const struct DBlock& diagblock)
{
	/**
	 Insert the DBlock diagblock in a cholmod_sparse matrix.
	 You must provide the location of the 6 columns of the block in the cholmod_sparse::i array
	 through a vector of indices index_chol_ix, and the nb of the topmost row (equivalently leftmost column) of the block.
	 */
	insertBlockColumnIndices((chol_int*)matrix->p+top_row_nb, index_chol_ix);
	insertDBlockRows((chol_int*)matrix->i, index_chol_ix, top_row_nb); // left_col_nb is also the top row nb on a diag block
	insertDBlockValues((double*)matrix->x, index_chol_ix, diagblock);
}

void StokesSolver::insertODBlock(cholmod_sparse *matrix, const vector<chol_int>& index_chol_ix,
								 chol_int top_row_nb, const struct ODBlock& offdiagblock)
{
	/**
	 Insert the ODBlock offdiagblock in a cholmod_sparse matrix.
	 You must provide the location of the 6 columns of the block in the cholmod_sparse::i array
	 through a vector of indices index_chol_ix, and the nb of of the topmost row of the block.
	 */
	insertODBlockRows((chol_int*)matrix->i, index_chol_ix, top_row_nb);
	insertODBlockValues((double*)matrix->x, index_chol_ix, offdiagblock);
}

void StokesSolver::completeResistanceMatrix()
{
	allocateResistanceMatrix();

	completeResistanceMatrix_MobileMobile();
	factorizeResistanceMatrix();
	completeResistanceMatrix_MobileFixed();
	completeResistanceMatrix_FixedFixed();
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

	// the vector index_chol_ix tracks the indices in the i and x arrays of the cholmod matrix for the 6 columns
	vector<chol_int> index_chol_ix(6);

	for (decltype(mobile_particle_nb) j=0; j<mobile_particle_nb; j++) {
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

		auto j6 = 6*j;
		auto od_nzero_nb = 5*(odbrows_table[j+1]-odbrows_table[j]);
		index_chol_ix[0] = 18*j+30*odbrows_table[j];
		for (unsigned col=1; col<6; col++) {
			index_chol_ix[col] = index_chol_ix[col-1]+dblocks_cntnonzero[col-1]+od_nzero_nb; // nb before previous + elements in previous
		}
		/********* 1: Insert the diagonal blocks elements *********/
		insertDBlock(chol_res_matrix, index_chol_ix, j6, dblocks[j]);
		for (unsigned col=0; col<6; col++) {
			index_chol_ix[col] += dblocks_cntnonzero[col];
		}
		/********  2  : off-diagonal blocks blocks elements ***********/
		for (chol_int k = odbrows_table[j]; k<odbrows_table[j+1]; k++) {
			insertODBlock(chol_res_matrix, index_chol_ix, odbrows[k], odblocks[k]);
			for (unsigned col=0; col<6; col++) {
				index_chol_ix[col] += 5;// 5 non-zero elements per columns in odblocks
			}
		}
	}
	// tell cholmod where the last column stops
	// nb of non-zero elements in column col_nb-1 is 1 (it is the last element of the diagonal)
	int nzero_nb_last_col = 1;
	((chol_int*)chol_res_matrix->p)[6*mobile_particle_nb] = ((chol_int*)chol_res_matrix->p)[6*mobile_particle_nb-1]+nzero_nb_last_col;
}

void StokesSolver::completeResistanceMatrix_FixedFixed()
{
	// this function is commented, but you are strongly advised to read
	// the description of storage in the header file first :)
	if (mobile_particle_nb == np) {
		return;
	}
	// the vector index_chol_ix tracks the indices in the i and x arrays of the cholmod matrix for the 6 columns
	vector<chol_int> index_chol_ix(6);

	for (decltype(fixed_particle_nb) j=0; j<fixed_particle_nb; j++) {
		/******* Initialize index_chol_ix for this column of blocks **************/
		// associated with particle j are 6 columns in the matrix:
		// { 6j, ... , 6j+5 }
		// the number of non-zero elements before column 6j is:
		// - 18*j from j diagonal blocks
		// - 30*odbrows_table_ff[j] from odbrows_table_ff[j] off-diagonal blocks
		//
		// the number of non-zero elements before column 6j+1 is:
		// - number of non-zero before column 6j + number of non-zero in column 6*j
		// (in 6j: dblocks_cntnonzero[0] elements in diagonal block, plus 5*(odbrows_table_ff[j+1]-odbrows_table_ff[j])
		//
		// for 6j+2 --> 6j+5: same idea

		auto od_nzero_nb = 5*(odbrows_table_ff[j+1]-odbrows_table_ff[j]);
		index_chol_ix[0] = 18*j+30*odbrows_table_ff[j];
		for (unsigned col=1; col<6; col++) {
			index_chol_ix[col] = index_chol_ix[col-1]+dblocks_cntnonzero[col-1]+od_nzero_nb; // nb before previous + elements in previous
		}

		/********* 1: Insert the diagonal blocks elements *********/
		insertDBlock(chol_res_matrix_ff, index_chol_ix, 6*j, dblocks_ff[j]);
		for (unsigned col=0; col<6; col++) {
			index_chol_ix[col] +=  dblocks_cntnonzero[col];
		}
		/********  2  : off-diagonal blocks blocks elements ***********/
		for (auto k = odbrows_table_ff[j]; k<odbrows_table_ff[j+1]; k++) {
			insertODBlock(chol_res_matrix_ff, index_chol_ix, odbrows_ff[k], odblocks_ff[k]);
			for (unsigned col=0; col<6; col++) {
				index_chol_ix[col] += 5;// 5 non-zero elements per columns in odblocks
			}
		}
	}
	// tell cholmod where the last column stops
	// nb of non-zero elements in column col_nb-1 is 1 (it is the last element of the diagonal)
	int nzero_nb_last_col = 1;
	((chol_int*)chol_res_matrix_ff->p)[6*fixed_particle_nb] = ((chol_int*)chol_res_matrix_ff->p)[6*fixed_particle_nb-1]+nzero_nb_last_col;
}

void StokesSolver::completeResistanceMatrix_MobileFixed()
{
	// this function is commented, but you are strongly advised to read
	// the description of storage in the header file first :)
	if (fixed_particle_nb == 0) {
		return;
	}
	// the vector index_chol_ix tracks the indices in the i and x arrays of the cholmod matrix for the 6 columns
	vector<chol_int> index_chol_ix(6);
	for (decltype(mobile_particle_nb) j=0; j<mobile_particle_nb; j++) {

		/******* Initialize index_chol_ix for this column of blocks **************/
		// associated with particle j are 6 columns in the matrix:
		// { 6j, ... , 6j+5 }
		// the number of non-zero elements before column 6j is:
		// - 30*odbrows_table_mf[j] from odbrows_table_mf[j] off-diagonal blocks
		//
		// the number of non-zero elements before column 6j+1 is:
		// - number of non-zero before column 6j + number of non-zero in column 6*j
		// (in 6j: 5*(odbrows_table_mf[j+1]-odbrows_table_mf[j])
		//
		// for 6j+2 --> 6j+5: same idea

		auto od_nzero_nb = 5*(odbrows_table_mf[j+1]-odbrows_table_mf[j]);
		index_chol_ix[0] = 30*odbrows_table_mf[j];
		for (unsigned col=1; col<6; col++) {
			index_chol_ix[col] = index_chol_ix[col-1]+od_nzero_nb; // nb before previous + elements in previous
		}
		insertBlockColumnIndices((chol_int*)chol_res_matrix_mf->p+6*j, index_chol_ix);

		for (auto k=odbrows_table_mf[j]; k<odbrows_table_mf[j+1]; k++) {
			insertODBlock(chol_res_matrix_mf, index_chol_ix, odbrows_mf[k], odblocks_mf[k]);
			for (unsigned col=0; col<6; col++) {
				index_chol_ix[col] += 5;// 5 non-zero elements per columns in odblocks
			}
		}
	}
	// tell cholmod where the last column stops
	// nb of non-zero elements in column col_nb-1 is 5*(odbrows_table_mf[col_nb]-odbrows_table_mf[col_nb-1]) (as it is an off diagonal matrix)
	auto nzero_nb_last_col = 5*(odbrows_table_mf[mobile_particle_nb]-odbrows_table_mf[mobile_particle_nb-1]);
	// so the last element in chol_res_matrix_mf->p must be
	((chol_int*)chol_res_matrix_mf->p)[6*mobile_particle_nb] = ((chol_int*)chol_res_matrix_mf->p)[6*mobile_particle_nb-1] + nzero_nb_last_col;
}

void StokesSolver::resetResistanceMatrix(chol_int nb_of_interactions_mm,
										 chol_int nb_of_interactions_mf,
										 chol_int nb_of_interactions_ff,
										 const vector<struct DBlock> &reset_resmat_dblocks,
										 bool matrix_pattern_changed)
{
	for (auto i=0u; i<dblocks.size(); i++) {
		dblocks[i] = reset_resmat_dblocks[i];
	}
	odbrows.clear();
	odblocks.resize(nb_of_interactions_mm);
	for (auto &b: odblocks) {
		resetODBlock(b);
	}
	mobile_matrix_done = false;

	// for the mixed problem
	odbrows_mf.clear();
	odblocks_mf.resize(nb_of_interactions_mf);
	for (auto &b: odblocks_mf) {
		resetODBlock(b);
	}

	// for the mixed problem
	for (auto i=0u; i<dblocks_ff.size(); i++) {
		dblocks_ff[i] = reset_resmat_dblocks[i+mobile_particle_nb];
	}
	odbrows_ff.clear();
	odblocks_ff.resize(nb_of_interactions_ff);
	for (auto &b: odblocks_ff) {
		resetODBlock(b);
	}
	current_column = 0;
	to_be_factorized = matrix_pattern_changed;
	if (to_be_factorized && chol_L != nullptr) {
		CHOL_FUNC(free_factor) (&chol_L, &chol_c);
	}
}

void StokesSolver::resetRHS()
{
	auto size = chol_rhs->nrow;
	for (decltype(size) i=0; i<size; i++) {
		((double*)chol_rhs->x)[i] = 0;
	}
}

void StokesSolver::resetRHStorque()
{
	auto size = chol_rhs->nrow;
	for (decltype(size) i=3; i<size; i+=6) {
		((double*)chol_rhs->x)[i  ] = 0;
		((double*)chol_rhs->x)[i+1] = 0;
		((double*)chol_rhs->x)[i+2] = 0;
	}
}

void StokesSolver::addToRHSForce(chol_int i, const vec3d& force_i)
{
	if (i < mobile_particle_nb) {
		auto i6 = 6*i;
		((double*)chol_rhs->x)[i6  ] += force_i.x;
		((double*)chol_rhs->x)[i6+1] += force_i.y;
		((double*)chol_rhs->x)[i6+2] += force_i.z;
	}
}

void StokesSolver::addToRHSTorque(chol_int i, const vec3d& torque_i)
{
	if (i < mobile_particle_nb) {
		auto i6_3 = 6*i+3;
		((double*)chol_rhs->x)[i6_3  ] += torque_i.x;
		((double*)chol_rhs->x)[i6_3+1] += torque_i.y;
		((double*)chol_rhs->x)[i6_3+2] += torque_i.z;
	}
}

void StokesSolver::addToRHS(double* rhs)
{
	auto size = chol_rhs->nrow;
	for (decltype(size) i=0; i<size; i++) {
		((double*)chol_rhs->x)[i] += rhs[i];
	}
}

void StokesSolver::addToRHS(const vector<double>& force)
{
	for (unsigned i=0; i<force.size(); i++) {
		((double*)chol_rhs->x)[i] += force[i];
	}
}

void StokesSolver::addToRHS(chol_int first_particle, const vector<double>& force)
{
	auto shift = 6*first_particle;
	for (unsigned i=0; i<force.size(); i++) {
		((double*)chol_rhs->x)[i+shift] += force[i];
	}
}

void StokesSolver::setRHS(const vector<vec3d>& force_and_torque)
{
	auto size = 2*mobile_particle_nb;
	// if (force_and_torque.size() != size) {
	// 	throw runtime_error("StokesSolver: setRHS with incompatible vector size\n");
	// }
	for (decltype(size) i=0; i<size; i++) {
		auto i3 = 3*i;
		((double*)chol_rhs->x)[i3  ] = force_and_torque[i].x;
		((double*)chol_rhs->x)[i3+1] = force_and_torque[i].y;
		((double*)chol_rhs->x)[i3+2] = force_and_torque[i].z;
	}
}

void StokesSolver::setRHS(const vector<vec3d>& force,
                          const vector<vec3d>& torque)
{
	auto size = mobile_particle_nb;
	// if (force_and_torque.size() != size) {
	// 	throw runtime_error("StokesSolver: setRHS with incompatible vector size\n");
	// }
	for (decltype(size) i=0; i<size; i++) {
		auto i6 = 6*i;
		((double*)chol_rhs->x)[i6  ] = force[i].x;
		((double*)chol_rhs->x)[i6+1] = force[i].y;
		((double*)chol_rhs->x)[i6+2] = force[i].z;
		((double*)chol_rhs->x)[i6+3] = torque[i].x;
		((double*)chol_rhs->x)[i6+4] = torque[i].y;
		((double*)chol_rhs->x)[i6+5] = torque[i].z;
	}
}

void StokesSolver::setRHSForce(chol_int i, const vec3d& force_i)
{
	if (i < mobile_particle_nb) {
		auto i6 = 6*i;
		((double*)chol_rhs->x)[i6  ] = force_i.x;
		((double*)chol_rhs->x)[i6+1] = force_i.y;
		((double*)chol_rhs->x)[i6+2] = force_i.z;
	}
}

void StokesSolver::setRHSTorque(chol_int i, const vec3d& torque_i)
{
	if (i < mobile_particle_nb) {
		auto i6_3 = 6*i+3;
		((double*)chol_rhs->x)[i6_3  ] = torque_i.x;
		((double*)chol_rhs->x)[i6_3+1] = torque_i.y;
		((double*)chol_rhs->x)[i6_3+2] = torque_i.z;
	}
}

// Computes X = L*RHS
void StokesSolver::compute_LTRHS(vector<vec3d> &F, vector<vec3d> &T)
{
	/*
	 Cholmod gives a factorizationof a permutated resistance
	 matrix Lc*Lc^T = P*RFU*P^T

	 That means P*L = Lc, with  L*L^T = RFU

	 So for a rhs Y:
	 X = L*Y = P^T*Lc*Y
		*/
	if (!chol_L->is_ll) {
		cerr << " The factorization is LDL^T. compute_LTRHS(F, T) only works for LL^T factorization." << endl;
	}
	double alpha[] = {1, 0};
	double beta[] = {0, 0};
	int transpose = 0;
	cholmod_factor* chol_L_copy = CHOL_FUNC(copy_factor) (chol_L, &chol_c);
	cholmod_sparse* chol_L_sparse = CHOL_FUNC(factor_to_sparse) (chol_L_copy, &chol_c);
	CHOL_FUNC(sdmult) (chol_L_sparse, transpose, alpha, beta, chol_rhs, chol_Psolution, &chol_c); // chol_Psolution = Lc*Y
	// chol_solution = P^T*chol_Psolution
	CHOL_FUNC(solve2) (CHOLMOD_Pt,                // solve chol_res_mat*chol_solution = chol_rhs
	                   chol_L,                   // Cholesky factor
	                   chol_Psolution,           // RHS
	                   NULL,                     // subpart of RHS, NULL means all
	                   &chol_solution,           // solution stored here
	                   NULL,                     // subpart of solution, NULL means all
	                   &chol_solveY_workspace,   // reusable workspace
	                   &chol_solveE_workspace,   // reusable workspace
	                   &chol_c);
	auto size = chol_solution->nrow/6;
	for (decltype(size) i=0; i<size; i++) {
		auto i6 = 6*i;
		F[i].x = ((double*)chol_solution->x)[i6  ];
		F[i].y = ((double*)chol_solution->x)[i6+1];
		F[i].z = ((double*)chol_solution->x)[i6+2];
		T[i].x = ((double*)chol_solution->x)[i6+3];
		T[i].y = ((double*)chol_solution->x)[i6+4];
		T[i].z = ((double*)chol_solution->x)[i6+5];
	}
	CHOL_FUNC(free_sparse) (&chol_L_sparse, &chol_c);
	CHOL_FUNC(free_factor) (&chol_L_copy, &chol_c);
}

void StokesSolver::solve(vector<vec3d> &velocity, vector<vec3d> &ang_velocity)
{
	CHOL_FUNC(solve2) (CHOLMOD_A,                // solve chol_res_mat*chol_solution = chol_rhs
	                   chol_L,                   // Cholesky factor
	                   chol_rhs,                 // RHS
	                   NULL,                     // subpart of RHS, NULL means all
	                   &chol_solution,           // solution stored here
	                   NULL,                     // subpart of solution, NULL means all
	                   &chol_solveY_workspace,   // reusable workspace
	                   &chol_solveE_workspace,   // reusable workspace
	                   &chol_c);
	auto size = chol_solution->nrow/6;
	if (size>velocity.size() || size>ang_velocity.size()) {
		// we don't try to resize things here, left to the caller
		// In particular this is sometimes called with size < velo.size(), so no resizing is safer.
		throw runtime_error(" StokesSolver::solve: allocated solution vector too small");
	}
	double *solx = ((double*)chol_solution->x);
	for (decltype(size) i=0; i<size; i++) {
		auto i6 = 6*i;
		velocity[i].x     = solx[i6  ];
		velocity[i].y     = solx[i6+1];
		velocity[i].z     = solx[i6+2];
		ang_velocity[i].x = solx[i6+3];
		ang_velocity[i].y = solx[i6+4];
		ang_velocity[i].z = solx[i6+5];
	}
#ifdef DEV
#ifdef USE_GPU
  cholmod_l_gpu_stats(&chol_c); // debug only
#endif
#endif
}

void StokesSolver::solve2D(vector<vec3d> &velocity, vector<vec3d> &ang_velocity)
{
	CHOL_FUNC(solve2) (CHOLMOD_A,                // solve chol_res_mat*chol_solution = chol_rhs
	                   chol_L,                   // Cholesky factor
	                   chol_rhs,                 // RHS
	                   NULL,                     // subpart of RHS, NULL means all
	                   &chol_solution,           // solution stored here
	                   NULL,                     // subpart of solution, NULL means all
	                   &chol_solveY_workspace,   // reusable workspace
	                   &chol_solveE_workspace,   // reusable workspace
	                   &chol_c);
	auto size = chol_solution->nrow/6;
	if (size>velocity.size() || size>ang_velocity.size()) {
		// we don't try to resize things here, left to the caller
		// In particular this is sometimes called with size < velo.size(), so no resizing is safer.
		throw runtime_error(" StokesSolver::solve: allocated solution vector too small");
	}
	double *solx = ((double*)chol_solution->x);
	for (decltype(size) i=0; i<size; i++) {
		auto i6 = 6*i;
		velocity[i].x     = solx[i6  ];
		velocity[i].y     = 0;
		velocity[i].z     = solx[i6+2];
		ang_velocity[i].x = 0;
		ang_velocity[i].y = solx[i6+4];
		ang_velocity[i].z = 0;
	}
#ifdef DEV
#ifdef USE_GPU
  cholmod_l_gpu_stats(&chol_c); // debug only
#endif
#endif
}

void StokesSolver::vec3dToDouble(double *a, const vector<vec3d>& b, const vector<vec3d>& c)
{
	for (unsigned i=0; i<b.size(); i++) {
		auto i6 = 6*i;
		a[i6  ] = b[i].x;
		a[i6+1] = b[i].y;
		a[i6+2] = b[i].z;
		a[i6+3] = c[i].x;
		a[i6+4] = c[i].y;
		a[i6+5] = c[i].z;
	}
}

void StokesSolver::doubleToVec3d(double *a, vector<vec3d>& b, vector<vec3d>& c)
{
	for (unsigned i=0; i<b.size(); i++) {
		auto i6 = 6*i;
		b[i].x = a[i6  ];
		b[i].y = a[i6+1];
		b[i].z = a[i6+2];
		c[i].x = a[i6+3];
		c[i].y = a[i6+4];
		c[i].z = a[i6+5];
	}
}

void StokesSolver::multiply_by_RFU_mm(vector<double>& velocity, vector<double>& force)
{
	double one[] = {1, 0};
	double zero[] = {0, 0};
	// chol_vel_mob->x = velocity.data();
	// see multiply_by_RFU_mf for the rationale about copying the data, not the pointer
	auto size = chol_vel_mob->nrow;
	for (decltype(size) i=0; i<size; i++) {
		((double*)chol_vel_mob->x)[i] = velocity[i];
	}
	CHOL_FUNC(sdmult) (chol_res_matrix, 1, one, zero, chol_vel_mob, chol_force_mob, &chol_c);
	for (unsigned i=0; i<force.size(); i++) {
		force[i] = ((double*)chol_force_mob->x)[i];
	}
}

void StokesSolver::multiply_by_RFU_mf(vector<double>& velocity, vector<double>& force)
{
	double one[] = {1, 0};
	double zero[] = {0, 0};
	// chol_vel_fix->x = velocity.data();
	// yes, it is evil!!! :). If you call this with a function having vector<double> velocity
	// as a local variable, when the function returns velocity is freed and so chol_vel_fix->x is.
	// The StokesSolver has no way to know about this, so on the next call to StokesSolver
	// using chol_vel_fix->x, chaos...
	// [ Note that this very function is fine as it will reassign a valid pointer to chol_vel_fix->x,
	// but other functions will not (that's how I noticed it was evil :)) ]
	auto size = chol_vel_fix->nrow;
	for (decltype(size) i=0; i<size; i++) {
		((double*)chol_vel_fix->x)[i] = velocity[i];
	}
	CHOL_FUNC(sdmult) (chol_res_matrix_mf, 1, one, zero, chol_vel_fix, chol_force_mob, &chol_c);
	for (unsigned i=0; i<force.size(); i++) {
		force[i] = ((double*)chol_force_mob->x)[i];
	}
}

void StokesSolver::multiply_by_RFU_fm(vector<double>& velocity, vector<double>& force)
{
	double one[] = {1, 0};
	double zero[] = {0, 0};
	// chol_vel_mob->x = velocity.data(); // see multiply_by_RFU_mf for the rationale about copying the data, not the pointer
	auto size = chol_vel_mob->nrow;
	for (decltype(size) i=0u; i<size; i++) {
		((double*)chol_vel_mob->x)[i] = velocity[i];
	}
	CHOL_FUNC(sdmult) (chol_res_matrix_mf, 0, one, zero, chol_vel_mob, chol_force_fix, &chol_c);
	for (unsigned i=0; i<force.size(); i++) {
		force[i] = ((double*)chol_force_fix->x)[i];
	}
}

void StokesSolver::multiply_by_RFU_fm(vector<vec3d>& velocity,
									  vector<vec3d>& ang_velocity,
									  vector<vec3d>& force,
									  vector<vec3d>& torque)
{
	double one[] = {1, 0};
	double zero[] = {0, 0};
	vec3dToDouble((double*)chol_vel_mob->x, velocity, ang_velocity);
	CHOL_FUNC(sdmult) (chol_res_matrix_mf, 0, one, zero, chol_vel_mob, chol_force_fix, &chol_c);
	doubleToVec3d((double*)chol_force_fix->x, force, torque);
}

void StokesSolver::multiply_by_RFU_ff(vector<double>& velocity, vector<double>& force)
{
	double one[] = {1, 0};
	double zero[] = {0, 0};
	// chol_vel_fix->x = velocity.data();
	// see multiply_by_RFU_mf for the rationale about copying the data, not the pointer
	auto size = chol_vel_fix->nrow;
	for (decltype(size) i=0; i<size; i++) {
		((double*)chol_vel_fix->x)[i] = velocity[i];
	}
	CHOL_FUNC(sdmult) (chol_res_matrix_ff, 1, one, zero, chol_vel_fix, chol_force_fix, &chol_c);
	for (unsigned i=0; i<force.size(); i++) {
		force[i] = ((double*)chol_force_fix->x)[i];
	}
}

void StokesSolver::multiply_by_RFU_ff(vector<vec3d>& velocity,
                                      vector<vec3d>& ang_velocity,
                                      vector<vec3d>& force,
                                      vector<vec3d>& torque)
{
	double one[] = {1, 0};
	double zero[] = {0, 0};
	vec3dToDouble((double*)chol_vel_fix->x, velocity, ang_velocity);
	CHOL_FUNC(sdmult) (chol_res_matrix_ff, 1, one, zero, chol_vel_fix, chol_force_fix, &chol_c);
	doubleToVec3d((double*)chol_force_fix->x, force, torque);
}

// testing function, don't use it in production code, very slow and unclean
void StokesSolver::multiplySolutionByResMat(double* vec)
{
	CHOL_FUNC(solve2) (CHOLMOD_A,                // solve chol_res_mat*chol_solution = chol_rhs
	               chol_L,                   // Cholesky factor
	               chol_rhs,                 // RHS
	               NULL,                     // subpart of RHS, NULL means all
	               &chol_solution,           // solution stored here
	               NULL,                     // subpart of solution, NULL means all
	               &chol_solveY_workspace,   // reusable workspace
	               &chol_solveE_workspace,   // reusable workspace
								 &chol_c);
	cholmod_dense* r;
	r = CHOL_FUNC(copy_dense) (chol_rhs, &chol_c);
	double one[] = {1, 0};
	double zero[] = {0, 0};
	CHOL_FUNC(sdmult) (chol_res_matrix, 1, one, zero, chol_solution, r, &chol_c);
	auto size = r->nrow;
	for (decltype(size) i=0; i<size; i++) {
		vec[i] = ((double*)r->x)[i];
	}
	CHOL_FUNC(free_dense) (&r, &chol_c);
}

void StokesSolver::solvingIsDone()
{
	CHOL_FUNC(free_sparse) (&chol_res_matrix_mf, &chol_c);
	CHOL_FUNC(free_sparse) (&chol_res_matrix_ff, &chol_c);
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
	odbrows_table_ff.resize(fixed_particle_nb+1);
	dblocks_ff.resize(fixed_particle_nb);
#ifdef USE_GPU
	if (cudaSetDevice(CUDA_DEVICE) == cudaSuccess) {
		cout << " Using CUDA Device : " << CUDA_DEVICE << endl;
		chol_c.useGPU = 1;
	} else {
		cout << " CUDA Device " << CUDA_DEVICE << " not found " << endl;
		chol_c.useGPU = 0;
	};
#endif
	auto size_mm = 6*mobile_particle_nb;
	auto size_ff = 6*fixed_particle_nb;
	chol_rhs       = CHOL_FUNC(allocate_dense) (size_mm, 1, size_mm, CHOLMOD_REAL, &chol_c);
	chol_vel_mob   = CHOL_FUNC(allocate_dense) (size_mm, 1, size_mm, CHOLMOD_REAL, &chol_c);
	chol_vel_fix   = CHOL_FUNC(allocate_dense) (size_ff, 1, size_ff, CHOLMOD_REAL, &chol_c);
	chol_force_mob = CHOL_FUNC(allocate_dense) (size_mm, 1, size_mm, CHOLMOD_REAL, &chol_c);
	chol_force_fix = CHOL_FUNC(allocate_dense) (size_ff, 1, size_ff, CHOLMOD_REAL, &chol_c);
	chol_Psolution = CHOL_FUNC(allocate_dense) (size_mm, 1, size_mm, CHOLMOD_REAL, &chol_c); // used for Brownian motion
	for (decltype(size_mm) i=0; i<size_mm; i++) {
		((double*)chol_rhs->x)[i] = 0;
	}
	chol_L = nullptr;
}

void StokesSolver::allocateResistanceMatrix()
{
	// CHOLMOD parameters
	int stype = -1; // 1 is symmetric, stored upper triangular (UT), -1 is LT
	int sorted = 1;		/* TRUE if columns sorted, FALSE otherwise */
	int packed = 1;		/* TRUE if matrix packed, FALSE otherwise */
	// allocate
	auto size_mm = 6*dblocks.size();
	auto nzmax = 18*dblocks.size(); // diagonal blocks
	nzmax += 30*odblocks.size();   // off-diagonal
	if (!chol_res_matrix) {
		chol_res_matrix = CHOL_FUNC(allocate_sparse) (size_mm, size_mm, nzmax, sorted, packed, stype, CHOLMOD_REAL, &chol_c);
	}
	if (size_mm != chol_res_matrix->nrow){
		CHOL_FUNC(free_sparse) (&chol_res_matrix, &chol_c);
		chol_res_matrix = CHOL_FUNC(allocate_sparse) (size_mm, size_mm, nzmax, sorted, packed, stype, CHOLMOD_REAL, &chol_c);
	}
	if (nzmax > chol_res_matrix->nzmax) {
		CHOL_FUNC(reallocate_sparse) (nzmax, chol_res_matrix, &chol_c);
	}
	auto col_nb = 6*mobile_particle_nb;
	auto row_nb = 6*fixed_particle_nb;
	nzmax = 30*(chol_int)odblocks_mf.size();  // off-diagonal
	int stype_mf = 0; // non-symmetric matrix
	chol_res_matrix_mf = CHOL_FUNC(allocate_sparse) (row_nb, col_nb, nzmax, sorted, packed, stype_mf, CHOLMOD_REAL, &chol_c);

	auto size_ff = 6*fixed_particle_nb;
	nzmax = 18*(chol_int)dblocks_ff.size(); // diagonal blocks
	nzmax += 30*(chol_int)odblocks_ff.size();  // off-diagonal
	chol_res_matrix_ff = CHOL_FUNC(allocate_sparse) (size_ff, size_ff, nzmax, sorted, packed, stype, CHOLMOD_REAL, &chol_c);
}

void StokesSolver::startNewColumn()
{
	if (current_column == mobile_particle_nb) {
		mobile_matrix_done = true;
		odbrows_table[mobile_particle_nb] = (chol_int)odbrows.size();
		odbrows_table_mf[mobile_particle_nb] = (chol_int)odbrows_mf.size();
	}
	if (mobile_matrix_done) {
		odbrows_table_ff[current_column-mobile_particle_nb] = (chol_int)odbrows_ff.size();
	} else {
		odbrows_table[current_column] = (chol_int)odbrows.size();
		odbrows_table_mf[current_column] = (chol_int)odbrows_mf.size();
	}
	current_column++;
}

void StokesSolver::matrixFillingDone()
{
	for (chol_int i=current_column; i<odbrows_table.size(); i++) {
		odbrows_table[i] = (chol_int)odbrows.size();
		odbrows_table_mf[i] = (chol_int)odbrows_mf.size();
	}
	mobile_matrix_done = true;
	for (chol_int i=current_column; i<odbrows_table_ff.size(); i++) {
		odbrows_table_ff[i] = (chol_int)odbrows_ff.size();
	}
	completeResistanceMatrix();
}

void StokesSolver::factorizeResistanceMatrix()
{
	/*reference code */
	chol_c.supernodal = CHOLMOD_SUPERNODAL;

	if (to_be_factorized) {
#ifdef USE_METIS
		chol_c.nmethods = 1;
		// it is the best one on sizes of a few 1000 particles, but may be worth rechecking at some point
		chol_c.method[0].ordering = CHOLMOD_METIS;
#endif
		chol_L = CHOL_FUNC(analyze) (chol_res_matrix, &chol_c);
	}
	CHOL_FUNC(factorize) (chol_res_matrix, chol_L, &chol_c);

	if (chol_c.status) {
		// This should not happen if the matrix is sym-pos-def.
		// If this is used, then it is most probably the sign
		// that the matrix you are passing is not what you think it is.
		cerr << " factorization failed. forcing simplicial algorithm... " << endl;
		chol_c.supernodal = CHOLMOD_SIMPLICIAL;
		if (to_be_factorized) {
			chol_L = CHOL_FUNC(analyze) (chol_res_matrix, &chol_c);
		}
		CHOL_FUNC(factorize) (chol_res_matrix, chol_L, &chol_c);
		cerr << " factorization status " << chol_c.status << " final_ll ( 0 is LDL, 1 is LL ) " << chol_c.final_ll <<endl;
		chol_c.supernodal = CHOLMOD_SUPERNODAL;
	}
}

// testing
void StokesSolver::printResistanceMatrix(ostream& out, string sparse_or_dense)
{
	if (sparse_or_dense == "sparse") {
		//		out << endl<< " chol res " << endl;
		auto size = chol_res_matrix->nrow;
		for (decltype(size) i = 0; i<size; i++) {
			for (auto k =((chol_int*)chol_res_matrix->p)[i]; k<((chol_int*)chol_res_matrix->p)[i+1]; k++) {
				out << i << " " << ((chol_int*)chol_res_matrix->i)[k] << " " << ((double*)chol_res_matrix->x)[k] << endl;
			}
		}
	}
	if (sparse_or_dense == "dense") {
		cholmod_dense* dense_res = CHOL_FUNC(sparse_to_dense) (chol_res_matrix, &chol_c);
		auto size = chol_res_matrix->nrow;
		for (decltype(size) i = 0; i<size; i++) {
			// if(i==0){
			// 	for (int j = 0; j < size/6; j++) {
			// 		out << j << "\t \t \t \t \t \t" ;
			// 	}
			// 	out << endl;
			// }
			for (decltype(size) j = 0; j<size; j++) {
				out << setprecision(3) << ((double*)dense_res->x)[i+j*size] << "\t" ;
			}
			out << endl;
		}
		out << endl;
		CHOL_FUNC(free_dense) (&dense_res, &chol_c);
	}
}

void StokesSolver::printFactor(ostream& out)
{
	cholmod_factor* chol_L_copy = CHOL_FUNC(copy_factor) (chol_L, &chol_c);
	cholmod_sparse* chol_L_sparse = CHOL_FUNC(transpose) (CHOL_FUNC(factor_to_sparse) (chol_L_copy, &chol_c), 1, &chol_c);
	cholmod_dense* chol_L_dense = CHOL_FUNC(sparse_to_dense) (chol_L_sparse, &chol_c);
	int transpose = 1;
	double alpha[] = {1, 0};
	double beta[] = {0, 0};
	cholmod_dense *dense_res = CHOL_FUNC(sparse_to_dense) (chol_res_matrix, &chol_c);
	CHOL_FUNC(sdmult) (chol_L_sparse, transpose, alpha, beta, chol_L_dense, dense_res, &chol_c);
	auto size = chol_res_matrix->nrow;
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
	for (decltype(size) i = 0; i<size; i++) {
		// if(i==0){
		// 	for (int j = 0; j < size/6; j++) {
		// 		out << j << "\t \t \t \t \t \t" ;
		// 	}
		// 	out << endl;
		// }
		for (decltype(size) j=0; j<size; j++) {
			out << setprecision(3) << ((double*)dense_res->x)[i+j*size] << "\t" ;
		}
		out << endl;
	}
	out << endl;
	CHOL_FUNC(free_sparse) (&chol_L_sparse, &chol_c);
	//		cholmod_free_dense(&chol_LT_dense, &chol_c);
	//		cholmod_free_sparse(&chol_PTL_sparse, &chol_c);
	CHOL_FUNC(free_factor) (&chol_L_copy, &chol_c);
	CHOL_FUNC(free_dense) (&dense_res, &chol_c);
}

// testing
void StokesSolver::printRHS()
{
	auto size = chol_rhs->nrow;
	for (decltype(size) i=0; i<size; i++) {
		cout << i << " (part " << " " << (i-i%6)/6 << " )  " << ((double*)chol_rhs->x)[i] << endl;
	}
}
