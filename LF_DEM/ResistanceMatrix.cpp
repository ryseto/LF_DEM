#include <stdexcept>
#include <iostream>
#include <iomanip>
#include "ResistanceMatrix.h"

void ResistanceMatrix::init(int linear_block_size, cholmod_common *chol_common)
{
	chol_c = chol_common;
	matrix_size = linear_block_size; // nb of blocks

	dblocks.resize(matrix_size);
	odbrows_table.resize(matrix_size + 1);

	db_layout = MatrixBlocks::getDBlockLayout();
	odb_layout = MatrixBlocks::getODBlockLayout();
	block_size = (int)db_layout.size();

	dblocks_cntnonzero.resize(block_size);
	for (decltype(block_size) i=0; i<block_size; i++) {
		dblocks_cntnonzero[i] = (int)db_layout[i].size();
	}
}

void ResistanceMatrix::reset(int nb_of_interactions,
														 const std::vector<struct MatrixBlocks::DBlock> &reset_resmat_dblocks)
{
	if (chol_res_matrix != NULL) {
		cholmod_free_sparse(&chol_res_matrix, chol_c);
	}
	for (auto i=0u; i<dblocks.size(); i++) {
		dblocks[i] = reset_resmat_dblocks[i];
	}
	odbrows.clear();
	odblocks.resize(nb_of_interactions);
	for (auto &b: odblocks) {
		resetODBlock(b);
	}
	current_column = 0;
}

void ResistanceMatrix::allocateCholmodMatrix()
{
	// CHOLMOD parameters
	int stype = -1; // 1 is symmetric, stored upper triangular (UT), -1 is LT
	int sorted = 0;		/* TRUE if columns sorted, FALSE otherwise */
	int packed = 1;		/* TRUE if matrix packed, FALSE otherwise */

	auto size = block_size*matrix_size;

	// nb of non-zero values
	int nz = 18*(int)dblocks.size(); // diagonal blocks
	nz += 30*(int)odblocks.size();   // off-diagonal
	chol_res_matrix = cholmod_allocate_sparse(size, size, nz, sorted, packed, stype, CHOLMOD_REAL, chol_c);
}


void ResistanceMatrix::setOffDiagBlock(int block_index, const struct MatrixBlocks::ODBlock& b)
{
	odbrows.push_back(6*block_index);
	auto i = odbrows.size()-1;
	odblocks[i] = b;
}


void ResistanceMatrix::insertDBlockValues(const std::vector<int> &index_chol_ix,
                                          const struct MatrixBlocks::DBlock& b)
{
	auto pcol0 = index_chol_ix[0];
	auto pcol1 = index_chol_ix[1];
	auto pcol2 = index_chol_ix[2];
	auto pcol3 = index_chol_ix[3];
	auto pcol4 = index_chol_ix[4];
	auto pcol5 = index_chol_ix[5];
	double *matrix_x = (double*)chol_res_matrix->x;
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

void ResistanceMatrix::insertODBlockValues(const std::vector<int>& index_chol_ix,
									   const struct MatrixBlocks::ODBlock& b)
{
	double *matrix_x = (double*)chol_res_matrix->x;

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

void ResistanceMatrix::insertBlockColumnIndices(int *matrix_p, const std::vector<int>& pvalues)
{
	/**
	 Insert the starting indices (pvalues) for the 6 columns corresponding to a column of blocks in a cholmod_sparse::p array.
	 You must give a pointer to cholmod_sparse::p[first_column] as matrix_p, *not* the bare cholmod_sparse::p
	 */
	for (int col=0; col<6; col++) {
		matrix_p[col] = pvalues[col];
	}
}

void ResistanceMatrix::insertDBlockRows(const std::vector<int>& index_values, int top_row_nb)
{
	/**
	 Insert row numbers for the diagonal block with top row top_row_nb in a cholmod_sparse::i array.
	 You need to provide the location of the 6 columns of the block in the cholmod_sparse::i array
	 through a std::vector of indices index_values.
	 */
	int* matrix_i = (int*)chol_res_matrix->i;
	for (int col=0; col<6; col++) {
		auto slim = dblocks_cntnonzero[col];
		auto index_start = index_values[col];
		for (decltype(slim) s=0; s<slim; s++) {
			matrix_i[index_start+s] = top_row_nb + db_layout[col][s];
		}
	}
}

void ResistanceMatrix::insertODBlockRows(const std::vector<int>& index_values, int top_row_nb)
{
	/**
	 Insert row numbers for an off-diagonal block with top row top_row_nb in a cholmod_sparse::i array.
	 You need to provide the location of the 6 columns of the block in the cholmod_sparse::i array
	 through a std::vector of indices index_values.
	 */
	int* matrix_i = (int*)chol_res_matrix->i;
	for (int col=0; col<6; col++) {
		auto index_start = index_values[col];
		auto &layout = odb_layout[col];
		for (int s=0; s<5; s++) {
			matrix_i[index_start+s] = top_row_nb + layout[s];
		}
	}
}

void ResistanceMatrix::insertDBlock(const std::vector<int>& index_chol_ix,
								int top_row_nb, const struct MatrixBlocks::DBlock& diagblock)
{
	/**
	 Insert the DBlock diagblock in a cholmod_sparse matrix.
	 You must provide the location of the 6 columns of the block in the cholmod_sparse::i array
	 through a vector of indices index_chol_ix, and the nb of the topmost row (equivalently leftmost column) of the block.
	 */
	insertBlockColumnIndices((int*)chol_res_matrix->p+top_row_nb, index_chol_ix);
	insertDBlockRows(index_chol_ix, top_row_nb); // left_col_nb is also the top row nb on a diag block
	insertDBlockValues(index_chol_ix, diagblock);
}

void ResistanceMatrix::insertODBlock(const std::vector<int>& index_chol_ix,
								 int top_row_nb, const struct MatrixBlocks::ODBlock& offdiagblock)
{
	/**
	 Insert the ODBlock offdiagblock in a cholmod_sparse matrix.
	 You must provide the location of the 6 columns of the block in the cholmod_sparse::i array
	 through a vector of indices index_chol_ix, and the nb of of the topmost row of the block.
	 */
	insertODBlockRows(index_chol_ix, top_row_nb);
	insertODBlockValues(index_chol_ix, offdiagblock);
}

void ResistanceMatrix::completeResistanceMatrix()
{
	// this function is commented, but you are strongly advised to read
	// the description of storage in the header file first :)

	// the vector index_chol_ix tracks the indices in the i and x arrays of the cholmod matrix for the 6 columns
	std::vector<int> index_chol_ix(6);

	for (decltype(matrix_size) j=0; j<matrix_size; j++) {
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
		for (int col=1; col<6; col++) {
			index_chol_ix[col] = index_chol_ix[col-1]+dblocks_cntnonzero[col-1]+od_nzero_nb; // nb before previous + elements in previous
		}
		/********* 1: Insert the diagonal blocks elements *********/
		insertDBlock(index_chol_ix, j6, dblocks[j]);
		for (int col=0; col<6; col++) {
			index_chol_ix[col] += dblocks_cntnonzero[col];
		}
		/********  2  : off-diagonal blocks blocks elements ***********/
		for (int k = odbrows_table[j]; k<odbrows_table[j+1]; k++) {
			insertODBlock(index_chol_ix, odbrows[k], odblocks[k]);
			for (int col=0; col<6; col++) {
				index_chol_ix[col] += 5;// 5 non-zero elements per columns in odblocks
			}
		}
	}
	// tell cholmod where the last column stops
	// nb of non-zero elements in column col_nb-1 is 1 (it is the last element of the diagonal)
	int last_column = block_size*matrix_size;
	int nzero_nb_last_col = 1;
	((int*)chol_res_matrix->p)[last_column] = ((int*)chol_res_matrix->p)[last_column-1]+nzero_nb_last_col;
}


void ResistanceMatrix::startNewColumn(int nb_of_odblocks)
{
	if (current_column >= odbrows_table.size()) {
		throw std::runtime_error(" Error: filling more columns than initially declared.");
	}
	odbrows_table[current_column] = (unsigned int)odbrows.size();
	current_column++;
}

void ResistanceMatrix::matrixFillingDone()
{
	allocateCholmodMatrix();
	for (unsigned int i=current_column; i<odbrows_table.size(); i++) {
		odbrows_table[i] = (unsigned int)odbrows.size();
	}
	completeResistanceMatrix();
}

// testing
void ResistanceMatrix::print(std::ostream& out, std::string sparse_or_dense)
{
	if (sparse_or_dense == "sparse") {
		//		out << std::endl<< " chol res " << std::endl;
		auto size = chol_res_matrix->nrow;
		for (decltype(size) i = 0; i<size; i++) {
			for (auto k =((int*)chol_res_matrix->p)[i]; k<((int*)chol_res_matrix->p)[i+1]; k++) {
				out << i << " " << ((int*)chol_res_matrix->i)[k] << " " << ((double*)chol_res_matrix->x)[k] << std::endl;
			}
		}
	}
	if (sparse_or_dense == "dense") {
		cholmod_dense* dense_res = cholmod_sparse_to_dense(chol_res_matrix, chol_c);
		auto size = chol_res_matrix->nrow;
		for (decltype(size) i = 0; i<size; i++) {
			// if(i==0){
			// 	for (int j = 0; j < size/6; j++) {
			// 		out << j << "\t \t \t \t \t \t" ;
			// 	}
			// 	out << std::endl;
			// }
			for (decltype(size) j = 0; j<size; j++) {
				out << std::setprecision(3) << ((double*)dense_res->x)[i+j*size] << "\t" ;
			}
			out << std::endl;
		}
		out << std::endl;
		cholmod_free_dense(&dense_res, chol_c);
	}
}
