#ifndef __LF_DEM__ResistanceMatrix__
#define __LF_DEM__ResistanceMatrix__

#include <vector>
#include "cholmod.h"
#include "MatrixBlocks.h"

class ResistanceMatrix{
private:
  cholmod_common* chol_c ;
  cholmod_sparse* chol_res_matrix;
  int matrix_size;
  int block_size;
  std::vector<struct MatrixBlocks::DBlock> dblocks;
	std::vector<struct MatrixBlocks::ODBlock> odblocks;
	std::vector<int> odbrows;
	std::vector<int> odbrows_table;
  std::vector<std::vector <int> > odb_layout;
  std::vector<std::vector <int> > db_layout;
  std::vector<int> dblocks_cntnonzero;
  int current_column;

public:
  ResistanceMatrix():chol_res_matrix(NULL) {}
  void init(int linear_block_size, cholmod_common* chol_c);
  void setOffDiagBlock(int block_index, const struct MatrixBlocks::ODBlock& b);
  void allocateCholmodMatrix();
  void reset(int nb_of_interactions,
  					 const std::vector<struct MatrixBlocks::DBlock> &reset_resmat_dblocks);
  void insertDBlockValues(const std::vector<int> &index_chol_ix,
                          const struct MatrixBlocks::DBlock& b);
  void insertODBlockValues(const std::vector<int>& index_chol_ix,
                           const struct MatrixBlocks::ODBlock& b);
  void insertBlockColumnIndices(int *matrix_p, const std::vector<int>& pvalues);
  void insertDBlockRows(const std::vector<int>& index_values,
                        int top_row_nb);
  void insertODBlockRows(const std::vector<int>& index_values,
                         int top_row_nb);
  void insertDBlock(const std::vector<int>& index_chol_ix,
  								int top_row_nb, const struct MatrixBlocks::DBlock& diagblock);
  void insertODBlock(const std::vector<int>& index_chol_ix,
  								 int top_row_nb, const struct MatrixBlocks::ODBlock& offdiagblock);
  void completeResistanceMatrix();
  void startNewColumn();
  void matrixFillingDone();
  void print(std::ostream& out, std::string sparse_or_dense);
};

#endif /* defined(__LF_DEM__StokesSolver__) */
