/**
 \class StokesSolver
 \brief Class dealing with the linear algebra.
 \author Ryohei Seto
 \author Romain Mari
 */

//
//  StokesSolver.h
//  LF_DEM
//
//  Created by Romain Mari on 01/29/13.
//

#ifndef __LF_DEM__StokesSolver__
#define __LF_DEM__StokesSolver__
#include <vector>
#include <array>
#include "vec3d.h"
#include "MatrixBlocks.h"
#include "cholmod.h"
// uncomment below to use long integer cholmod (necessary for GPU)
#ifndef USE_GPU
#define CHOL_FUNC(NAME) cholmod_ ## NAME
typedef int chol_int;
#else
#define CHOL_FUNC(NAME) cholmod_l_ ## NAME
typedef long chol_int;
#endif

class StokesSolver{
	/*
	  This class provides solver for the resistance problem, which is solving for the Velocities
		in an equation of motion of the type:

	  ResistanceMatrix*Velocities = SomeForces

	  Here, Velocities are translational and rotational velocities (U and W) (a 6N vector, FT/UW case).
	  It is stored as (U_1, W_1, U_2, W_2, ...., U_N, W_N).
	  Accordingly, Forces are forces and torques ( F and T, stored as (F_1, T_1,  ..., F_N, T_N) ).
	  This solver handles the case for which the ResistanceMatrix is made of two-body short range interactions, i.e. is sparse.

		This solver can also handle the case of a "mixed problem", i. e. a problem in which some velocities are prescribed,
		and only the remaining velocities are unknowns of the resistance problem.

		Internally, it mainly acts as a wrapper around the Cholmod solver of the SuiteSparse library by Tim Davis.

	  ============================
	  ============================
	  The ResistanceMatrix has the following structure (all blocks are 6x6 if UW/FT version):

                         ================== 6N ==================
                         .........        ..........
                    ||  |  6x6   .        . 6x6    .            | ||
                    ||  |ODBlock .   0    . DBlock .            | ||
                    ||  |        .        .        .            | ||
                    ||  |...........................            | ||
                    ||  |        .        .                     | ||
                    ||  |     0  . DBlock .                     | ||
                    6N  |        .        .                     | 6N
 ResistanceMatrix = ||  |..................                     | ||
                    ||  |        .           .                  | ||
                    ||  |ODBlock .             .                | ||
                    ||  |        .               .              | ||
                    ||  |.........                 .            | ||
                    ||  |                            .          | ||
                    ||  |                              .........| ||
                    ||  |                              .        | ||
                    ||  |                              . DBlock | ||
                    ||  |                              .        | ||
                                                       .........
                        =================== 6N =================

	 This matrix is symmetric, so we only store its lower part. This implies DBlocks themselves are symmetric.
	 ODBlocks however are antisymmetric.

	 The non-zero 6x6 blocks are of two kinds:
	 - DBlock's ("diagonal") contain the contributions to the force/torque acting on particle i
	   coming from particle i's velocity/angular velocity itself: they are the Stokes drag terms ("self"-terms)
	   on the diagonal, plus parts of the interaction terms coming from lubrication or contact dashpots.
	 - ODBlock's ("off-diagonal") contain the contributions to the force (torque) acting on particle i
	   coming from other particles velocities and angular velocities.

	  ============================
		How to use
		============================

	  With cholmod library, the storage is first done in the "natural" block structures presented above.
	  At the end of filling, the storage is converted to the compressed-column form used by cholmod.

		For technical details on the storage of matrix elements in the DBlock's and ODBlock's, #
		see the MatrixBlocks.h header.

		Solving the resistance problem is done in two steps:
		1. Building the resistance matrix
		2. Build the rhs (the forces)
		3. Actually solving the linear algebra problem
		(Steps 1 and 2 can be interchanged, but in any case resetResistanceMatrix must be called before anything else).

		Detail of step 1, building the resistance matrix and the RHS
			a. start a new instance of resistance matrix with the method resetResistanceMatrix.
			b. filling it. The resistance matrix is filled columns by columns, from column 0 up, IN ORDER.
			    Declare the beginning of a new column with startNewColumn, and add the terms linked
			    to a pair interaction with addResistanceBlocks.
					NOTE: within a column, you do not need to fill the rows in order.
			c. call matrixFillingDone. This exports the temporary data storage into cholmod format,
			   and factorizes the matrix. Therefore, after calling this method, you cannot modify
			   the resistance matrix anymore.

		Step 2 is done with the *RHS* methods.
		Step 3 is done with the solve methods. After solving, solvingIsDone must be called.
		Once solvingIsDone has been called, the cholmod matrix object has been freed,
		so you cannot use the solver anymore without starting the process all over from step 1.
	 */

private:
	chol_int np;
	chol_int mobile_particle_nb;
	chol_int fixed_particle_nb;
	bool mobile_matrix_done;
	bool to_be_factorized;

	// Cholmod variables
	cholmod_factor* chol_L ;
	cholmod_common chol_c ;
	cholmod_dense* chol_rhs;
	cholmod_sparse* chol_res_matrix;
	cholmod_sparse* chol_res_matrix_mf;
	cholmod_sparse* chol_res_matrix_ff;
	cholmod_dense* chol_solution;
	cholmod_dense* chol_solveE_workspace;
	cholmod_dense* chol_solveY_workspace;
	cholmod_dense* chol_vel_mob;
	cholmod_dense* chol_vel_fix;
	cholmod_dense* chol_force_mob;
	cholmod_dense* chol_force_fix;
	// cholmod_dense* chol_PTsolution;
	cholmod_dense* chol_Psolution;
	// resistance matrix building
	chol_int dblocks_size;
	chol_int current_column;
	std::vector<struct DBlock> dblocks;
	std::vector<struct ODBlock> odblocks;
	std::vector<chol_int> odbrows;
	std::vector<chol_int> odbrows_table;
	std::vector<struct ODBlock> odblocks_mf;
	std::vector<chol_int> odbrows_mf;
	std::vector<chol_int> odbrows_table_mf;
	std::vector<struct DBlock> dblocks_ff;
	std::vector<struct ODBlock> odblocks_ff;
	std::vector<chol_int> odbrows_ff;
	std::vector<chol_int> odbrows_table_ff;
	std::vector<std::vector <chol_int> > odb_layout;
	std::vector<std::vector <chol_int> > db_layout;
	std::vector<chol_int> dblocks_cntnonzero;

	void factorizeResistanceMatrix();
	void allocateResistanceMatrix();
	void allocateRessources();
	void completeResistanceMatrix_MobileMobile();
	void completeResistanceMatrix_MobileFixed();
	void completeResistanceMatrix_FixedFixed();
	void insertODBlock(cholmod_sparse *matrix,
					   const std::vector<chol_int> &index_chol_ix,
					   chol_int top_row_nb,
					   const struct ODBlock &offdiagblock);
	void insertDBlock(cholmod_sparse *matrix,
					  const std::vector<chol_int> &index_chol_ix,
					  chol_int top_row_nb,
					  const struct DBlock &diagblock);
	void insertODBlockRows(chol_int *matrix_i,
						   const std::vector<chol_int> &index_values,
						   chol_int top_row_nb);
	void insertDBlockRows(chol_int *matrix_i,
						  const std::vector<chol_int> &index_values,
						  chol_int top_row_nb);
	void insertBlockColumnIndices(chol_int *matrix_p,
								  const std::vector<chol_int> &pvalues);
	void insertODBlockValues(double *matrix_x,
							 const std::vector<chol_int>& index_chol_ix,
							 const struct ODBlock& b);
	void insertDBlockValues(double *matrix_x,
							const std::vector<chol_int>& index_chol_ix,
							const struct DBlock& b);
	void setOffDiagBlock(chol_int jj, const struct ODBlock& b);
	void addToDiagBlock(chol_int ii, const struct DBlock &b);
	/*
	 completeResistanceMatrix() :
	 - transforms temporary arrays/vectors used to build resistance
	 matrix into Cholmod objects really used by
	 the solvers
	 - also performs Cholesky factorization
	 - must be called before solving the linear algebra problem
	 - must be called after all terms are added
	 */
	void completeResistanceMatrix();

public:
	StokesSolver();
	~StokesSolver();
	void init(chol_int np_total, chol_int np_mobile);
	void printResistanceMatrix(std::ostream&, std::string);
	void printFactor(std::ostream&);
	void printRHS();
	void convertDirectToIterative();
	// R_FU filling methods
	/* resetResistanceMatrix(string solver_type, chol_int nb_of_interactions) :
	 - initialize arrays/vectors used for building
	 - to be called before adding elements
	 - nb_of_interactions is the number of odblocks in the matrix
	 */
	void resetResistanceMatrix(chol_int nb_of_interactions_mm,
							   chol_int nb_of_interactions_mf,
							   chol_int nb_of_interactions_ff,
							   const std::vector<struct DBlock>& reset_resmat_dblocks,
							   bool matrix_pattern_changed);
	void addResistanceBlocks(chol_int i,
							 chol_int j,
							 const std::pair<struct DBlock, struct DBlock> &DiagBlocks_i_and_j,
							 const struct ODBlock& ODBlock_ij);
	/*
	 This must be called with order (ii < jj),
	 because we have to fill according to the lower-triangular
	 storage.
	 */
	void startNewColumn();
	void matrixFillingDone();
	/*
	   Right-hand vector access methods
	*/
	void resetRHS();
	void resetRHStorque();
	void addToRHS(double*);
	void addToRHS(const std::vector<double>&);
	void addToRHS(chol_int, const std::vector<double>&);
	void addToRHSForce(chol_int, const vec3d&);
	void addToRHSTorque(chol_int, const vec3d&);
	void setRHS(const std::vector<vec3d>&);
	void setRHS(const std::vector<vec3d>&, const std::vector<vec3d>& torque);
	void setRHSForce(chol_int, const vec3d&);
	void setRHSTorque(chol_int, const vec3d&);
	/*
	 solve(vec3d* velocity, vec3d* ang_velocity) :
	 - once the resistance matrix and the RHS vector are built
	 ( completeResistanceMatrix() must have been called )
	 - solves Resistance * velocity = RHS, and stores it in velocity array
	 */
	 void solve(vec3d* velocity, vec3d* ang_velocity);
	 void solve(std::vector<vec3d> &velocity, std::vector<vec3d> &ang_velocity);
	/*
	 compute_LTRHS(double* X) :
	 - once the resistance matrix and the RHS vector are built
	 ( completeResistanceMatrix() must have been called )
	 - solves X = L^T * RHS, and stores it in X array,
	 where L^t is the transpose of the Cholesky factor ( ResistanceMatrix = L L^t )
	 - works only for direct solver, as we need the Cholesky factor
	 */
	 void compute_LTRHS(std::vector<vec3d>&, std::vector<vec3d>&);
	 /*
	 solvingIsDone(bool free_Cholesky_factor) :
	 - deletes resistance matrix and some other arrays
	 (preconditionner, Cholesky factors, ...) used for solving
	 - should be called once every solve() call WITH THE SAME RESISTANCE MATRIX
	 is done. If matrix changes, solve() must be followed by
	 solvingIsDone()
	 */
	void solvingIsDone();
	void doubleToVec3d(double *a, std::vector<vec3d>& b, std::vector<vec3d>& c);
	void vec3dToDouble(double *a, const std::vector<vec3d>& b, const std::vector<vec3d>& c);

	void multiply_by_RFU_mm(std::vector<double>& velocity, std::vector<double>& force);
	void multiply_by_RFU_mf(std::vector<double>& velocity, std::vector<double>& force);
	void multiply_by_RFU_fm(std::vector<double>& velocity, std::vector<double>& force);
	void multiply_by_RFU_fm(std::vector<vec3d>& velocity,
							std::vector<vec3d>& ang_velocity,
							std::vector<vec3d>& force,
							std::vector<vec3d>& torque);
	void multiply_by_RFU_ff(std::vector<double>& velocity, std::vector<double>& force);
	void multiply_by_RFU_ff(std::vector<vec3d>& velocity,
							std::vector<vec3d>& ang_velocity,
							std::vector<vec3d>& force,
							std::vector<vec3d>& torque);
	// testing functions
	void multiplyByResMat(double *vec);
	void multiplySolutionByResMat(double *vec);
};

#endif /* defined(__LF_DEM__StokesSolver__) */
