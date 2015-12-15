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
//#define CHOLMOD_EXTRA
#include <vector>
#include <array>
#include "vec3d.h"
#include "cholmod.h"

struct ODBlock {
	std::array<double, 5> col0;
	std::array<double, 3> col1;
	std::array<double, 1> col2;
	std::array<double, 5> col3;
	std::array<double, 3> col4;
	std::array<double, 1> col5;
	int bla;
};

inline void resetODBlock(struct ODBlock &b)
{
	b.col0.fill(0);
	b.col1.fill(0);
	b.col2.fill(0);
	b.col3.fill(0);
	b.col4.fill(0);
	b.col5.fill(0);
}

struct DBlock {
	std::array<double, 5> col0;
	std::array<double, 3> col1;
	std::array<double, 1> col2;
	std::array<double, 3> col3;
	std::array<double, 2> col4;
	std::array<double, 1> col5;
};

inline void resetDBlock(struct DBlock& b)
{
	b.col0.fill(0);
	b.col1.fill(0);
	b.col2.fill(0);
	b.col3.fill(0);
	b.col4.fill(0);
	b.col5.fill(0);
}

class StokesSolver{
	/*
	  This class provides solver for the "Stokes" equation, which is an equation of motion of the type:

	  ResistanceMatrix*Velocity = SomeForces

	  Here, Velocity can be either translational velocity (U) alone (a 3N vector, F/U case),
	  or more generally translational and rotational velocities (U and W) (a 6N vector, FT/UW case).
	  In the latter case it is stored as (U_1, W_1, U_2, W_2, ...., U_N, W_N).
	  Accordingly, Forces are either forces alone (F) or torques and forces ( F and T, stored as (F_1, T_1,  ..., F_N, T_N) ).
	  This solver handles the case for which the ResistanceMatrix is made of two-body short range interactions, i.e. is sparse.

	  The descriptions given here are a priori written for the FT/UW case. Adaptation to the F/U case is straightforward.

	  ============================
	  ============================
	  The ResistanceMatrix has the following structure (all blocks are 6x6 if UW/FT version, 3x3 if U/F version):

	                    ================== 6N ==================
	                     .........        ..........
                    ||  |  6x6   .        . 6x6    .            | ||
                    ||  | dblock .   0    .odblock .            | ||
                    ||  |        .        .        .            | ||
                    ||  |...........................            | ||
                    ||  |        .        .                     | ||
                    ||  |     0  . dblock .                     | ||
                    6N  |        .        .                     | 6N
 ResistanceMatrix = ||  |..................                     | ||
                    ||  |        .           .                  | ||
                    ||  |odblock .             .                | ||
                    ||  |        .               .              | ||
                    ||  |.........                 .            | ||
										||  |                            .          | ||
										||  | 						       .........| ||
                    ||  |                              .        | ||
                    ||  |                              . dblock | ||
                    ||  |                              .        | ||
					                                   .........
	                    =================== 6N =================

	 This matrix is symmetric, so we only store its lower part. This implies dblocks themselves are symmetric.
	 odblocks however are antisymmetric.

	 The non-zero 6x6 blocks are of two kinds:
	 - dblocks ("diagonal") contain the contributions to the force/torque acting on particle i
	   coming from particle i's velocity/angular velocity itself: they are the Stokes drag terms ("self"-terms)
	   on the diagonal, plus parts of the interaction terms coming from lubrication or contact dashpots.
	 - odblocks ("off-diagonal") contain the contributions to the force (torque) acting on particle i
	   coming from other particles velocities and angular velocities.

	  ============================
	  ============================
	  Technical details on the storage of matrix elements:

	  Elements are stored differently depending on which library is used.

	  * CHOLMOD storage:
	  With cholmod library, the storage is first done in the "natural" block structures presented above.
	  At the end of filling, the storage is converted to the compressed-column form used by cholmod.

	  Natural form is done as follows
	  - dblocks: Due to symmetries (dblocks are symmetric, and their "B" subblocks are antisymmetric)
	             each dblock has 18 independent elements (6 diagonal terms and 12 off diagonal).
	             There are np such blocks.
	             They are all stored in 'double *dblocks', which is allocated with size 6*np.
				 Labeling in block associated with particle i is the following:
				 (subblocks names according to "Jeffrey" (Brady) ):


				 | dblocks[i].col0[0]    .                   .                 .                     .                    .               	|
"A11"(FU)| dblocks[i].col0[1]  dblocks[i].col1[0]     .                .                     .                    .               	|
				 | dblocks[i].col0[2]  dblocks[i].col1[1]  dblocks[i].col2[0]  .                     .                    .               	|
				 | 0                      .                   .                dblocks[i].col3[0]    .                    .               	|
"B11"(TU)| dblocks[i].col0[3]  0                     .                 dblocks[i].col3[1]	 dblocks[i].col4[0]     .                 |
				 | dblocks[i].col0[4]  dblocks[i].col1[2]  0                   dblocks[i].col3[2]  dblocks[i].col4[1]  	odblocks[t].col5[0] |
				 																																	"C11" subblock (TW)

	  - odblocks: 24 independent elements per block.
	              Organization is much closer to compressed-column form.
				  All the odblocks values are stored in a vector of struct ODBlock.
					A ODBlock contains 6 arrays called col0, col1, etc, each of them for a column in the block.
				  The corresponding locations in the ResistanceMatrix are stored in a vector* called odrows
                  and an array called odrows_table.

				  This works as follows:
				  Particle i interacts with particle j (i<j), we have 1 associated odblock in
				  the ResistanceMatrix. Say j appears as the m^th interaction involving i. Then:
				  odbrows[odbrows_table[i]+m] = j, and the elements are accessed via:

				  (only stored elements are shown)                                "\tilde B21"
					| odblocks[t].col0[0]    .                   .                 	0                     .                    .               	 	|
"A21"			| odblocks[t].col0[1]  odblocks[t].col1[0]     .                odblocks[t].col3[0]  0                     .                	|
					| odblocks[t].col0[2]  odblocks[t].col1[1]  odblocks[t].col2[0] odblocks[t].col3[1]  odblocks[t].col4[0]  0                 	|
					| 0                      .                   .                 	odblocks[t].col3[2]   .                    .               		|
"B21"			| odblocks[t].col0[3]  0                     .                	odblocks[t].col3[3]	 odblocks[t].col4[1]   .                 	|
					| odblocks[t].col0[4]  odblocks[t].col1[2]  0                   odblocks[t].col3[4]  odblocks[t].col4[2]  odblocks[t].col5[0] |
		                                                                          "C21"
	  with t=odbrows_table[i]+m

				Thanks to the symmetries of the resistance matrix, we can recover all the elements as:
	             (all elements)                                                   "\tilde B21"
					| odblocks[t].col0[0]  odblocks[t].col0[1]  odblocks[t].col0[2]  0                   -odblocks[t].col3[0] -odblocks[t].col3[1] |
"A21"			| odblocks[t].col0[1]  odblocks[t].col1[0]  odblocks[t].col1[1]  odblocks[t].col3[0]  0                   -odblocks[t].col4[0] |
					| odblocks[t].col0[2]  odblocks[t].col1[1]  odblocks[t].col2[0]	 odblocks[t].col3[1]  odblocks[t].col4[0]  0                   |
					| 0                   -odblocks[t].col0[3] -odblocks[t].col0[4]  odblocks[t].col3[2] 	odblocks[t].col3[3]  odblocks[t].col3[4] |
"B21"			| odblocks[t].col0[3]  0                   -odblocks[t].col1[2]  odblocks[t].col3[3]  odblocks[t].col4[1]  odblocks[t].col4[2] |
					| odblocks[t].col0[4]  odblocks[t].col1[2]  0                    odblocks[t].col3[4]  odblocks[t].col4[2]  odblocks[t].col5[0] |
		                                                                          "C21"
	  with t=odbrows_table[i]+m
	 */

private:
	int np;
	int mobile_particle_nb;

	bool mobile_matrix_done;

	// Cholmod variables
	cholmod_factor* chol_L ;
	cholmod_common chol_c ;
	cholmod_dense* chol_rhs;
	cholmod_sparse* chol_res_matrix;
	cholmod_sparse* chol_res_matrix_mf;
	cholmod_sparse* chol_res_matrix_FF;
	cholmod_dense* chol_solution;
	// cholmod_dense* chol_PTsolution;
	cholmod_dense* chol_Psolution;
	int stype;
	int sorted;
	int packed;
	int xtype;
	bool chol_L_to_be_freed;
	// resistance matrix building
	int odblocks_nb;
	int dblocks_size;
	std::vector<struct DBlock> dblocks;
	std::vector<struct ODBlock> odblocks;
	std::vector<int> odbrows;
	std::vector<int> odbrows_table;
	std::vector<struct ODBlock> odblocks_mf;
	std::vector<int> odbrows_mf;
	std::vector<int> odbrows_table_mf;
	std::vector<struct DBlock> dblocks_ff;
	std::vector<struct ODBlock> odblocks_ff;
	std::vector<int> odbrows_ff;
	std::vector<int> odbrows_table_ff;
	std::vector<std::vector <int> > odb_layout;
	std::vector<std::vector <int> > db_layout;
	std::vector<int> dblocks_cntnonzero;

	void factorizeResistanceMatrix();
    /*
     setColumn(const vec3d &nvec, int jj, double scaledXA, double scaledYB, double scaledYBtilde, double scaledYC) :
	 - appends alpha * |nvec><nvec| and corresponding indices
	 [ 3*jj, 3*jj+1, 3*jj+2 ] to column storage vectors
	 odblocks and odFrows_table
	 - this must be called with order, according to LT filling
	 */
	void setColumn(const vec3d& nvec, int jj,
			   double scaledXA, double scaledYA,
			   double scaledYB, double scaledYBtilde, double scaledYC);

	struct ODBlock buildODBlock(const vec3d &nvec,
				 							 double scaledXA,
				 							 double scaledYA, double scaledYB,
				 							 double scaledYBtilde, double scaledYC);
	void addToDBlock(struct DBlock &b, const vec3d& nvec,
												double scaledXA, double scaledYA,
			 								  double scaledYB, double scaledYC);

	void allocateResistanceMatrix();
	void allocateRessources();
	void setDiagBlockPreconditioner();
	void setIncCholPreconditioner();
	void setSpInvPreconditioner();
	void setSolverType(std::string);
	void fillCholmodFromDBlock(double *chol_x,
							   const std::vector<int>& index_chol_ix,
							   const struct DBlock& b);
	void fillCholmodFromODBlock(double *chol_x,
								const std::vector<int>& index_chol_ix,
								const struct ODBlock& b);
	void completeResistanceMatrix_MobileMobile();
	void completeResistanceMatrix_MobileFixed();
	void completeResistanceMatrix_FixedFixed();
	void fillCholmodDiag(cholmod_sparse *matrix,
									const struct DBlock &diagblock,
									int left_col_nb,
									std::vector<int> &index_chol_ix);
	void setMatrixPforBlockColumn(int *matrix_p, const std::vector<int> &pvalues);
	void setMatrixIforDiagBlock(int *matrix_i, const std::vector<int> &index_values, int top_row_nb);


public:
  ~StokesSolver();
	void init(int np);
	void printResistanceMatrix(std::ostream&, std::string);
	void printFactor(std::ostream&);
	void printRHS();
	void convertDirectToIterative();
    // R_FU filling methods
    /* resetResistanceMatrix(string solver_type, int nb_of_interactions) :
	 - initialize arrays/vectors used for building
	 - to be called before adding elements
	 - nb_of_interactions is the number of odblocks in the matrix
     */
    void resetResistanceMatrix(int nb_of_interactions,
							   const std::vector<struct DBlock>& reset_resmat_dblocks);
    /* addToDiag(int ii, double FUvalue, TWvalue) :
	 - adds FUvalue to diagonal elements to diagonal elements of FU matrix for particle ii
	 - adds TWvalue to diagonal elements to diagonal elements of TW matrix for particle ii
	 */
	//void addToDiag(int ii, double FUvalue, double TWvalue);
    /* addToDiagBlock(const vec3d &nvec, int ii, double scaledXA, double scaledYA, double scaledYB, double scaledYC);
	  Adds to block (ii, ii):
	 - scaledXA * |nvec><nvec| + scaledYA(1-|nvec><nvec|) on FU part
	 - scaledYB * e_ijk nvec_k on TU part
	 - scaledYC *(1 - |nvec><nvec|) on TW part
	 */
    void addToDiagBlock(const vec3d& nvec, int ii,
						double scaledXA, double scaledYA,
						double scaledYB, double scaledYC);
    /*
     setOffDiagBlock(const vec3d &nvec, int ii, int jj, double scaledXA, double scaledYB, double scaledYC) :
	 Sets (ii,jj) block with:
     -  scaledXA * |nvec><nvec| for FU part
	 -  scaledYB * e_ijk nvec_ij for TU part ( scaledYB is scaledYB_12(lambda) in Jeffrey & Onishi's notations)
	 -  scaledYBtilde * e_ijk nvec_ij    ( scaledYBtilde is scaledYB_12(1/lambda) in Jeffrey & Onishi's notations)
	 -  scaledYC *(1 - |nvec><nvec|) for TW part

	 This must be called with order (ii < jj),
	 because we have to fill according to the lower-triangular
	 storage.
	 */
	void setOffDiagBlock(const vec3d& nvec, int jj,
						 double scaledXA, double scaledYA, double scaledYB,
						 double scaledYBtilde, double scaledYC);
	/*
	 doneBlocks(int i) :
	 - to be called when all terms involving particle i have been added,
	 ie blocks in row i and column i are completed
	 */
	void doneBlocks(int i);
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

	/*
	   Right-hand vector access methods
	*/
    void resetRHS();
	void resetRHStorque();
	void addToRHS(double*);
    void addToRHSForce(int, double*);
    void addToRHSForce(int, const vec3d&);
    void addToRHSTorque(int, double*);
    void addToRHSTorque(int, const vec3d&);
    void setRHS(double*);
    void setRHSForce(int, const vec3d&);
    void setRHSTorque(int, const vec3d&);
    void getRHS(double*);
    /*
	 solve(vec3d* velocity, vec3d* ang_velocity) :
	 - once the resistance matrix and the RHS vector are built
	 ( completeResistanceMatrix() must have been called )
	 - solves Resistance * velocity = RHS, and stores it in velocity array
	 */
    void solve(vec3d* velocity, vec3d* ang_velocity);
    void solve(double* velocity);
    /*
	 solve_LT(double* X) :
	 - once the resistance matrix and the RHS vector are built
	 ( completeResistanceMatrix() must have been called )
	 - solves L^t * X = RHS, and stores it in X array,
	 where L^t is the transpose of the Cholesky factor ( ResistanceMatrix = L L^t )
	 - works only for direct solver, as we need the Cholesky factor
	 */
    void solve_LT(double* X);
    void solve_LT(vec3d* X, vec3d* ang_X);
    /*
	 compute_LTRHS(double* X) :
	 - once the resistance matrix and the RHS vector are built
	 ( completeResistanceMatrix() must have been called )
	 - solves X = L^T * RHS, and stores it in X array,
	 where L^t is the transpose of the Cholesky factor ( ResistanceMatrix = L L^t )
	 - works only for direct solver, as we need the Cholesky factor
	 */
    void compute_LTRHS(double* X);
    /*
	 solvingIsDone(bool free_Cholesky_factor) :
	 - deletes resistance matrix and some other arrays
	 (preconditionner, Cholesky factors, ...) used for solving
	 - should be called once every solve() call WITH THE SAME RESISTANCE MATRIX
	 is done. If matrix changes, solve() must be followed by
	 solvingIsDone()
	 */
	void solvingIsDone();
	void multiplyByResMat(double *vec);
	void multiplySolutionByResMat(double *vec);
};

#endif /* defined(__LF_DEM__StokesSolver__) */
