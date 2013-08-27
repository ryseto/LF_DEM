//
//  StokesSolver.h
//  LF_DEM
//
//  Created by Romain Mari on 01/29/13.
//

#ifndef __LF_DEM__StokesSolver__
#define __LF_DEM__StokesSolver__
//#define CHOLMOD_EXTRA
//#define TRILINOS
#include <vector>
#include "vec3d.h"
#include "cholmod.h"
#ifdef CHOLMOD_EXTRA
extern "C" {
#include "cholmod_extra.h"
}
#endif
#ifdef TRILINOS
#include "Epetra_SerialComm.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"
#include "Teuchos_RCP.hpp" // ref counting pointer
//#include "Ifpack2_Preconditioner.hpp" // incomplete Cholesky preconditioner
//#include "Ifpack_CrsIct.h"
#include "Ifpack_IC.h"
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;
using Teuchos::parameterList;

typedef double                                SCAL;
typedef Epetra_MultiVector                     VEC;
typedef Epetra_Operator                        MAT;
#endif
using namespace std;
//class System;

class StokesSolver{

	/*
	 
	  This class provides solver for the "Stokes" equation, which is an equation of motion of the type:
	  
	  ResistanceMatrix*Velocity = SomeForces

	  Here, Velocity can be either translational velocity (U) alone (a 3N vector, F/U case), 
	  or more generally translational and rotational velocities (U and W) (a 6N vector, FT/UW case). 
	  In the latter case it is stored as (U_1, W_1, U_2, W_2, ...., U_N, W_N).
	  Accordingly, Forces are either forces alone (F) or torques and forces ( F and T, stored as (F_1, T_1,  ..., F_N, T_N) ).
	  This solver handles the case for which the ResistanceMatrix is made of two-body short range interactions, ie is sparse.
	  
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

			     	 | 18*i      .        .        .         .         .     |
"A11" subblock (FU)	 | 18*i+1   18*i+6    .        .         .         .     |
                     | 18*i+2   18*i+7   18*i+10   .         .         .     |
					 | 18*i+3    .        .       18*i+12    .         .     |  
"B11" subblock (TU)	 | 18*i+4   18*i+8    .       18*i+13   18*i+15    .     |
				     | 18*i+5   18*i+9   18*i+11  18*i+14   18*i+16  18*i+17 |
                                                     "C11" subblock (TW)

	  - odblocks: 24 independent elements per block.
	              Organization is much closer to compressed-column form.
				  All the odblocks values are stored columnwise in 6 vectors called odblocks[0-5]
				  The corresponding locations in the ResistanceMatrix are stored in a vector* called odrows 
                  and an array called odrows_table.

				  This works as follows:
				  Particle i interacts with particle j (i<j), we have 1 associated odblock in
				  the ResistanceMatrix. Say j appears as the m^th interaction involving i. Then: 
				  odbrows[odbrows_table[i]+m] = j, and the elements are accessed via:

				  (only stored elements are shown)                                "\tilde B21"
		| odblocks[0][6*t  ]    .                   .                 odblocks[3][6*t  ]    .                   .                |
"A21"	| odblocks[0][6*t+1]  odblocks[1][4*t  ]    .                 odblocks[3][6*t+1]  odblocks[4][4*t  ]    .                |
		| odblocks[0][6*t+2]  odblocks[1][4*t+1]  odblocks[2][2*t  ]  odblocks[3][6*t+2]  odblocks[4][4*t+1]  odblocks[5][2*t  ] |
		| odblocks[0][6*t+3]    .                   .                 odblocks[3][6*t+3]    .                   .                |
"B21"	| odblocks[0][6*t+4]  odblocks[1][4*t+2]    .                 odblocks[3][6*t+4]  odblocks[4][4*t+2]    .                | 
		| odblocks[0][6*t+5]  odblocks[1][4*t+3]  odblocks[2][2*t+1]  odblocks[3][6*t+5]  odblocks[4][4*t+3]  odblocks[5][2*t+1] |
		                                                                          "C21"
	  with t=odbrows_table[i]+m
	  
	             (all elements)                                                   "\tilde B21"
		| odblocks[0][6*t  ]  odblocks[0][6*t+1]  odblocks[0][6*t+2]  odblocks[3][6*t  ] -odblocks[3][6*t+1] -odblocks[3][6*t+2] |
"A21"	| odblocks[0][6*t+1]  odblocks[1][4*t  ]  odblocks[1][4*t+1]  odblocks[3][6*t+1]  odblocks[4][4*t  ] -odblocks[4][4*t+1] |
		| odblocks[0][6*t+2]  odblocks[1][4*t+1]  odblocks[2][2*t  ]  odblocks[3][6*t+2]  odblocks[4][4*t+1]  odblocks[5][2*t  ] |
		| odblocks[0][6*t+3] -odblocks[0][6*t+4] -odblocks[0][6*t+5]  odblocks[3][6*t+3]  odblocks[3][6*t+4]  odblocks[3][6*t+5] |
"B21"	| odblocks[0][6*t+4]  odblocks[1][4*t+2] -odblocks[1][4*t+3]  odblocks[3][6*t+4]  odblocks[4][4*t+2]  odblocks[4][4*t+3] | 
		| odblocks[0][6*t+5]  odblocks[1][4*t+3]  odblocks[2][2*t+1]  odblocks[3][6*t+5]  odblocks[4][4*t+3]  odblocks[5][2*t+1] |
		                                                                          "C21"
	  with t=odbrows_table[i]+m


	  * TRILINOS storage:
	  TO BE DONE. TRILINOS NOT YET ADAPTED FOR FT/UW VERSION.
				  
	 */

private:
	int np;
	//	int np3;
	int np6;
	
    int res_matrix_linear_size;
	int odblocks_nb;
	int dblocks_size;
	bool brownian;
	
	bool _iterative;
	bool _direct;
	
	//	bool FTcoupling;

	// Cholmod variables
    cholmod_factor *chol_L ;
    cholmod_common chol_c ;
    cholmod_dense *chol_rhs;
    cholmod_sparse *chol_res_matrix;
	
    cholmod_dense *chol_solution;
    cholmod_dense *chol_PTsolution;
    cholmod_dense *chol_v_nonBrownian;
    cholmod_dense *chol_v_Brownian_init;
    cholmod_dense *chol_v_Brownian_mid;
    cholmod_dense *chol_brownian_rhs;
	bool chol_init;
	
    int stype;
    int sorted;
    int packed;
    int xtype;
	bool chol_L_to_be_freed;
	
    // resistance matrix building
    double *dblocks;
    vector <double> *odblocks;
    vector <int> odbrows;
	int *odbrows_table;
	int *current_index_positions;

    void factorizeResistanceMatrix();
	
    
#ifdef TRILINOS
    int MyPID;
    /* #ifdef EPETRA_MPI */
    /* 	// Initialize MPI */
    /* 	MPI_Init(&argc,&argv); */
    /* 	Epetra_MpiComm Comm(MPI_COMM_WORLD); */
    /* 	MyPID = Comm.MyPID(); */
    /* #else */
    Epetra_SerialComm Comm;
    //#endif
    RCP < Epetra_Map > Map;
    RCP < Epetra_Vector > tril_solution;
    RCP < Epetra_Vector > tril_rhs;
    Epetra_CrsMatrix *tril_res_matrix;
    Epetra_CrsMatrix *tril_l_precond;
	
	RCP < Belos::LinearProblem < SCAL, VEC, MAT > > tril_stokes_equation;
    RCP < Belos::SolverManager < SCAL, VEC, MAT > > tril_solver;
    Belos::SolverFactory<SCAL, VEC, MAT> tril_factory;
	
#endif
    // resistance matrix building for Trilinos
    int** columns;  // diagonal block stored first, then off-diag columns, with no particular order
    int* columns_nb; // nb of non-zero elements in each row
    int columns_max_nb; // max nb of non-zero elements per row
    double **values; // values corresponding to 'columns' array coordinates
    /*
     setRow(const vec3d &nvec, int ii, int jj, double alpha) :
	 - TRILINOS ONLY
	 - appends alpha * |nvec><nvec| and corresponding indices
	 to blocks [ 3*ii, 3*ii+1, 3*ii+2 ][ 3*jj, 3*jj+1, 3*jj+2 ]
	 AND symmetric [ 3*jj, 3*jj+1, 3*jj+2 ][ 3*ii, 3*ii+1, 3*ii+2 ].
	 
	 */
    void setRow(const vec3d &nvec, int ii, int jj, double scaledXA, double scaledYBtilde, double scaledYB, double scaledYC);
	
    /*
     setColumn(const vec3d &nvec, int jj, double scaledXA, double scaledYB, double scaledYBtilde, double scaledYC) :
	 - CHOLMOD ONLY
	 - appends alpha * |nvec><nvec| and corresponding indices
	 [ 3*jj, 3*jj+1, 3*jj+2 ] to column storage vectors
	 odblocks and odFrows_table
	 - this must be called with order, according to LT filling
	 */
    void setColumn(const vec3d &nvec, int ii, int jj, double scaledXA, double scaledYB, double scaledYBtilde, double scaledYC);
    void allocateResistanceMatrix();
    void completeResistanceMatrix_cholmod();
    void completeResistanceMatrix_trilinos();
    void allocateRessources();
    void setDiagBlockPreconditioner();
    void setIncCholPreconditioner();
    void setSpInvPreconditioner();
	void setSolverType(string);
	
public:
    ~StokesSolver();
	void init(int np, bool is_brownian);
    void initialize();
	bool direct(){
		return _direct;
	}
	bool iterative(){
		return _iterative;
	}
	void printResistanceMatrix(ofstream &, string);
	void printRHS();
	void convertDirectToIterative();
    // R_FU filling methods
    /* resetResistanceMatrix(string solver_type, int nb_of_interactions) :
	 - initialize arrays/vectors used for building
	 - to be called before adding elements
	 - possible solver_type: "direct" or "iterative"
	 - nb_of_interactions is the number of odblocks in the matrix
     */
    void resetResistanceMatrix(string solver_type, int nb_of_interactions);
	
    /* addToDiag(int ii, double FUvalue, TWvalue) :
	 - adds FUvalue to diagonal elements to diagonal elements of FU matrix for particle ii 
	 - adds TWvalue to diagonal elements to diagonal elements of TW matrix for particle ii 
	 */
    void addToDiag(int ii, double FUvalue, double TWvalue);
	
    /* addToDiagBlock(const vec3d &nvec, int ii, double scaledXA, double scaledYB, double scaledYC);
	  Adds to block (ii, ii):
	 - scaledXA * |nvec><nvec| on FU part
	 - scaledYB * e_ijk nvec_ij on TU part (NOT IMPLEMENTED YET)
	 - scaledYC *(1 - |nvec><nvec|) on TW part (NOT IMPLEMENTED YET)
	 */
    void addToDiagBlock(const vec3d &nvec, int ii, double scaledXA, double scaledYB, double scaledYC);
	
    /*
     setOffDiagBlock(const vec3d &nvec, int ii, int jj, double scaledXA, double scaledYB, double scaledYC) :
	 Sets (ii,jj) block with:
     -  scaledXA * |nvec><nvec| for FU part
	 -  scaledYB * e_ijk nvec_ij for TU part ( scaledYB is scaledYB_12(lambda) in Jeffrey & Onishi's notations)
	 -  scaledYBtilde * e_ijk nvec_ij    ( scaledYBtilde is scaledYB_12(1/lambda) in Jeffrey & Onishi's notations)
	 -  scaledYC *(1 - |nvec><nvec|) for TW part

     If the solver is Trilinos, it also sets the symmetric block (jj, ii)

	 This must be called with order (ii < jj) if CHOLMOD is used,
	 because we have to fill according to the lower-triangular
	 storage.
	 */
    void setOffDiagBlock(const vec3d &nvec, int ii, int jj, double scaledXA, double scaledYB, double scaledYBtilde, double scaledYC);
	
	/*
	 doneBlocks(int i) :
	 - to be called when all terms involving particle i have been added,
	 ie blocks in row i and column i are completed
	 */
	void doneBlocks(int i);
	
	/*
	 completeResistanceMatrix() :
	 - transforms temporary arrays/vectors used to build resistance
	 matrix into Cholmod or Epetra objects really used by
	 the solvers
	 - also builds preconditionner when needed (if TRILINOS)
	 - also performs Cholesky factorization (if CHOLMOD)
	 - must be called before solving the linear algebra problem
	 - must be called after all terms are added
	 */
    void completeResistanceMatrix();
    
	
    void resetRHS();
    void addToRHSForce(int, double *);
    void addToRHSForce(int, const vec3d &);
    void addToRHSTorque(int, double *);
    void addToRHSTorque(int, const vec3d &);
    void addToRHS(double*);
    void setRHS(double*);
    void getRHS(double*);
	
    /*
	 solve(double* velocity) :
	 - once the resistance matrix and the RHS vector are built
	 ( completeResistanceMatrix() must have been called )
	 - solves Resistance * velocity = RHS, and stores it in velocity array
	 */
    void solve(double* velocity);
	
    /*
	 solve_CholTrans(double* velocity) :
	 - once the resistance matrix and the RHS vector are built
	 ( completeResistanceMatrix() must have been called )
	 - solves L^t * velocity = RHS, and stores it in velocity array,
	 where L^t is the transpose of the Cholesky factor ( ResistanceMatrix = L L^t )
	 - works only for direct solver, as we need the Cholesky factor
	 */
    void solve_CholTrans(double* velocity);
	
    /*
	 solvingIsDone(bool free_Cholesky_factor) :
	 - deletes resistance matrix and some other arrays
	 (preconditionner, Cholesky factors, ...) used for solving
	 - should be called once every solve() call WITH THE SAME RESISTANCE MATRIX
	 is done. If matrix changes, solve() must be followed by
	 solvingIsDone()
	 */
    void solvingIsDone();
	
	
	
};
#endif /* defined(__LF_DEM__StokesSolver__) */
