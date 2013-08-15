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

	  Here, Velocity can be either translational velocity (U) alone (a 3N vector), or more generally translational and rotational velocities (U and W) (a 6N vector).
	  Accordingly, Forces are either forces alone (F) or torques and forces (F and T).
	  This solver handles the case for which the ResistanceMatrix is made of two-body short range interactions, ie is sparse.

	  The ResistanceMatrix has the following structure (all blocks are 3x3):
	  (if only force-velocity resistance matrix is considered, then ResistanceMatrix is limited tothe top left corner)

	                    =========== 3Nx3N FU submatrix =========||========== 3Nx3N FW submatrix =========
	                     .........        ..........            :.........        ..........
                    ||  |  3x3   .        .        .            :        .        .        .           | ||
                    ||  | dblock .   0    .odblock .            :sdblock .  0     .odblock .           | ||
                    ||  |  (FU)  .        .        .            : (FW)   .        .        .           | ||
                    ||  |...........................            :...........................           | ||
                    ||  |        .        .                     :       .        .                     | ||
                    ||  |     0  . dblock .                     :  0    .sdblock .                     | ||
                    FU  |        .        .                     :       .        .                     | 3Nx3N FW submatrix
                    ||  |        ..........                     :       ..........                     | ||
                    ||  |                    .                  :                   .                  | || 
                    ||  |                      .                :                     .                | ||     
                    ||  |                        .              :                       .              | ||    
                    ||  |                          .            :                         .            | ||
					||  |                            .          :                           .          | ||
					||  | 						       .........:                             .........| ||
                    ||  |                              .        :                             .        | ||
                    ||  |                              . dblock :                             .sdblock | ||
                    ||  |                              .        :                             .        | ||
 ResistanceMatrix =     |---------------------------------------:--------------------------------------| =
					||  |        .		  .        .            :        .        .        .           | ||
					||  |sdblock .  0     .odblock .            : dblock .  0     .odblock .           | ||
					||  | (TU)   .		  .        .	        :  (WT)  .        .        .           | ||
					||  |...........................  		    :...........................           | ||
                    ||  |        .        .                     :        .        .                    | ||
                    ||  |    0   .sdblock .                     : 0      . dblock .                    | ||
                    ||  |        .        .                     :        .        .                    | ||
                    ||  |        ..........                     :        ..........                    | ||
                    TU  |                    .                  :                   .                  | 3Nx3N TW submatrix
                    ||  |                      .                :                     .                | ||                     .
                    ||  |                        .              :                       .              | ||    
                    ||  |                          .            :                         .            | ||
					||  |                            .          :                           .          | ||
					||  | 				   		       .........:                             .........| ||
                    ||  |                              .        :                             .        | ||
                    ||  |                              .sdblock :                             . dblock | ||
                    ||  |                              .        :                             .        | ||
					||  |                              .........:                             .........| || 
	                    ========== 3Nx3N TU submatrix ==========||========== 3Nx3N TW submatrix =========					  

	 This matrix is symmetric, so we only store its upper or lower part.

	 The non-zero 3x3 blocks are of three kinds:
	 - dblocks ("diagonal") contain the contributions to the force (resp torque) acting on particle i 
	   coming from particle i's velocity (angular velocity) itself: they are the Stokes drag terms ("self"-terms) 
	   on the diagonal, plus parts of the interaction terms coming from lubrication or contact dashpots. 
	 - sdblocks ("secondary-diagonal") contain the contributions to the force (torque) acting on particle i 
	   coming from particle i's angular velocity (velocity) itself: they are only interaction terms.
	 - odblocks ("off-diagonal") contain the contributions to the force (torque) acting on particle i 
	   coming from other particles velocities and angular velocities.




	 */

private:
	int np;
	int np3;
	//    int dof;
    int res_matrix_linear_size;
	int dblocks_nb;
	int sdblocks_nb;
	int dblocks_element_nb;
	bool brownian;
	
	bool _iterative;
	bool _direct;
	
	bool FTcoupling;

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
    vector <int> rows;
    double *dblocks;
    vector <double> *odblocks;
    int *ploc;
	
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
     appendToRow(const vec3d &nvec, int ii, int jj, double alpha) :
	 - TRILINOS ONLY
	 - appends alpha * |nvec><nvec| and corresponding indices
	 to blocks [ 3*ii, 3*ii+1, 3*ii+2 ][ 3*jj, 3*jj+1, 3*jj+2 ]
	 AND symmetric [ 3*jj, 3*jj+1, 3*jj+2 ][ 3*ii, 3*ii+1, 3*ii+2 ].
	 
	 */
    void appendToRow(const vec3d &nvec, int ii, int jj, double alpha);
	
    /*
     appendToColumn(const vec3d &nvec, int jj, double alpha) :
	 - CHOLMOD ONLY
	 - appends alpha * |nvec><nvec| and corresponding indices
	 [ 3*jj, 3*jj+1, 3*jj+2 ] to column storage vectors
	 odblocks and ploc
	 - this must be called with order, according to LT filling
	 */
    void appendToColumn(const vec3d &nvec, int ii, int jj, double alpha);
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
	void printResistanceMatrix();
	void convertDirectToIterative();
    // R_FU filling methods
    /* resetResistanceMatrix() :
	 - initialize arrays/vectors used for building
	 - to be called before adding elements
	 - possible arguments: "direct" or "iterative"
     */
    void resetResistanceMatrix(string solver_type);
	
    /* addToDiag(int ii, double alpha) :
	 - adds alpha to diagonal elements 3*ii, 3*ii+1, and 3*ii+2
	 */
    void addToDiag(int ii, double alpha);
	
    /* addToDiagBlock(const vec3d &nvec, int ii, double alpha);
	 - adds alpha * |nvec><nvec| to diagonal block
	 [ 3*ii, 3*ii+1, 3*ii+2 ][ 3*ii, 3*ii+1, 3*ii+2 ]
	 */
    void addToDiagBlock(const vec3d &nvec, int ii, double alpha);
	
    /*
     appendToOffDiagBlock(const vec3d &nvec, int ii, int jj, double alpha) :
     - appends alpha * |nvec><nvec| and corresponding indices
	 to block [ 3*jj, 3*jj+1, 3*jj+2 ][ 3*ii, 3*ii+1, 3*ii+2 ]
	 
     - if the solver is Trilinos, it also appends it to the
	 symmetric block [ 3*ii, 3*ii+1, 3*ii+2 ][ 3*jj, 3*jj+1, 3*jj+2 ]
	 
	 - this must be called with order (ii < jj) if CHOLMOD is used,
	 because we have to fill according to the lower-triangular
	 storage.
	 */
    void appendToOffDiagBlock(const vec3d &nvec, int ii, int jj, double alpha);
	
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
    void addToRHS(int, const double);
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
