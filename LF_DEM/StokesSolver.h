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

using namespace std;

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



//class System;

class StokesSolver{

 private:
 
	int np, np3;
    int dof;
    int linalg_size;
	bool brownian;

	bool _iterative;
	bool _direct;

	// Cholmod variables
    cholmod_sparse *chol_rfu_matrix;
    cholmod_dense *chol_solution;
    cholmod_dense *chol_v_nonBrownian;
    cholmod_dense *chol_v_Brownian_init;
    cholmod_dense *chol_v_Brownian_mid;
    cholmod_dense *chol_brownian_rhs;
    cholmod_dense *chol_rhs;
    int stype;
    int sorted;
    int packed;
    int xtype;
	bool chol_L_to_be_freed;
	
    // resistance matrix building
    vector <int> rows;
    double *diag_values;
    vector <double> *off_diag_values;
    int *ploc;

    void factorizeRFU();

    
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
    Epetra_CrsMatrix *tril_rfu_matrix;
    Epetra_CrsMatrix *tril_l_precond;

	RCP < Belos::LinearProblem < SCAL, VEC, MAT > > tril_stokes_equation;
    RCP < Belos::SolverManager < SCAL, VEC, MAT > > tril_solver;
    Belos::SolverFactory<SCAL, VEC, MAT> tril_factory;

    // resistance matrix building
    int** columns;  // diagonal block stored first, then off-diag columns, with no particular order
    int* columns_nb; // nb of non-zero elements in each row
    int columns_max_nb; // max nb of non-zero elements per row
    double **values; // values corresponding to 'columns' array coordinates


#endif
    

    /* 
     appendToRow(const vec3d &nvec, int ii, int jj, double alpha) :
	  - TRILINOS ONLY
      - appends alpha * |nvec><nvec| and corresponding indices 
        to blocks [ 3*ii, 3*ii+1, 3*ii+2 ][ 3*jj, 3*jj+1, 3*jj+2 ] 
	AND symmetric [ 3*jj, 3*jj+1, 3*jj+2 ][ 3*ii, 3*ii+1, 3*ii+2 ].
	  
    */
    void appendToRow_RFU(const vec3d &nvec, int ii, int jj, double alpha);

    /* 
     appendToColumn(const vec3d &nvec, int jj, double alpha) :
	  - CHOLMOD ONLY
      - appends alpha * |nvec><nvec| and corresponding indices 
        [ 3*jj, 3*jj+1, 3*jj+2 ] to column storage vectors 
	off_diag_values and ploc 
      - this must be called with order, according to LT filling
    */
    void appendToColumn_RFU(const vec3d &nvec, int ii, int jj, double alpha); 

    void allocate_RFU();
    void complete_RFU_cholmod();
    void complete_RFU_trilinos();
    
    void allocateRessources();

    void setDiagBlockPreconditioner();
    void setIncCholPreconditioner();
    void setSpInvPreconditioner();

	void setSolverType(string);

	void print_RFU();
 public:

    StokesSolver(int np, bool is_brownian);
    ~StokesSolver();
    void initialize();

	bool direct(){
		return _direct;
	}
	bool iterative(){
		return _iterative;
	}

	void convertDirectToIterative();
    // R_FU filling methods
    /* prepareNewBuild_RFU() :
        - initialize arrays/vectors used for building
        - to be called before adding elements
		- possible arguments: "direct" or "iterative"
     */
    void prepareNewBuild_RFU(string solver_type);

    /* addToDiag_RFU(int ii, double alpha) :
       - adds alpha to diagonal elements 3*ii, 3*ii+1, and 3*ii+2
    */
    void addToDiag_RFU(int ii, double alpha);

    /* addToDiagBlock_RFU(const vec3d &nvec, int ii, double alpha);
       - adds alpha * |nvec><nvec| to diagonal block 
         [ 3*ii, 3*ii+1, 3*ii+2 ][ 3*ii, 3*ii+1, 3*ii+2 ] 
    */
    void addToDiagBlock_RFU(const vec3d &nvec, int ii, double alpha);

    /* 
     appendToOffDiagBlock_RFU(const vec3d &nvec, int ii, int jj, double alpha) :
     - appends alpha * |nvec><nvec| and corresponding indices 
       to block [ 3*jj, 3*jj+1, 3*jj+2 ][ 3*ii, 3*ii+1, 3*ii+2 ]

     - if the solver is Trilinos, it also appends it to the 
       symmetric block [ 3*ii, 3*ii+1, 3*ii+2 ][ 3*jj, 3*jj+1, 3*jj+2 ]
       
      - this must be called with order (ii < jj) if CHOLMOD is used, 
        because we have to fill according to the lower-triangular
	storage.
    */
    void appendToOffDiagBlock_RFU(const vec3d &nvec, int ii, int jj, double alpha); 

	/* 
	 doneBlocks(int i) :
	  - to be called when all terms involving particle i have been added,
	    ie blocks in row i and column i are completed
	*/
	void doneBlocks(int i);

	/* 
	 complete_RFU() :
      - transforms temporary arrays/vectors used to build resistance
        matrix into Cholmod or Epetra objects really used by 
	the solvers
      - also builds preconditionner when needed (if TRILINOS)
      - also performs Cholesky factorization (if CHOLMOD)
      - must be called before solving the linear algebra problem
      - must be called after all terms are added
    */
    void complete_RFU();
    

    void resetRHS();
    void addToRHS(int, const double);
    void setRHS(double*);

    /* 
    solve(double* velocity) :
      - once the RFU matrix and the RHS vector are built 
        ( complete_RFU() must have been called )
      - solves RFU * velocity = RHS, and stores it in velocity array
    */
    void solve(double* velocity);

    /* 
    solve_CholTrans(double* velocity) :
      - once the RFU matrix and the RHS vector are built 
        ( complete_RFU() must have been called )
      - solves L^t * velocity = RHS, and stores it in velocity array,
	    where L^t is the transpose of the Cholesky factor ( RFU = L L^t )
	  - works only for direct solver, as we need the Cholesky factor
    */
    void solve_CholTrans(double* velocity);

    /* 
    solvingIsDone(bool free_Cholesky_factor) :
      - deletes RFU matrix and some other arrays 
        (preconditionner, Cholesky factors, ...) used for solving
      - should be called once every call to solve() WITH THE SAME RFU
        is done. If RFU changes, solve() must be followed by 
        solvingIsDone()
    */
    void solvingIsDone();



    cholmod_factor *chol_L ;
    cholmod_common chol_c ;

};
#endif /* defined(__LF_DEM__StokesSolver__) */
