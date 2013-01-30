//
//  StokesSolver.h
//  LF_DEM
//
//  Created by Romain Mari on 01/29/13.
//

#ifndef __LF_DEM__StokesSolver__
#define __LF_DEM__StokesSolver__

//#define CHOLMOD
#define TRILINOS

#include "vec3d.h"

#ifdef CHOLMOD
#include "cholmod.h"
#endif

#ifdef TRILINOS
#include "Epetra_SerialComm.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"
#include "Teuchos_RCP.hpp" // ref counting pointer
#include "Ifpack2_Preconditioner.hpp" // incomplete Cholesky preconditioner
#endif

#ifdef TRILINOS
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;
using Teuchos::parameterList;

typedef double                                SCAL;
//typedef Teuchos::ScalarTraits<SCAL>          SCT;
//typedef SCT::magnitudeType                    MT;
typedef Epetra_MultiVector                     VEC;
typedef Epetra_Operator                        MAT;
//typedef Belos::MultiVecTraits<SCAL,VEC>      MVT;
//typedef Belos::OperatorTraits<SCAL,VEC,MAT>  OPT;

#endif


using namespace std;
//class System;

class StokesSolver{

 private:
    //    System *sys;

    int np, np3;
    int dof;
    int linalg_size;
#ifdef CHOLMOD
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
    
    // resistance matrix building
    vector <int> rows;
    double *diag_values;
    vector <double> *off_diag_values;
    int *ploc;
    int last_updated_colblock;

    factorizeRFU();
#endif
    
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
    RCP < Epetra_MultiVector > tril_solution;
    RCP < Epetra_MultiVector > tril_rhs;
    Epetra_CrsMatrix *tril_rfu_matrix;
    Epetra_CrsMatrix *tril_l_precond;
    //	RCP < ParameterList > params;
    //	RCP < Belos::LinearProblem < SCAL, VEC, MAT > > tril_stokes_equation;
    RCP < Belos::SolverManager < SCAL, VEC, MAT > > tril_solver;
    Belos::SolverFactory<SCAL, VEC, MAT> tril_factory;

    // resistance matrix building
    int** columns;  // diagonal block stored first, then off-diag columns, with no particular order
    int* columns_nb;
    int columns_max_nb;
    double **values;
#endif
    

#ifdef TRILINOS
    /* 
     appendToRow(const vec3d &nvec, int ii, int jj, double alpha) :
      - appends alpha * |nvec><nvec| and corresponding indices 
        to blocks [ 3*ii, 3*ii+1, 3*ii+2 ][ 3*jj, 3*jj+1, 3*jj+2 ] 
	AND symmetric [ 3*jj, 3*jj+1, 3*jj+2 ][ 3*ii, 3*ii+1, 3*ii+2 ].
    */
    void appendToRow_RFU(const vec3d &nvec, int ii, int jj, double alpha);
#endif
#ifdef CHOLMOD
    /* 
     appendToColumn(const vec3d &nvec, int jj, double alpha) :
      - appends alpha * |nvec><nvec| and corresponding indices 
        [ 3*jj, 3*jj+1, 3*jj+2 ] to column storage vectors 
	off_diag_values and ploc 
      - this must be called with order, according to LT filling
    */
    void appendToColumn_RFU(const vec3d &nvec, int ii, int jj, double alpha); 
#endif

    void allocate_RFU();
    
    void allocateRessources();

    void buildDiagBlockPreconditioner();
    void buildIncCholPreconditioner();
    
 public:

    StokesSolver(int);
    ~StokesSolver();
    void initialize();



    // R_FU filling methods
    /* prepareNewBuild_RFU() :
        - initialize arrays/vectors used for building
        - to be called before adding elements
     */
    void prepareNewBuild_RFU();

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
    complete_RFU() :
      - transforms temporary arrays/vectors used to build resistance
        matrix into Cholmod or Epetra objects really used by 
	the solvers
      - also builds preconditionner when needed (TRILINOS)
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
    solvingIsDone() :
      - deletes RFU matrix and some other arrays 
        (preconditionner, Cholesky factors, ...) used for solving
      - should be called once every call to solve() WITH THE SAME RFU
        is done. If RFU changes, solve() must be followed by 
        solvingIsDone()
    */
    void solvingIsDone();

#ifdef CHOLMOD
    cholmod_factor *chol_L ;
    cholmod_common chol_c ;
#endif

};
#endif /* defined(__LF_DEM__StokesSolver__) */
