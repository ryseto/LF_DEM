#include "fstream"
#include "StokesSolver.h"
#ifdef TRILINOS
#include <BelosCGIteration.hpp>
#endif
using namespace std;
#define DELETE(x) if(x){delete [] x; x = NULL;}

/******************************************************
 *                                                     *
 *                   Public Methods                    *
 *                                                     *
 ******************************************************/

StokesSolver::~StokesSolver(){
    if (!dblocks) {
		delete [] dblocks;
	}
	if (!odblocks) {
		delete [] odblocks;
	}
	if (!odbrows_table) {
		delete [] odbrows_table;
	}
    if (!current_index_positions) {
		delete [] current_index_positions;
	}
	if(!chol_solution) {
		cholmod_free_dense(&chol_solution, &chol_c);
	}
	if(!chol_rhs) {
		cholmod_free_dense(&chol_rhs, &chol_c);
	}
	if (brownian) {
		cholmod_free_dense(&chol_brownian_rhs, &chol_c);
	}
	if(!chol_res_matrix) {
		cholmod_free_sparse(&chol_res_matrix, &chol_c);
	}
	if(chol_init) {
		cholmod_finish(&chol_c);
	}
	
#ifdef TRILINOS
	for (int i=0; i<res_matrix_linear_size; i++) {
		delete [] columns[i];
	}
	delete [] columns;
	delete [] columns_nb;
	for (int i=0; i<res_matrix_linear_size; i++){
		delete [] values[i];
	}
	delete [] values;
#endif
}

void
StokesSolver::init(int n, bool is_brownian){
	np = n;
    np6 = 6*np;
	brownian = is_brownian;
	// initializing values that can be changed later
	_direct = true;
	_iterative = false;
	chol_init = false;
	//	FTcoupling = false;
}

void
StokesSolver::initialize(){
	// CHOLMOD parameters
	stype = -1; // 1 is symmetric, stored upper triangular (UT), -1 is LT
	sorted = 0;		/* TRUE if columns sorted, FALSE otherwise*/
	packed = 1;		/* TRUE if matrix packed, FALSE otherwise */
	xtype = CHOLMOD_REAL;
#ifdef TRILINOS
	// TRILINOS init and parameters
	// initialize solver
	RCP<ParameterList> solverParams = parameterList();
	// parameters to be tuned (and understood!)
	int blocksize = 10;
	int maxiters = 400;
	double tol = 1.e-6;
	solverParams->set("Block Size", blocksize);              // Blocksize to be used by iterative solver
	solverParams->set("Maximum Iterations", maxiters);       // Maximum number of iterations allowed
	solverParams->set("Convergence Tolerance", tol);         // Relative convergence tolerance requested
	solverParams->set("Verbosity", Belos::Errors + Belos::Warnings);
	tril_solver = tril_factory.create ("CG", solverParams);
	// initialize empty linear problem
    tril_stokes_equation = rcp(new Belos::LinearProblem <SCAL, VEC, MAT> ());
#endif
	// resistance matrix characteristics (see header for matrix description)
	res_matrix_linear_size = np6;
	dblocks_size = 18*np;
	allocateRessources();
	chol_L_to_be_freed = false;
}

/************* Matrix filling methods **********************/

// Diagonal Terms, FT/UW version
void
StokesSolver::addToDiag(int ii, double FUvalue, double TWvalue){
	if (direct()) {
		int ii18 = 18*ii;
		dblocks[ii18   ] += FUvalue;
		dblocks[ii18+6 ] += FUvalue;
		dblocks[ii18+10] += FUvalue;
		dblocks[ii18+12] += TWvalue;
		dblocks[ii18+15] += TWvalue;
		dblocks[ii18+17] += TWvalue;
	}
#ifdef TRILINOS
	if (iterative()) {
		cerr << " Error : StokesSolver::addToDiag(const vec3d &nvec, int ii, double FUvalue, double TWvalue) not implemented for TRILINOS yet ! " << endl;
		exit(1);
		int ii3 = 3*ii;
		for (int j = 0; j < 3; j ++) {
			values[ii3+j][j] += FUvalue;
		}
	}
#endif
}

// Diagonal Blocks Terms, FT/UW version
void

StokesSolver::addToDiagBlock(const vec3d &nvec, int ii, double scaledXA, double scaledYA, double scaledYB, double scaledYC){
	double n0n0 = nvec.x*nvec.x;
	double n0n1 = nvec.x*nvec.y;
	double n0n2 = nvec.x*nvec.z;
	double n1n1 = nvec.y*nvec.y;
	double n1n2 = nvec.y*nvec.z;
	double n2n2 = nvec.z*nvec.z;
	double one_n0n0 = 1-n0n0;
	double one_n1n1 = 1-n1n1;
	double one_n2n2 = 1-n2n2;
	if (direct()) {
		int ii18 = 18*ii;
		// (*,0)
		dblocks[ii18   ] += scaledXA*n0n0 + scaledYA*one_n0n0;        // 00 element of the dblock
		dblocks[ii18+1 ] += (scaledXA-scaledYA)*n0n1;        // 10
		dblocks[ii18+2 ] += (scaledXA-scaledYA)*n0n2;        // 20
		dblocks[ii18+3 ] += 0;                    // 30 
		dblocks[ii18+4 ] += -scaledYB*nvec.z;     // 40
		dblocks[ii18+5 ] += +scaledYB*nvec.y;      // 50
		// (*,1)
		dblocks[ii18+6 ] += scaledXA*n1n1 + scaledYA*one_n1n1;        // 11
		dblocks[ii18+7 ] += (scaledXA-scaledYA)*n1n2;        // 21
		dblocks[ii18+8 ] += 0;                    // 41
		dblocks[ii18+9 ] += -scaledYB*nvec.x;     // 51
		// (*,2)
		dblocks[ii18+10] += scaledXA*n2n2 + scaledYA*one_n2n2;        // 22
		dblocks[ii18+11] += 0;                    // 32
		// (*,3)
		dblocks[ii18+12] += scaledYC*one_n0n0;    // 33
		dblocks[ii18+13] += -scaledYC*n0n1;       // 43
		dblocks[ii18+14] += -scaledYC*n0n2;       // 53
		// (*,4)
		dblocks[ii18+15] += scaledYC*one_n1n1;    // 44
		dblocks[ii18+16] += -scaledYC*n1n2;       // 54
		// (*,5)
		dblocks[ii18+17] += scaledYC*one_n2n2;    // 55
	}
#ifdef TRILINOS
	if (iterative()) {
		int iidof = dof*ii;
		cerr << " Error : StokesSolver::addToDiagBlock(const vec3d &nvec, int ii, double scaledXA, double scaledYB, double scaledYC) not implemented for TRILINOS yet ! " << endl;
		exit(1);
		values[iidof  ][0] += scaledXA*n0n0; // 00
		values[iidof  ][1] += scaledXA*n0n1; // 01
		values[iidof  ][2] += scaledXA*n0n2; // 02
		values[iidof+1][0] += scaledXA*n0n1; // 10
		values[iidof+1][1] += scaledXA*n1n1; // 11
		values[iidof+1][2] += scaledXA*n1n2; // 12
		values[iidof+2][0] += scaledXA*n0n2; // 20
		values[iidof+2][1] += scaledXA*n1n2; // 21
		values[iidof+2][2] += scaledXA*n2n2; // 22
	}
#endif
}


// Off-Diagonal Blocks Terms, FT/UW version
void
StokesSolver::setOffDiagBlock(const vec3d &nvec, int ii, int jj,
							  double scaledXA, double scaledYA, double scaledYB, double scaledYBtilde, double scaledYC){
	if (direct()) {
		setColumn(nvec, jj, scaledXA, scaledYA, scaledYB, scaledYBtilde, scaledYC);
	}
#ifdef TRILINOS
	if (iterative()) {
		setRow(nvec, ii, jj, scaledXA, scaledYA, scaledYB, scaledYBtilde, scaledYC);
	}
#endif
	return;
	ii = 0; // prevents gcc warning when compiled without TRILINOS
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
void
StokesSolver::completeResistanceMatrix_cholmod(){

	// this function is commented, but you are strongly advised to read 
	// the description of storage in the header file first :)

    // declare the last 2 values of odbrows_table
    odbrows_table[np-1] = (unsigned int)odbrows.size();
    odbrows_table[np] = (unsigned int)odbrows.size();
    allocateResistanceMatrix();
    // fill
	for (int j=0; j<np; j++) {
		// associated with particle j are 6 columns in the matrix:
		// { 6j, ... , 6j+5 }
		int j6 = 6*j;
		int j18 = 18*j;
		int j21 = 21*j;
		// the number of non-zero elements before column 6j is:
		// - 21*j from j diagonal blocks
		// - 36*odbrows_table[j] from odbrows_table[j] off-diagonal blocks
		//
		// the number of non-zero elements before column 6j+1 is:
		// - number of non-zero before column 6j + number of non-zero in column 6*j
		// (in 6j: 6 elements in diagonal block, plus 6*(odbrows_table[j+1]-odbrows_table[j])
		//
		// for 6j+2 --> 6j+5: same idea
		int od_nzero_nb = 6*(odbrows_table[j+1]-odbrows_table[j]);
		((int*)chol_res_matrix->p)[j6  ] = j21   + 36*odbrows_table[j];
		((int*)chol_res_matrix->p)[j6+1] = ((int*)chol_res_matrix->p)[j6] + 6 + od_nzero_nb;
		((int*)chol_res_matrix->p)[j6+2] = ((int*)chol_res_matrix->p)[j6+1] + 5 + od_nzero_nb;
		((int*)chol_res_matrix->p)[j6+3] = ((int*)chol_res_matrix->p)[j6+2] + 4 + od_nzero_nb;
		((int*)chol_res_matrix->p)[j6+4] = ((int*)chol_res_matrix->p)[j6+3] + 3 + od_nzero_nb;
		((int*)chol_res_matrix->p)[j6+5] = ((int*)chol_res_matrix->p)[j6+4] + 2 + od_nzero_nb;	
		int pj6   = ((int*)chol_res_matrix->p)[j6];
		int pj6_1 = ((int*)chol_res_matrix->p)[j6+1];
		int pj6_2 = ((int*)chol_res_matrix->p)[j6+2];
		int pj6_3 = ((int*)chol_res_matrix->p)[j6+3];
		int pj6_4 = ((int*)chol_res_matrix->p)[j6+4];
		int pj6_5 = ((int*)chol_res_matrix->p)[j6+5];
		// diagonal block row indices (21)
		((int*)chol_res_matrix->i)[pj6  ] = j6;   // column j6
		((int*)chol_res_matrix->i)[pj6+1] = j6+1;
		((int*)chol_res_matrix->i)[pj6+2] = j6+2;
		((int*)chol_res_matrix->i)[pj6+3] = j6+3;
		((int*)chol_res_matrix->i)[pj6+4] = j6+4;
		((int*)chol_res_matrix->i)[pj6+5] = j6+5;
		((int*)chol_res_matrix->i)[pj6_1  ] = j6+1;    // column j6+1
		((int*)chol_res_matrix->i)[pj6_1+1] = j6+2;
		((int*)chol_res_matrix->i)[pj6_1+2] = j6+3;
		((int*)chol_res_matrix->i)[pj6_1+3] = j6+4;
		((int*)chol_res_matrix->i)[pj6_1+4] = j6+5;
		((int*)chol_res_matrix->i)[pj6_2  ] = j6+2;    // column j6+2
		((int*)chol_res_matrix->i)[pj6_2+1] = j6+3;
		((int*)chol_res_matrix->i)[pj6_2+2] = j6+4;
		((int*)chol_res_matrix->i)[pj6_2+3] = j6+5;
		((int*)chol_res_matrix->i)[pj6_3  ] = j6+3;    // column j6+3
		((int*)chol_res_matrix->i)[pj6_3+1] = j6+4;
		((int*)chol_res_matrix->i)[pj6_3+2] = j6+5;
		((int*)chol_res_matrix->i)[pj6_4  ] = j6+4;    // column j6+4
		((int*)chol_res_matrix->i)[pj6_4+1] = j6+5;
		((int*)chol_res_matrix->i)[pj6_5  ] = j6+5;    // column j6+5
		// diagonal blocks row values (21)
		((double*)chol_res_matrix->x)[pj6  ] = dblocks[j18];   // column j6
		((double*)chol_res_matrix->x)[pj6+1] = dblocks[j18+1];
		((double*)chol_res_matrix->x)[pj6+2] = dblocks[j18+2];
		((double*)chol_res_matrix->x)[pj6+3] = dblocks[j18+3]; 
		((double*)chol_res_matrix->x)[pj6+4] = dblocks[j18+4];
		((double*)chol_res_matrix->x)[pj6+5] = dblocks[j18+5];
		((double*)chol_res_matrix->x)[pj6_1  ] = dblocks[j18+6];   // column j6+1
		((double*)chol_res_matrix->x)[pj6_1+1] = dblocks[j18+7];
		((double*)chol_res_matrix->x)[pj6_1+2] = -dblocks[j18+4];  // anti-symmetry
		((double*)chol_res_matrix->x)[pj6_1+3] = dblocks[j18+8];
		((double*)chol_res_matrix->x)[pj6_1+4] = dblocks[j18+9];
		((double*)chol_res_matrix->x)[pj6_2  ] = dblocks[j18+10];   // column j6+2
		((double*)chol_res_matrix->x)[pj6_2+1] = -dblocks[j18+5];   // anti-symmetry
		((double*)chol_res_matrix->x)[pj6_2+2] = -dblocks[j18+9];   // anti-symmetry
		((double*)chol_res_matrix->x)[pj6_2+3] = dblocks[j18+11];
		((double*)chol_res_matrix->x)[pj6_3  ] = dblocks[j18+12];   // column j6+3
		((double*)chol_res_matrix->x)[pj6_3+1] = dblocks[j18+13];
		((double*)chol_res_matrix->x)[pj6_3+2] = dblocks[j18+14];
		((double*)chol_res_matrix->x)[pj6_4  ] = dblocks[j18+15];   // column j6+4
		((double*)chol_res_matrix->x)[pj6_4+1] = dblocks[j18+16];
		((double*)chol_res_matrix->x)[pj6_5  ] = dblocks[j18+17];   // column j6+5
		/*****  2  : off-diagonal blocks row indices and values ***********/
		// 36 non-zero elements per block
		for(int k = odbrows_table[j]; k < odbrows_table[j+1]; k++){
			int u = 6*(k-odbrows_table[j]); 
			// we are filling the "k-odbFrows_table[j]"th off-diag block of the column.
			// For column j6, for exemple, the indices of the non-zero values are:
			// pj6 for all non-zero elements before column j6,
			// + 6 for the diagonal block of column j6
			// + u (=6*(k-odbFrows_table[j])) for the off-diag blocks of j6
			// + index inside the current block
			for(int s=0; s<6;s++){
				((int*)chol_res_matrix->i)[pj6+6 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_1+5 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_2+4 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_3+3 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_4+2 + u +s] = odbrows[k]+s;
				((int*)chol_res_matrix->i)[pj6_5+1 + u +s] = odbrows[k]+s;
			}
			int k6 = 6*k;
			int k4 = 4*k;
			int k2 = 2*k;
			((double*)chol_res_matrix->x)[pj6+6 + u   ]   = odblocks[0][k6  ]; // A   // column j6
			((double*)chol_res_matrix->x)[pj6+6 + u +1]   = odblocks[0][k6+1];
			((double*)chol_res_matrix->x)[pj6+6 + u +2]   = odblocks[0][k6+2];
			((double*)chol_res_matrix->x)[pj6+6 + u +3]   = odblocks[0][k6+3]; // B
			((double*)chol_res_matrix->x)[pj6+6 + u +4]   = odblocks[0][k6+4];
			((double*)chol_res_matrix->x)[pj6+6 + u +5]   = odblocks[0][k6+5];
			//
			((double*)chol_res_matrix->x)[pj6_1+5 + u   ] = odblocks[0][k6+1]; // A  // column j6+1
			((double*)chol_res_matrix->x)[pj6_1+5 + u +1] = odblocks[1][k4  ];
			((double*)chol_res_matrix->x)[pj6_1+5 + u +2] = odblocks[1][k4+1];
			((double*)chol_res_matrix->x)[pj6_1+5 + u +3] = -odblocks[0][k6+4]; // B
			((double*)chol_res_matrix->x)[pj6_1+5 + u +4] = odblocks[1][k4+2];
			((double*)chol_res_matrix->x)[pj6_1+5 + u +5] = odblocks[1][k4+3];
			//
			((double*)chol_res_matrix->x)[pj6_2+4 + u   ] = odblocks[0][k6+2]; // A // column j6+2
			((double*)chol_res_matrix->x)[pj6_2+4 + u +1] = odblocks[1][k4+1];
			((double*)chol_res_matrix->x)[pj6_2+4 + u +2] = odblocks[2][k2  ];
			((double*)chol_res_matrix->x)[pj6_2+4 + u +3] = -odblocks[0][k6+5]; // B
			((double*)chol_res_matrix->x)[pj6_2+4 + u +4] = -odblocks[1][k4+3];
			((double*)chol_res_matrix->x)[pj6_2+4 + u +5] = odblocks[2][k2+1];
			//
			((double*)chol_res_matrix->x)[pj6_3+3 + u   ] = odblocks[3][k6  ]; // Btilde   // column j6+3
			((double*)chol_res_matrix->x)[pj6_3+3 + u +1] = odblocks[3][k6+1];
			((double*)chol_res_matrix->x)[pj6_3+3 + u +2] = odblocks[3][k6+2];
			((double*)chol_res_matrix->x)[pj6_3+3 + u +3] = odblocks[3][k6+3]; // C
			((double*)chol_res_matrix->x)[pj6_3+3 + u +4] = odblocks[3][k6+4];
			((double*)chol_res_matrix->x)[pj6_3+3 + u +5] = odblocks[3][k6+5];
			
			((double*)chol_res_matrix->x)[pj6_4+2 + u   ] = -odblocks[3][k6+1]; // Btilde // column j6+4
			((double*)chol_res_matrix->x)[pj6_4+2 + u +1] = odblocks[4][k4  ];
			((double*)chol_res_matrix->x)[pj6_4+2 + u +2] = odblocks[4][k4+1];
			((double*)chol_res_matrix->x)[pj6_4+2 + u +3] = odblocks[3][k6+4]; // C
			((double*)chol_res_matrix->x)[pj6_4+2 + u +4] = odblocks[4][k4+2];
			((double*)chol_res_matrix->x)[pj6_4+2 + u +5] = odblocks[4][k4+3];
			//
			((double*)chol_res_matrix->x)[pj6_5+1 + u   ] = -odblocks[3][k6+2]; // Btilde // column j6+5
			((double*)chol_res_matrix->x)[pj6_5+1 + u +1] = -odblocks[4][k4+1];
			((double*)chol_res_matrix->x)[pj6_5+1 + u +2] = odblocks[5][k2  ];
			((double*)chol_res_matrix->x)[pj6_5+1 + u +3] = odblocks[3][k6+5]; // C
			((double*)chol_res_matrix->x)[pj6_5+1 + u +4] = odblocks[4][k4+3];
			((double*)chol_res_matrix->x)[pj6_5+1 + u +5] = odblocks[5][k2+1];
		}
	}
    ((int*)chol_res_matrix->p)[np6] = ((int*)chol_res_matrix->p)[np6-1]+1;
	factorizeResistanceMatrix();
}


/*************** Epetra_CrsMatrix Filling *************
 Epetra_CrsMatrix we are using are defined in row major order.
 
 Epetra_CrsMatrix must be stored completely.
 
 Epetra_CrsMatrix elements are not accessed directly for filling.
 Instead we use user friendly methods, that take one row at a time.
 
 *****************************************************/

void
StokesSolver::completeResistanceMatrix_trilinos(){
#ifdef TRILINOS
    for (int i = 0; i < res_matrix_linear_size; i++) {
		tril_res_matrix->InsertGlobalValues(i, columns_nb[i] , values[i], columns[i]);
    }
    // FillComplete matrix before building the preconditioner
    tril_res_matrix->FillComplete();
	tril_stokes_equation->setOperator(rcp(tril_res_matrix, false));
	//setDiagBlockPreconditioner();
	//	setIncCholPreconditioner();
	//	setSpInvPreconditioner();
#endif
}

void
StokesSolver::completeResistanceMatrix(){
	if (direct()) {
		completeResistanceMatrix_cholmod();
	}
#ifdef TRILINOS
	if (iterative()) {
		completeResistanceMatrix_trilinos();
	}
#endif
}

void
StokesSolver::resetResistanceMatrix(string solver_type, int nb_of_interactions, double *reset_resmat_dblocks){
	setSolverType(solver_type);
	odblocks_nb = nb_of_interactions;
	if (direct()) {
		for (int k=0; k<dblocks_size; k++) {
			dblocks[k] = reset_resmat_dblocks[k];
		}
		odbrows.clear();
		odblocks[0].resize(6*odblocks_nb);
		odblocks[1].resize(4*odblocks_nb);
		odblocks[2].resize(2*odblocks_nb);
		odblocks[3].resize(6*odblocks_nb);
		odblocks[4].resize(4*odblocks_nb);
		odblocks[5].resize(2*odblocks_nb);
		for (int k=0; k<6; k++) {
			for (unsigned int i=0; i<odblocks[k].size(); i++) {
				odblocks[k][i] = 0;
			}
			current_index_positions[k] = 0;
		}
		odbrows_table[0] = 0;
	}
#ifdef TRILINOS
	/* reset_resmat_dblocks is not adopted to reset.
	 */
	if (iterative()) {
		tril_res_matrix = new Epetra_CrsMatrix(Copy, *Map, 20*dof+dof );
		tril_l_precond = new Epetra_CrsMatrix(Copy, *Map, 3);
		tril_res_matrix->PutScalar(0);
		for (int i=0; i<res_matrix_linear_size; i++){
			for (int j=0; j<columns_max_nb; j++){
				columns[i][j] = -1;
				values[i][j] = 0;
			}
		}
		// declare the diagonal blocks
		for (int i=0; i<np; i++){
			int idof = dof*i;
			for (int j=0; j<dof; j++){
				columns[idof  ][j] = idof+j;
				columns[idof+1][j] = idof+j;
				columns[idof+2][j] = idof+j;
			}
			columns_nb[idof  ] = dof;
			columns_nb[idof+1] = dof;
			columns_nb[idof+2] = dof;
		}
	}
#endif
}

void
StokesSolver::resetRHS(){
	if (direct()) {
		for (int i=0; i<res_matrix_linear_size; i++){
			((double*)chol_rhs->x)[i] = 0;
		}
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_rhs->PutScalar(0);
	}
#endif
}

void
StokesSolver::addToRHSForce(int i, double *force_i){
	int i6 = 6*i;
	if (direct()) {
		for (int u=0; u<3; u++) {
			((double*)chol_rhs->x)[i6+u] += force_i[u];
		}
	}
#ifdef TRILINOS
	if (iterative()) {
		for(int u=0; u<3;u++){
			tril_rhs->SumIntoGlobalValue(i6+u, 0, force_i[u]);
		}
	}
#endif
}

void
StokesSolver::addToRHSForce(int i, const vec3d &force_i){
	int i6 = 6*i;
	if (direct()) {
		((double*)chol_rhs->x)[i6  ] += force_i.x;
		((double*)chol_rhs->x)[i6+1] += force_i.y;
		((double*)chol_rhs->x)[i6+2] += force_i.z;
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_rhs->SumIntoGlobalValue(i6  , 0, force_i.x);
		tril_rhs->SumIntoGlobalValue(i6+1, 0, force_i.y);
		tril_rhs->SumIntoGlobalValue(i6+2, 0, force_i.z);
	}
#endif
}


void
StokesSolver::addToRHSTorque(int i, double *torque_i){
	int i6_3 = 6*i+3;
	if (direct()) {
		for(int u=0; u<3;u++){
			((double*)chol_rhs->x)[i6_3+u] += torque_i[u];
		}
	}
#ifdef TRILINOS
	if (iterative()) {
		for(int u=0; u<3;u++){
			tril_rhs->SumIntoGlobalValue(i6_3+u, 0, torque_i[u]);
		}
	}
#endif
}

void
StokesSolver::addToRHSTorque(int i, const vec3d &torque_i){
	int i6_3 = 6*i+3;
	if (direct()) {
		((double*)chol_rhs->x)[i6_3  ] += torque_i.x;
		((double*)chol_rhs->x)[i6_3+1] += torque_i.y;
		((double*)chol_rhs->x)[i6_3+2] += torque_i.z;
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_rhs->SumIntoGlobalValue(i6_3  , 0, torque_i.x);
		tril_rhs->SumIntoGlobalValue(i6_3+1, 0, torque_i.y);
		tril_rhs->SumIntoGlobalValue(i6_3+2, 0, torque_i.z);
	}
#endif
}

void
StokesSolver::addToRHS(double *rhs){
	if (direct()) {
		for (int i=0; i<res_matrix_linear_size; i++) {
			((double*)chol_rhs->x)[i] += rhs[i];
		}
	}
	
#ifdef TRILINOS
	if (iterative()) {
		for (int i=0; i<res_matrix_linear_size; i++) {
			tril_rhs->SumIntoGlobalValue(i, 0, rhs[i]);
		}
	}
#endif
}

void
StokesSolver::setRHS(double* rhs){
	if (direct()) {
		for (int i=0; i<res_matrix_linear_size; i++) {
			((double*)chol_rhs->x)[i] = rhs[i];
		}
	}
	if (iterative()) {
		cerr << " Error : StokesSolver::setRHS(double* rhs) not implemented for TRILINOS yet ! " << endl;
		exit(1);
	}
}

void
StokesSolver::setRHSForce(int i, const vec3d &force_i){
	int i6 = 6*i;
	if (direct()) {
		((double*)chol_rhs->x)[i6  ] = force_i.x;
		((double*)chol_rhs->x)[i6+1] = force_i.y;
		((double*)chol_rhs->x)[i6+2] = force_i.z;
	}
#ifdef TRILINOS
	if (iterative()) {
		cerr << " Error : StokesSolver:: setRHSForce(int i, const vec3d &force_i) not implemented for TRILINOS yet ! " << endl;
		exit(1);
	}
#endif
}

void
StokesSolver::setRHSTorque(int i, const vec3d &torque_i){
	int i6_3 = 6*i+3;
	if (direct()) {
		((double*)chol_rhs->x)[i6_3  ] = torque_i.x;
		((double*)chol_rhs->x)[i6_3+1] = torque_i.y;
		((double*)chol_rhs->x)[i6_3+2] = torque_i.z;
	}
#ifdef TRILINOS
	if (iterative()) {
		cerr << " Error : StokesSolver:: setRHSTorque(int i, const vec3d &torque_i) not implemented for TRILINOS yet ! " << endl;
		exit(1);		
	}
#endif
}

void
StokesSolver::getRHS(double* rhs){
	if (direct()) {
		for (int i=0; i<res_matrix_linear_size; i++) {
			rhs[i] = ((double*)chol_rhs->x)[i];
		}
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_rhs->ExtractCopy(rhs);
	}
#endif
}



// Computes X = L*RHS
void
StokesSolver::compute_LTRHS(double* X){
	if (direct()) {
		/*
		  Cholmod gives a factorizationof a permutated resistance
		  matrix Lc*Lc^T = P*RFU*P^T
		  
		  That means P*L = Lc, with  L*L^T = RFU
		  
		  So for a rhs Y:
		  X = L*Y = P^T*Lc*Y
		*/
		if(!chol_L->is_ll){
			cerr << " The factorization is LDL^T. compute_LTRHS(double* X) only works for LL^T factorization." << endl;
		}

		double alpha [2] = {1,0};
		double beta [2] = {0,0};
		int transpose = 0;

		cholmod_sparse* chol_L_sparse = cholmod_factor_to_sparse(cholmod_copy_factor(chol_L, &chol_c), &chol_c);
		cholmod_sdmult(chol_L_sparse, transpose, alpha, beta, chol_rhs, chol_Psolution, &chol_c) ; // chol_Psolution = Lc*Y
		chol_solution = cholmod_solve(CHOLMOD_Pt, chol_L, chol_Psolution, &chol_c) ; // chol_solution = P^T*chol_Psolution

		for (int i=0; i<res_matrix_linear_size; i++) {
			X[i] = ((double*)chol_solution->x)[i];
		}
		cholmod_free_sparse(&chol_L_sparse, &chol_c);
		cholmod_free_dense(&chol_solution, &chol_c);
	}
#ifdef TRILINOS
	if (iterative()) {
		cerr << " StokesSolver::compute_LTRHS(double* velocity) not implemented for iterative solver." << endl;
		exit(1);
	}
#endif
}


// Finds solutions to L^T X = RHS
void
StokesSolver::solve_LT(double* X){
	if (direct()) {
		chol_PTsolution = cholmod_solve(CHOLMOD_Lt, chol_L, chol_rhs, &chol_c) ;
		chol_solution = cholmod_solve(CHOLMOD_Pt, chol_L, chol_PTsolution, &chol_c) ;
		for (int i=0; i<res_matrix_linear_size; i++) {
			X[i] = ((double*)chol_solution->x)[i];
		}
		cholmod_free_dense(&chol_solution, &chol_c);
		cholmod_free_dense(&chol_PTsolution, &chol_c);
	}
#ifdef TRILINOS
	if (iterative()) {
		cerr << " StokesSolver::solve_CholTrans(double* velocity) not implemented for iterative solver." << endl;
		exit(1);
	}
#endif
}

void
StokesSolver::solve_LT(vec3d* X, vec3d* ang_X){
	if (direct()) {
		chol_PTsolution = cholmod_solve (CHOLMOD_Lt, chol_L, chol_rhs, &chol_c) ;
		chol_solution = cholmod_solve (CHOLMOD_Pt, chol_L, chol_PTsolution, &chol_c) ;
		for (int i=0; i<np; i++) {
			int i6 = 6*i;
			X[i].x = ((double*)chol_solution->x)[i6  ];
			X[i].y = ((double*)chol_solution->x)[i6+1];
			X[i].z = ((double*)chol_solution->x)[i6+2];
			ang_X[i].x = ((double*)chol_solution->x)[i6+3];
			ang_X[i].y = ((double*)chol_solution->x)[i6+4];
			ang_X[i].z = ((double*)chol_solution->x)[i6+5];
		}				
		cholmod_free_dense(&chol_solution, &chol_c);
		cholmod_free_dense(&chol_PTsolution, &chol_c);
	}
#ifdef TRILINOS
	if (iterative()) {
		cerr << " StokesSolver::solve_CholTrans(double* velocity) not implemented for iterative solver." << endl;
		exit(1);
	}
#endif
}

void
StokesSolver::solve(vec3d* velocity, vec3d* ang_velocity){
	if (direct()) {
		chol_solution = cholmod_solve (CHOLMOD_A, chol_L, chol_rhs, &chol_c) ;
		for (int i=0; i<np; i++) {
			int i6 = 6*i;
			velocity[i].x = ((double*)chol_solution->x)[i6  ];
			velocity[i].y = ((double*)chol_solution->x)[i6+1];
			velocity[i].z = ((double*)chol_solution->x)[i6+2];
			ang_velocity[i].x = ((double*)chol_solution->x)[i6+3];
			ang_velocity[i].y = ((double*)chol_solution->x)[i6+4];
			ang_velocity[i].z = ((double*)chol_solution->x)[i6+5];
		}				
		cholmod_free_dense(&chol_solution, &chol_c);
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_stokes_equation->setLHS(tril_solution);
		tril_stokes_equation->setRHS(tril_rhs);
		bool set_success = tril_stokes_equation->setProblem();
		if (!set_success) {
			cerr << "ERROR: StokesSolver::solve : Belos::LinearProblem failed to set up correctly" << endl;
			exit(1);
		}
		tril_solver->setProblem (tril_stokes_equation);
		Belos::ReturnType ret = tril_solver->solve();
		if (ret != Belos::Converged) {
			cerr << " Warning: StokesSolver::solve : Belos::Solver did not converge" << endl;
		}
		tril_solver->getNumIters();
		//		cout << " iterations " << tril_solver->getNumIters() << endl;
		tril_solution->ExtractCopy(&velocity);
	}
#endif
}
void
StokesSolver::solve(double* velocity){
	if (direct()) {
		chol_solution = cholmod_solve (CHOLMOD_A, chol_L, chol_rhs, &chol_c) ;
		for (int i=0; i<res_matrix_linear_size; i++) {
			velocity[i] = ((double*)chol_solution->x)[i];
		}				
		cholmod_free_dense(&chol_solution, &chol_c);
	}
#ifdef TRILINOS
	if (iterative()) {
		tril_stokes_equation->setLHS(tril_solution);
		tril_stokes_equation->setRHS(tril_rhs);
		bool set_success = tril_stokes_equation->setProblem();
		if (!set_success) {
			cerr << "ERROR: StokesSolver::solve : Belos::LinearProblem failed to set up correctly" << endl;
			exit(1);
		}
		tril_solver->setProblem (tril_stokes_equation);
		Belos::ReturnType ret = tril_solver->solve();
		if (ret != Belos::Converged) {
			cerr << " Warning: StokesSolver::solve : Belos::Solver did not converge" << endl;
		}
		tril_solver->getNumIters();
		//		cout << " iterations " << tril_solver->getNumIters() << endl;
		tril_solution->ExtractCopy(&velocity);
	}
#endif
}

void
StokesSolver::convertDirectToIterative(){
#ifdef TRILINOS
	// don't free the Cholesky factor, but rememver to do it when solvingIsDone
	cholmod_free_sparse(&chol_res_matrix, &chol_c);
	cholmod_free_dense(&chol_solution, &chol_c);
	chol_L_to_be_freed = true;
	// convert RHS
	for (int i=0; i<res_matrix_linear_size;i++) {
		tril_rhs->ReplaceGlobalValue(i, 0, ((double*)chol_rhs->x)[i]);
	}
	setSolverType("iterative");
#else
	cerr << " Error: StokesSolver::convertDirectToIterative() : no iterative solver. Compile withe Trilinos. " << endl;
	exit(1);
#endif
}

void
StokesSolver::solvingIsDone(){
	if (direct()) {
		cholmod_free_factor(&chol_L, &chol_c);
		cholmod_free_sparse(&chol_res_matrix, &chol_c);
		//	cholmod_free_dense(&chol_rhs, &chol_c);
	}
#ifdef TRILINOS
	if (iterative()) {
		delete tril_res_matrix;
		delete tril_l_precond;
		if (chol_L_to_be_freed) {
			cholmod_free_factor(&chol_L, &chol_c);
		}
	}
#endif
}

/******************************************************
 *                                                     *
 *                  Private Methods                    *
 *                                                     *
 ******************************************************/

void
StokesSolver::allocateRessources(){
#ifdef TRILINOS
    int maxnum_interactionpair_per_particle = 20;
    columns_max_nb = 6*maxnum_interactionpair_per_particle;
    int numlhs = 1;
    int numrhs = 1;
    Map = rcp(new Epetra_Map(res_matrix_linear_size, 0, Comm));
    tril_solution = rcp(new Epetra_Vector(*Map, numlhs));
    tril_rhs = rcp(new Epetra_Vector(*Map, numrhs));
    //	tril_res_matrix = rcp( new MAT(res_matrix_linear_size) );
    columns = new int* [res_matrix_linear_size];
    for (int i=0; i<res_matrix_linear_size; i++) {
		columns[i] = new int [columns_max_nb];
		for (int j=0; j<columns_max_nb; j++) {
			columns[i][j] = -1;
		}
    }
    values = new double* [res_matrix_linear_size];
    for (int i=0; i<res_matrix_linear_size; i++) {
		values[i] = new double [columns_max_nb];
		for (int j=0; j<columns_max_nb; j++) {
			values[i][j] = 0.;
		}
    }
    columns_nb = new int [res_matrix_linear_size];
    for (int i=0; i<res_matrix_linear_size; i++) {
		columns_nb[i] = 0;
	}
#endif

    dblocks = new double [dblocks_size];
    odblocks = new vector <double> [6];
	current_index_positions = new int [6];
    odbrows_table = new int [np+1];

    cholmod_start (&chol_c);
	chol_init = true;

    chol_rhs = cholmod_allocate_dense(np6, 1, np6, xtype, &chol_c);
	chol_Psolution  = cholmod_allocate_dense(np6, 1, np6, xtype, &chol_c); // used for Brownian motion
	for (int i=0; i<np6; i++) {
		((double*)chol_rhs->x)[i] = 0;
	}
	
    chol_L = NULL;
}

// only needed for Cholmod
void
StokesSolver::allocateResistanceMatrix(){
	// allocate
	int nzmax; // non-zero values
	nzmax = 21*np; // diagonal blocks
	nzmax += 36*odblocks_nb;  // off-diagonal
	chol_res_matrix = cholmod_allocate_sparse(np6, np6, nzmax, sorted, packed, stype,xtype, &chol_c);
}

void
StokesSolver::doneBlocks(int i){
	if (direct()) {
		odbrows_table[i+1] = (unsigned int)odbrows.size();
	}
}

// odblocks fillings, for FT/UW version
void
StokesSolver::setColumn(const vec3d &nvec, int jj, double scaledXA, double scaledYA, double scaledYB, double scaledYBtilde, double scaledYC){
	double n0n0 = nvec.x*nvec.x;
	double n0n1 = nvec.x*nvec.y;
	double n0n2 = nvec.x*nvec.z;
	double n1n1 = nvec.y*nvec.y;
	double n1n2 = nvec.y*nvec.z;
	double n2n2 = nvec.z*nvec.z;
	double one_n0n0 = 1-n0n0;
	double one_n1n1 = 1-n1n1;
	double one_n2n2 = 1-n2n2;
	
	odbrows.push_back(6*jj);
	
	int i0 = current_index_positions[0];
	int i1 = current_index_positions[1];
	int i2 = current_index_positions[2];
	int i3 = current_index_positions[3];
	int i4 = current_index_positions[4];
	int i5 = current_index_positions[5];
	
	odblocks[0][i0  ] = scaledXA*n0n0 + scaledYA*one_n0n0; // column 0
	odblocks[0][i0+1] = (scaledXA-scaledYA)*n0n1;
	odblocks[0][i0+2] = (scaledXA-scaledYA)*n0n2;
	odblocks[0][i0+3] = 0;
	odblocks[0][i0+4] = -scaledYB*nvec.z;
	odblocks[0][i0+5] = scaledYB*nvec.y;
	odblocks[1][i1  ] = scaledXA*n1n1 + scaledYA*one_n1n1; // column 1
	odblocks[1][i1+1] = (scaledXA-scaledYA)*n1n2;
	odblocks[1][i1+2] = 0;
	//	cout << " scaledYB	" << scaledYB << " scaledYBtilde	" << scaledYBtilde << endl;
	odblocks[1][i1+3] = -scaledYB*nvec.x;
	odblocks[2][i2  ] = scaledXA*n2n2 + scaledYA*one_n2n2; // column 2
	odblocks[2][i2+1] = 0;
	odblocks[3][i3  ] = 0; // column 3
	odblocks[3][i3+1] = scaledYBtilde*nvec.z;
	odblocks[3][i3+2] = -scaledYBtilde*nvec.y;
	odblocks[3][i3+3] = scaledYC*one_n0n0;
	odblocks[3][i3+4] = -scaledYC*n0n1;
	odblocks[3][i3+5] = -scaledYC*n0n2;
	odblocks[4][i4  ] = 0; // column 4
	odblocks[4][i4+1] = scaledYBtilde*nvec.x;
	odblocks[4][i4+2] = scaledYC*one_n1n1;
	odblocks[4][i4+3] = -scaledYC*n1n2;
	odblocks[5][i5  ] = 0; // column 5
	odblocks[5][i5+1] = scaledYC*one_n2n2;
	
	current_index_positions[0] += 6;
	current_index_positions[1] += 4;
	current_index_positions[2] += 2;
	current_index_positions[3] += 6;
	current_index_positions[4] += 4;
	current_index_positions[5] += 2;
}

#ifdef TRILINOS
void
StokesSolver::setRow(const vec3d &nvec, int ii, int jj, double scaledXA, double scaledYA, double scaledYB, double scaledYBtilde, double scaledYC){
	cerr << " Error : StokesSolver::addToDiag(const vec3d &nvec, int ii, double FUvalue, double TWvalue) not implemented for TRILINOS yet ! " << endl;
	exit(1);

    int ii3 = 3*ii;
    int ii3_1 = ii3+1;
    int ii3_2 = ii3+2;
    int jj3 = 3*jj;
    int jj3_1 = jj3+1;
    int jj3_2 = jj3+2;
    
    double scaledXA_n0 = scaledXA*nvec.x;
    double scaledXA_n1 = scaledXA*nvec.y;
    double scaledXA_n2 = scaledXA*nvec.z;
    double scaledXA_n1n0 = scaledXA_n0*nvec.y;
    double scaledXA_n2n1 = scaledXA_n1*nvec.z;
    double scaledXA_n0n2 = scaledXA_n2*nvec.x;
	
    // declare ii and jj new columns, and update column nb
    int last_col_nb_ii = columns_nb[ii3];
    int last_col_nb_jj = columns_nb[jj3];
	
    columns[ii3  ][last_col_nb_ii  ] = jj3  ;
    columns[ii3  ][last_col_nb_ii+1] = jj3_1;
    columns[ii3  ][last_col_nb_ii+2] = jj3_2;
    columns[ii3_1][last_col_nb_ii  ] = jj3  ;
    columns[ii3_1][last_col_nb_ii+1] = jj3_1;
    columns[ii3_1][last_col_nb_ii+2] = jj3_2;
    columns[ii3_2][last_col_nb_ii  ] = jj3  ;
    columns[ii3_2][last_col_nb_ii+1] = jj3_1;
    columns[ii3_2][last_col_nb_ii+2] = jj3_2;
    
    columns[jj3  ][last_col_nb_jj  ] = ii3  ;
    columns[jj3  ][last_col_nb_jj+1] = ii3_1;
    columns[jj3  ][last_col_nb_jj+2] = ii3_2;
    columns[jj3_1][last_col_nb_jj  ] = ii3  ;
    columns[jj3_1][last_col_nb_jj+1] = ii3_1;
    columns[jj3_1][last_col_nb_jj+2] = ii3_2;
    columns[jj3_2][last_col_nb_jj  ] = ii3  ;
    columns[jj3_2][last_col_nb_jj+1] = ii3_1;
    columns[jj3_2][last_col_nb_jj+2] = ii3_2;
    
    columns_nb[ii3]   += 3;
    columns_nb[ii3_1] += 3;
    columns_nb[ii3_2] += 3;
    columns_nb[jj3]   += 3;
    columns_nb[jj3_1] += 3;
    columns_nb[jj3_2] += 3;
    
    // set values
    values[ii3  ][last_col_nb_ii  ] = scaledXA_n0*nvec.x; // 00
    values[ii3  ][last_col_nb_ii+1] = scaledXA_n1n0;      // 01
    values[ii3  ][last_col_nb_ii+2] = scaledXA_n0n2;      // 02
    values[ii3_1][last_col_nb_ii  ] = scaledXA_n1n0;      // 10
    values[ii3_1][last_col_nb_ii+1] = scaledXA_n1*nvec.y; // 11
    values[ii3_1][last_col_nb_ii+2] = scaledXA_n2n1;      // 12
    values[ii3_2][last_col_nb_ii  ] = scaledXA_n0n2;      // 20
    values[ii3_2][last_col_nb_ii+1] = scaledXA_n2n1;      // 21
    values[ii3_2][last_col_nb_ii+2] = scaledXA_n2*nvec.z; // 22
    
    values[jj3  ][last_col_nb_jj  ] = scaledXA_n0*nvec.x; // 00
    values[jj3  ][last_col_nb_jj+1] = scaledXA_n1n0;      // 01
    values[jj3  ][last_col_nb_jj+2] = scaledXA_n0n2;      // 02
    values[jj3_1][last_col_nb_jj  ] = scaledXA_n1n0;      // 10
    values[jj3_1][last_col_nb_jj+1] = scaledXA_n1*nvec.y; // 11
    values[jj3_1][last_col_nb_jj+2] = scaledXA_n2n1;      // 12
    values[jj3_2][last_col_nb_jj  ] = scaledXA_n0n2;      // 20
    values[jj3_2][last_col_nb_jj+1] = scaledXA_n2n1;      // 21
    values[jj3_2][last_col_nb_jj+2] = scaledXA_n2*nvec.z; // 22
}
#endif

void
StokesSolver::factorizeResistanceMatrix(){
	/*debug
	chol_c.nmethods = 1;
	//   	chol_c.method[0].ordering = CHOLMOD_NATURAL;
	//   	chol_c.method[0].ordering = CHOLMOD_GIVEN;
	int *perm = new int [np];
	int *fset = new int [np];
	for(int i = 0;i<np-1;i++){
		perm[i] = i+1;
		fset[i] = i;
	}
	perm[np-1]=0;
	fset[np-1]=np-1;

	chol_L = cholmod_analyze_p(chol_res_matrix, perm, fset, np, &chol_c);
	double beta [2] = {0,0};
	cholmod_factorize_p(chol_res_matrix, beta, fset, np, chol_L, &chol_c);
	 end debug*/


	/*reference code */
	//	chol_c.nmethods = 1;
	//   	chol_c.method[0].ordering = CHOLMOD_NATURAL;// force natural ordering (=no ordering) at the moment
	//	chol_c.postorder = 0 ;

	chol_c.supernodal = CHOLMOD_SUPERNODAL;
	chol_L = cholmod_analyze(chol_res_matrix, &chol_c);
	cholmod_factorize(chol_res_matrix, chol_L, &chol_c);
	//	cout << chol_L->ordering << endl;

    if (chol_c.status) {
		// Cholesky decomposition has failed: usually because matrix is incorrectly found to be positive-definite
		// It is very often enough to force another preconditioner to solve the problem.
		cerr << " factorization failed. forcing simplicial algorithm... " << endl;
		chol_c.supernodal = CHOLMOD_SIMPLICIAL;
		chol_L = cholmod_analyze (chol_res_matrix, &chol_c);
		cholmod_factorize (chol_res_matrix, chol_L, &chol_c);
		cerr << " factorization status " << chol_c.status << " final_ll ( 0 is LDL, 1 is LL ) " <<  chol_c.final_ll <<endl;
		chol_c.supernodal = CHOLMOD_SUPERNODAL;
    }
	
}

#ifdef TRILINOS
/*************  Preconditioners *********************/

/*
 buildDiagBlockPreconditioner() :
 A block diagonal (left-)preconditioner, almost similar to the
 one described in Amit Kumar's PhD Thesis: it is zero everywhere
 except along the diagonal where diagonal 3x3 block are the
 ones of R_FU:
 
        ..........
       |          .                             |
       | R_FU(i,j).     0                       |
       |          .                             |
       | .....................                  |
       |          .          .                  |
       |     0    . R_FU(i,j).                  |
       |          .          .                  |
       |          ............                  |
 P =   |                       .                |
       |                         .              |
       |                           .            |
       |                              ..........|
       |                             .          |
       |                             . R_FU(i,j)|
       |                             .          |
	                                  ...........
 
 This method stores P^{-1} in tril_l_precond.
 */
void
StokesSolver::setDiagBlockPreconditioner(){

    double a00, a01, a02, a11, a12, a22;
    double det, idet;
    double *precond_row = new double [3];
    int *indices = new int [3];
    for (int i = 0; i < np; i++) {
		int i3 = 3*i;
		
		indices[0] = i3;
		indices[1] = i3+1;
		indices[2] = i3+2;
		
	    // +2.5*r in the diagonal ? --> Amit Kumar, PhD Thesis
		a00 = values[i3][0];
		a01 = values[i3][1];
		a02 = values[i3][2];
		a11 = values[i3+1][1];
		a12 = values[i3+1][2];
		a22 = values[i3+2][2];
		
		det = a00*(a22*a11-a12*a12)+a01*(-a01*a22+2*a12*a02)-a02*a02*a11;
		idet = 1/det;
		
		// row i3
		precond_row[0] = idet*(a11*a22-a12*a12);
		precond_row[1] = idet*(a02*a12-a01*a22);
		precond_row[2] = idet*(a01*a12-a02*a11);
		
		tril_l_precond->InsertGlobalValues(i3, 3, precond_row, indices);
		
		// row i3+1
		precond_row[0] = precond_row[1]; // symmetric matrix!
		precond_row[1] = idet*(a00*a22-a02*a02);
		precond_row[2] = idet*(a02*a01-a00*a12);
		
		tril_l_precond->InsertGlobalValues(i3+1, 3, precond_row, indices);
		
		// row i3+2
		precond_row[1] = precond_row[2]; // symmetric matrix!
		precond_row[0] = idet*(a01*a12-a02*a11);
		precond_row[2] = idet*(a00*a11-a01*a01);
		
		tril_l_precond->InsertGlobalValues(i3+2, 3, precond_row, indices);
    }
	
    tril_l_precond->FillComplete();
	
	// give it to the LinearProblem
    tril_stokes_equation->setLeftPrec(rcp(tril_l_precond, false));
	
    delete [] precond_row;
    delete [] indices;
}

/*
 setIncCholPreconditioner() :
 A incomplete Cholesky factorization (left-)preconditioner.
 It uses Ifpack routines.
 
 Right now it seems that the routine sends back a diagonal preconditionner,
 with values being the inverse of the ones on R_FU diagonal.
 I (Romain) don't understand this behavior so far.
 */
void
StokesSolver::setIncCholPreconditioner(){
	//  parameters to be tuned
    int fill_level = 0;
	//    double drop_tolerance = 1.;
	
	RCP <Ifpack_IC> tril_ICT_precond = rcp(new Ifpack_IC(tril_res_matrix));
	
	ParameterList precondParams;
	//	precondParams.set("fact: drop tolerance", drop_tolerance);
	precondParams.set("fact: ict level-of-fill", fill_level);
	
	tril_ICT_precond->SetParameters(precondParams);
	tril_ICT_precond->Initialize();
	tril_ICT_precond->Compute();
	
	
	/*****	 TESTING *****
	 cout << " non zeros : " << tril_ICT_precond->NumGlobalNonzeros() << " " << tril_ICT_precond->IsInitialized() << " " << tril_ICT_precond->IsComputed() << endl;
	 int nb;
	 double *values = new double [res_matrix_linear_size];
	 int *indices = new int [res_matrix_linear_size];
	 
	 Epetra_CrsMatrix precU = tril_ICT_precond->U();
	 Epetra_Vector precD = tril_ICT_precond->D();
	 precD.ExtractCopy(values);
	 
	 cout << " Diagonal " << endl;
	 for(int i=0; i<res_matrix_linear_size; i++)
	 cout << i << " " << values[i] << endl;
	 
	 cout << " Upper " << endl;
	 for(int i=0; i<res_matrix_linear_size; i++){
	 precU.ExtractGlobalRowCopy(i, res_matrix_linear_size, nb, values, indices);
	 for(int j=0; j<nb; j++)
	 cout << i << " " << indices[j] << " " << values[j] << endl;
	 }
	 
	 cout << " Original Matrix Diagonal " << endl;
	 for(int i=0; i<res_matrix_linear_size; i++){
	 
	 tril_res_matrix->ExtractGlobalRowCopy(i, res_matrix_linear_size, nb, values, indices);
	 for(int j=0; j<nb; j++){
	 if(indices[j] == i )
	 cout << i << " " << indices[j] << " " << 1./values[j] << endl;
	 }
	 }
	 
	 delete [] values;
	 delete [] indices;
	 ***** END TESTING ********/
	
	// template conversion, to make Ifpack preconditioner compatible with belos
	RCP<Belos::EpetraPrecOp> belos_ICT_precond = rcp (new Belos::EpetraPrecOp(tril_ICT_precond));
	
	tril_stokes_equation->setLeftPrec(belos_ICT_precond);
}
#endif

#ifdef TRILINOS
void
StokesSolver::matrixChol2Tril(const cholmod_sparse *C, Epetra_CrsMatrix* &T){
	vector <double> row_values;
	vector <int> row_indices;
	for(int i=0; i< res_matrix_linear_size;i++){
		int nz = C->
		C->x
	}
}
#endif

#ifdef CHOLMOD_EXTRA
/*
 setSpInvPreconditioner() :
 A sparse inverse (left-)preconditioner.
 cholmod-extra routine by Jaakko Luttinen
 
 */
void
StokesSolver::setSpInvPreconditioner(){
	cholmod_sparse *sparse_inv = cholmod_spinv(chol_L, &chol_c);
	cholmod_free_sparse(&sparse_inv, &chol_c);
}
#endif

void
StokesSolver::setSolverType(string solver_type){

	if (solver_type == "direct") {
		_direct = true;
		_iterative = false;
	} else {
		if (solver_type == "iterative") {
			cerr << " Error : StokesSolver::setSolverType(string solver_type) : extended Linear algebra not implemented with TRILINOS yet. "<< endl;
			exit(1);
#ifdef TRILINOS
			_direct = false;
			_iterative = true;
#else
			cerr << " Error : StokesSolver::setSolverType(string solver_type) : 'iterative' solver asked, but needs to be compiled with Trilinos library."<< endl;
			exit(1);
#endif
		} else {
			cerr << " Error : StokesSolver::setSolverType(string solver_type) : Unknown solver type '" << solver_type << "'"<< endl;
			exit(1);
		}
	}
}

// testing
void
StokesSolver::printResistanceMatrix(ostream &out, string sparse_or_dense){
	if (direct()) {
		if(sparse_or_dense=="sparse"){
			//		out << endl<< " chol res " << endl;
			for (int i = 0; i < res_matrix_linear_size; i++) {
				for (int k =((int*)chol_res_matrix->p)[i] ; k < ((int*)chol_res_matrix->p)[i+1]; k++) {
					out << i << " " << ((int*)chol_res_matrix->i)[k] << " " << ((double*)chol_res_matrix->x)[k] << endl;
				}
			}
		}
		if(sparse_or_dense=="dense"){
			cholmod_dense *dense_res = cholmod_sparse_to_dense(chol_res_matrix,&chol_c); 
			for (int i = 0; i < res_matrix_linear_size; i++) {
				// if(i==0){
				// 	for (int j = 0; j < res_matrix_linear_size/6; j++) {
				// 		out << j << "\t \t \t \t \t \t" ;
				// 	}
				// 	out << endl;
				// }
				for (int j = 0; j < res_matrix_linear_size; j++) {
					out <<setprecision(3) <<  ((double*)dense_res->x)[i+j*res_matrix_linear_size] << "\t" ;
				}
				out << endl;
			}
			out << endl;
			cholmod_free_dense(&dense_res, &chol_c);
		}

	}
	//	exit(1);
#ifdef TRILINOS
	if (iterative()) {
		int int_nb = 100;
		double *val = new double [int_nb];
		int *ind = new int [int_nb];
		int nz;
		cout << endl<< " tril res " << endl;
		for (int i = 0; i < res_matrix_linear_size; i++) {
			tril_res_matrix->ExtractGlobalRowCopy(i, int_nb, nz, val, ind);
			//	   cout << i << " " << nz << endl;
			for (int j = 0; j < nz; j++) {
				cout << i << " " << ind[j] << " " << val[j] << endl;
			}
		}
		// cout << "precond " << endl;
		// for(int i = 0; i < res_matrix_linear_size; i++){
		//   tril_l_precond->ExtractGlobalRowCopy(i, int_nb, nz, val, ind);
		//   cout << " line " << i << " " << nz << endl;
		//    for(int j = 0; j < nz; j++){
		//      cout << i << " " << ind[j] << " " << val[j] << endl;
		//    }
		// }
		delete [] val;
		delete [] ind;
	}
#endif
}

void
StokesSolver::printFactor(ostream &out){
	if (direct()) {
		cholmod_factor* chol_L_copy = cholmod_copy_factor(chol_L, &chol_c);
		cholmod_sparse* chol_L_sparse = cholmod_transpose(cholmod_factor_to_sparse(chol_L_copy, &chol_c), 1,  &chol_c);
		cholmod_dense* chol_L_dense = cholmod_sparse_to_dense(chol_L_sparse, &chol_c);
		//		cholmod_sparse* chol_PTL_sparse = cholmod_dense_to_sparse(cholmod_solve(CHOLMOD_Pt, chol_L, chol_L_dense, &chol_c), 1, &chol_c) ; // chol_solution = P^T*chol_Psolution
		
		//		cholmod_dense* chol_LT_dense = cholmod_sparse_to_dense(cholmod_transpose(chol_L_sparse, 1, &chol_c), &chol_c);


		int transpose = 1;
		double alpha [2] = {1,0};
		double beta [2] = {0,0};

		cholmod_dense *dense_res = cholmod_sparse_to_dense(chol_res_matrix,&chol_c); 
		cholmod_sdmult(chol_L_sparse, transpose, alpha, beta, chol_L_dense, dense_res, &chol_c);

		// out << " Cholesky factor" << endl;
		// for (int i = 0; i < res_matrix_linear_size; i++) {
			
		// 		// if(i==0){
		// 		// 	for (int j = 0; j < res_matrix_linear_size/6; j++) {
		// 		// 		out << j << "\t \t \t \t \t \t" ;
		// 		// 	}
		// 		// 	out << endl;
		// 		// }
		// 	for (int j = 0; j < res_matrix_linear_size; j++) {
		// 		out <<setprecision(3) <<  ((double*)chol_L_dense->x)[i+j*res_matrix_linear_size] << "\t" ;
		// 	}
		// 	out << endl;
		// }
		// out << endl;

		//		out << " Cholesky squared  " << endl;
		for (int i = 0; i < res_matrix_linear_size; i++) {
			
				// if(i==0){
				// 	for (int j = 0; j < res_matrix_linear_size/6; j++) {
				// 		out << j << "\t \t \t \t \t \t" ;
				// 	}
				// 	out << endl;
				// }
			for (int j = 0; j < res_matrix_linear_size; j++) {
				out <<setprecision(3) <<  ((double*)dense_res->x)[i+j*res_matrix_linear_size] << "\t" ;
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

}

// testing
void
StokesSolver::printRHS(){
	if (direct()) {
		for (int i = 0; i < res_matrix_linear_size; i++) {
			cout << i << " (part " << " " << (i-i%6)/6 << " )  " << ((double*)chol_rhs->x)[i] <<  endl;
		}
	}
}
