#include "PairwiseResistanceVelocitySolver.h"
#include "StdInteractionManager.h"
#include "DimerManager.h"
#include "ForceComponent.h"

namespace Dynamics {

PairwiseResistanceVelocitySolver::PairwiseResistanceVelocitySolver(double stokes_coeff) :
sd_coeff(stokes_coeff),
pairwise_resistance_changed(true)
{}

PairwiseResistanceVelocitySolver::PairwiseResistanceVelocitySolver(double stokes_coeff, 
													 			   const std::vector<double> &particle_radius) : 
PairwiseResistanceVelocitySolver(stokes_coeff)
{
	radius = particle_radius;
	np = particle_radius.size();
	np_mobile = np;
	initStokesSolver();
}

template<class VelAssignor>
PairwiseResistanceVelocitySolver::PairwiseResistanceVelocitySolver(double stokes_coeff, 
													 const std::vector<double> &particle_radius,
													 const VelAssignor &va) : 
PairwiseResistanceVelocitySolver(stokes_coeff)
{
	radius = particle_radius;
	np = particle_radius.size();
	velo_assignor = std::unique_ptr<VelAssignor>(new VelAssignor (va));
	np_mobile = np - velo_assignor->np_fixed;
	initStokesSolver();
}

void PairwiseResistanceVelocitySolver::initStokesSolver()
{
	if (velo_assignor) {
		nb_blocks_ff.resize(velo_assignor->np_fixed, 0);
		nb_blocks_mf.resize(np_mobile, 0);
	} else {
		nb_blocks_ff.resize(0);
		nb_blocks_mf.resize(0);
	}
	nb_blocks_mm.resize(np_mobile, 0);
	stokes_solver.init(np, np_mobile);
}

void PairwiseResistanceVelocitySolver::declareResistance(unsigned p0, unsigned p1)
{
	if (p1 < np_mobile) { // i and j mobile
		nb_blocks_mm[p0]++;
	} else if (p0 >= np_mobile) { // i and j fixed
		nb_blocks_ff[p0-np_mobile]++;
	} else {
		nb_blocks_mf[p0]++;
	}
	pairwise_resistance_changed = true;
}

void PairwiseResistanceVelocitySolver::eraseResistance(unsigned p0, unsigned p1)
{
	if (p1 < np_mobile) { // i and j mobile
		nb_blocks_mm[p0]--;
	} else if (p0 >= np_mobile) { // i and j fixed
		nb_blocks_ff[p0-np_mobile]--;
	} else {
		nb_blocks_mf[p0]--;
	}
	pairwise_resistance_changed = true;
}

void PairwiseResistanceVelocitySolver::buildDiagonalBlocks()
{
	unsigned int size_mm = 0;
	for (auto bnb: nb_blocks_mm) {
		size_mm += bnb;
	}
	unsigned int size_mf = 0;
	for (auto bnb: nb_blocks_mf) {
		size_mf += bnb;
	}
	unsigned int size_ff = 0;
	for (auto bnb: nb_blocks_ff) {
		size_ff += bnb;
	}

	// create a new resistance matrix in stokes_solver
	/* [note]
	 * The resistance matrix is reset with resistance_matrix_dblock
	 */
	std::vector <struct DBlock> resistance_matrix_dblock (np);
	for (unsigned i=0; i<np; i++) {
		double FUvalue = sd_coeff*radius[i];
		double TWvalue = sd_coeff*pow(radius[i], 3)*4.0/3;
		resistance_matrix_dblock[i].col0[0] = FUvalue;
		resistance_matrix_dblock[i].col1[0] = FUvalue;
		resistance_matrix_dblock[i].col2[0] = FUvalue;
		resistance_matrix_dblock[i].col3[0] = TWvalue;
		resistance_matrix_dblock[i].col4[0] = TWvalue;
		resistance_matrix_dblock[i].col5[0] = TWvalue;
	}
	stokes_solver.resetResistanceMatrix(size_mm, size_mf, size_ff,
										resistance_matrix_dblock, pairwise_resistance_changed);
	pairwise_resistance_changed = false;
}

void PairwiseResistanceVelocitySolver::buildResistanceMatrix(const Interactions::StdInteractionManager &interaction_manager)
{
	/**
	 \brief Builds the resistance matrix

	 The built matrix \f$R_{\mathrm{FU}}\f$ (in Bossis and Brady \cite
	 brady_stokesian_1988 notations) contains pairwise resistances,
	 coming from lubrication or contact dashpots.
	 */
	buildDiagonalBlocks();
	for (unsigned i=0; i<np-1; i ++) {
		stokes_solver.startNewColumn();
		for (auto it : interaction_manager.interactions_pp[i]) {
			auto j = it->partner(i);
			if (j > i) {
				if (it->hasPairwiseResistance()) {
					stokes_solver.addResistanceBlocks(i, j,
													  it->RFU_DBlocks(),
													  it->RFU_ODBlock());
				}
			}
		}
	}
	stokes_solver.matrixFillingDone();
}

void PairwiseResistanceVelocitySolver::buildResistanceMatrix(const Interactions::StdInteractionManager &interaction_manager,
															 const Interactions::Dimer::DimerManager &dimer_manager)
{
	/**
	 \brief Builds the resistance matrix

	 The built matrix \f$R_{\mathrm{FU}}\f$ (in Bossis and Brady \cite
	 brady_stokesian_1988 notations) contains pairwise resistances,
	 coming from lubrication or contact dashpots.
	 */
	buildDiagonalBlocks();
	for (unsigned i=0; i<np-1; i ++) {
		stokes_solver.startNewColumn();
		for (auto it : interaction_manager.interactions_pp[i]) {
			auto j = it->partner(i);
			if (j > i) {
				if (it->hasPairwiseResistance()) {
					stokes_solver.addResistanceBlocks(i, j,
													  it->RFU_DBlocks(),
													  it->RFU_ODBlock());
				}
			}
		}
		for (auto it : dimer_manager.interactions_pp[i]) {
			auto j = it->partner(i);
			if (j > i) {
				stokes_solver.addResistanceBlocks(i, j,
												  it->RFU_DBlocks(),
												  it->RFU_ODBlock());
			}
		}
	}
	stokes_solver.matrixFillingDone();
}

void PairwiseResistanceVelocitySolver::setSolverRHS(const ForceComponent &fc)
{
	if (fc.has_torque) {
		for (unsigned i=0; i<np; i++) {
			stokes_solver.setRHSForce(i, fc.force[i]);
			stokes_solver.setRHSTorque(i, fc.torque[i]);
		}
	} else {
		for (unsigned i=0; i<np; i++) {
			stokes_solver.setRHSForce(i, fc.force[i]);
		}
		stokes_solver.resetRHStorque();
	}
}

void PairwiseResistanceVelocitySolver::setSolverRHS(std::vector<vec3d> &force,
													std::vector<vec3d> &torque)
{
	for (unsigned i=0; i<np; i++) {
		stokes_solver.setRHSForce(i, force[i]);
		stokes_solver.setRHSTorque(i, torque[i]);
	}
}

void PairwiseResistanceVelocitySolver::resetRHS()
{
	stokes_solver.resetRHS();
}

void PairwiseResistanceVelocitySolver::compute_LTRHS(std::vector<vec3d> &force,
													std::vector<vec3d> &torque)
{
	stokes_solver.compute_LTRHS(force, torque);
}

void PairwiseResistanceVelocitySolver::addToSolverRHS(const ForceComponent &fc)
{
	if (fc.has_torque) {
		for (unsigned i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, fc.force[i]);
			stokes_solver.addToRHSTorque(i, fc.torque[i]);
		}
	} else {
		for (unsigned i=0; i<np; i++) {
			stokes_solver.addToRHSForce(i, fc.force[i]);
		}
	}
}

void PairwiseResistanceVelocitySolver::solve(std::vector<vec3d> &na_velocity, 
											 std::vector<vec3d> &na_ang_velocity)
{
	stokes_solver.solve(na_velocity, na_ang_velocity);
}

void PairwiseResistanceVelocitySolver::solvingIsDone()
{
	stokes_solver.solvingIsDone();
}


}