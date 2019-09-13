#ifndef __LF_DEM__PairwiseResistanceVelocitySolver__
#define __LF_DEM__PairwiseResistanceVelocitySolver__

#include <vector>
#include <memory>
#include <map>

#include "VelocityAssignor.h"
#include "ParticleConfig.h"
#include "StokesSolver.h"

namespace Interactions {
	class StdInteractionManager;
namespace Dimer {
	class DimerManager;
}
}
struct ForceComponent;

namespace Dynamics {

class PairwiseResistanceVelocitySolver {
public:
	PairwiseResistanceVelocitySolver(double stokes_coeff, 
							         const std::vector<double> &particle_radius);
	template<class VelAssignor>
	PairwiseResistanceVelocitySolver(double stokes_coeff, 
									 const std::vector<double> &particle_radius,
									 const VelAssignor &va);

	void declareResistance(unsigned p0, unsigned p1);
	void eraseResistance(unsigned p0, unsigned p1);
	void buildResistanceMatrix(const Interactions::StdInteractionManager &interaction_manager);
	void buildResistanceMatrix(const Interactions::StdInteractionManager &interaction_manager,
							   const Interactions::Dimer::DimerManager &dimer_manager);
	void setSolverRHS(const ForceComponent &fc);
	void addToSolverRHS(const ForceComponent &fc);

	void setSolverRHS(std::vector<vec3d> &force,
					  std::vector<vec3d> &torque);
	void resetRHS();
	void compute_LTRHS(std::vector<vec3d> &force,
					   std::vector<vec3d> &torque);
	void solve(std::vector<vec3d> &na_velocity, 
				std::vector<vec3d> &na_ang_velocity); // get V

	void computeVelocity(std::vector<vec3d> &vel, std::vector<vec3d> &ang_vel,
						 const std::vector<vec3d> &force, const std::vector<vec3d> &torque);
	void solvingIsDone();

	std::unique_ptr<VelocityAssignor> velo_assignor;
private:
	PairwiseResistanceVelocitySolver(double stokes_coeff);
	void buildDiagonalBlocks();
	void initStokesSolver();

	unsigned np_mobile;
	unsigned np;

	double sd_coeff;
	std::vector<double> radius;

	std::vector<unsigned int> nb_blocks_mm;
	std::vector<unsigned int> nb_blocks_mf;
	std::vector<unsigned int> nb_blocks_ff;
	bool pairwise_resistance_changed;

	StokesSolver stokes_solver;
};

}
#endif