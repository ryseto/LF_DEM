#ifndef __LF_DEM__PairwiseResistanceVelocitySolver__
#define __LF_DEM__PairwiseResistanceVelocitySolver__

#include <vector>
#include <memory>

#include "VelocityAssignor.h"
#include "ForceComponent.h"
#include "StokesSolver.h"

namespace Interactions {
	class InteractionSet;
}

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
	void buildResistanceMatrix(const Interactions::InteractionSet &interactions);
	void setSolverRHS(const ForceComponent &fc);
	void addToSolverRHS(const ForceComponent &fc);

	void setSolverRHS(vector<vec3d> &force,
					  vector<vec3d> &torque);
	void resetRHS();
	void compute_LTRHS(vector<vec3d> &force,
					   vector<vec3d> &torque);
	void solve(vector<vec3d> &na_velocity, 
				vector<vec3d> &na_ang_velocity); // get V

	void computeVelocity(std::vector<vec3d> &vel, std::vector<vec3d> &ang_vel,
						 const std::vector<vec3d> &force, const std::vector<vec3d> &torque);
private:
	PairwiseResistanceVelocitySolver(double stokes_coeff);
	unsigned np_mobile;
	double sd_coeff;
	std::vector<double> radius;

	std::vector<unsigned int> nb_blocks_mm;
	std::vector<unsigned int> nb_blocks_mf;
	std::vector<unsigned int> nb_blocks_ff;
	bool pairwise_resistance_changed;

	StokesSolver stokes_solver;
	std::unique_ptr<VelocityAssignor> velo_assignor;
};

}
#endif