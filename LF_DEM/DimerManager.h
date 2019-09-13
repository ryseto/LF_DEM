#ifndef __LF_DEM__DimerManager__
#define __LF_DEM__DimerManager__

#include <memory>
#include <vector>
#include <map>
#include "InteractionManager.h"
#include "InteractionManagerOutput.h"
#include "Dimer.h"

class ParticleConfig;
class ParticleVelocity;

namespace Parameters {
	struct ParameterSet;
}

namespace Geometry {
	class PairwiseConfig;
}

namespace Dynamics {
	class PairwiseResistanceVelocitySolver;
}

struct ForceComponent;
namespace Interactions
{
namespace Dimer
{

class DimerManager : public InteractionManager<Dimer> {
public:
	DimerManager(unsigned np,
				 std::shared_ptr<PairManager> &pairmanager,
				 ParticleConfig *config,
				 ParticleVelocity *background_vel,
				 std::shared_ptr<Geometry::PairwiseConfig> pd,
				 Parameters::ParameterSet *params,
				 std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver,
				 const std::vector <struct DimerState>& ds);
	DimerManager(unsigned np,
				 std::shared_ptr<PairManager> &pairmanager,
				 ParticleConfig *config,
				 ParticleVelocity *background_vel,
				 std::shared_ptr<Geometry::PairwiseConfig> pd,
				 Parameters::ParameterSet *params,
				 std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver,
				 const std::vector <struct UnloadedDimerState>& uds);
	void updateInteractions(double dt, ParticleVelocity *vel);

	double getMaxRelativeVelocity(ParticleVelocity *vel);

	void declareForceComponents(std::map<std::string, ForceComponent> &force_components);
	void setForceToParticle(const std::string &component, std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void addUpStress(std::vector<Sym2Tensor> &cstress_XF, ParticleVelocity *vel);
	void addUpStressDashpot(std::vector<Sym2Tensor> &cstress_XF, ParticleVelocity *vel);
	void addUpStressSpring(std::vector<Sym2Tensor> &cstress_XF);

	void setDimers(const std::vector <struct DimerState>& ds);
	void setDimers(const std::vector <struct UnloadedDimerState>& ds);

	std::vector <struct DimerState> getDimers() const;
	void saveState();
	void restoreState();

	
private:
	ParticleConfig *conf;
	ParticleVelocity *velinf;
	std::shared_ptr<Geometry::PairwiseConfig> pdist;
	Parameters::ParameterSet *p;
	std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> solver;

	void createNewInteraction(unsigned i, unsigned j);
	void setDashpotForceToParticle(std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void setSpringForceToParticle(std::vector<vec3d> &force, std::vector<vec3d> &torque);

	friend void output(const DimerManager &manager, ParticleVelocity *vel, const struct Parameters::ParameterSet &p, OutputData &outdata_int);
};

double maxRangeStdInteraction(const Parameters::ParameterSet &params, const std::vector<double> &radii);

bool hasPairwiseResistanceStdInteraction(const Parameters::ParameterSet &p);

} // namespace Dimer
} // namespace Interactions
#endif