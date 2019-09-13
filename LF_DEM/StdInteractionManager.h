#ifndef __LF_DEM__StdInteractionManager__
#define __LF_DEM__StdInteractionManager__

#include <memory>
#include <vector>
#include <map>
#include "InteractionManager.h"
#include "InteractionManagerOutput.h"
#include "StdInteraction.h"

class ParticleConfig;
class ParticleVelocity;
class ParticleVelocityGrad;


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

class StdInteractionManager : public InteractionManager<StdInteraction> {
public:
	StdInteractionManager(unsigned np,
						  std::shared_ptr<PairManager> &pairmanager,
						  ParticleConfig *config,
						  ParticleVelocity *background_vel,
						  ParticleVelocityGrad *background_velgrad,
						  std::shared_ptr<Geometry::PairwiseConfig> pd,
						  Parameters::ParameterSet *params,
						  std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver);
	StdInteractionManager(unsigned np,
						  std::shared_ptr<PairManager> &pairmanager,
						  ParticleConfig *config,
						  ParticleVelocity *background_vel,
						  ParticleVelocityGrad *background_velgrad,
						  std::shared_ptr<Geometry::PairwiseConfig> pd,
						  Parameters::ParameterSet *params,
						  std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver,
						  const std::vector <struct contact_state>& cs);
	virtual void checkNewInteractions();
	void updateInteractions(double dt, ParticleVelocity *vel);
	void updateInteractions(); // not for time evolution

	double getMaxRelativeVelocity(ParticleVelocity *vel);

	void declareForceComponents(std::map<std::string, ForceComponent> &force_components);
	void setForceToParticle(const std::string &component, std::vector<vec3d> &force, std::vector<vec3d> &torque);
	void addUpInteractionStressME(std::vector<Sym2Tensor> &stress_comp);
	void addUpInteractionStressGU(std::vector<Sym2Tensor> &stress_comp, ParticleVelocity *vel);
	void addUpContactStressXF(std::vector<Sym2Tensor> &cstress_XF, ParticleVelocity *vel);
	void addUpContactSpringStressXF(std::vector<Sym2Tensor> &cstress_XF);
	void addUpContactDashpotStressXF(std::vector<Sym2Tensor> &cstress_XF, ParticleVelocity *vel);
	void addUpRepulsiveStressXF(std::vector<Sym2Tensor> &rstress_XF);

	void setContacts(const std::vector <struct contact_state>& cs);
	std::vector <struct contact_state> getContacts() const;
	void saveState();
	void restoreState();

	
private:
	ParticleConfig *conf;
	ParticleVelocity *velinf;
	ParticleVelocityGrad *Einf;
	std::shared_ptr<Geometry::PairwiseConfig> pdist;
	Parameters::ParameterSet *p;
	std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> solver;

	double calcInteractionRange(unsigned i, unsigned j);
	void createNewInteraction(unsigned i, unsigned j, double scaled_interaction_range);
	
	void setDashpotForceToParticle(std::vector<vec3d> &force,
								   std::vector<vec3d> &torque);
	void setHydroForceToParticle_squeeze(std::vector<vec3d> &force,
										 std::vector<vec3d> &torque);
	void setHydroForceToParticle_squeeze_tangential(std::vector<vec3d> &force,
													std::vector<vec3d> &torque);
	void setContactForceToParticle(std::vector<vec3d> &force,
								   std::vector<vec3d> &torque);
	void setRepulsiveForceToParticle(std::vector<vec3d> &force,
									 std::vector<vec3d> &torque);
	void setTActAdhesionForceToParticle(std::vector<vec3d> &force,
								        std::vector<vec3d> &torque);
	void checkInputParams(std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver);

	friend void output(const StdInteractionManager &manager, ParticleVelocity *vel, const struct Parameters::ParameterSet &p, OutputData &outdata_int);
};

double maxRangeStdInteraction(const Parameters::ParameterSet &params, const std::vector<double> &radii);

bool hasPairwiseResistanceStdInteraction(const Parameters::ParameterSet &p);

}
#endif