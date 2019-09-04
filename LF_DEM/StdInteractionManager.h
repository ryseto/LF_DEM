#ifndef __LF_DEM__StdInteractionManager__
#define __LF_DEM__StdInteractionManager__

#include <memory>
#include <vector>
#include <map>
#include "InteractionManager.h"
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

namespace Interactions
{

class StdInteractionManager : public InteractionManager<StdInteraction> {
public:
	StdInteractionManager(unsigned np,
						  ParticleConfig *config,
						  ParticleVelocity *background_vel,
						  ParticleVelocityGrad *background_velgrad,
						  std::shared_ptr<Geometry::PairwiseConfig> pd,
						  std::shared_ptr<Parameters::ParameterSet> params,
						  std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver);
	StdInteractionManager(unsigned np,
						  ParticleConfig *config,
						  ParticleVelocity *background_vel,
						  ParticleVelocityGrad *background_velgrad,
						  std::shared_ptr<Geometry::PairwiseConfig> pd,
						  std::shared_ptr<Parameters::ParameterSet> params,
						  std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver,
						  const std::vector <struct contact_state>& cs);
	void checkNewInteractions();
	void updateInteractions(double dt, ParticleVelocity *vel);
	double getMaxRelativeVelocity(ParticleVelocity *vel);

	void declareForceComponents(std::map<std::string, ForceComponent> &force_components);
	void setForceToParticle(const std::string &component, std::vector<vec3d> &force, std::vector<vec3d> &torque);

	void setContacts(const std::vector <struct contact_state>& cs);

	std::vector <struct contact_state> getContacts() const;

	
private:
	ParticleConfig *conf;
	ParticleVelocity *velinf;
	ParticleVelocityGrad *Einf;
	std::shared_ptr<Geometry::PairwiseConfig> pdist;
	std::shared_ptr<Parameters::ParameterSet> p;
	std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> solver;

	double calcInteractionRange(unsigned i, unsigned j);
	void createNewInteraction(unsigned i, unsigned j, double scaled_interaction_range);
	void checkInputParams(std::shared_ptr<Parameters::ParameterSet> params,
						  std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver);
	
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

};

double maxRangeStdInteraction(const Parameters::ParameterSet &params, const std::vector<double> &radii);

}
#endif