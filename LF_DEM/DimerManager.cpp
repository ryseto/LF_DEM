#include "DimerManager.h"
#include "Dimer.h"
#include "PairwiseConfig.h"
#include "ParameterSet.h"

namespace Interactions {
namespace Dimer {

DimerManager::DimerManager(unsigned np,
						   std::shared_ptr<PairManager> &pairmanager,
						   ParticleConfig *config,
						   ParticleVelocity *background_vel,
						   std::shared_ptr<Geometry::PairwiseConfig> pd,
						   Parameters::ParameterSet *params,
						   std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver,
						   const std::vector <struct DimerState>& ds) :
InteractionManager<Dimer>(np, pairmanager),
conf(config),
velinf(background_vel),
pdist(pd),
p(params),
solver(vel_solver)
{
	setDimers(ds);
}

DimerManager::DimerManager(unsigned np,
						   std::shared_ptr<PairManager> &pairmanager,
						   ParticleConfig *config,
						   ParticleVelocity *background_vel,
						   std::shared_ptr<Geometry::PairwiseConfig> pd,
						   Parameters::ParameterSet *params,
						   std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver,
						   const std::vector <struct UnloadedDimerState>& uds) :
InteractionManager<Dimer>(np, pairmanager),
conf(config),
velinf(background_vel),
pdist(pd),
p(params),
solver(vel_solver)
{
	setDimers(uds);
}


void DimerManager::setDimers(const std::vector <struct DimerState> &dimer_states)
{
	DimerParams params = p->dimer;
	for (auto ds: dimer_states) {
		if (ds.p0 > interactions_pp.size() || ds.p0 > interactions_pp.size()) {
			throw std::runtime_error(" Inconsistency between configuration and dimer files.");
		}
		struct PairId pairid = {ds.p0, ds.p1, conf->radius[ds.p0], conf->radius[ds.p1]};
		addInteraction(ds.p0, ds.p1, std::make_shared<Dimer>(pairid, pdist->getSeparation(ds.p0, ds.p1), ds, params));
	}
}

void DimerManager::setDimers(const std::vector <struct UnloadedDimerState> &dimer_states)
{
	DimerParams params = p->dimer;
	for (auto ds: dimer_states) {
		if (ds.p0 > interactions_pp.size() || ds.p0 > interactions_pp.size()) {
			throw std::runtime_error(" Inconsistency between configuration and dimer files.");
		}
		struct PairId pairid = {ds.p0, ds.p1, conf->radius[ds.p0], conf->radius[ds.p1]};
		addInteraction(ds.p0, ds.p1, std::make_shared<Dimer>(pairid, pdist->getSeparation(ds.p0, ds.p1), ds, params));
	}
}

std::vector <struct DimerState> DimerManager::getDimers() const
{
	/**
		\brief Get the list of dimers with their state variables.
	 */
	std::vector <struct DimerState> ds;
	for (const auto &dimer: interactions) {
		ds.push_back(dimer->getState());
	}
	return ds;
}

void DimerManager::updateInteractions(double dt, ParticleVelocity *vel)
{
	/**
	 \brief Updates the state of active interactions.

	 To be called after particle moved.
	 Note that this routine does not look for new interactions (this is done by checkNewInteraction), it only updates already known active interactions.
	 It however desactivate interactions removes interactions that became inactive (ie when the distance between particles gets larger than the interaction range).

	 */
	unsigned i, j;
	struct PairVelocity pvel;

	pdist->setVelocityState(vel);
	for (auto &dimer: interactions) {
		std::tie(i, j) = dimer->get_par_num();
		pdist->getVelocities(i, j, pvel);
		dimer->applyTimeStep(dt, pdist->getSeparation(i, j), pvel);
	}
}

void DimerManager::saveState()
{
	for (auto &inter: interactions) {
		inter->saveState();
	}
}

void DimerManager::restoreState()
{
	for (auto &inter: interactions) {
		inter->restoreState();
	}
}

double DimerManager::getMaxRelativeVelocity(ParticleVelocity *vel)
{
	/**
	 */
	unsigned i, j;
	struct PairVelocity pvel;
	pdist->setVelocityState(vel);

	double max_velocity = 0;
	for (const auto &dimer: interactions) {
		std::tie(i, j) = dimer->get_par_num();
		pdist->getVelocities(i, j, pvel);

		double velocity = dimer->getContactVelocity(pvel).norm();
		if (velocity > max_velocity) {
			max_velocity = velocity;
		}
		velocity = dimer->getRotationDiffVelocity(pvel).norm();
		if (velocity > max_velocity) {
			max_velocity = velocity;
		}
	}
	return max_velocity;
}

void DimerManager::declareForceComponents(std::map<std::string, ForceComponent> &force_components)
{
	/******* Contact force, spring part ***********/
	unsigned np = conf->position.size();
	force_components["dimer_spring"] = ForceComponent(np, RateDependence::independent, true);
	declared_forces.push_back("dimer_spring");
	force_components["dimer_dashpot"] = ForceComponent(np, RateDependence::proportional, true);
	declared_forces.push_back("dimer_dashpot");	
}

void DimerManager::setForceToParticle(const std::string &component, std::vector<vec3d> &force, std::vector<vec3d> &torque)
{
	if (component == "dimer_spring") {
		setSpringForceToParticle(force, torque);
	} else if (component == "dimer_dashpot") {
		setDashpotForceToParticle(force, torque);
	} else {
		throw std::runtime_error(" DimerManager::Unknown force component \""+component+"\"");
	}
}

void DimerManager::setDashpotForceToParticle(std::vector<vec3d> &force,
											 std::vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	vec3d f0, f1, t0, t1;
	unsigned int i, j;
	struct PairVelocity pvel;

	pdist->setVelocityState(velinf);
	for (const auto &dimer: interactions) {
		std::tie(i, j) = dimer->get_par_num();
		pdist->getVelocities(i, j, pvel);
		std::tie(f0, f1, t0, t1) = dimer->getForceTorqueDashpot(pvel);
		force[i] += f0;
		force[j] += f1;
		torque[i] += t0;
		torque[j] += t1;
	}
}

void DimerManager::setSpringForceToParticle(std::vector<vec3d> &force,
											 std::vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	vec3d f0, f1, t0, t1;
	unsigned int i, j;

	for (const auto &dimer: interactions) {
		std::tie(i, j) = dimer->get_par_num();
		std::tie(f0, f1, t0, t1) = dimer->getForceTorqueSpring();
		force[i] += f0;
		force[j] += f1;
		torque[i] += t0;
		torque[j] += t1;
	}
}

void DimerManager::addUpStress(std::vector<Sym2Tensor> &cstress_XF, ParticleVelocity *vel)
{
	struct PairVelocity pvel;
	pdist->setVelocityState(vel);
	for (const auto &dimer: interactions) {
		unsigned int i, j;
		std::tie(i, j) = dimer->get_par_num();
		pdist->getVelocities(i, j, pvel);
		auto stress = dimer->getStress(pvel);
		cstress_XF[i] += stress;
		cstress_XF[j] += stress;
	}
}

void DimerManager::addUpStressDashpot(std::vector<Sym2Tensor> &cstress_XF, ParticleVelocity *vel)
{
	struct PairVelocity pvel;
	pdist->setVelocityState(vel);
	for (const auto &dimer: interactions) {
		unsigned int i, j;
		std::tie(i, j) = dimer->get_par_num();
		pdist->getVelocities(i, j, pvel);
		auto stress = dimer->getDashpotStress(pvel);
		cstress_XF[i] += stress;
		cstress_XF[j] += stress;
	}
}

void DimerManager::addUpStressSpring(std::vector<Sym2Tensor> &cstress_XF)
{
	for (const auto &dimer: interactions) {
		unsigned int i, j;
		std::tie(i, j) = dimer->get_par_num();
		auto stress = dimer->getSpringStress();
		cstress_XF[i] += stress;
		cstress_XF[j] += stress;
	}
}

bool hasPairwiseResistanceDimer(const Parameters::ParameterSet &p) {
	return true;
}

} // namespace Dimer
} // namespace Interactions
