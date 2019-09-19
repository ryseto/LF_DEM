#include "StdInteractionManager.h"
#include "StdInteraction.h"
#include "PairwiseConfig.h"
#include "ParameterSet.h"

namespace Interactions {

StdInteractionManager::StdInteractionManager(unsigned np,
											 std::shared_ptr<PairManager> &pairmanager,
											 ParticleConfig *config,
											 ParticleVelocity *background_vel,
											 ParticleVelocityGrad *background_velgrad,
											 std::shared_ptr<Geometry::PairwiseConfig> pd,
											 Parameters::ParameterSet *params,
											 std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver) :
InteractionManager<StdInteraction>(np, pairmanager),
conf(config),
velinf(background_vel),
Einf(background_velgrad),
pdist(pd),
p(params),
solver(vel_solver)
{
	checkInputParams(vel_solver);
	checkNewInteractions();
}

StdInteractionManager::StdInteractionManager(unsigned np,
											 std::shared_ptr<PairManager> &pairmanager,
											 ParticleConfig *config,
											 ParticleVelocity *background_vel,
											 ParticleVelocityGrad *background_velgrad,
											 std::shared_ptr<Geometry::PairwiseConfig> pd,
											 Parameters::ParameterSet *params,
											 std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver,
											 const std::vector <struct contact_state>& cs) :
InteractionManager<StdInteraction>(np, pairmanager),
conf(config),
velinf(background_vel),
Einf(background_velgrad),
pdist(pd),
p(params),
solver(vel_solver)
{
	checkInputParams(vel_solver);
	checkNewInteractions();
	setContacts(cs);
}

void StdInteractionManager::checkInputParams(std::shared_ptr<Dynamics::PairwiseResistanceVelocitySolver> vel_solver)
{
	if (p->brownian > 0 && p->lub.model != "none" && (p->contact.relaxation_time > 0 || p->contact.relaxation_time_tan >0) ) {
		p->contact.relaxation_time = -1;
		p->contact.relaxation_time_tan = -1;
		std::cerr << " WARNING: overriding user input of contact relaxation times in the case of Brownian simulation with lubrication." << std::endl;
	}

	if (!vel_solver && (has_lubrication(p->lub) || has_dashpot(p->contact))) {
		throw std::runtime_error(" InteractionManager:: Error: Dashpot and lubrication require a pairwise velocity solver.");
	}

	if (p->contact.critical_load > 0) {
		if (p->contact.friction_model != FrictionModel::criticalload) {
			p->contact.friction_model = FrictionModel::criticalload;
			std::cerr << "Warning : critical load simulation -> switched to friction_model=2" << std::endl;
		}
	}
}

double StdInteractionManager::calcInteractionRange(unsigned i, unsigned j)
{
	auto range = calcContactRange(conf->radius[i], conf->radius[j]);
	if (has_lubrication(p->lub)) {
		auto r = Lub::calcLubricationRange(p->lub.max_gap, conf->radius[i], conf->radius[j]);
		if (r > range) {
			range = r;
		}
	}
	if (has_repulsion(p->repulsion)) {
		auto r = calcRepulsiveForceRange(p->repulsion, conf->radius[i], conf->radius[j]);
		if (r > range) {
			range = r;
		}
	}
	return range;
}



void StdInteractionManager::createNewInteraction(unsigned i, unsigned j, double scaled_interaction_range)
{
	// new interaction
	Interactions::StdInteractionParams params;
	params.add(p->contact);
	if (has_lubrication(p->lub)) {
		params.add(p->lub);
		params.lubp->max_gap = 2*(Lub::calcLubricationRange(p->lub.max_gap, conf->radius[i], conf->radius[j])/(conf->radius[i]+conf->radius[j]) - 1);
	}
	if (has_repulsion(p->repulsion)) {
		params.add(p->repulsion);
		params.repp->max_gap = calcRepulsiveForceRange(p->repulsion, conf->radius[i], conf->radius[j]);
	}

	struct PairId pairid;
	pairid.p0 = i;
	pairid.p1 = j;
	pairid.a0 = conf->radius[i];
	pairid.a1 = conf->radius[j];
	addInteraction(i, j, std::make_shared<StdInteraction>(pairid, pdist->getSeparation(i, j), scaled_interaction_range, 
															std::move(params), solver.get()));
}

void StdInteractionManager::checkNewInteractions()
{
	/**
	 \brief Checks if there are new pairs of interacting particles. If so, creates and sets up the corresponding Interaction objects.

	 To be called after particle moved.
	 */
	for (unsigned i=0; i<conf->position.size()-1; i++) {
		for (unsigned j : pdist->neighborhood(i)) {
			if (j > i) {
				if (!pair_manager->areInteracting(i, j)) {
					double sq_dist = pdist->getSeparation(i, j).sq_norm();
					double scaled_interaction_range = calcInteractionRange(i, j);
					double sq_dist_lim = scaled_interaction_range*scaled_interaction_range;
					if (sq_dist < sq_dist_lim) {
						createNewInteraction(i, j, scaled_interaction_range);
					}
				}
			}
		}
	}
}

void StdInteractionManager::updateInteractions(double dt, ParticleVelocity *vel)
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
	for (unsigned k=0; k<interactions.size(); k++) {
		bool deactivated = false;
		std::tie(i, j) = interactions[k]->get_par_num();
		pdist->getVelocities(i, j, pvel);
		interactions[k]->updateState(pvel,
									 pdist->getSeparation(i, j), 
									 dt,
									 deactivated);
		if (deactivated) {
			removeInteraction(k);
		}
	}
}

void StdInteractionManager::updateInteractions()
{
	/**
	 \brief Updates the state of active interactions in absence of time evolution.
	 */
	unsigned i, j;
	struct PairVelocity pvel;
	for (unsigned k=0; k<interactions.size(); k++) {
		bool deactivated = false;
		std::tie(i, j) = interactions[k]->get_par_num();
		interactions[k]->updateState(pvel,
									 pdist->getSeparation(i, j), 
									 0,
									 deactivated);
		if (deactivated) {
			removeInteraction(k);
		}
	}
}

void StdInteractionManager::saveState()
{
	for (auto &inter: interactions) {
		inter->saveState();
	}
}

void StdInteractionManager::restoreState()
{
	for (auto &inter: interactions) {
		inter->restoreState();
	}
}

double StdInteractionManager::getMaxRelativeVelocity(ParticleVelocity *vel)
{
	/**
	 */
	unsigned i, j;
	struct PairVelocity pvel;
	pdist->setVelocityState(vel);

	double max_velocity = 0;
	bool friction = has_friction(p->contact.friction_model);
	bool rolling = has_rolling_friction(p->contact);

	for (const auto &inter: interactions) {
		std::tie(i, j) = inter->get_par_num();
		pdist->getVelocities(i, j, pvel);

		double velocity = std::abs(inter->getNormalVelocity(pvel));
		if (velocity > max_velocity) {
			max_velocity = velocity;
		}
		if (inter->contact) {
			if (friction) {
				velocity = inter->contact->getSlidingVelocity(pvel).norm();
				if (velocity > max_velocity) {
					max_velocity = velocity;
				}
			}
			if (rolling) {
				velocity = inter->contact->getRollingVelocity(pvel).norm();
				if (velocity > max_velocity) {
					max_velocity = velocity;
				}
			}
		}
	}
	return max_velocity;
}

void StdInteractionManager::declareForceComponents(std::map<std::string, ForceComponent> &force_components)
{
	// Only declare in force components the forces on the rhs of
	// R_FU*(U-U_inf) = RHS
	// These forces will be used to compute the na_velo_components,
	// a set of components that must add up exactly to the total non-affine velocity

	/******* Contact force, spring part ***********/
	bool torque = true;

	/******* Contact force, spring part ***********/
	unsigned np = conf->position.size();
	if (Interactions::has_friction(p->contact.friction_model)) {
		force_components["contact"] = ForceComponent(np, RateDependence::independent, torque);
		declared_forces.push_back("contact");
	} else {
		force_components["contact"] = ForceComponent(np, RateDependence::independent, !torque);
		declared_forces.push_back("contact");
	}
	
	if (Interactions::has_dashpot(p->contact)) {
		/******* Contact force, dashpot part, U_inf only, i.e. R_FU^{dashpot}.U_inf ***********/
		force_components["dashpot"] = ForceComponent(np, RateDependence::proportional, torque);
		declared_forces.push_back("dashpot");
	}

	/*********** Hydro force, i.e.  R_FE:E_inf *****************/
	if (p->lub.model == "normal") {
		force_components["hydro"] = ForceComponent(np, RateDependence::proportional, !torque);
		declared_forces.push_back("hydro");
	}
	if (p->lub.model == "tangential") {
		force_components["hydro"] = ForceComponent(np, RateDependence::proportional, torque);
		declared_forces.push_back("hydro");
	}
	
	if (Interactions::has_repulsion(p->repulsion)) {
		force_components["repulsion"] = ForceComponent(np, RateDependence::independent, !torque);
		declared_forces.push_back("repulsion");
	}

	// if (Interactions::has_TActAdhesion(p->TA_adhesion)) {
	// 	force_components["delayed_adhesion"] = ForceComponent(np, RateDependence::independent, !torque);
	// 	declared_forces.push_back("delayed_adhesion");
	// }
}

void StdInteractionManager::setForceToParticle(const std::string &component, std::vector<vec3d> &force, std::vector<vec3d> &torque)
{
	if (component == "contact") {
		setContactForceToParticle(force, torque);
	} else if (component == "dashpot") {
		setDashpotForceToParticle(force, torque);
	} else if (component == "hydro") {
		if (p->lub.model == "normal") {
			setHydroForceToParticle_squeeze(force, torque);
		}
		if (p->lub.model == "tangential") {
			setHydroForceToParticle_squeeze_tangential(force, torque);
		}
	} else if (component == "repulsion") {
		setRepulsiveForceToParticle(force, torque);
	} else if (component == "delayed_adhesion") {
		setTActAdhesionForceToParticle(force, torque);
	} else {
		throw std::runtime_error(" StdInteractionManager::Unknown force component.");
	}
}



void StdInteractionManager::setDashpotForceToParticle(std::vector<vec3d> &force,
											   		  std::vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	vec3d GEi, GEj, HEi, HEj;
	unsigned int i, j;
	struct PairVelocity pvel;

	pdist->setVelocityState(velinf);
	for (const auto &inter: interactions) {
		if (inter->contact && inter->contact->dashpot) {
			std::tie(i, j) = inter->get_par_num();
			pdist->getVelocities(i, j, pvel);
			std::tie(GEi, GEj, HEi, HEj) = inter->contact->dashpot->getForcesTorques(pvel);
			force[i] += GEi;
			force[j] += GEj;
			torque[i] += HEi;
			torque[j] += HEj;
		}
	}
}

void StdInteractionManager::setHydroForceToParticle_squeeze(std::vector<vec3d> &force,
											 				std::vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	vec3d GEi, GEj;
	unsigned int i, j;
	for (const auto &inter: interactions) {
		if (inter->lubrication) {
			std::tie(i, j) = inter->get_par_num();
			std::tie(GEi, GEj) = inter->lubrication->calcGE(0.5*(Einf->E[i] + Einf->E[j])); // G*E_\infty term
			force[i] += GEi;
			force[j] += GEj;
		}
	}
}

void StdInteractionManager::setHydroForceToParticle_squeeze_tangential(std::vector<vec3d> &force,
																	   std::vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	vec3d GEi, GEj, HEi, HEj;
	unsigned int i, j;
	for (const auto &inter: interactions) {
		if (inter->lubrication) {
			std::tie(i, j) = inter->get_par_num();
			std::tie(GEi, GEj, HEi, HEj) = inter->lubrication->calcGEHE(0.5*(Einf->E[i] + Einf->E[j])); // G*E_\infty term, no gamma dot
			force[i] += GEi;
			force[j] += GEj;
			torque[i] += HEi;
			torque[j] += HEj;
		}
	}
}

void StdInteractionManager::setContactForceToParticle(std::vector<vec3d> &force,
									   			std::vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	for (const auto &inter: interactions) {
		if (inter->contact) {
			inter->contact->addUpForceTorque(force, torque);
		}
	}
}

void StdInteractionManager::setRepulsiveForceToParticle(std::vector<vec3d> &force,
										 		 std::vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	for (const auto &inter: interactions) {
		inter->repulsion->addUpForce(force);
	}
}

void StdInteractionManager::setTActAdhesionForceToParticle(std::vector<vec3d> &force,
										            std::vector<vec3d> &torque)
{
	for (auto &f: force) {
		f.reset();
	}
	for (auto &t: torque) {
		t.reset();
	}
	unsigned int i, j;
	for (const auto &inter: interactions) {
		std::tie(i, j) = inter->get_par_num();
		inter->delayed_adhesion->addUpForce(force[i], force[j]);
	}
}


void StdInteractionManager::setContacts(const std::vector <struct contact_state>& cs)
{
	/**
		\brief Set a list of contacts with their state variables.

		Used to restart the simulation from a given state.
	 */
	for (const auto& c : cs) {
		bool existing_interaction = false;
		for (auto &inter: interactions) {
			unsigned int p0, p1;
			std::tie(p0, p1) = inter->get_par_num();
			if (p0 == c.p0 && p1 == c.p1) {
				inter->contact->setState(c);
				existing_interaction = true;
			}
		}
		if (!existing_interaction) {
			throw std::runtime_error(" StdInteractionManager:: inconsistent initial interactions.");			
		}
	}
}

std::vector <struct contact_state> StdInteractionManager::getContacts() const
{
	/**
		\brief Get the list of contacts with their state variables.

		Used to output a configuration including contact info. Useful if you want to restart from exact same configuration.
	 */
	std::vector <struct contact_state> cs;
	for (const auto &inter: interactions) {
		if (inter->contact) {
			cs.push_back(inter->contact->getState());
		}
	}
	return cs;
}

void StdInteractionManager::addUpInteractionStressME(std::vector<Sym2Tensor> &stress_comp)
{
	for (const auto &inter: interactions) {
		if (inter->lubrication) {
			unsigned int i, j;
			std::tie(i, j) = inter->get_par_num();
			inter->lubrication->addMEStresslet(0.5*(Einf->E[i] + Einf->E[j]),
											 stress_comp[i],
											 stress_comp[j]); // R_SE:Einf-R_SU*v
		}
	}
}

void StdInteractionManager::addUpInteractionStressGU(std::vector<Sym2Tensor> &stress_comp,
													 ParticleVelocity *vel)
{
	struct PairVelocity pvel;
	pdist->setVelocityState(vel);
	for (const auto &inter: interactions) {
		if (inter->lubrication) {
			unsigned int i, j;
			std::tie(i, j) = inter->get_par_num();
			pdist->getVelocities(i, j, pvel);
			inter->lubrication->addGUStresslet(pvel, stress_comp[i], stress_comp[j]);
		}
	}
}

void StdInteractionManager::addUpContactStressXF(std::vector<Sym2Tensor> &cstress_XF, ParticleVelocity *vel)
{
	struct PairVelocity pvel;
	pdist->setVelocityState(vel);
	for (auto &inter: interactions) {
		if (inter->contact) {
			unsigned int i, j;
			std::tie(i, j) = inter->get_par_num();
			pdist->getVelocities(i, j, pvel);
			inter->contact->addUpStress(cstress_XF[i], cstress_XF[j], pvel); // - rF_cont
		}
	}
}

void StdInteractionManager::addUpContactSpringStressXF(std::vector<Sym2Tensor> &cstress_XF)
{
	for (const auto &inter: interactions) {
		if (inter->contact) {
			unsigned int i, j;
			std::tie(i, j) = inter->get_par_num();
			inter->contact->addUpStressSpring(cstress_XF[i], cstress_XF[j]); // - rF_cont
		}
	}
}

void StdInteractionManager::addUpContactDashpotStressXF(std::vector<Sym2Tensor> &cstress_XF, ParticleVelocity *vel)
{
	struct PairVelocity pvel;
	pdist->setVelocityState(vel);
	for (const auto &inter: interactions) {
		if (inter->contact) {
			if (inter->contact->dashpot) {
				unsigned int i, j;
				std::tie(i, j) = inter->get_par_num();
				pdist->getVelocities(i, j, pvel);
				inter->contact->addUpStressDashpot(cstress_XF[i], cstress_XF[j], pvel); // - rF_cont
			}
		}
	}
}

void StdInteractionManager::addUpRepulsiveStressXF(std::vector<Sym2Tensor> &rstress_XF)
{
	for (auto &inter: interactions) {
		if (inter->repulsion) {
			unsigned int i, j;
			std::tie(i, j) = inter->get_par_num();
			inter->repulsion->addUpStressXF(rstress_XF[i], rstress_XF[j]); // - rF_rep
		}
	}
}

// void StdInteractionManager::addUpDelayedAdhesionStressXF(std::vector<Sym2Tensor> &rstress_XF)
// {
// 	for (auto &inter: interactions) {
// 		if (inter->delayed_adhesion) {
// 			unsigned int i, j;
// 			std::tie(i, j) = inter.get_par_num();
// 			inter->delayed_adhesion->addUpStressXF(rstress_XF[i], rstress_XF[j], inter->rvec); // - rF_rep
// 		}
// 	}
// }


double maxRangeStdInteraction(const Parameters::ParameterSet &params, const std::vector<double> &radii)
{
	double max_range = 0;
	for (auto r1: radii) {
		for (auto r2: radii) {
			auto range = calcContactRange(r1, r2);
			if (has_lubrication(params.lub)) {
				auto r = Lub::calcLubricationRange(params.lub.max_gap, r1, r2);
				if (r > range) {
					range = r;
				}
			}
			if (has_repulsion(params.repulsion)) {
				auto r = calcRepulsiveForceRange(params.repulsion, r1, r2);
				if (r > range) {
					range = r;
				}
			}
			if (range > max_range) {
				max_range = range;
			}
		}
	}
	return max_range;
}

bool hasPairwiseResistanceStdInteraction(const Parameters::ParameterSet &p) {
	return Interactions::has_lubrication(p.lub) || Interactions::has_dashpot(p.contact);

}

} // namespace Interactions