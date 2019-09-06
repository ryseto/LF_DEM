#ifndef __LF_DEM__InteractionManager__
#define __LF_DEM__InteractionManager__

#include <memory>
#include <vector>
#include <string>

#include "ForceComponent.h"

class StdInteractionManagerOutput;

namespace Interactions
{

template <class InteractionT>
class InteractionManager {
public:
	InteractionManager(unsigned np) :
		interactions_pp(np),
		particle_partners(np) {};
	virtual void checkNewInteractions() = 0;
	// virtual void updateInteractions(double dt) = 0;
	virtual void declareForceComponents(std::map<std::string, ForceComponent> &force_components) = 0;
	virtual void setForceToParticle(const std::string &component, std::vector<vec3d> &force, std::vector<vec3d> &torque) = 0;

	std::vector< std::string > declared_forces;
	std::vector< std::shared_ptr<InteractionT> > interactions;
	std::vector< std::vector< std::shared_ptr<InteractionT> > > interactions_pp; // per particles

	bool areInteracting(unsigned i, unsigned j) const;

	std::size_t size() const {return interactions.size();};

	typename std::vector< std::shared_ptr<InteractionT> >::iterator begin() {return interactions.begin();}
	typename std::vector< std::shared_ptr<InteractionT> >::iterator end() {return interactions.end();}
	typename std::vector< std::shared_ptr<InteractionT> >::const_iterator begin() const {return interactions.begin();}
	typename std::vector< std::shared_ptr<InteractionT> >::const_iterator end() const {return interactions.end();}

protected:
	void removeInteraction(unsigned inter_idx);
	void addInteraction(unsigned i, unsigned j, std::shared_ptr<InteractionT> inter);

private:
	std::vector< std::vector<unsigned> >particle_partners;
	void removePartners(unsigned i, unsigned j);
	void removeInteractionPP(unsigned i, std::shared_ptr<InteractionT> dead_inter);

	friend StdInteractionManagerOutput;
};

template <class InteractionT>
void InteractionManager<InteractionT>::addInteraction(unsigned i,
									    unsigned j,
									    std::shared_ptr<InteractionT> inter)
{
	if (!areInteracting(i, j)) {
		interactions.push_back(inter);

		// tell i and j their new partner and interaction
		particle_partners[i].push_back(j);
		particle_partners[j].push_back(i);
		interactions_pp[i].push_back(inter);
		interactions_pp[j].push_back(inter);
	}
}

template <class InteractionT>
void InteractionManager<InteractionT>::removeInteraction(unsigned inter_idx)
{
	unsigned i, j;
	std::tie(i, j) = interactions[inter_idx]->get_par_num();

	removePartners(i, j);
	removeInteractionPP(i, interactions[inter_idx]);
	interactions[inter_idx] = interactions[interactions.size()-1];
	interactions.pop_back();
}

template <class InteractionT>
void InteractionManager<InteractionT>::removePartners(unsigned i, unsigned j)
{
	auto &neighi = particle_partners[i];
	auto &neighj = particle_partners[j];
	auto l = neighi.back();
	if (l != j) {
		for (unsigned k=0; k<neighi.size(); k++) {
			if (neighi[k] == j) {
				neighi[k] = l;
				break;
			}
		}
	}
	neighi.pop_back();
	l = neighj.back();
	if (l != i) {
		for (unsigned k=0; k<neighj.size(); k++) {
			if (neighj[k] == i) {
				neighj[k] = l;
				break;
			}
		}
	}
	neighj.pop_back();
}

template <class InteractionT>
bool InteractionManager<InteractionT>::areInteracting(unsigned i, unsigned j) const
{
	for (auto k : particle_partners[i]) {
		if (j == k) {
			return true;
		}
	}
	return false;
}

template <class InteractionT>
void InteractionManager<InteractionT>::removeInteractionPP(unsigned i, std::shared_ptr<InteractionT> dead_inter)
{
	auto last_inter = interactions_pp[i].back();
	if (last_inter != dead_inter) {
		for (auto &inter: interactions_pp[i]) {
			if (inter == dead_inter) {
				inter = last_inter;
				break;
			}
		}
	}
	interactions_pp[i].pop_back();
}
}
#endif