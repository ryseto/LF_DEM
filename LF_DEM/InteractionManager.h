#ifndef __LF_DEM__InteractionManager__
#define __LF_DEM__InteractionManager__

#include <memory>
#include <vector>
#include <string>

#include "ForceComponent.h"
#include "PairManager.h"

class StdInteractionManagerOutput;

namespace Interactions
{

template <class InteractionT>
class InteractionManager {
public:
	InteractionManager(unsigned np, std::shared_ptr<PairManager> &pairmanager) :
		interactions_pp(np),
		pair_manager(pairmanager) {};
	virtual void declareForceComponents(std::map<std::string, ForceComponent> &force_components) = 0;
	virtual void setForceToParticle(const std::string &component, std::vector<vec3d> &force, std::vector<vec3d> &torque) = 0;

	std::vector< std::string > declared_forces;
	std::vector< std::shared_ptr<InteractionT> > interactions;
	std::vector< std::vector< std::shared_ptr<InteractionT> > > interactions_pp; // per particles

	std::size_t size() const {return interactions.size();};

	typename std::vector< std::shared_ptr<InteractionT> >::iterator begin() {return interactions.begin();}
	typename std::vector< std::shared_ptr<InteractionT> >::iterator end() {return interactions.end();}
	typename std::vector< std::shared_ptr<InteractionT> >::const_iterator begin() const {return interactions.begin();}
	typename std::vector< std::shared_ptr<InteractionT> >::const_iterator end() const {return interactions.end();}

protected:
	void removeInteraction(unsigned inter_idx);
	void addInteraction(unsigned i, unsigned j, std::shared_ptr<InteractionT> inter);
	std::shared_ptr<PairManager> pair_manager;

private:
	void removeInteractionPP(unsigned i, std::shared_ptr<InteractionT> dead_inter);
	friend StdInteractionManagerOutput;
};

template <class InteractionT>
void InteractionManager<InteractionT>::addInteraction(unsigned i,
                                                      unsigned j,
                                                      std::shared_ptr<InteractionT> inter)
{
    std::cout << "Testttttttttttttt = (0) " << std::endl;
	if (!pair_manager->areInteracting(i, j)) {
		interactions.push_back(inter);

		// tell i and j their new partner and interaction
		pair_manager->addPartners(i, j);
		interactions_pp[i].push_back(inter);
		interactions_pp[j].push_back(inter);
	}
}

template <class InteractionT>
void InteractionManager<InteractionT>::removeInteraction(unsigned inter_idx)
{
	unsigned i, j;
	std::tie(i, j) = interactions[inter_idx]->get_par_num();

	pair_manager->removePartners(i, j);
	removeInteractionPP(i, interactions[inter_idx]);
	interactions[inter_idx] = interactions[interactions.size()-1];
	interactions.pop_back();
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
