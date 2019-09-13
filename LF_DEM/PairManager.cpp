#include "PairManager.h"

namespace Interactions {

void PairManager::addPartners(unsigned i, unsigned j)
{
	particle_partners[i].push_back(j);
	particle_partners[j].push_back(i);
}

void PairManager::removePartners(unsigned i, unsigned j)
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

bool PairManager::areInteracting(unsigned i, unsigned j) const
{
	for (auto k : particle_partners[i]) {
		if (j == k) {
			return true;
		}
	}
	return false;
}

} // namespace Interactions
