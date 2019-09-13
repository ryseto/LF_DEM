#ifndef __LF_DEM__PairManager__
#define __LF_DEM__PairManager__

#include <vector>

namespace Interactions {

class PairManager {
public:
	PairManager(unsigned np) :
		particle_partners(np) {};
	bool areInteracting(unsigned i, unsigned j) const;
	void removePartners(unsigned i, unsigned j);
	void addPartners(unsigned i, unsigned j);

private:
	std::vector< std::vector<unsigned> >particle_partners;
};

} // namespace Interactions
#endif