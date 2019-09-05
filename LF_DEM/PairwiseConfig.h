#ifndef __LF_DEM__PairwiseDistance__
#define __LF_DEM__PairwiseDistance__

#include <memory>
#include <vector>
#include "ParticleConfig.h"
#include "PairVelocity.h"
#include "Box3d.h"
#include "vec3d.h"

namespace BC {
	class LeesEdwardsBC;
	class KraynikReineltBC;
}

namespace Boxing {
	// class BoxSet;
	class SimpleShearBoxSet;
	class ExtensionalShearBoxSet;
}
namespace Geometry {

class PairwiseConfig {
public:
	virtual void updateAfterParticleMove() = 0;
	virtual void updateAfterParticleMove(unsigned i) = 0;
	virtual void updateAfterDeformation() = 0;

	virtual vec3d getSeparation(unsigned i, unsigned j) const = 0;
	virtual const std::vector <int>& neighborhood(unsigned i) const = 0;

	virtual void getVelocities(unsigned i, unsigned j, struct PairVelocity &k) const = 0;
	void setVelocityState(ParticleVelocity *vel);
protected:
	PairwiseConfig(std::shared_ptr<ParticleConfig> config); 
	std::shared_ptr<ParticleConfig> conf;
	ParticleVelocity *velo;
};

class LeesEdwardsPairwiseConfig : public PairwiseConfig {
public:
	LeesEdwardsPairwiseConfig(std::shared_ptr<ParticleConfig> config,
							  std::shared_ptr<BC::LeesEdwardsBC> lebc,
							  double max_interaction_range);
	~LeesEdwardsPairwiseConfig();
	void updateAfterParticleMove();
	void updateAfterParticleMove(unsigned i);
	void updateAfterDeformation();

	vec3d getSeparation(unsigned i, unsigned j) const;
	const std::vector <int>& neighborhood(unsigned i) const;
	void getVelocities(unsigned i, unsigned j, struct PairVelocity &k) const;

private:
	std::shared_ptr<BC::LeesEdwardsBC> lees;
	Boxing::SimpleShearBoxSet *boxset;
	vec3d getVelOffset(unsigned i, unsigned j) const;
};

class KraynikReineltPairwiseConfig : public PairwiseConfig {
public:
	KraynikReineltPairwiseConfig(std::shared_ptr<ParticleConfig> config, 
								 std::shared_ptr<BC::KraynikReineltBC> krbc,
								 double max_interaction_range);
	~KraynikReineltPairwiseConfig();
	void updateAfterParticleMove();
	void updateAfterParticleMove(unsigned i);
	void updateAfterDeformation();

	vec3d getSeparation(unsigned i, unsigned j) const;
	const std::vector <int>& neighborhood(unsigned i) const;
	void getVelocities(unsigned i, unsigned j, struct PairVelocity &k) const;
	Boxing::ExtensionalShearBoxSet *boxset;

private:
	std::shared_ptr<BC::KraynikReineltBC> kr;
	struct box3d container_ext_flow;
	bool twodimension;

	vec3d getVelOffset(unsigned i, unsigned j) const;
	void periodizeExtFlow(unsigned i, bool &pd_transport) const;
	vec3d pdShiftExtFlow(vec3d pos_diff, unsigned i, unsigned j) const;
	void periodizeDiffExtFlow(vec3d& pos_diff, unsigned i, unsigned j) const;
};


// class OpenSystemPairwiseDistance : public PairwiseDistance {
// // TBD
// };

} // namespace Geometry
#endif