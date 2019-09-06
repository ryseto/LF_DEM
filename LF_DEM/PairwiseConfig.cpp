#include "PairwiseConfig.h"
#include "SimpleShearBoxSet.h"
#include "ExtensionalShearBoxSet.h"
#include "LeesEdwards.h"
#include "KraynikReinelt.h"

namespace Geometry {

PairwiseConfig::PairwiseConfig(std::shared_ptr<ParticleConfig> config) :
conf(config),
velo(nullptr)
// backgroung_velo(velinf)
{}

void PairwiseConfig::setVelocityState(ParticleVelocity *vel)
{
	velo = vel;
}

LeesEdwardsPairwiseConfig::LeesEdwardsPairwiseConfig(std::shared_ptr<ParticleConfig> config,
													 std::shared_ptr<BC::LeesEdwardsBC> lebc,
													 double max_interaction_range) :
PairwiseConfig(config),
lees(lebc)
{
	boxset = new Boxing::SimpleShearBoxSet(max_interaction_range, 
										   conf->position.size(), 
										   *lees);
	// boxset = std::make_unique<Boxing::SimpleShearBoxSet>(max_interaction_range, 
	// 													 conf->position.size(), 
	// 													 *lees);
	updateAfterParticleMove();
	updateAfterDeformation();
}

LeesEdwardsPairwiseConfig::~LeesEdwardsPairwiseConfig()
{
	delete boxset;
}

void LeesEdwardsPairwiseConfig::updateAfterParticleMove()
{
	for (unsigned i=0; i<conf->position.size(); i++) {
		updateAfterParticleMove(i);
	}
} 

void LeesEdwardsPairwiseConfig::updateAfterParticleMove(unsigned i)
{
	boxset->box(i, conf->position[i]);
} 

void LeesEdwardsPairwiseConfig::updateAfterDeformation()
{
	boxset->update(*lees);
} 

const std::vector <int>& LeesEdwardsPairwiseConfig::neighborhood(unsigned i) const
{
	return boxset->neighborhood(i);
}

vec3d LeesEdwardsPairwiseConfig::getSeparation(unsigned i, unsigned j) const
{
	vec3d rvec = conf->position[j]-conf->position[i];
	lees->periodizeDiff(rvec);
	return rvec;
}

vec3d LeesEdwardsPairwiseConfig::getVelOffset(unsigned i, unsigned j) const
{
	vec3d rvec = conf->position[j]-conf->position[i];
	vec3d rvec_cp = rvec;
	lees->periodizeDiff(rvec);
	return lees->getVelDifference(rvec - rvec_cp);
}

void LeesEdwardsPairwiseConfig::getVelocities(unsigned i, unsigned j, struct PairVelocity &k) const
{
	vec3d pbc_vel_offset;
	if (velo->type == VelocityType::total) {
		pbc_vel_offset = getVelOffset(i, j);
	}
	k.U[0] = velo->vel[i];
	k.U[1] = velo->vel[j] + pbc_vel_offset;
	k.O[0] = velo->ang_vel[i];
	k.O[1] = velo->ang_vel[j];
}

KraynikReineltPairwiseConfig::KraynikReineltPairwiseConfig(std::shared_ptr<ParticleConfig> config,
															std::shared_ptr<BC::KraynikReineltBC> krbc,
															double max_interaction_range) :
PairwiseConfig(config),
kr(krbc),
twodimension(krbc->getContainer().ly == 0)
{
	// extensional flow
	auto container = kr->getContainer();
	double dl = max_interaction_range;
	int num_x = (int)(container.lx/dl);
	int num_z = (int)(container.lz/dl);
	container_ext_flow = { 3*dl*(num_x+1)+2*dl,
							container.ly,
							2*dl*(num_z+1)+2*dl };
	vec3d box_origin(dl*(num_x+2), 0, dl*(num_z+2));
	boxset = new Boxing::ExtensionalShearBoxSet(max_interaction_range, 
												box_origin, conf->position.size(), 
												container_ext_flow, 
												container.ly != 0);
	// boxset = std::make_unique<Boxing::ExtensionalShearBoxSet>(max_interaction_range, 
	// 														  box_origin, conf->position.size(), 
	// 														  container_ext_flow, 
	// 														  container.ly != 0);

	krbc->setBoxSet(boxset, &container_ext_flow);
	updateAfterParticleMove();
	updateAfterDeformation();
}

KraynikReineltPairwiseConfig::~KraynikReineltPairwiseConfig()
{
	delete boxset;
}

void KraynikReineltPairwiseConfig::updateAfterParticleMove()
{
	for (unsigned i=0; i<conf->position.size(); i++) {
		updateAfterParticleMove(i);
	}
} 

void KraynikReineltPairwiseConfig::updateAfterParticleMove(unsigned i)
{
	boxset->box(i, conf->position[i]);
} 

void KraynikReineltPairwiseConfig::updateAfterDeformation()
{
	boxset->updateExtFlow(container_ext_flow, kr->ext_ax);	
} 

const std::vector <int>& KraynikReineltPairwiseConfig::neighborhood(unsigned i) const
{
	return boxset->neighborhood(i);
}

vec3d KraynikReineltPairwiseConfig::getSeparation(unsigned i, unsigned j) const
{
	vec3d rvec = conf->position[j]-conf->position[i];
	periodizeDiffExtFlow(rvec, i, j);
	return rvec;
}

vec3d KraynikReineltPairwiseConfig::getVelOffset(unsigned i, unsigned j) const
{
	vec3d rvec = conf->position[j]-conf->position[i];
	return kr->getVelDifference(pdShiftExtFlow(rvec, i, j));
}

vec3d KraynikReineltPairwiseConfig::pdShiftExtFlow(vec3d pos_diff, unsigned i, unsigned j) const
{
	auto container = kr->getContainer();
	auto pd_shift = boxset->periodicDiffShift(i, j, kr->ext_ax);
	if (twodimension == false) {
		if (pos_diff.y > 0.5*container.ly) {
			pd_shift.y -= container.ly;
		} else if (pos_diff.y < -0.5*container.ly) {
			pd_shift.y += container.ly;
		}
	}
	return pd_shift;
}

void KraynikReineltPairwiseConfig::periodizeDiffExtFlow(vec3d& pos_diff, unsigned i, unsigned j) const
{
	auto pd_shift = pdShiftExtFlow(pos_diff, i, j);
	if (pd_shift.is_nonzero()) {
		pos_diff += pd_shift;
	}
}

void KraynikReineltPairwiseConfig::getVelocities(unsigned i, unsigned j, struct PairVelocity &k) const
{
	vec3d pbc_vel_offset;
	if (velo->type == VelocityType::total) {
		pbc_vel_offset = getVelOffset(i, j);
	}
	k.U[0] = velo->vel[i];
	k.U[1] = velo->vel[j] + pbc_vel_offset;
	k.O[0] = velo->ang_vel[i];
	k.O[1] = velo->ang_vel[j];
}

} // namespace Geometry