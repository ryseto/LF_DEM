#include "Dimer.h"

namespace Interactions {

namespace Dimer {

Spring::Spring(SpringState state, double stiffness) :
k(stiffness),
relaxed_length(state.relaxed_length),
saved_relaxed_length(state.relaxed_length),
stretch(state.stretch),
saved_stretch(state.stretch)
{}

void Spring::saveState()
{
	saved_stretch = stretch;
	saved_relaxed_length = relaxed_length;
}

void Spring::restoreState()
{
	stretch = saved_stretch;
	relaxed_length = saved_relaxed_length;
}

vec3d Spring::getForce() const
{
	return k*stretch;
}

vec3d Spring::getStretch() const
{
	return stretch;
}

void Spring::setStretch(const vec3d &new_stretch)
{
	stretch = new_stretch;
}


double Spring::getRelaxedLength() const
{
	return relaxed_length;
}

void Spring::setRelaxedLength(double l)
{
	relaxed_length = l;
}

void Spring::incrementStretch(const vec3d &delta_stretch)
{
	stretch += delta_stretch;
}

Dashpot::Dashpot(double resistance, Dimer *parent) : 
res(resistance),
rotres(resistance),
dimer(parent)
{}


struct ODBlock Dashpot::RFU_ODBlock() const
{
	ResistanceBlocks::ODBlockBuilder builder;

	builder.addIdentity(ResistanceBlocks::Label::ForceVel, -res);
	builder.addVectorProduct(ResistanceBlocks::Label::ForceAngVel, res*dimer->rvec/2.);
	builder.addVectorProduct(ResistanceBlocks::Label::TorqueVel, -res*dimer->rvec/2.);
	builder.addComplementaryDyadic(ResistanceBlocks::Label::TorqueAngVel, dimer->rvec/2., -res);

	// rotational velocity damping
	builder.addIdentity(ResistanceBlocks::Label::TorqueAngVel, rotres);

	return builder.block;
}

std::pair<struct DBlock, struct DBlock> Dashpot::RFU_DBlocks() const
{
	ResistanceBlocks::DBlockBuilder builder;

	builder.addIdentity(ResistanceBlocks::Label::ForceVel, res);
	builder.addComplementaryDyadic(ResistanceBlocks::Label::TorqueAngVel, dimer->rvec/2., res);
	builder.addVectorProduct(ResistanceBlocks::Label::TorqueVel, res*dimer->rvec/2.);

	// rotational velocity damping
	builder.addIdentity(ResistanceBlocks::Label::TorqueAngVel, -rotres);

	return std::make_pair(builder.block, builder.block);
}

std::pair<vec3d, vec3d> Dashpot::getForceTorque(const struct PairVelocity &vel) const
{
	vec3d force = -res*dimer->getContactVelocity(vel);
	vec3d torque = cross(dimer->rvec/2, force);

	torque += -rotres*dimer->getRotationDiffVelocity(vel);

	return std::make_pair(force, torque);
}

void Dimer::checkHomegeneity()
{
	if (a0 != a1) {
		throw std::runtime_error("Only homogeneous dimers are implemented");
	}
}

Dimer::Dimer(const PairId &pairid, vec3d sep, const struct UnloadedDimerState &uds, DimerParams p) :
PairwiseInteraction(pairid, sep),
dashpot(p.resistance, this),
sliding_spring({getSeparation() - uds.relaxed_length*getUnitSeparation(), uds.relaxed_length}, p.stiffness),
rotation_spring({vec3d(), 0}, p.stiffness)
{
	checkHomegeneity();
}


Dimer::Dimer(const PairId &pairid, vec3d sep, const struct DimerState &ds, DimerParams p) :
PairwiseInteraction(pairid, sep),
dashpot(p.resistance, this),
sliding_spring(ds.sliding_st, p.stiffness),
rotation_spring(ds.rotation_st, p.stiffness)
{
	checkHomegeneity();
}

struct DimerState Dimer::getState() const
{
	struct DimerState s;
	s.p0 = p0;
	s.p1 = p1;
	s.sliding_st.stretch = sliding_spring.getStretch();
	s.sliding_st.relaxed_length = rotation_spring.getRelaxedLength();
	s.rotation_st.stretch = sliding_spring.getStretch();
	s.rotation_st.relaxed_length = rotation_spring.getRelaxedLength();
	return s;
}

void Dimer::saveState()
{
	sliding_spring.saveState();
	rotation_spring.saveState();
}

void Dimer::restoreState()
{
	sliding_spring.restoreState();
	rotation_spring.restoreState();
}


vec3d Dimer::getContactVelocity(const struct PairVelocity &vel)
{
	return vel.U[1] - vel.U[0] - cross(vel.O[0] + vel.O[1], rvec/2.);
}

vec3d Dimer::getRotationDiffVelocity(const struct PairVelocity &vel)
{
	return vel.O[1]-vel.O[0];
}

void Dimer::applyTimeStep(double dt, const struct PairVelocity &vel)
{
	sliding_spring.incrementStretch(getContactVelocity(vel)*dt);
	rotation_spring.incrementStretch(getRotationDiffVelocity(vel)*dt);
}

std::pair<vec3d, vec3d> Dimer::getForceTorque(const struct PairVelocity &vel) const
{
	std::pair<vec3d, vec3d> fts = getForceTorqueSpring();
	std::pair<vec3d, vec3d> ftd = getForceTorqueDashpot(vel);

	return std::make_pair(fts.first + ftd.first, fts.second + ftd.second);
}


std::pair<vec3d, vec3d> Dimer::getForceTorqueSpring() const
{
	//Sliding spring
	vec3d force = sliding_spring.getForce();
	vec3d torque = cross(rvec/2, force);

	//Rotation spring
	torque += rotation_spring.getForce();

	return std::make_pair(force, torque);
}

std::pair<vec3d, vec3d> Dimer::getForceTorqueDashpot(const struct PairVelocity &vel) const
{
	return dashpot.getForceTorque(vel);
}

Sym2Tensor Dimer::getStress(const struct PairVelocity &vel) const
{
	auto ft = getForceTorque(vel);
	return outer_sym(rvec, ft.first);
}

Sym2Tensor Dimer::getDashpotStress(const struct PairVelocity &vel) const
{
	auto ft = getForceTorqueDashpot(vel);
	return outer_sym(rvec, ft.first);
}

Sym2Tensor Dimer::getSpringStress() const
{
	auto ft = getForceTorqueSpring();
	return outer_sym(rvec, ft.first);
}

} // namespace Dimer

} // namespace Interactions
