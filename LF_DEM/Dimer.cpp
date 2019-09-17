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
	// !!! Block (p1, p0)
	// Force and torque on p1 from motion of p0
	ResistanceBlocks::ODBlockBuilder builder;

	builder.addIdentity(ResistanceBlocks::Label::ForceVel, -res);
	builder.addVectorProduct(ResistanceBlocks::Label::ForceAngVel, res*dimer->rvec/2.);
	builder.addVectorProduct(ResistanceBlocks::Label::TorqueVel, res*dimer->rvec/2.);
	builder.addComplementaryDyadic(ResistanceBlocks::Label::TorqueAngVel, dimer->rvec/2., -res);

	// rotational velocity damping
	builder.addIdentity(ResistanceBlocks::Label::TorqueAngVel, rotres);

	return builder.block;
}

std::pair<struct DBlock, struct DBlock> Dashpot::RFU_DBlocks() const
{
	ResistanceBlocks::DBlockBuilder builder0, builder1;

	builder0.addIdentity(ResistanceBlocks::Label::ForceVel, res);
	builder0.addComplementaryDyadic(ResistanceBlocks::Label::TorqueAngVel, dimer->rvec/2., res);
	builder0.addVectorProduct(ResistanceBlocks::Label::TorqueVel, res*dimer->rvec/2.);

	// rotational velocity damping
	builder0.addIdentity(ResistanceBlocks::Label::TorqueAngVel, rotres);

	builder1.addIdentity(ResistanceBlocks::Label::ForceVel, res);
	builder1.addComplementaryDyadic(ResistanceBlocks::Label::TorqueAngVel, dimer->rvec/2., res);
	builder1.addVectorProduct(ResistanceBlocks::Label::TorqueVel, -res*dimer->rvec/2.);

	// rotational velocity damping
	builder1.addIdentity(ResistanceBlocks::Label::TorqueAngVel, rotres);

	return std::make_pair(builder0.block, builder1.block);
}

std::tuple<vec3d, vec3d, vec3d, vec3d> Dashpot::getForceTorque(const struct PairVelocity &vel) const
{
	vec3d force_sliding = res*dimer->getContactVelocity(vel);
	vec3d torque_sliding = cross(dimer->rvec/2, force_sliding);

	vec3d torque_rotation = rotres*dimer->getRotationDiffVelocity(vel);

	return std::make_tuple(force_sliding, -force_sliding, torque_sliding + torque_rotation, torque_sliding - torque_rotation);
}

void Dimer::checkHomegeneity()
{
	if (a0 != a1) {
		throw std::runtime_error("Only homogeneous dimers are implemented");
	}
}

Dimer::Dimer(const PairId &pairid, vec3d sep, const struct UnloadedDimerState &uds, DimerParams p) :
PairwiseInteraction(pairid, sep),
dashpot(p.relaxation_time*p.stiffness, this),
sliding_spring({getSeparation() - uds.relaxed_length*getUnitSeparation(), uds.relaxed_length}, p.stiffness),
rotation_spring({vec3d(), 0}, p.stiffness)
{
	checkHomegeneity();
}


Dimer::Dimer(const PairId &pairid, vec3d sep, const struct DimerState &ds, DimerParams p) :
PairwiseInteraction(pairid, sep),
dashpot(p.relaxation_time*p.stiffness, this),
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

void Dimer::applyTimeStep(double dt, vec3d sep, const struct PairVelocity &vel)
{
	setSeparation(sep);
	sliding_spring.incrementStretch(getContactVelocity(vel)*dt);
	rotation_spring.incrementStretch(getRotationDiffVelocity(vel)*dt);
}

std::tuple<vec3d, vec3d, vec3d, vec3d> Dimer::getForceTorque(const struct PairVelocity &vel) const
{
	auto fts = getForceTorqueSpring();
	auto ftd = getForceTorqueDashpot(vel);

	return std::make_tuple(std::get<0>(fts) + std::get<0>(ftd), 
						   std::get<1>(fts) + std::get<1>(ftd), 
						   std::get<2>(fts) + std::get<2>(ftd), 
						   std::get<3>(fts) + std::get<3>(ftd));
}


std::tuple<vec3d, vec3d, vec3d, vec3d> Dimer::getForceTorqueSpring() const
{
	//Sliding spring
	vec3d force_sliding = sliding_spring.getForce();
	vec3d torque_sliding = cross(rvec/2, force_sliding);

	//Rotation spring
	vec3d torque_rotation = rotation_spring.getForce();

	return std::make_tuple(force_sliding, -force_sliding, torque_sliding+torque_rotation, torque_sliding-torque_rotation);
}

std::tuple<vec3d, vec3d, vec3d, vec3d> Dimer::getForceTorqueDashpot(const struct PairVelocity &vel) const
{
	return dashpot.getForceTorque(vel);
}

Sym2Tensor Dimer::getStress(const struct PairVelocity &vel) const
{
	auto ft = getForceTorque(vel);
	return outer_sym(rvec, std::get<0>(ft));
}

Sym2Tensor Dimer::getDashpotStress(const struct PairVelocity &vel) const
{
	auto ft = getForceTorqueDashpot(vel);
	return outer_sym(rvec, std::get<0>(ft));
}

Sym2Tensor Dimer::getSpringStress() const
{
	auto ft = getForceTorqueSpring();
	return outer_sym(rvec, std::get<0>(ft));
}

struct ODBlock Dimer::RFU_ODBlock() const
{
	return dashpot.RFU_ODBlock();
}

std::pair<struct DBlock, struct DBlock> Dimer::RFU_DBlocks() const
{
	return dashpot.RFU_DBlocks();
}

} // namespace Dimer

} // namespace Interactions
