#include "Dimer.h"

namespace Interactions {

namespace Dimer {

Spring::Spring(double stiffness) :
k(stiffness),
stretch(0),
saved_stretch(0)
{}

void Spring::saveState()
{
	saved_stretch = stretch;
}

void Spring::restoreState()
{
	stretch = saved_stretch;
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


struct DimerState Dimer::getState() const
{
	struct DimerState s;
	s.p0 = p0;
	s.p1 = p1;
	s.sliding_stretch = sliding_spring.getStretch();
	s.rotation_stretch = rotation_spring.getStretch();
	return s;
}

void Dimer::setState(const struct DimerState &ds)
{
	sliding_spring.setStretch(ds.sliding_stretch);
	rotation_spring.setStretch(ds.rotation_stretch);
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
	//Sliding spring
	vec3d force = sliding_spring.getForce();
	vec3d torque = cross(rvec/2, force);

	//Rotation spring
	torque += rotation_spring.getForce();

	//Dashpot
	auto dashpot_ft = dashpot.getForceTorque(vel);
	force += dashpot_ft.first;
	torque += dashpot_ft.second;

	return std::make_pair(force, torque);
}

Sym2Tensor Dimer::getStress(const struct PairVelocity &vel) const
{
	auto ft = getForceTorque(vel);
	return outer_sym(rvec, ft.first);
}

namespace io 
{

std::vector <struct DimerState> readStatesBStream(std::istream &input)
{
	unsigned ndimer;
	input.read((char*)&ndimer, sizeof(decltype(ndimer)));
	std::vector <struct DimerState> states;
	for (unsigned i=0; i<ndimer; i++) {
		unsigned p0, p1;
		double dt_x, dt_y, dt_z, dr_x, dr_y, dr_z;
		input.read((char*)&p0, sizeof(unsigned));
		input.read((char*)&p1, sizeof(unsigned));
		input.read((char*)&dt_x, sizeof(decltype(dt_x)));
		input.read((char*)&dt_y, sizeof(decltype(dt_y)));
		input.read((char*)&dt_z, sizeof(decltype(dt_z)));
		input.read((char*)&dr_x, sizeof(decltype(dr_x)));
		input.read((char*)&dr_y, sizeof(decltype(dr_y)));
		input.read((char*)&dr_z, sizeof(decltype(dr_z)));
		struct DimerState s;
		s.p0 = (int)p0;
		s.p1 = (int)p1;
		s.sliding_stretch = vec3d(dt_x, dt_y, dt_z);
		s.rotation_stretch = vec3d(dr_x, dr_y, dr_z);
		states.push_back(s);
	}
	return states;
}


void writeStatesBStream(std::ostream &conf_export, const std::vector <struct DimerState> &ds)
{
	unsigned ndimer = ds.size();
	conf_export.write((char*)&ndimer, sizeof(unsigned int));
	for (unsigned i=0; i<ds.size(); i++) {
		conf_export.write((char*)&(ds[i].p0), sizeof(unsigned int));
		conf_export.write((char*)&(ds[i].p1), sizeof(unsigned int));
		conf_export.write((char*)&(ds[i].sliding_stretch.x), sizeof(double));
		conf_export.write((char*)&(ds[i].sliding_stretch.y), sizeof(double));
		conf_export.write((char*)&(ds[i].sliding_stretch.z), sizeof(double));
		conf_export.write((char*)&(ds[i].rotation_stretch.x), sizeof(double));
		conf_export.write((char*)&(ds[i].rotation_stretch.y), sizeof(double));
		conf_export.write((char*)&(ds[i].rotation_stretch.z), sizeof(double));
	}
}

} //namespace io

} // namespace Dimer

} // namespace Interactions
