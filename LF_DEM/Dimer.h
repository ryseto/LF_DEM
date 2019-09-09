#ifndef __LF_DEM__Dimer__
#define __LF_DEM__Dimer__

#include <vector>
#include "vec3d.h"
#include "MatrixBlocks.h"
#include "PairwiseInteraction.h"
#include "Sym2Tensor.h"
#include "PairVelocity.h"


namespace Interactions {

namespace Dimer {

class Dimer;

struct DimerState {
	unsigned int p0;
	unsigned int p1;
	vec3d sliding_stretch;
	vec3d rotation_stretch;
};

namespace io 
{
std::vector <struct DimerState> readStatesBStream(std::istream &input);
void writeStatesBStream(std::ostream &conf_export, const std::vector <struct DimerState> &ds);
}

class Spring {
public:
	Spring(double stiffness);

	void saveState();
	void restoreState();
	vec3d getForce() const;
	vec3d getStretch() const;
	void setStretch(const vec3d &new_stretch);
	void incrementStretch(const vec3d &delta_stretch);

private:
	double k;
	vec3d stretch;
	vec3d saved_stretch;
};

class Dashpot {
public:
	Dashpot(double resistance, Dimer *dimer);
	std::pair<vec3d, vec3d> getForceTorque(const struct PairVelocity &vel) const;
	struct ODBlock RFU_ODBlock() const;
	std::pair<struct DBlock, struct DBlock> RFU_DBlocks() const;
private:
	double res;
	double rotres;
	Dimer *dimer;
};

class Dimer : public PairwiseInteraction {
private:
	Spring sliding_spring;
	Spring rotation_spring;
	Dashpot dashpot;

protected:
	vec3d getContactVelocity(const struct PairVelocity &vel);
	vec3d getRotationDiffVelocity(const struct PairVelocity &vel);

	friend Dashpot;
	friend Spring;

public:
	Dimer(double rest_separation, double stiffness, double resistance);

	struct DimerState getState() const;
	void setState(const struct DimerState& ds);

	void saveState();
	void restoreState();
	void applyTimeStep(double dt, const struct PairVelocity &vel);

	std::pair<vec3d, vec3d> getForceTorque(const struct PairVelocity &vel) const;
	Sym2Tensor getStress(const struct PairVelocity &vel) const;
};

} // namespace Dimer

} // namespace Interactions
#endif /* defined(__LF_DEM__Dimer__) */

