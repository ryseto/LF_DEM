#ifndef __LF_DEM__Dimer__
#define __LF_DEM__Dimer__

#include <vector>
#include "vec3d.h"
#include "MatrixBlocks.h"
#include "PairwiseInteraction.h"
#include "Sym2Tensor.h"
#include "PairVelocity.h"
#include "DimerParams.h"
#include "DimerState.h"


namespace Interactions {

namespace Dimer {

class Dimer;
class DimerManager;

namespace io 
{
std::vector <struct DimerState> readStatesBStream(std::istream &input);
std::vector <struct UnloadedDimerState> readTxtDimer(const std::string& filename);
void writeStatesBStream(std::ostream &conf_export, const std::vector <struct DimerState> &ds);
}

class Spring {
public:
	Spring(SpringState state, double stiffness);

	void saveState();
	void restoreState();
	vec3d getForce() const;
	vec3d getStretch() const;
	double getRelaxedLength() const;
	void setStretch(const vec3d &new_stretch);
	void setRelaxedLength(double l);
	void incrementStretch(const vec3d &delta_stretch);

private:
	double k;
	double relaxed_length;
	double saved_relaxed_length;
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
	Dashpot dashpot;
	Spring sliding_spring;
	Spring rotation_spring;
	void checkHomegeneity();

protected:
	vec3d getContactVelocity(const struct PairVelocity &vel);
	vec3d getRotationDiffVelocity(const struct PairVelocity &vel);

	friend Dashpot;
	friend Spring;
	friend DimerManager;

public:
	Dimer(const PairId &pairid, vec3d sep, const struct UnloadedDimerState &uds, DimerParams p);
	Dimer(const PairId &pairid, vec3d sep, const struct DimerState &uds, DimerParams p);

	struct DimerState getState() const;
	// void setSpringState(const struct SpringState& ds);

	void saveState();
	void restoreState();
	void applyTimeStep(double dt, const struct PairVelocity &vel);

	std::pair<vec3d, vec3d> getForceTorque(const struct PairVelocity &vel) const;
	std::pair<vec3d, vec3d> getForceTorqueDashpot(const struct PairVelocity &vel) const;
	std::pair<vec3d, vec3d> getForceTorqueSpring() const;

	Sym2Tensor getStress(const struct PairVelocity &vel) const;
	Sym2Tensor getDashpotStress(const struct PairVelocity &vel) const;
	Sym2Tensor getSpringStress() const;

	double getRelaxedLength() const {return sliding_spring.getRelaxedLength();};

};

} // namespace Dimer

} // namespace Interactions
#endif /* defined(__LF_DEM__Dimer__) */

