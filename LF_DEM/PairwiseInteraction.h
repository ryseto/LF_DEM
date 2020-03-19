//
//  PairwiseInteraction.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/10/12.
//  Copyright (c) 2019 Ryohei Seto and Romain Mari. All rights reserved.
//

#ifndef __LF_DEM__PairwiseInteraction__
#define __LF_DEM__PairwiseInteraction__

#include "vec3d.h"
#include "Sym2Tensor.h"
#include "PairVelocity.h"


namespace Geometry {
class PairwiseConfig;
}

namespace Interactions
{
class InteractionSet;
class Contact;
class ContactDashpot;
namespace Lub {
class Lubrication;
}
class PotentialForce;
class RepulsiveForce;
class vanDerWaalsForce;
namespace TActAdhesion {
class TimeActivatedAdhesion;
}
namespace ActAdhesion {
class ActivatedAdhesion;
}

struct PairId {
	unsigned p0; 
	unsigned p1;
	double a0;
	double a1; 
};

class PairwiseInteraction {

protected:
	unsigned p0;
	unsigned p1;
	double a0;
	double a1;

	//======= relative position/velocity data  =========//
	double reduced_gap; // gap between particles (dimensionless gap = s - 2, s = 2r/(a1+a2) )
	double r; // center-center distance
	vec3d rvec; // vector center to center
	vec3d nvec; // normal vector


	void setSeparation(const vec3d &sep);

	double gap2separation(double gap);  // gap is in units of (a0+a1)/2, separation is in std unit length
	double separation2gap(double sep);  // gap is in units of (a0+a1)/2, separation is in std unit length

	friend Contact;
	friend ContactDashpot;
	friend Lub::Lubrication;
	friend PotentialForce;
	friend RepulsiveForce;
	friend vanDerWaalsForce;
	friend TActAdhesion::TimeActivatedAdhesion;
	friend ActAdhesion::ActivatedAdhesion;
	friend InteractionSet;

public:
	PairwiseInteraction(const PairId &data, vec3d sep);

	virtual void saveState() = 0;
	virtual void restoreState() = 0;


	unsigned partner(unsigned i) const {return (i == p0 ? p1 : p0);}
	std::pair<unsigned, unsigned> get_par_num() const {return std::make_pair(p0, p1);}
	unsigned get_p0() const {return p0;}
	unsigned get_p1() const {return p1;}
	double getReducedGap() const {return reduced_gap;}
	double getGap() const {return r-a0-a1;}
	vec3d getSeparation() const {return rvec;}
	vec3d getUnitSeparation() const {return nvec;}
	double getNormalVelocity(const struct PairVelocity &vel) const;

};

} // namespace Interactions

#endif /* defined(__LF_DEM__PairwiseInteraction__) */
