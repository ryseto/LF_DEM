//
//  BoxSet.h
//  LF_DEM
//
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class SimpleShearBoxSet
 \brief Set of Box objects making a partition of the simulation box for Lees-Edwards.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__SimpleShearBoxSet__
#define __LF_DEM__SimpleShearBoxSet__
#include "BoxSet.h"

namespace BC {
	class LeesEdwardsBC;
}

namespace Boxing {


class SimpleShearBoxSet : public BoxSet {
private:
	std::set <Box*> BulkBoxes;
	std::set <Box*> TopBoxes;
	std::set <Box*> BottomBoxes;
	std::set <Box*> TopBottomBoxes;
	/*****
	 WhichBox(vec3d pos)
	 returns a pointer on the box containg position pos
	 *****/
	
	void updateNeighbors(const BC::LeesEdwardsBC &lees);
	// init methods
	void positionBoxes();
	void assignNeighbors(const BC::LeesEdwardsBC &lees);
	void assignNeighborsBulk(const BC::LeesEdwardsBC &lees);
	void assignNeighborsTop(const BC::LeesEdwardsBC &lees);
	void assignNeighborsBottom(const BC::LeesEdwardsBC &lees);
	void assignNeighborsTopBottom(const BC::LeesEdwardsBC &lees);
	std::vector<vec3d> top_probing_positions;
	std::vector<vec3d> bottom_probing_positions;
public:
	SimpleShearBoxSet(double interaction_dist, unsigned np, const BC::LeesEdwardsBC &lees);
	/*****
	 update()
	 
	 To be called at each time step.
	 It updates the neighborhood relations betwenn boxes.
	 Those relations change at each time step for boxes on top or bottom
	 *****/
	void update(const BC::LeesEdwardsBC &lees);
	/*****
	 is_boxed()
	 
	 Can be called before calling an other method of BoxSet.
	 If false, than calls to other method may usually be avoided.
	 They can be performed anyway though, and they are normally safe (and useless)
	 
	 is_boxed() tells if the boxing is effective.
	 If the system size is small, the neighborhood of a box may contain
	 the whole system. If this happens, the boxing consists of only one box (as it
	 is useless, if not ill defined, to do something else), and is_boxed() returns false.
	 
	 *****/
	
	void box(int i, vec3d position);
	
};

}

#endif /* defined(__LF_DEM__SimpleShearBoxSet__) */
