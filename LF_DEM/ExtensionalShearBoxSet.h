//
//  BoxSet.h
//  LF_DEM
//
//  Copyright (c) 2013-2019 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class ExtensionalShearBoxSet
 \brief Set of Box objects making a partition of the simulation box for Kraynik-Reinelt boundary conditions.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__ExtensionalShearBoxSet__
#define __LF_DEM__ExtensionalShearBoxSet__

#include "BoxSet.h"
#include "KraynikReinelt.h"

namespace Boxing {

class ExtensionalShearBoxSet : public BoxSet {
private:
	int nb_movingneighbor_part;
	/*****
	 WhichBox(vec3d pos)
	 returns a pointer on the box containg position pos
	 *****/
	void updateNeighborsExtFlow(const struct Geometry::box3d &sysbox,
							    const struct BC::ExtFlowAxis &ax);
	// init methods
	void positionBoxes();
	void assignNeighborsStaticExtFlow(Geometry::box3d sysbox);
	void assignNeighobrsDynamicExtFlow(Box& bx, const vec3d &pos, const Geometry::box3d &sysbox);
	unsigned int labelRight(int unsigned label_) {
		return label_+yz_box_nb;
	}
	unsigned int labelLeft(int unsigned label_) {
		return label_-yz_box_nb;
	}
	unsigned int labelUp(int unsigned label_) {
		return label_+1;
	}
	unsigned int labelDown(int unsigned label_) {
		return label_-1;
	}
	std::vector<int> next_labels;
	std::vector < std::vector<int> > next_of_periodicimage;
	vec3d origin_ext_flow;
	std::vector<vec3d> box_corners;
	bool twodimension;
	
public:
	ExtensionalShearBoxSet(double interaction_dist,
						   const vec3d origin,
						   unsigned np, 
						   Geometry::box3d sysbox,
						   bool twodimension_);
	/*****
	 update()
	 
	 To be called at each time step.
	 It updates the neighborhood relations betwenn boxes.
	 Those relations change at each time step for boxes on top or bottom
	 *****/
	void updateExtFlow(const struct Geometry::box3d &sysbox,
					   const struct BC::ExtFlowAxis &ax);
	void box(int i, vec3d position);
	void yaplotBox(std::ofstream &fout_boxing);
	void checkNeighobrType(std::set<unsigned int> &box_labels);
	int boxType(int i)
	{
		return boxMap[i]->type;
	}
	int boxNType(int i)
	{
		return boxMap[i]->type_neighborhood;
	}
	vec3d periodicDiffShift(int i, int j, const struct BC::ExtFlowAxis &ax);
};

}

#endif /* defined(__LF_DEM__ExtensionalShearBoxSet__) */
