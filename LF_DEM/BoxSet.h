//
//  BoxSet.h
//  LF_DEM
//
//  Copyright (c) 2013-2015 Ryohei Seto and Romain Mari. All rights reserved.
//

/**
 \class BoxSet
 \brief Set of Box objects making a partition of the simulation box.
 \author Ryohei Seto
 \author Romain Mari
 */

#ifndef __LF_DEM__BoxSet__
#define __LF_DEM__BoxSet__
#include <set>
#include <vector>
#include <fstream>
#include "vec3d.h"
// #include "Geometry.h"
#include "Box.h"

namespace Boxing {

class BoxSet {
public:
	virtual void box(int i, vec3d position) = 0;
	/*****
	 neighborhood_begin(int i) and neighborhood_end(int i)
	 gives iterators to beginning and ending point of the container including
	 all particles in the box containing particle i and in the adjacent boxes.
	 *****/
	const std::vector <int>& neighborhood(int i) const;

	~BoxSet() {};
	void printBoxNetwork() const;
	void printBoxContainers() const;
	void printNeighborhoodContainers() const;
	void printBoxMap() const;

protected:
	double box_xsize;
	double box_ysize;
	double box_zsize;
	std::size_t x_box_nb;
	std::size_t y_box_nb;
	std::size_t z_box_nb;
	std::size_t yz_box_nb;
	std::size_t box_nb;
	std::set <Box*> Boxes;
	std::vector <Box*> box_labels;
	std::vector <Box*> boxMap;

	bool _is_boxed;
	bool is_boxed() const;

	Box* whichBox(const vec3d&) const;
	Box* whichBox(unsigned int box_label) const;
	unsigned int whichBoxLabel(const vec3d&) const;
	void addToBox(int i, Box* b);
	void allocateBoxes();

};

// class OpenSystemBoxSet : public BoxSet {
// // TBD
// };

}

#endif /* defined(__LF_DEM__BoxSet__) */
